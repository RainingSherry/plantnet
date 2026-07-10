from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import time
import traceback
from pathlib import Path

import numpy as np

from .config import MODEL_RUNNER, PAPER_ROOT, PROJECT_ROOT, PYTHON_BIN, REMOTE_RESULT_ROOT, resolved_run_config, sha256_payload
from .io_utils import git_revision, require_disk_space, sha256_file, utc_now, verify_checksum_manifest, write_checksum_manifest, write_json
from .remote_store import RemoteStore
from .run_scVICAR import (
    canonicalize, cli_args, code_digests, run_monitored, runtime_identity,
    validate_raw_outputs,
)


STRESS_DATASETS = ("Blood_BoneMarrow", "Human_Pancreas_3", "PRJNA895163")
STRESS_VARIANTS = ("fixed", "topology_full")
CONTAMINATION_RATIOS = (0.0, 0.25, 0.5, 0.75, 1.0)
ESTIMATORS = ("current", "uniform_sample", "full")
STRESS_VERSION = "stress_v1"


def stress_code_digests() -> dict[str, str]:
    result = code_digests()
    for path in (
        Path(__file__).resolve(),
        PROJECT_ROOT / "papers/scVICAR/code/stress_orchestrate.py",
        PROJECT_ROOT / "papers/scVICAR/code/stress_analysis.py",
    ):
        result[str(path.relative_to(PROJECT_ROOT))] = sha256_file(path)
    return dict(sorted(result.items()))


def ratio_slug(value: float) -> str:
    return f"{int(round(100 * value)):03d}pct"


def save_topology_bundle(raw_dir: Path, run_dir: Path) -> None:
    names = {
        "indices": "neighbor_indices.npy",
        "base_probs": "neighbor_base_probs.npy",
        "similarity": "neighbor_similarity.npy",
        "distance": "neighbor_distance.npy",
        "edge_reliability": "edge_reliability.npy",
        "node_gate": "node_gate.npy",
        "labels": "labels.npy",
        "predicted": "eval_kmeans_known_k.npy",
        "mapped": "eval_kmeans_known_k_mapped.npy",
    }
    arrays = {key: np.load(raw_dir / name) for key, name in names.items() if (raw_dir / name).is_file()}
    np.savez_compressed(run_dir / "topology_diagnostics.npz", **arrays)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one frozen scVICAR graph-stress task")
    parser.add_argument("--dataset", required=True, choices=STRESS_DATASETS)
    parser.add_argument("--data-path", required=True, type=Path)
    parser.add_argument("--variant", required=True, choices=STRESS_VARIANTS)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--contamination", required=True, type=float, choices=CONTAMINATION_RATIOS)
    parser.add_argument("--estimator", default="current", choices=ESTIMATORS)
    parser.add_argument("--gpu", required=True, type=int, choices=range(1, 7))
    parser.add_argument("--staging-root", type=Path, default=PAPER_ROOT / ".staging" / "stress")
    parser.add_argument("--epochs", type=int)
    parser.add_argument("--upload", action="store_true")
    parser.add_argument("--cleanup-after-upload", action="store_true")
    args = parser.parse_args()

    if not args.data_path.is_file():
        raise FileNotFoundError(args.data_path)
    args.staging_root.mkdir(parents=True, exist_ok=True)
    require_disk_space(args.staging_root, 5.0)
    config = resolved_run_config(args.dataset, args.variant, args.seed, args.epochs)
    config.update(
        {
            "execution_mode": "stress",
            "stress_version": STRESS_VERSION,
            "stress_bad_edge_ratio": float(args.contamination),
            "neighbor_estimator": args.estimator,
        }
    )
    data_digest = sha256_file(args.data_path)
    sources = stress_code_digests()
    runtime = runtime_identity()
    freeze = json.loads((PAPER_ROOT / "experiments/protocol_v1/source_freeze.json").read_text(encoding="utf-8"))
    stress_freeze = json.loads((PAPER_ROOT / "experiments/stress_v1/source_freeze.json").read_text(encoding="utf-8"))
    if stress_freeze["parent_primary_freeze_hash"] != freeze["freeze_hash"]:
        raise RuntimeError("Stress freeze does not reference the active primary freeze")
    if stress_freeze["code_sha256"] != sources or stress_freeze["runtime"] != runtime:
        raise RuntimeError("Stress code/runtime changed after stress_v1 freeze")
    if freeze["dataset_sha256"][args.dataset] != data_digest:
        raise RuntimeError("Stress task data does not match the frozen primary dataset")
    payload = {
        "config": config, "data_sha256": data_digest, "code_sha256": sources,
        "stress_freeze_hash": stress_freeze["freeze_hash"], "runtime": runtime,
    }
    config_hash = sha256_payload(payload)
    run_id = (
        f"{args.dataset}--{args.variant}--contam{ratio_slug(args.contamination)}--"
        f"{args.estimator}--seed{args.seed}--{config_hash}"
    )
    run_dir = args.staging_root / run_id
    raw_dir = run_dir / "raw"
    remote_dir = (
        f"{REMOTE_RESULT_ROOT}/runs/{STRESS_VERSION}/{args.dataset}/{args.variant}/"
        f"contamination_{ratio_slug(args.contamination)}/{args.estimator}/seed_{args.seed}/{config_hash}"
    )
    store = RemoteStore() if args.upload else None
    if store and store.exists(f"{remote_dir}/COMPLETED"):
        result = store.run(f"cd {shlex.quote(remote_dir)} && sha256sum -c SHA256SUMS", check=False)
        if result.returncode != 0:
            raise RuntimeError(f"Remote stress checksum verification failed: {remote_dir}")
        if args.cleanup_after_upload and (run_dir / "COMPLETED").is_file():
            verify_checksum_manifest(run_dir)
            shutil.rmtree(run_dir)
        print(f"SKIP verified remote complete: {remote_dir}")
        return 0
    if (run_dir / "COMPLETED").is_file():
        verify_checksum_manifest(run_dir)
        if store:
            store.upload_directory_atomic(run_dir, remote_dir)
            if args.cleanup_after_upload:
                shutil.rmtree(run_dir)
        print(remote_dir if store else run_dir)
        return 0
    run_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    runner_config = {key: value for key, value in config.items() if key != "stress_version"}
    command = cli_args(runner_config, args.data_path.resolve(), raw_dir, args.gpu)
    metadata = {
        "run_id": run_id,
        "dataset": args.dataset,
        "variant": args.variant,
        "seed": args.seed,
        "gpu": args.gpu,
        "config_hash": config_hash,
        "protocol": config,
        "source_freeze_hash": freeze["freeze_hash"],
        "stress_freeze_hash": stress_freeze["freeze_hash"],
        "data": {"sha256": data_digest, "size_bytes": args.data_path.stat().st_size},
        "code_sha256": sources,
        "git": git_revision(PROJECT_ROOT),
        "runtime_identity": runtime,
        "remote_dir": remote_dir,
        "label_usage": "cross-label edge injection only when contamination > 0",
        "started_utc": utc_now(),
        "status": "running",
    }
    write_json(run_dir / "run_metadata.json", metadata)
    write_json(run_dir / "resolved_config.json", payload)
    env = os.environ.copy()
    env.update(
        {
            "CUDA_VISIBLE_DEVICES": str(args.gpu),
            "MPLCONFIGDIR": "/tmp/scvicar-matplotlib",
            "NUMBA_CACHE_DIR": "/tmp/scvicar-numba",
            "OMP_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
            "PYTORCH_CUDA_ALLOC_CONF": "expandable_segments:True",
            "PYTHONHASHSEED": str(args.seed), "CUBLAS_WORKSPACE_CONFIG": ":4096:8",
        }
    )
    metadata["environment"] = {key: env[key] for key in (
        "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
        "OPENBLAS_NUM_THREADS", "PYTHONHASHSEED", "CUBLAS_WORKSPACE_CONFIG",
    )}
    write_json(run_dir / "run_metadata.json", metadata)
    start = time.monotonic()
    local_complete = False
    try:
        returncode, resources = run_monitored(command, run_dir / "run.log", env, args.gpu)
        if returncode:
            raise RuntimeError(f"Stress model runner exited with code {returncode}")
        validation = validate_raw_outputs(raw_dir)
        canonicalize(raw_dir, run_dir)
        save_topology_bundle(raw_dir, run_dir)
        shutil.rmtree(raw_dir)
        metadata.update(status="complete", completed_utc=utc_now(), runtime_seconds=time.monotonic() - start, validation=validation, resources=resources)
        write_json(run_dir / "run_metadata.json", metadata)
        (run_dir / "COMPLETED").write_text(metadata["completed_utc"] + "\n", encoding="utf-8")
        write_checksum_manifest(run_dir)
        local_complete = True
        if store:
            store.upload_directory_atomic(run_dir, remote_dir)
            if args.cleanup_after_upload:
                shutil.rmtree(run_dir)
        print(remote_dir if store else run_dir)
        return 0
    except Exception as exc:
        if local_complete:
            write_json(
                run_dir.with_name(run_dir.name + "--UPLOAD_FAILED.json"),
                {"run_id": run_id, "failed_utc": utc_now(), "error": str(exc),
                 "traceback": traceback.format_exc(), "local_result_preserved_and_checksummed": True},
            )
            raise
        metadata.update(status="failed", failed_utc=utc_now(), runtime_seconds=time.monotonic() - start, error=str(exc), traceback=traceback.format_exc())
        write_json(run_dir / "run_metadata.json", metadata)
        (run_dir / "FAILED").write_text(str(exc) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
