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

from .config import PAPER_ROOT, PROJECT_ROOT, REMOTE_RESULT_ROOT, resolved_run_config, sha256_payload
from .io_utils import (
    git_revision, require_disk_space, sha256_file, utc_now, verify_checksum_manifest,
    write_checksum_manifest, write_json,
)
from .prepare_full_label_sensitivity import VERSION
from .remote_store import RemoteStore
from .run_scVICAR import (
    MODEL_RUNNER, PYTHON_BIN, canonicalize, cli_args, code_digests, run_monitored,
    runtime_identity, validate_raw_outputs,
)


SENSITIVITY_VARIANTS = ("nomix", "fixed", "topology_full")


def sensitivity_code_digests() -> dict[str, str]:
    result = code_digests()
    for path in (
        Path(__file__).resolve(),
        PROJECT_ROOT / "papers/scVICAR/code/full_label_sensitivity_orchestrate.py",
        PROJECT_ROOT / "papers/scVICAR/code/freeze_full_label_sensitivity.py",
        PROJECT_ROOT / "papers/scVICAR/code/prepare_full_label_sensitivity.py",
    ):
        result[path.relative_to(PROJECT_ROOT).as_posix()] = sha256_file(path)
    return dict(sorted(result.items()))


def manifest_lookup(path: Path) -> dict[str, dict]:
    rows = json.loads(path.read_text(encoding="utf-8"))
    lookup = {row["name"]: row for row in rows}
    if len(rows) != 6 or len(lookup) != 6:
        raise ValueError("Full-label sensitivity requires exactly six unique datasets")
    return lookup


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one frozen full-label sensitivity task")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--data-path", required=True, type=Path)
    parser.add_argument("--variant", required=True, choices=SENSITIVITY_VARIANTS)
    parser.add_argument("--seed", required=True, type=int, choices=(42, 2024, 3407))
    parser.add_argument("--gpu", required=True, type=int, choices=range(1, 7))
    parser.add_argument(
        "--manifest", type=Path,
        default=PAPER_ROOT / f"manifests/{VERSION}/dataset_manifest.json",
    )
    parser.add_argument(
        "--freeze", type=Path,
        default=PAPER_ROOT / f"experiments/{VERSION}/source_freeze.json",
    )
    parser.add_argument(
        "--staging-root", type=Path,
        default=PAPER_ROOT / f".staging/{VERSION}",
    )
    parser.add_argument("--upload", action="store_true")
    parser.add_argument("--cleanup-after-upload", action="store_true")
    args = parser.parse_args()
    if not args.data_path.is_file():
        raise FileNotFoundError(args.data_path)
    args.staging_root.mkdir(parents=True, exist_ok=True)
    require_disk_space(args.staging_root, 5.0)
    datasets = manifest_lookup(args.manifest)
    if args.dataset not in datasets:
        raise KeyError(args.dataset)
    dataset = datasets[args.dataset]
    data_digest = sha256_file(args.data_path)
    if data_digest != dataset["sha256"]:
        raise ValueError("Full-label dataset digest mismatch")
    config = resolved_run_config(args.dataset, args.variant, args.seed)
    config.update({
        "protocol_version": VERSION,
        "n_clusters": int(dataset["n_labels"]),
        "execution_mode": "formal_sensitivity",
    })
    sources = sensitivity_code_digests()
    runtime = runtime_identity()
    freeze = json.loads(args.freeze.read_text(encoding="utf-8"))
    if freeze["version"] != VERSION or freeze["code_sha256"] != sources or freeze["runtime"] != runtime:
        raise RuntimeError("Full-label sensitivity source/runtime changed after freeze")
    if freeze["dataset_sha256"].get(args.dataset) != data_digest:
        raise RuntimeError("Full-label sensitivity data changed after freeze")
    payload = {
        "config": config, "data_sha256": data_digest, "code_sha256": sources,
        "sensitivity_freeze_hash": freeze["freeze_hash"], "runtime": runtime,
    }
    config_hash = sha256_payload(payload)
    run_id = f"{args.dataset}--{args.variant}--seed{args.seed}--{config_hash}"
    run_dir = args.staging_root / run_id
    raw_dir = run_dir / "raw"
    remote_dir = (
        f"{REMOTE_RESULT_ROOT}/runs/{VERSION}/{args.dataset}/{args.variant}/"
        f"seed_{args.seed}/{config_hash}"
    )
    store = RemoteStore() if args.upload else None
    if store and store.exists(f"{remote_dir}/COMPLETED"):
        result = store.run(f"cd {shlex.quote(remote_dir)} && sha256sum -c SHA256SUMS", check=False)
        if result.returncode != 0:
            raise RuntimeError(f"Remote sensitivity checksum failed: {remote_dir}")
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
    if run_dir.exists() and any(run_dir.iterdir()):
        run_dir = args.staging_root / f"{run_id}--retry-{utc_now().replace(':', '')}"
        raw_dir = run_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)
    command_config = {key: value for key, value in config.items() if key != "protocol_version"}
    command = cli_args(command_config, args.data_path.resolve(), raw_dir, args.gpu)
    env = os.environ.copy()
    env.update({
        "CUDA_VISIBLE_DEVICES": str(args.gpu), "MPLCONFIGDIR": "/tmp/scvicar-matplotlib",
        "NUMBA_CACHE_DIR": "/tmp/scvicar-numba", "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1",
        "PYTHONHASHSEED": str(args.seed), "CUBLAS_WORKSPACE_CONFIG": ":4096:8",
    })
    metadata = {
        "run_id": run_id, "dataset": args.dataset, "variant": args.variant,
        "seed": args.seed, "gpu": args.gpu, "config_hash": config_hash,
        "protocol": config, "sensitivity_freeze_hash": freeze["freeze_hash"],
        "parent_primary_freeze_hash": freeze["parent_primary_freeze_hash"],
        "data": {"sha256": data_digest, "size_bytes": args.data_path.stat().st_size},
        "code_sha256": sources, "runtime_identity": runtime,
        "git": git_revision(PROJECT_ROOT), "remote_dir": remote_dir,
        "label_usage": "post-hoc known-K evaluation only; no clean graph/training/model-selection use",
        "environment": {key: env[key] for key in (
            "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "OPENBLAS_NUM_THREADS", "PYTHONHASHSEED", "CUBLAS_WORKSPACE_CONFIG",
        )},
        "started_utc": utc_now(), "status": "running",
    }
    write_json(run_dir / "run_metadata.json", metadata)
    write_json(run_dir / "resolved_config.json", payload)
    start = time.monotonic()
    local_complete = False
    try:
        returncode, resources = run_monitored(command, run_dir / "run.log", env, args.gpu)
        if returncode:
            raise RuntimeError(f"Sensitivity model runner exited with code {returncode}")
        validation = validate_raw_outputs(raw_dir)
        if validation["n_cells"] != int(dataset["n_cells"]):
            raise ValueError("Sensitivity runner changed the full cell set")
        if validation["n_labels"] != int(dataset["n_labels"]):
            raise ValueError("Sensitivity runner changed the full label set")
        canonicalize(raw_dir, run_dir)
        shutil.rmtree(raw_dir)
        metadata.update(
            status="complete", completed_utc=utc_now(), runtime_seconds=time.monotonic() - start,
            validation=validation, resources=resources,
        )
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
                {"run_id": run_id, "error": str(exc), "failed_utc": utc_now(),
                 "traceback": traceback.format_exc(), "local_complete_preserved": True},
            )
            raise
        metadata.update(status="failed", failed_utc=utc_now(), error=str(exc), traceback=traceback.format_exc())
        write_json(run_dir / "run_metadata.json", metadata)
        (run_dir / "FAILED").write_text(str(exc) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
