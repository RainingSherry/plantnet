from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy as np

from .config import MODEL_RUNNER, PAPER_ROOT, PROJECT_ROOT, PYTHON_BIN, REMOTE_RESULT_ROOT, resolved_run_config, sha256_payload
from .io_utils import (
    git_revision,
    require_disk_space,
    sha256_file,
    utc_now,
    verify_checksum_manifest,
    write_checksum_manifest,
    write_json,
)
from .remote_store import RemoteStore


BOOL_ARGS = {"scale_input"}


def code_digests() -> dict[str, str]:
    paths = [
        Path(__file__).resolve(),
        PROJECT_ROOT / "papers/scVICAR/code/config.py",
        PROJECT_ROOT / "papers/scVICAR/code/orchestrate.py",
        PROJECT_ROOT / "papers/scVICAR/code/io_utils.py",
        PROJECT_ROOT / "papers/scVICAR/code/remote_store.py",
        MODEL_RUNNER,
        PROJECT_ROOT / "methods/DeepLearning/scMAE_family.py",
        PROJECT_ROOT / "methods/shared_utils.py",
        PROJECT_ROOT / "experimental_retired_models/NeighborMix_scMAE/model.py",
        PROJECT_ROOT / "experimental_retired_models/RG_NeighborMix_scMAE/model.py",
        PROJECT_ROOT / "experimental_retired_models/RG_NeighborMix_scMAE/mixing.py",
        PROJECT_ROOT / "experimental_retired_models/RG_NeighborMix_scMAE/neighbor_graph.py",
        PROJECT_ROOT / "experimental_retired_models/RG_NeighborMix_scMAE/diagnostics.py",
    ]
    return {str(path.relative_to(PROJECT_ROOT)): sha256_file(path) for path in paths}


def runtime_identity() -> dict:
    query = (
        "import importlib.metadata as m,json;"
        "names=['numpy','scipy','scikit-learn','torch','scanpy','anndata'];"
        "available={d.metadata['Name'].lower():d.version for d in m.distributions()};"
        "print(json.dumps({n:available.get(n.lower()) for n in names},sort_keys=True))"
    )
    packages = json.loads(subprocess.run(
        [str(PYTHON_BIN), "-c", query], text=True, stdout=subprocess.PIPE,
        stderr=subprocess.PIPE, check=True,
    ).stdout)
    return {
        "python": str(PYTHON_BIN),
        "python_sha256": sha256_file(PYTHON_BIN.resolve()),
        "packages": packages,
    }


def gpu_memory_mib(physical_gpu: int) -> float | None:
    result = subprocess.run(
        [
            "nvidia-smi", f"--id={physical_gpu}",
            "--query-gpu=memory.used", "--format=csv,noheader,nounits",
        ],
        text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, check=False,
    )
    if result.returncode != 0:
        return None
    try:
        return float(result.stdout.strip().splitlines()[0])
    except (IndexError, ValueError):
        return None


def run_monitored(command: list[str], log_path: Path, env: dict[str, str], gpu: int) -> tuple[int, dict]:
    baseline = gpu_memory_mib(gpu)
    peak = baseline
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command, cwd=PROJECT_ROOT, env=env, stdout=log, stderr=subprocess.STDOUT
        )
        while process.poll() is None:
            observed = gpu_memory_mib(gpu)
            if observed is not None:
                peak = observed if peak is None else max(peak, observed)
            time.sleep(0.5)
        returncode = process.wait()
    return returncode, {
        "device_memory_start_mib": baseline,
        "device_memory_peak_mib": peak,
        "device_memory_peak_delta_mib": (
            peak - baseline if peak is not None and baseline is not None else None
        ),
        "measurement": "nvidia-smi device memory sampled every 0.5 s on the isolated physical GPU",
    }


def cli_args(config: dict, data_path: Path, save_dir: Path, gpu: int) -> list[str]:
    args = [
        str(PYTHON_BIN), str(MODEL_RUNNER),
        "--data_path", str(data_path),
        "--save_dir", str(save_dir),
        "--gpu", str(gpu),
        "--no_save_h5ad",
    ]
    excluded = {"protocol_version", "execution_mode"}
    for key, value in config.items():
        if key in excluded:
            continue
        args.extend([f"--{key}", str(value).lower() if key in BOOL_ARGS else str(value)])
    return args


def required_outputs(raw_dir: Path) -> list[Path]:
    names = [
        "args.json", "dataset_profile.json", "preprocess_config.json",
        "training_history.json", "metrics.json", "eval_fixed.csv",
        "eval_metrics.json", "summary.json", "embedding_final.npy",
        "labels.npy", "gene_names.npy", "eval_kmeans_known_k.npy",
    ]
    return [raw_dir / name for name in names]


def validate_raw_outputs(raw_dir: Path) -> dict:
    missing = [str(path) for path in required_outputs(raw_dir) if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing formal outputs: " + ", ".join(missing))
    embedding = np.load(raw_dir / "embedding_final.npy", mmap_mode="r")
    labels = np.load(raw_dir / "labels.npy", mmap_mode="r")
    pred = np.load(raw_dir / "eval_kmeans_known_k.npy", mmap_mode="r")
    if embedding.ndim != 2 or labels.shape != (embedding.shape[0],) or pred.shape != labels.shape:
        raise ValueError("Embedding/label/prediction shapes are inconsistent")
    if not np.isfinite(np.asarray(embedding)).all():
        raise ValueError("Embedding contains NaN or Inf")
    metrics = json.loads((raw_dir / "metrics.json").read_text(encoding="utf-8"))
    known = metrics.get("kmeans_known_k", {})
    if not known.get("uses_known_k", False):
        raise ValueError("Formal primary metric must use the fixed known-K protocol")
    return {
        "n_cells": int(embedding.shape[0]),
        "embedding_dim": int(embedding.shape[1]),
        "n_labels": int(np.unique(labels).size),
        "n_pred_clusters": int(np.unique(pred).size),
        "metrics": known,
    }


def canonicalize(raw_dir: Path, run_dir: Path) -> None:
    embedding = np.load(raw_dir / "embedding_final.npy").astype(np.float32, copy=False)
    labels = np.load(raw_dir / "labels.npy").astype(np.int32, copy=False)
    pred = np.load(raw_dir / "eval_kmeans_known_k.npy").astype(np.int32, copy=False)
    mapped = np.load(raw_dir / "eval_kmeans_known_k_mapped.npy").astype(np.int32, copy=False)
    np.savez_compressed(run_dir / "embedding_float32.npz", embedding=embedding)
    np.savez_compressed(run_dir / "clusters.npz", labels=labels, predicted=pred, mapped=mapped)
    keep_names = [
        "args.json", "dataset_profile.json", "preprocess_config.json", "training_history.json",
        "metrics.json", "eval_fixed.csv", "eval_metrics.json", "summary.json",
        "neighbor_graph_profile.json", "train_knn_diagnostics.json", "edge_weight_summary.json",
        "neighbor_reliability_summary.json", "gate_summary.json", "pseudo_perturbation_summary.json",
        "contrast_diagnostics.json", "embedding_geometry_summary.csv", "per_cell_type_metrics.csv",
        "rare_cell_effect_summary.csv", "confusion_matrix_raw.csv", "confusion_matrix_mapped.csv",
    ]
    for name in keep_names:
        source = raw_dir / name
        if source.is_file():
            shutil.copy2(source, run_dir / name)


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one preregistered matched-backbone scVICAR task")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--data-path", required=True, type=Path)
    parser.add_argument("--variant", required=True)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--gpu", required=True, type=int, choices=range(1, 7))
    parser.add_argument("--staging-root", type=Path, default=PAPER_ROOT / ".staging" / "runs")
    parser.add_argument("--epochs", type=int)
    parser.add_argument("--upload", action="store_true")
    parser.add_argument("--cleanup-after-upload", action="store_true")
    parser.add_argument("--force-local", action="store_true")
    args = parser.parse_args()

    if not PYTHON_BIN.is_file():
        raise FileNotFoundError(PYTHON_BIN)
    if not MODEL_RUNNER.is_file():
        raise FileNotFoundError(MODEL_RUNNER)
    if not args.data_path.is_file():
        raise FileNotFoundError(args.data_path)
    args.staging_root.mkdir(parents=True, exist_ok=True)
    require_disk_space(args.staging_root, 5.0)

    config = resolved_run_config(args.dataset, args.variant, args.seed, args.epochs)
    data_digest = sha256_file(args.data_path)
    source_digests = code_digests()
    runtime = runtime_identity()
    freeze_path = PAPER_ROOT / "experiments" / "protocol_v1" / "source_freeze.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8")) if freeze_path.is_file() else None
    if freeze is not None and config["execution_mode"] == "formal":
        if freeze.get("code_sha256") != source_digests:
            raise RuntimeError("Current source digests do not match the frozen protocol_v1 source identity")
        if freeze.get("runtime") != runtime:
            raise RuntimeError("Current Python runtime does not match the frozen protocol_v1 environment identity")
        expected_data = freeze.get("dataset_sha256", {}).get(args.dataset)
        if expected_data != data_digest:
            raise RuntimeError(f"Dataset digest does not match frozen protocol_v1 for {args.dataset}")
    hash_payload = {
        "config": config,
        "data_sha256": data_digest,
        "runner": str(MODEL_RUNNER.relative_to(PROJECT_ROOT)),
        "code_sha256": source_digests,
        "source_freeze_hash": freeze.get("freeze_hash") if freeze else None,
        "runtime": runtime,
    }
    config_hash = sha256_payload(hash_payload)
    run_dir = args.staging_root / args.dataset / args.variant / f"seed_{args.seed}" / config_hash
    raw_dir = run_dir / "raw"
    remote_dir = f"{REMOTE_RESULT_ROOT}/runs/{config['protocol_version']}/{args.dataset}/{args.variant}/seed_{args.seed}/{config_hash}"
    store = RemoteStore() if args.upload else None
    if store and store.exists(f"{remote_dir}/COMPLETED") and not args.force_local:
        result = store.run(
            f"cd {shlex.quote(remote_dir)} && sha256sum -c SHA256SUMS", check=False
        )
        if result.returncode != 0:
            raise RuntimeError(f"Remote checksum verification failed: {remote_dir}")
        print(f"SKIP verified remote complete: {remote_dir}")
        return 0
    if (run_dir / "COMPLETED").is_file() and not args.force_local:
        verify_checksum_manifest(run_dir)
        if store:
            store.upload_directory_atomic(run_dir, remote_dir)
            if args.cleanup_after_upload:
                shutil.rmtree(run_dir)
        print(remote_dir if store else run_dir)
        return 0

    run_dir.mkdir(parents=True, exist_ok=True)
    raw_dir.mkdir(parents=True, exist_ok=True)
    command = cli_args(config, args.data_path.resolve(), raw_dir, args.gpu)
    env = os.environ.copy()
    env.update(
        {
            "CUDA_VISIBLE_DEVICES": str(args.gpu),
            "MPLCONFIGDIR": "/tmp/scvicar-matplotlib",
            "NUMBA_CACHE_DIR": "/tmp/scvicar-numba",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "PYTHONHASHSEED": str(args.seed),
            "CUBLAS_WORKSPACE_CONFIG": ":4096:8",
        }
    )
    metadata = {
        "run_id": f"{args.dataset}--{args.variant}--seed{args.seed}--{config_hash}",
        "dataset": args.dataset,
        "variant": args.variant,
        "seed": args.seed,
        "gpu": args.gpu,
        "config_hash": config_hash,
        "protocol": config,
        "data": {"sha256": data_digest, "size_bytes": args.data_path.stat().st_size},
        "git": git_revision(PROJECT_ROOT),
        "code_sha256": source_digests,
        "runtime_identity": runtime,
        "python": {
            "launcher": str(PYTHON_BIN),
            "version": subprocess.run(
                [str(PYTHON_BIN), "--version"], text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=True,
            ).stdout.strip(),
        },
        "environment": {
            key: env[key] for key in (
                "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "PYTHONHASHSEED", "CUBLAS_WORKSPACE_CONFIG",
            )
        },
        "remote_dir": remote_dir,
        "started_utc": utc_now(),
        "status": "running",
    }
    write_json(run_dir / "run_metadata.json", metadata)
    write_json(run_dir / "resolved_config.json", hash_payload)
    start = time.monotonic()
    log_path = run_dir / "run.log"
    try:
        returncode, resources = run_monitored(command, log_path, env, args.gpu)
        if returncode != 0:
            raise RuntimeError(f"Model runner exited with code {returncode}")
        validation = validate_raw_outputs(raw_dir)
        canonicalize(raw_dir, run_dir)
        shutil.rmtree(raw_dir)
        metadata.update(
            {
                "status": "complete",
                "completed_utc": utc_now(),
                "runtime_seconds": time.monotonic() - start,
                "resources": resources,
                "validation": validation,
            }
        )
        write_json(run_dir / "run_metadata.json", metadata)
        (run_dir / "COMPLETED").write_text(metadata["completed_utc"] + "\n", encoding="utf-8")
        write_checksum_manifest(run_dir)
        if store:
            store.upload_directory_atomic(run_dir, remote_dir)
            if args.cleanup_after_upload:
                shutil.rmtree(run_dir)
        print(remote_dir if store else run_dir)
        return 0
    except Exception as exc:
        metadata.update(
            {
                "status": "failed",
                "failed_utc": utc_now(),
                "runtime_seconds": time.monotonic() - start,
                "error": str(exc),
                "traceback": traceback.format_exc(),
            }
        )
        write_json(run_dir / "run_metadata.json", metadata)
        (run_dir / "FAILED").write_text(str(exc) + "\n", encoding="utf-8")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
