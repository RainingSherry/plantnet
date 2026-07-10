from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import time
import traceback
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/scvicar-matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/scvicar-numba")

import numpy as np
import anndata as ad
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.preprocessing import LabelEncoder

from methods.DeepLearning import scMAE_family as family

from .config import DATASETS, PAPER_ROOT, PROJECT_ROOT, REMOTE_RESULT_ROOT, SEEDS, sha256_payload
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


BASELINE_VERSION = "baselines_v1"


@dataclass(frozen=True)
class BaselineSpec:
    python: str
    script: str
    gpu: bool
    args: tuple[str, ...] = ()


BASELINES = {
    "pca_kmeans": BaselineSpec(
        "/data/luolie/conda/base/bin/python", "methods/Traditional/PCA_KMeans/run.py", False,
        ("--label_key", "resolved_label", "--n_top_genes", "1000", "--pca_dim", "128"),
    ),
    "scmae": BaselineSpec(
        "/data/luolie/conda/envs/scclubench-sccdcg-h100/bin/python",
        "methods/DeepLearning/scMAE/run.py",
        True,
        (
            "--label_key", "resolved_label", "--input_mode", "auto",
            "--n_top_genes", "1000", "--hidden_size", "128",
            "--mask_prob", "0.4", "--epochs", "80", "--no_save_h5ad",
        ),
    ),
    "scvi": BaselineSpec(
        "/data/luolie/conda/base/bin/python", "methods/DeepLearning/scVI/run.py", True,
        ("--n_top_genes", "2000", "--epochs", "200"),
    ),
    "scdcc": BaselineSpec(
        "/data/luolie/conda/base/bin/python", "methods/DeepLearning/scDCC/run.py", True,
        ("--pretrain_epochs", "400", "--epochs", "200"),
    ),
    "scdeepcluster": BaselineSpec(
        "/data/luolie/conda/envs/scclubench-tf212-cpu/bin/python",
        "methods/DeepLearning/scDeepCluster/run.py", False,
        ("--pretrain_epochs", "400", "--maxiter", "20000", "--no_cuda"),
    ),
    "scrcl": BaselineSpec(
        "/data/luolie/conda/envs/plantnet-scrcl/bin/python", "methods/DeepLearning/scRCL/run.py", True,
        ("--label_key", "resolved_label", "--hvg", "2000", "--epochs", "200", "--loss_chunk_size", "1024"),
    ),
}


def baseline_code_digests(script: Path) -> dict[str, str]:
    """Hash the adapter, method package, and shared benchmark utilities."""
    paths = set(script.parent.rglob("*.py"))
    paths.update(
        path for path in (
            PROJECT_ROOT / "methods/utils.py",
            PROJECT_ROOT / "methods/evaluation.py",
            PROJECT_ROOT / "methods/preprocess.py",
            PROJECT_ROOT / "methods/shared_utils.py",
            PROJECT_ROOT / "methods/DeepLearning/scMAE_family.py",
            PROJECT_ROOT / "papers/scVICAR/code/config.py",
            PROJECT_ROOT / "papers/scVICAR/code/baseline_orchestrate.py",
            PROJECT_ROOT / "papers/scVICAR/code/io_utils.py",
            PROJECT_ROOT / "papers/scVICAR/code/remote_store.py",
            Path(__file__).resolve(),
        ) if path.is_file()
    )
    return {
        path.relative_to(PROJECT_ROOT).as_posix(): sha256_file(path)
        for path in sorted(paths)
        if "__pycache__" not in path.parts
    }


def baseline_identity(method: str) -> dict:
    spec = BASELINES[method]
    script = PROJECT_ROOT / spec.script
    python = Path(spec.python)
    lock = subprocess.run(
        [str(python), "-m", "pip", "freeze", "--all"], text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True,
    ).stdout.splitlines()
    distributions = sorted(line.strip() for line in lock if line.strip())
    return {
        "python": spec.python,
        "python_sha256": sha256_file(python.resolve()),
        "distribution_lock": distributions,
        "distribution_lock_hash": sha256_payload(distributions),
        "script": spec.script,
        "gpu": spec.gpu,
        "runner_args": list(spec.args),
        "code_sha256": baseline_code_digests(script),
    }


def canonicalize_baseline(
    raw: Path, output: Path, data_path: Path, dataset: str, method: str, seed: int
) -> dict:
    embedding_path = raw / "embedding_final.npy"
    labels_path = raw / "labels.npy"
    cell_ids_path = raw / "cell_ids.npy"
    if not embedding_path.is_file() or not labels_path.is_file() or not cell_ids_path.is_file():
        raise FileNotFoundError("Baseline must produce embedding_final.npy, labels.npy, and cell_ids.npy")
    embedding = np.load(embedding_path).astype(np.float32)
    upstream_labels = np.load(labels_path).astype(np.int64)
    upstream_cell_ids = np.load(cell_ids_path).astype(str)
    canonical = ad.read_h5ad(data_path, backed="r")
    expected_cells = int(canonical.n_obs)
    if "resolved_label" not in canonical.obs:
        raise KeyError("Canonical dataset is missing resolved_label")
    labels = LabelEncoder().fit_transform(canonical.obs["resolved_label"].astype(str)).astype(np.int64)
    canonical_cell_ids = np.asarray(canonical.obs_names, dtype=str)
    canonical.file.close()
    if embedding.ndim != 2 or upstream_labels.shape != (embedding.shape[0],) or not np.isfinite(embedding).all():
        raise ValueError("Invalid baseline embedding/labels")
    if embedding.shape[0] != expected_cells:
        raise ValueError(f"Baseline changed the frozen cell set: {embedding.shape[0]} != {expected_cells}")
    if upstream_cell_ids.shape != canonical_cell_ids.shape or not np.array_equal(upstream_cell_ids, canonical_cell_ids):
        raise ValueError("Baseline cell IDs do not exactly match the frozen canonical cell order")
    if adjusted_rand_score(labels, upstream_labels) != 1.0:
        raise ValueError("Baseline output labels do not match the frozen canonical cell order")
    k = DATASETS[dataset].expected_clusters
    pred = KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(embedding)
    metrics, mapped = family.compute_kmeans_metrics(labels, pred)
    np.savez_compressed(output / "embedding_float32.npz", embedding=embedding)
    np.savez_compressed(
        output / "clusters.npz", labels=labels, predicted=pred.astype(np.int32),
        mapped=mapped.astype(np.int32), cell_ids=canonical_cell_ids,
    )
    write_json(output / "metrics.json", {"kmeans_known_k": metrics})
    for name in ("args.json", "run_config.json", "config.json", "preprocess_config.json", "training_history.json"):
        if (raw / name).is_file():
            shutil.copy2(raw / name, output / f"upstream_{name}")
    return {
        "dataset": dataset, "method": method, "seed": seed, "n_cells": int(embedding.shape[0]),
        "expected_cells": int(expected_cells), "cell_set_complete": bool(embedding.shape[0] == expected_cells),
        "embedding_dim": int(embedding.shape[1]), "metrics": metrics,
        "cell_id_order_exact": True,
        "evaluation": "post-hoc known-K KMeans; per-cell labels are excluded from optimization",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Run one immutable scVICAR external baseline")
    parser.add_argument("--dataset", required=True, choices=sorted(DATASETS))
    parser.add_argument("--data-path", required=True, type=Path)
    parser.add_argument("--method", required=True, choices=sorted(BASELINES))
    parser.add_argument("--seed", required=True, type=int, choices=SEEDS)
    parser.add_argument("--gpu", type=int, default=2, choices=range(1, 7))
    parser.add_argument("--staging-root", type=Path, default=PAPER_ROOT / ".staging/baselines")
    parser.add_argument("--upload", action="store_true")
    parser.add_argument("--cleanup-after-upload", action="store_true")
    args = parser.parse_args()
    args.staging_root.mkdir(parents=True, exist_ok=True)
    require_disk_space(args.staging_root, 5.0)
    spec = BASELINES[args.method]
    python = Path(spec.python); script = PROJECT_ROOT / spec.script
    if not python.is_file() or not script.is_file() or not args.data_path.is_file():
        raise FileNotFoundError((python, script, args.data_path))
    freeze = json.loads((PAPER_ROOT / "experiments/protocol_v1/source_freeze.json").read_text())
    dataset_manifest = PAPER_ROOT / "manifests/dataset_upload/dataset_manifest.json"
    if sha256_file(dataset_manifest) != freeze["dataset_manifest_sha256"]:
        raise RuntimeError("Canonical dataset manifest changed after the primary protocol freeze")
    data_digest = sha256_file(args.data_path)
    if freeze["dataset_sha256"][args.dataset] != data_digest:
        raise RuntimeError("Baseline data does not match frozen canonical dataset")
    identity = baseline_identity(args.method)
    frozen_dataset = freeze["protocol"]["datasets"][args.dataset]
    known_k = int(frozen_dataset["expected_clusters"])
    if known_k != DATASETS[args.dataset].expected_clusters:
        raise RuntimeError("Live dataset K disagrees with the primary protocol freeze")
    if args.seed not in freeze["protocol"]["seeds"]:
        raise RuntimeError("Seed is absent from the primary protocol freeze")
    known_k_training = args.method in {"scdcc", "scdeepcluster", "scrcl"}
    protocol = {
        "baseline_version": BASELINE_VERSION, "dataset": args.dataset, "method": args.method,
        "seed": args.seed, "known_k": known_k,
        "known_k_training": known_k_training,
        "per_cell_labels_in_optimization": False,
        **identity,
    }
    baseline_freeze_path = PAPER_ROOT / "experiments/baselines_v1/source_freeze.json"
    if not baseline_freeze_path.is_file():
        raise FileNotFoundError(f"Freeze external baselines before execution: {baseline_freeze_path}")
    baseline_freeze = json.loads(baseline_freeze_path.read_text(encoding="utf-8"))
    if baseline_freeze["primary_freeze_hash"] != freeze["freeze_hash"]:
        raise RuntimeError("Baseline freeze does not reference the active primary freeze")
    if baseline_freeze["methods"][args.method] != identity:
        raise RuntimeError(f"Baseline source changed after freeze: {args.method}")
    config_hash = sha256_payload({
        "protocol": protocol, "data_sha256": data_digest,
        "baseline_freeze_hash": baseline_freeze["freeze_hash"],
    })
    run_id = f"{args.dataset}--{args.method}--seed{args.seed}--{config_hash}"
    run_dir = args.staging_root / run_id; raw = run_dir / "raw"
    remote = f"{REMOTE_RESULT_ROOT}/runs/{BASELINE_VERSION}/{args.dataset}/{args.method}/seed_{args.seed}/{config_hash}"
    store = RemoteStore() if args.upload else None
    if store and store.exists(f"{remote}/COMPLETED"):
        result = store.run(
            f"cd {shlex.quote(remote)} && sha256sum -c SHA256SUMS", check=False
        )
        if result.returncode != 0:
            raise RuntimeError(f"Remote baseline checksum verification failed: {remote}")
        if args.cleanup_after_upload and (run_dir / "COMPLETED").is_file():
            verify_checksum_manifest(run_dir)
            shutil.rmtree(run_dir)
        print(f"SKIP verified remote complete: {remote}"); return 0
    if (run_dir / "COMPLETED").is_file():
        verify_checksum_manifest(run_dir)
        if store:
            store.upload_directory_atomic(run_dir, remote)
            if args.cleanup_after_upload:
                shutil.rmtree(run_dir)
        print(remote if store else run_dir)
        return 0
    if run_dir.exists() and any(run_dir.iterdir()):
        suffix = utc_now().replace(":", "").replace("+", "_")
        run_dir = args.staging_root / f"{run_id}--retry-{suffix}"
        raw = run_dir / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    command = [str(python), str(script), "--data_path", str(args.data_path.resolve()), "--save_dir", str(raw),
               "--n_clusters", str(known_k), "--seed", str(args.seed), *spec.args]
    if spec.gpu:
        command.extend(["--gpu", "0"])
    env = os.environ.copy()
    if spec.gpu:
        env["CUDA_VISIBLE_DEVICES"] = str(args.gpu)
    for key, value in {
        "MPLCONFIGDIR": "/tmp/scvicar-matplotlib",
        "NUMBA_CACHE_DIR": "/tmp/scvicar-numba",
        "OMP_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "PYTHONHASHSEED": str(args.seed),
        "CUBLAS_WORKSPACE_CONFIG": ":4096:8",
        "TF_DETERMINISTIC_OPS": "1",
        "TF_ENABLE_ONEDNN_OPTS": "0",
    }.items():
        env[key] = value
    metadata = {
        "run_id": run_id, "dataset": args.dataset, "method": args.method, "seed": args.seed,
        "variant": args.method, "execution_mode": "formal",
        "config_hash": config_hash, "protocol": protocol, "source_freeze_hash": freeze["freeze_hash"],
        "baseline_freeze_hash": baseline_freeze["freeze_hash"],
        "code_sha256": identity["code_sha256"],
        "known_k_training": known_k_training,
        "per_cell_labels_in_optimization": False,
        "data": {"sha256": data_digest, "size_bytes": args.data_path.stat().st_size},
        "git": git_revision(PROJECT_ROOT),
        "python": {
            "launcher": str(python),
            "version": subprocess.run(
                [str(python), "--version"], text=True, stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT, check=True,
            ).stdout.strip(),
        },
        "environment": {
            key: env[key] for key in (
                "CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS", "PYTHONHASHSEED", "CUBLAS_WORKSPACE_CONFIG",
                "TF_DETERMINISTIC_OPS", "TF_ENABLE_ONEDNN_OPTS",
            ) if key in env
        },
        "remote_dir": remote, "started_utc": utc_now(), "status": "running",
    }
    write_json(run_dir / "run_metadata.json", metadata); write_json(run_dir / "resolved_config.json", protocol)
    start = time.monotonic()
    local_complete = False
    try:
        with (run_dir / "run.log").open("w", encoding="utf-8") as log:
            result = subprocess.run(command, cwd=PROJECT_ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError(f"Baseline exited with code {result.returncode}")
        validation = canonicalize_baseline(
            raw, run_dir, args.data_path, args.dataset, args.method, args.seed
        )
        shutil.rmtree(raw)
        metadata.update(status="complete", completed_utc=utc_now(), runtime_seconds=time.monotonic() - start, validation=validation)
        write_json(run_dir / "run_metadata.json", metadata)
        (run_dir / "COMPLETED").write_text(metadata["completed_utc"] + "\n", encoding="utf-8")
        write_checksum_manifest(run_dir)
        local_complete = True
        if store:
            store.upload_directory_atomic(run_dir, remote)
            if args.cleanup_after_upload: shutil.rmtree(run_dir)
        print(remote if store else run_dir); return 0
    except Exception as exc:
        if local_complete:
            write_json(
                run_dir.with_name(run_dir.name + "--UPLOAD_FAILED.json"),
                {
                    "run_id": run_id, "remote_dir": remote,
                    "failed_utc": utc_now(), "error": str(exc),
                    "traceback": traceback.format_exc(),
                    "local_result_preserved_and_checksummed": True,
                },
            )
            raise
        metadata.update(status="failed", failed_utc=utc_now(), error=str(exc), traceback=traceback.format_exc())
        write_json(run_dir / "run_metadata.json", metadata); (run_dir / "FAILED").write_text(str(exc) + "\n")
        raise


if __name__ == "__main__":
    raise SystemExit(main())
