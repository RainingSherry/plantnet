#!/usr/bin/env python3
"""Benchmark-tuned Optuna search for NeighborMix, RG-NeighborMix, and RG+CL."""

from __future__ import annotations

import argparse
import csv
import fcntl
import hashlib
import json
import os
import queue
import shlex
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

from run_scmae_all_methods_suite import (
    DATA_ROOT_DEFAULT,
    DEFAULT_SEEDS,
    METRIC_NAMES,
    PROJECT_ROOT,
    TIMEOUT_SECONDS,
    DatasetSpec,
    discover_datasets,
    git_info,
    normalize_metrics,
    planned_prepared_dataset,
    prepare_dataset,
    read_json_or_none,
    validate_gpus,
    write_json,
)


PYTHON_BIN_DEFAULT = Path("/data/luolie/conda/envs/scclubench-sccdcg-h100/bin/python")
DEFAULT_DATASETS = [
    "Macosko",
    "Melanoma_5K",
    "Pollen",
    "Quake_Smart-seq2_Lung",
    "Tosches",
    "Wang",
    "worm_neuron_cell",
    "hrvatin_geo",
]
FULL16_DATASETS = [
    "Bach",
    "Baron",
    "Guo",
    "Limb_Muscle",
    "Macosko",
    "Melanoma_5K",
    "Pollen",
    "Quake_10x_Spleen",
    "Quake_Smart-seq2_Lung",
    "Shekhar",
    "Tosches",
    "Wang",
    "Young",
    "hrvatin",
    "hrvatin_geo",
    "worm_neuron_cell",
]
DEFAULT_GPUS = [1, 2, 3, 4, 5, 6]
MODEL_KEYS = ["neighbormix_scmae", "rg_neighbormix_scmae", "rg_neighbormix_scmae_contrast_safe"]
FINAL_ARTIFACT_STAGES = {"stage4", "stage7"}
FULL16_STAGES = {"stage5", "stage6", "stage7"}
FULL16_EXPECTED_ROWS = len(FULL16_DATASETS) * len(DEFAULT_SEEDS)
CONVERTED_CACHE_DIRS = [
    PROJECT_ROOT / "result" / "scmae_retired_methods_20260706_full" / "converted_data",
    PROJECT_ROOT / "result" / "scmae_all_methods_20260705_full" / "converted_data",
]


@dataclass(frozen=True)
class ModelFamily:
    key: str
    runner: Path
    method_name: str
    variant_name: str
    default_params: dict[str, Any]
    extra_args: list[str] = field(default_factory=lambda: ["--no_save_h5ad"])


@dataclass(frozen=True)
class TrialSpec:
    stage: str
    model: ModelFamily
    trial_id: str
    params: dict[str, Any]
    optuna_number: int | None = None


@dataclass(frozen=True)
class RunTask:
    trial: TrialSpec
    dataset: DatasetSpec
    seed: int


MODEL_FAMILIES = {
    "neighbormix_scmae": ModelFamily(
        key="neighbormix_scmae",
        runner=PROJECT_ROOT / "experimental_retired_models" / "NeighborMix_scMAE" / "run.py",
        method_name="NeighborMix_scMAE",
        variant_name="nm_scmae_mid",
        default_params={
            "n_top_genes": 1000,
            "hidden_size": 128,
            "dropout": 0.0,
            "masked_data_weight": 0.75,
            "mask_loss_weight": 0.7,
            "batch_size": 256,
            "lr": 1e-3,
            "mask_ratio": 0.4,
            "use_pseudo": True,
            "pseudo_weight": 0.3,
            "alpha": 0.9,
            "neighbor_k": 5,
            "mix_neighbors": 4,
            "tau": 0.2,
            "knn_pca_dim": 50,
        },
    ),
    "rg_neighbormix_scmae": ModelFamily(
        key="rg_neighbormix_scmae",
        runner=PROJECT_ROOT / "experimental_retired_models" / "RG_NeighborMix_scMAE" / "run.py",
        method_name="RG_NeighborMix_scMAE",
        variant_name="rg_nm_v1_reliability",
        default_params={
            "n_top_genes": 1000,
            "hidden_size": 128,
            "dropout": 0.0,
            "masked_data_weight": 0.75,
            "mask_loss_weight": 0.7,
            "batch_size": 256,
            "lr": 1e-3,
            "mask_ratio": 0.4,
            "mix_mode": "reliability",
            "gate_mode": "topology",
            "edge_reliability_mode": "sim_mutual_snn_distance",
            "neighbor_k": 10,
            "mix_neighbors": 4,
            "tau": 0.2,
            "knn_pca_dim": 50,
            "pseudo_weight": 0.3,
            "gate_max": 0.15,
            "gate_min": 0.0,
            "contrast_weight": 0.0,
        },
    ),
    "rg_neighbormix_scmae_contrast_safe": ModelFamily(
        key="rg_neighbormix_scmae_contrast_safe",
        runner=PROJECT_ROOT / "experimental_retired_models" / "RG_NeighborMix_scMAE" / "run.py",
        method_name="RG_NeighborMix_scMAE",
        variant_name="rg_neighbormix_scmae_contrast_safe",
        default_params={
            "n_top_genes": 1000,
            "hidden_size": 128,
            "dropout": 0.0,
            "masked_data_weight": 0.75,
            "mask_loss_weight": 0.7,
            "batch_size": 256,
            "lr": 1e-3,
            "mask_ratio": 0.4,
            "mix_mode": "reliability",
            "gate_mode": "topology",
            "edge_reliability_mode": "sim_mutual_snn_distance",
            "neighbor_k": 10,
            "mix_neighbors": 4,
            "tau": 0.2,
            "knn_pca_dim": 50,
            "pseudo_weight": 0.3,
            "gate_max": 0.15,
            "gate_min": 0.0,
            "contrast_weight": 0.01,
            "contrast_temperature": 0.5,
            "contrast_start_epoch": 20,
            "contrast_neighbor_positive_weight": 0.0,
            "contrast_hard_negative_weight": 0.0,
            "contrast_projection_dim": 64,
            "contrast_min_negatives": 16,
            "contrast_partition_mode": "kmeans",
        },
    ),
}


def make_out_root(out_root: Path | None) -> Path:
    if out_root is not None:
        return out_root.resolve()
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return PROJECT_ROOT / "result" / f"scmae_nm_rg_contrast_param_search_20260707_{stamp}"


def locked_prepare_dataset(dataset: DatasetSpec, converted_dir: Path) -> DatasetSpec:
    converted_dir.mkdir(parents=True, exist_ok=True)
    lock_path = converted_dir / f".{dataset.name}.prepare.lock"
    with lock_path.open("w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        try:
            cached = reuse_cached_prepared_dataset(dataset, converted_dir)
            if cached is not None:
                return cached
            return prepare_dataset(dataset, converted_dir)
        finally:
            fcntl.flock(lock_handle, fcntl.LOCK_UN)


def reuse_cached_prepared_dataset(dataset: DatasetSpec, converted_dir: Path) -> DatasetSpec | None:
    out_name = "hrvatin_geo.h5ad" if dataset.kind == "hrvatin_geo" else f"{dataset.name}.h5ad"
    meta_name = "hrvatin_geo.meta.json" if dataset.kind == "hrvatin_geo" else f"{dataset.name}.meta.json"
    out_path = converted_dir / out_name
    meta_path = converted_dir / meta_name
    if out_path.exists() and meta_path.exists():
        meta = read_json_or_none(meta_path) or {}
        return DatasetSpec(dataset.name, dataset.source_path, dataset.kind, out_path, int(meta["n_clusters"]))

    for cache_dir in CONVERTED_CACHE_DIRS:
        cache_h5ad = cache_dir / out_name
        cache_meta = cache_dir / meta_name
        if not (cache_h5ad.exists() and cache_meta.exists()):
            continue
        if not out_path.exists():
            out_path.symlink_to(cache_h5ad.resolve())
        if not meta_path.exists():
            meta_path.symlink_to(cache_meta.resolve())
        meta = read_json_or_none(meta_path) or {}
        if "n_clusters" not in meta:
            continue
        return DatasetSpec(dataset.name, dataset.source_path, dataset.kind, out_path, int(meta["n_clusters"]))
    return None


def parameter_signature(params: dict[str, Any]) -> str:
    payload = json.dumps(params, sort_keys=True, separators=(",", ":"))
    return hashlib.sha1(payload.encode("utf-8")).hexdigest()[:12]


def trial_out_dir(out_root: Path, trial: TrialSpec, dataset: DatasetSpec, seed: int) -> Path:
    return out_root / "runs" / trial.stage / trial.model.key / trial.trial_id / dataset.name / f"seed_{seed}"


def bool_string(value: bool) -> str:
    return "true" if bool(value) else "false"


def command_for(python_bin: Path, task: RunTask, out_dir: Path, gpu: int, epochs: int) -> list[str]:
    dataset = task.dataset
    if dataset.prepared_path is None or dataset.n_clusters is None:
        raise RuntimeError(f"Dataset was not prepared: {dataset.name}")
    model = task.trial.model
    cmd = [
        str(python_bin),
        str(model.runner),
        "--data_path",
        str(dataset.prepared_path),
        "--save_dir",
        str(out_dir),
        "--dataset_name",
        dataset.name,
        "--method_name",
        model.method_name,
        "--variant_name",
        model.variant_name,
        "--n_clusters",
        str(dataset.n_clusters),
        "--seed",
        str(task.seed),
        "--epochs",
        str(epochs),
        "--gpu",
        str(gpu),
    ]
    for key, value in sorted(task.trial.params.items()):
        if value is None:
            continue
        cmd.append(f"--{key}")
        cmd.append(bool_string(value) if isinstance(value, bool) else str(value))
    if task.trial.stage not in FINAL_ARTIFACT_STAGES:
        cmd.append("--lightweight_outputs")
    cmd.extend(model.extra_args)
    return cmd


def run_environment(python_bin: Path, gpu: int) -> dict[str, str]:
    env = os.environ.copy()
    env["PYTHONPATH"] = str(PROJECT_ROOT) + ":" + str(PROJECT_ROOT / "methods") + ":" + env.get("PYTHONPATH", "")
    env["PYTHONNOUSERSITE"] = "1"
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    env.setdefault("OMP_NUM_THREADS", "2")
    env.setdefault("OPENBLAS_NUM_THREADS", "2")
    env.setdefault("MKL_NUM_THREADS", "2")
    env.setdefault("NUMEXPR_NUM_THREADS", "2")
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    return env


def write_run_environment(out_dir: Path, python_bin: Path, gpu: int, trial: TrialSpec) -> None:
    payload: dict[str, Any] = {
        "python_executable": str(python_bin),
        "runtime_env": "scclubench-sccdcg-h100",
        "gpu": gpu,
        "model": trial.model.key,
        "stage": trial.stage,
        "trial_id": trial.trial_id,
        "parameter_signature": parameter_signature(trial.params),
    }
    try:
        probe = subprocess.run(
            [str(python_bin), "-c", "import sys, torch; print(sys.version.split()[0]); print(torch.__version__)"],
            text=True,
            capture_output=True,
            timeout=20,
        )
        if probe.returncode == 0:
            lines = [line.strip() for line in probe.stdout.splitlines() if line.strip()]
            if lines:
                payload["python_version"] = lines[0]
            if len(lines) > 1:
                payload["torch_version"] = lines[1]
    except Exception as exc:
        payload["probe_error"] = repr(exc)
    write_json(out_dir / "environment.json", payload)


def required_artifacts_missing(out_dir: Path, lightweight: bool) -> list[str]:
    required = ["metrics.json", "args.json"]
    if not lightweight:
        required.extend(["embedding_final.npy", "labels.npy"])
    return [name for name in required if not (out_dir / name).exists()]


def load_metrics(out_dir: Path) -> dict[str, Any]:
    return normalize_metrics(read_json_or_none(out_dir / "metrics.json") or {})


def load_contrast_valid_fraction(out_dir: Path) -> float | None:
    diag = read_json_or_none(out_dir / "contrast_diagnostics.json") or {}
    if not diag.get("contrast_enabled", False):
        return None
    value = diag.get("mean_valid_batch_fraction")
    try:
        return float(value)
    except Exception:
        return None


def status_row(task: RunTask, out_dir: Path) -> dict[str, Any]:
    status = read_json_or_none(out_dir / "status.json") or {}
    metrics = load_metrics(out_dir)
    row = {
        "stage": task.trial.stage,
        "trial_id": task.trial.trial_id,
        "optuna_number": task.trial.optuna_number if task.trial.optuna_number is not None else "",
        "dataset": task.dataset.name,
        "model": task.trial.model.key,
        "seed": task.seed,
        "status": status.get("status", "pending"),
        **{metric: metrics.get(metric, "") for metric in METRIC_NAMES},
        "runtime_seconds": status.get("runtime_seconds", ""),
        "gpu": status.get("gpu", ""),
        "command": status.get("command", ""),
        "save_dir": str(out_dir),
        "error": status.get("error", ""),
        "parameter_signature": parameter_signature(task.trial.params),
        "optuna_value": status.get("optuna_value", ""),
        "contrast_valid_batch_fraction": status.get("contrast_valid_batch_fraction", ""),
        "commit_sha": status.get("commit_sha", ""),
        "branch": status.get("branch", ""),
    }
    return row


def run_task(
    task: RunTask,
    out_root: Path,
    python_bin: Path,
    gpu_queue: queue.Queue[int],
    epochs: int,
    commit_sha: str,
    branch: str,
    resume: bool,
    rerun_failed: bool,
    dry_run: bool,
) -> dict[str, Any]:
    out_dir = trial_out_dir(out_root, task.trial, task.dataset, task.seed)
    status_path = out_dir / "status.json"
    if resume and status_path.exists():
        previous = read_json_or_none(status_path) or {}
        previous_status = str(previous.get("status", "unknown"))
        if previous_status == "success" or (previous_status in {"failed", "timeout", "incomplete", "contrast_invalid"} and not rerun_failed):
            return status_row(task, out_dir)

    out_dir.mkdir(parents=True, exist_ok=True)
    gpu = gpu_queue.get()
    try:
        cmd = command_for(python_bin, task, out_dir, gpu, epochs)
        write_json(
            out_dir / "trial_params.json",
            {
                "stage": task.trial.stage,
                "trial_id": task.trial.trial_id,
                "model": task.trial.model.key,
                "variant_name": task.trial.model.variant_name,
                "params": task.trial.params,
                "parameter_signature": parameter_signature(task.trial.params),
                "command": cmd,
            },
        )
        write_run_environment(out_dir, python_bin, gpu, task.trial)
        (out_dir / "command.txt").write_text(shlex.join(cmd) + "\n", encoding="utf-8")
        start = datetime.now().isoformat()
        if dry_run:
            status = {
                "status": "dry_run",
                "runtime_seconds": 0.0,
                "gpu": gpu,
                "command": shlex.join(cmd),
                "error": "",
                "python_executable": str(python_bin),
                "commit_sha": commit_sha,
                "branch": branch,
                "start_time": start,
                "end_time": datetime.now().isoformat(),
            }
            write_json(status_path, status)
            return status_row(task, out_dir)

        t0 = time.time()
        try:
            proc = subprocess.run(
                cmd,
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                timeout=TIMEOUT_SECONDS,
                env=run_environment(python_bin, gpu),
            )
            elapsed = round(time.time() - t0, 1)
            (out_dir / "run.log").write_text(
                f"Command: {shlex.join(cmd)}\nExit code: {proc.returncode}\n\n"
                f"=== STDOUT ===\n{proc.stdout}\n\n=== STDERR ===\n{proc.stderr}\n",
                encoding="utf-8",
            )
            missing = required_artifacts_missing(out_dir, lightweight=task.trial.stage not in FINAL_ARTIFACT_STAGES)
            if proc.returncode != 0:
                state = "failed"
                error = f"Exit code {proc.returncode}"
            elif missing:
                state = "incomplete"
                error = "Missing required artifacts: " + ", ".join(missing)
            else:
                state = "success"
                error = ""
        except subprocess.TimeoutExpired as exc:
            elapsed = round(time.time() - t0, 1)
            stdout = exc.stdout.decode() if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            stderr = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
            (out_dir / "run.log").write_text(
                f"Command timed out after {TIMEOUT_SECONDS}s\nCommand: {shlex.join(cmd)}\n"
                f"\n=== STDOUT ===\n{stdout}\n\n=== STDERR ===\n{stderr}\n",
                encoding="utf-8",
            )
            state = "timeout"
            error = f"Execution exceeded {TIMEOUT_SECONDS}s"

        contrast_valid = load_contrast_valid_fraction(out_dir)
        if task.trial.model.key == "rg_neighbormix_scmae_contrast_safe" and state == "success" and contrast_valid is not None and contrast_valid < 0.3:
            state = "contrast_invalid"
            error = f"Contrast valid row fraction below hard threshold: {contrast_valid:.4f}"

        status = {
            "status": state,
            "runtime_seconds": elapsed,
            "gpu": gpu,
            "command": shlex.join(cmd),
            "error": error,
            "python_executable": str(python_bin),
            "commit_sha": commit_sha,
            "branch": branch,
            "start_time": start,
            "end_time": datetime.now().isoformat(),
            "contrast_valid_batch_fraction": "" if contrast_valid is None else contrast_valid,
        }
        write_json(status_path, status)
        return status_row(task, out_dir)
    finally:
        gpu_queue.put(gpu)


def run_trial(
    out_root: Path,
    trial: TrialSpec,
    datasets: list[DatasetSpec],
    seeds: list[int],
    python_bin: Path,
    gpus: list[int],
    max_workers: int,
    epochs: int,
    commit_sha: str,
    branch: str,
    resume: bool,
    rerun_failed: bool,
    dry_run: bool,
) -> list[dict[str, Any]]:
    tasks = [RunTask(trial, dataset, seed) for dataset in datasets for seed in seeds]
    gpu_queue: queue.Queue[int] = queue.Queue()
    for gpu in gpus:
        gpu_queue.put(gpu)
    rows: list[dict[str, Any]] = []
    worker_count = max(1, min(int(max_workers), len(gpus), len(tasks)))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = [
            executor.submit(
                run_task,
                task,
                out_root,
                python_bin,
                gpu_queue,
                epochs,
                commit_sha,
                branch,
                resume,
                rerun_failed,
                dry_run,
            )
            for task in tasks
        ]
        for future in as_completed(futures):
            rows.append(future.result())
    return sorted(rows, key=lambda row: (row["dataset"], int(row["seed"])))


def rows_from_status_files(out_root: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for status_path in sorted((out_root / "runs").glob("*/*/trial_*/*/seed_*/status.json")):
        parts = status_path.relative_to(out_root / "runs").parts
        if len(parts) < 6:
            continue
        stage, model_key, trial_id, dataset_name, seed_dir, _ = parts
        seed = int(seed_dir.replace("seed_", ""))
        model = MODEL_FAMILIES.get(model_key)
        if model is None:
            continue
        params = read_json_or_none(status_path.parent / "trial_params.json") or {}
        trial = TrialSpec(stage=stage, model=model, trial_id=trial_id, params=params.get("params", {}))
        rows.append(status_row(RunTask(trial, DatasetSpec(dataset_name, Path(""), "h5ad"), seed), status_path.parent))
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def write_master_and_summaries(out_root: Path) -> list[dict[str, Any]]:
    import statistics
    from collections import defaultdict

    rows = rows_from_status_files(out_root)
    write_csv(out_root / "run_master.csv", rows)
    write_json(out_root / "run_master.json", rows)

    grouped_dataset: dict[tuple[str, str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    grouped_config: dict[tuple[str, str, str, str], list[dict[str, Any]]] = defaultdict(list)
    grouped_model: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped_dataset[(row["stage"], row["model"], row["trial_id"], row["parameter_signature"], row["dataset"])].append(row)
        grouped_config[(row["stage"], row["model"], row["trial_id"], row["parameter_signature"])].append(row)
        if row.get("status") == "success":
            grouped_model[row["model"]].append(row)

    dataset_rows: list[dict[str, Any]] = []
    for key, group in sorted(grouped_dataset.items()):
        out = {
            "stage": key[0],
            "model": key[1],
            "trial_id": key[2],
            "parameter_signature": key[3],
            "dataset": key[4],
            "n_rows": len(group),
            "n_success": sum(1 for row in group if row.get("status") == "success"),
        }
        for metric in METRIC_NAMES:
            vals = [float(row[metric]) for row in group if row.get("status") == "success" and row.get(metric) not in {"", None}]
            if vals:
                out[f"{metric}_mean"] = statistics.mean(vals)
                out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        dataset_rows.append(out)
    write_csv(out_root / "config_dataset_summary.csv", dataset_rows)

    config_rows: list[dict[str, Any]] = []
    for key, group in sorted(grouped_config.items()):
        success = [row for row in group if row.get("status") == "success"]
        contrast_vals = [
            float(row["contrast_valid_batch_fraction"])
            for row in group
            if row.get("contrast_valid_batch_fraction") not in {"", None}
        ]
        out = {
            "stage": key[0],
            "model": key[1],
            "trial_id": key[2],
            "parameter_signature": key[3],
            "n_rows": len(group),
            "n_success": len(success),
            "fail_rate": 1.0 - len(success) / max(1, len(group)),
            "mean_runtime_seconds": statistics.mean(
                [float(row["runtime_seconds"]) for row in group if row.get("runtime_seconds") not in {"", None}]
            )
            if any(row.get("runtime_seconds") not in {"", None} for row in group)
            else "",
            "mean_contrast_valid_batch_fraction": statistics.mean(contrast_vals) if contrast_vals else "",
        }
        for metric in METRIC_NAMES:
            vals = [float(row[metric]) for row in success if row.get(metric) not in {"", None}]
            if vals:
                out[f"{metric}_mean"] = statistics.mean(vals)
                out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        config_rows.append(out)
    config_rows.sort(key=lambda row: (row["model"], row.get("ari_mean", -1) if row.get("ari_mean") != "" else -1), reverse=True)
    write_csv(out_root / "global_config_summary.csv", config_rows)
    write_csv(out_root / "candidate_funnel.csv", config_rows)
    write_csv(
        out_root / "runtime_summary.csv",
        [
            {
                "stage": row["stage"],
                "model": row["model"],
                "trial_id": row["trial_id"],
                "parameter_signature": row["parameter_signature"],
                "mean_runtime_seconds": row.get("mean_runtime_seconds", ""),
                "fail_rate": row.get("fail_rate", ""),
                "n_success": row.get("n_success", ""),
                "n_rows": row.get("n_rows", ""),
            }
            for row in config_rows
        ],
    )

    model_rows: list[dict[str, Any]] = []
    for model, group in sorted(grouped_model.items()):
        out = {"model": model, "n_success": len(group)}
        for metric in METRIC_NAMES:
            vals = [float(row[metric]) for row in group if row.get(metric) not in {"", None}]
            if vals:
                out[f"{metric}_mean"] = statistics.mean(vals)
                out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        model_rows.append(out)
    write_csv(out_root / "summary_by_model.csv", model_rows)

    grouped_model_dataset: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        if row.get("status") == "success":
            grouped_model_dataset[(row["model"], row["dataset"])].append(row)
    model_dataset_rows: list[dict[str, Any]] = []
    for (model, dataset), group in sorted(grouped_model_dataset.items()):
        out = {"model": model, "dataset": dataset, "n_success": len(group)}
        for metric in METRIC_NAMES:
            vals = [float(row[metric]) for row in group if row.get(metric) not in {"", None}]
            if vals:
                out[f"{metric}_mean"] = statistics.mean(vals)
                out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        model_dataset_rows.append(out)
    write_csv(out_root / "summary_by_model_dataset.csv", model_dataset_rows)

    active_param_rows: list[dict[str, Any]] = []
    seen_params: set[tuple[str, str, str, str]] = set()
    for params_path in sorted((out_root / "runs").glob("*/*/trial_*/*/seed_*/trial_params.json")):
        payload = read_json_or_none(params_path) or {}
        params = payload.get("params", {})
        if not isinstance(params, dict):
            continue
        rel = params_path.relative_to(out_root / "runs").parts
        if len(rel) < 6:
            continue
        stage, model, trial_id = rel[0], rel[1], rel[2]
        signature = str(payload.get("parameter_signature") or parameter_signature(params))
        key = (stage, model, trial_id, signature)
        if key in seen_params:
            continue
        seen_params.add(key)
        row = {"stage": stage, "model": model, "trial_id": trial_id, "parameter_signature": signature}
        row.update(params)
        active_param_rows.append(row)
    write_csv(out_root / "active_param_summary.csv", active_param_rows)
    return rows


def stage_priority(stage: str) -> int:
    return {
        "stage7": 70,
        "stage6": 60,
        "stage5": 50,
        "stage4": 40,
        "stage3": 30,
        "stage2": 20,
        "stage0": 10,
        "smoke": 0,
    }.get(stage, -1)


def dedupe_rows_for_full16(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    wanted = set(FULL16_DATASETS)
    best: dict[tuple[str, str, str, str], dict[str, Any]] = {}
    for row in rows:
        if row.get("model") not in MODEL_KEYS or row.get("dataset") not in wanted:
            continue
        if row.get("stage") in {"smoke", "stage0"}:
            continue
        signature = str(row.get("parameter_signature") or "")
        if not signature:
            continue
        key = (row["model"], signature, row["dataset"], str(row["seed"]))
        current = best.get(key)
        if current is None:
            best[key] = row
            continue
        row_rank = (row.get("status") == "success", stage_priority(str(row.get("stage", ""))))
        current_rank = (current.get("status") == "success", stage_priority(str(current.get("stage", ""))))
        if row_rank > current_rank:
            best[key] = row
    return sorted(best.values(), key=lambda row: (row["model"], row["parameter_signature"], row["dataset"], int(row["seed"])))


def params_for_signature(out_root: Path, model_key: str, signature: str) -> dict[str, Any]:
    for params_path in sorted((out_root / "runs").glob(f"*/*/{'trial_*'}/*/seed_*/trial_params.json")):
        payload = read_json_or_none(params_path) or {}
        if payload.get("model") != model_key:
            continue
        params = payload.get("params")
        if isinstance(params, dict) and str(payload.get("parameter_signature") or parameter_signature(params)) == signature:
            return params
    return {}


def write_full16_summaries(out_root: Path) -> list[dict[str, Any]]:
    import statistics
    from collections import defaultdict

    rows = dedupe_rows_for_full16(rows_from_status_files(out_root))
    grouped_dataset: dict[tuple[str, str, str], list[dict[str, Any]]] = defaultdict(list)
    grouped_config: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped_dataset[(row["model"], row["parameter_signature"], row["dataset"])].append(row)
        grouped_config[(row["model"], row["parameter_signature"])].append(row)

    dataset_rows: list[dict[str, Any]] = []
    for (model, signature, dataset), group in sorted(grouped_dataset.items()):
        stages = ",".join(sorted({str(row.get("stage", "")) for row in group}, key=stage_priority))
        trials = ",".join(sorted({str(row.get("trial_id", "")) for row in group}))
        out: dict[str, Any] = {
            "model": model,
            "parameter_signature": signature,
            "dataset": dataset,
            "source_stages": stages,
            "source_trials": trials,
            "n_rows": len(group),
            "n_success": sum(1 for row in group if row.get("status") == "success"),
        }
        for metric in METRIC_NAMES:
            vals = [float(row[metric]) for row in group if row.get("status") == "success" and row.get(metric) not in {"", None}]
            if vals:
                out[f"{metric}_mean"] = statistics.mean(vals)
                out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        dataset_rows.append(out)
    write_csv(out_root / "summary_by_model_dataset_full16.csv", dataset_rows)

    config_rows: list[dict[str, Any]] = []
    for (model, signature), group in sorted(grouped_config.items()):
        success = [row for row in group if row.get("status") == "success"]
        contrast_vals = [
            float(row["contrast_valid_batch_fraction"])
            for row in group
            if row.get("contrast_valid_batch_fraction") not in {"", None}
        ]
        out = {
            "model": model,
            "parameter_signature": signature,
            "source_stages": ",".join(sorted({str(row.get("stage", "")) for row in group}, key=stage_priority)),
            "source_trials": ",".join(sorted({str(row.get("trial_id", "")) for row in group})),
            "n_rows": len(group),
            "n_success": len(success),
            "fail_rate": 1.0 - len(success) / max(1, len(group)),
            "mean_runtime_seconds": statistics.mean(
                [float(row["runtime_seconds"]) for row in group if row.get("runtime_seconds") not in {"", None}]
            )
            if any(row.get("runtime_seconds") not in {"", None} for row in group)
            else "",
            "mean_contrast_valid_batch_fraction": statistics.mean(contrast_vals) if contrast_vals else "",
            "full16_complete": len(group) >= FULL16_EXPECTED_ROWS and len(success) >= FULL16_EXPECTED_ROWS,
        }
        for metric in METRIC_NAMES:
            vals = [float(row[metric]) for row in success if row.get(metric) not in {"", None}]
            if vals:
                out[f"{metric}_mean"] = statistics.mean(vals)
                out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
        config_rows.append(out)
    config_rows.sort(key=lambda row: (row["model"], row.get("full16_complete") is True, float(row.get("ari_mean") or -1)), reverse=True)
    write_csv(out_root / "global_config_summary_full16.csv", config_rows)

    model_rows: list[dict[str, Any]] = []
    for model in MODEL_KEYS:
        complete = [
            row
            for row in config_rows
            if row.get("model") == model and row.get("full16_complete") is True and row.get("ari_mean") not in {"", None}
        ]
        complete.sort(key=lambda row: float(row.get("ari_mean") or -1), reverse=True)
        if complete:
            model_rows.append(complete[0])
    write_csv(out_root / "summary_by_model_full16.csv", model_rows)
    write_selected_params_full16(out_root, config_rows)
    return rows


def read_top_full16_configs(out_root: Path, model_key: str, top_n: int, complete_only: bool = True) -> list[dict[str, Any]]:
    path = out_root / "global_config_summary_full16.csv"
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("model") != model_key or row.get("ari_mean") in {"", None}:
                continue
            if complete_only and str(row.get("full16_complete", "")).lower() != "true":
                continue
            rows.append(row)
    rows.sort(key=lambda row: float(row.get("ari_mean") or -1), reverse=True)
    configs: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        signature = row["parameter_signature"]
        if signature in seen:
            continue
        seen.add(signature)
        params = params_for_signature(out_root, model_key, signature)
        if params:
            configs.append(params)
        if len(configs) >= top_n:
            break
    return configs


def write_selected_params_full16(out_root: Path, config_rows: list[dict[str, Any]] | None = None) -> None:
    if config_rows is None:
        path = out_root / "global_config_summary_full16.csv"
        if not path.exists():
            return
        with path.open(newline="", encoding="utf-8") as handle:
            config_rows = list(csv.DictReader(handle))
    selected: dict[str, dict[str, Any]] = {}
    for model in MODEL_KEYS:
        candidates = [
            row
            for row in config_rows
            if row.get("model") == model
            and str(row.get("full16_complete", "")).lower() == "true"
            and row.get("ari_mean") not in {"", None}
        ]
        candidates.sort(key=lambda row: float(row.get("ari_mean") or -1), reverse=True)
        if not candidates:
            continue
        best = candidates[0]
        selected[model] = {
            "summary": best,
            "params": params_for_signature(out_root, model, str(best.get("parameter_signature") or "")),
        }
    write_json(out_root / "selected_params_full16.json", selected)
    lines = ["# Full-16 Benchmark-Tuned Selection Report", ""]
    for model, payload in selected.items():
        summary = payload["summary"]
        lines.append(
            f"- {model}: signature={summary.get('parameter_signature')} mean_ARI={summary.get('ari_mean')} "
            f"mean_ACC={summary.get('acc_mean')} mean_NMI={summary.get('nmi_mean')} "
            f"n_success={summary.get('n_success')}/{summary.get('n_rows')}"
        )
    (out_root / "selection_report_full16.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def canonical_lookup(rows: list[dict[str, Any]]) -> dict[tuple[str, str, str], dict[str, float]]:
    lookup: dict[tuple[str, str, str], dict[str, float]] = {}
    for row in rows:
        if row.get("stage") != "stage0" or row.get("status") != "success":
            continue
        metrics = {}
        for metric in METRIC_NAMES:
            if row.get(metric) not in {"", None}:
                metrics[metric] = float(row[metric])
        lookup[(row["model"], row["dataset"], str(row["seed"]))] = metrics
    return lookup


def objective_from_rows(model_key: str, rows: list[dict[str, Any]], canonical: dict[tuple[str, str, str], dict[str, float]]) -> tuple[float, dict[str, Any]]:
    import statistics

    expected = len(rows)
    success = [row for row in rows if row.get("status") == "success"]
    ari_vals = [float(row["ari"]) for row in success if row.get("ari") not in {"", None}]
    mean_ari = statistics.mean(ari_vals) if ari_vals else 0.0
    fail_rate = 1.0 - len(success) / max(1, expected)
    catastrophic = 0
    compared = 0
    for row in success:
        ref = canonical.get((model_key, row["dataset"], str(row["seed"])))
        if ref and "ari" in ref and row.get("ari") not in {"", None}:
            compared += 1
            catastrophic += int(float(row["ari"]) - float(ref["ari"]) < -0.10)
    catastrophic_rate = catastrophic / max(1, compared) if compared else 0.0
    runtimes = [float(row["runtime_seconds"]) for row in rows if row.get("runtime_seconds") not in {"", None}]
    mean_runtime = statistics.mean(runtimes) if runtimes else 0.0
    runtime_penalty = max(0.0, (mean_runtime / 3600.0) - 1.0)
    contrast_vals = [
        float(row["contrast_valid_batch_fraction"])
        for row in rows
        if row.get("contrast_valid_batch_fraction") not in {"", None}
    ]
    mean_contrast_valid = statistics.mean(contrast_vals) if contrast_vals else None
    contrast_penalty = 0.0
    hard_contrast_invalid = False
    if model_key == "rg_neighbormix_scmae_contrast_safe" and mean_contrast_valid is not None:
        hard_contrast_invalid = mean_contrast_valid < 0.3
        contrast_penalty = max(0.0, 0.7 - mean_contrast_valid) * 0.05
    score = mean_ari - 0.20 * fail_rate - 0.10 * catastrophic_rate - 0.02 * runtime_penalty - contrast_penalty
    if hard_contrast_invalid:
        score -= 0.50
    details = {
        "mean_ari": mean_ari,
        "success_count": len(success),
        "expected_count": expected,
        "fail_rate": fail_rate,
        "catastrophic_rate": catastrophic_rate,
        "mean_runtime_seconds": mean_runtime,
        "runtime_penalty": runtime_penalty,
        "mean_contrast_valid_batch_fraction": "" if mean_contrast_valid is None else mean_contrast_valid,
        "hard_contrast_invalid": hard_contrast_invalid,
        "objective": score,
    }
    return score, details


def common_train_params(trial, base: dict[str, Any] | None = None, local: bool = False) -> dict[str, Any]:
    if local and base:
        hidden_choices = sorted({64, 128, 256, 512, int(base.get("hidden_size", 128))})
        batch_choices = sorted({128, 256, 512, int(base.get("batch_size", 256))})
        n_top_choices = sorted({1000, 1500, 2000, 3000, int(base.get("n_top_genes", 1000))})
        lr_low = max(1e-4, float(base.get("lr", 1e-3)) / 2)
        lr_high = min(5e-3, float(base.get("lr", 1e-3)) * 2)
        return {
            "hidden_size": trial.suggest_categorical("hidden_size", hidden_choices),
            "dropout": trial.suggest_float("dropout", max(0.0, float(base.get("dropout", 0.0)) - 0.08), min(0.35, float(base.get("dropout", 0.0)) + 0.08)),
            "batch_size": trial.suggest_categorical("batch_size", batch_choices),
            "lr": trial.suggest_float("lr", lr_low, lr_high, log=True),
            "mask_ratio": trial.suggest_float("mask_ratio", max(0.15, float(base.get("mask_ratio", 0.4)) - 0.12), min(0.65, float(base.get("mask_ratio", 0.4)) + 0.12)),
            "masked_data_weight": trial.suggest_float("masked_data_weight", max(0.45, float(base.get("masked_data_weight", 0.75)) - 0.12), min(0.95, float(base.get("masked_data_weight", 0.75)) + 0.12)),
            "mask_loss_weight": trial.suggest_float("mask_loss_weight", max(0.35, float(base.get("mask_loss_weight", 0.7)) - 0.15), min(0.95, float(base.get("mask_loss_weight", 0.7)) + 0.15)),
            "n_top_genes": trial.suggest_categorical("n_top_genes", n_top_choices),
        }
    return {
        "hidden_size": trial.suggest_categorical("hidden_size", [64, 128, 256, 512]),
        "dropout": trial.suggest_float("dropout", 0.0, 0.3),
        "batch_size": trial.suggest_categorical("batch_size", [128, 256, 512]),
        "lr": trial.suggest_float("lr", 3e-4, 3e-3, log=True),
        "mask_ratio": trial.suggest_float("mask_ratio", 0.2, 0.6),
        "masked_data_weight": trial.suggest_float("masked_data_weight", 0.5, 0.9),
        "mask_loss_weight": trial.suggest_float("mask_loss_weight", 0.4, 0.9),
        "n_top_genes": trial.suggest_categorical("n_top_genes", [1000, 1500, 2000, 3000]),
    }


def local_int_choices(base_value: int, allowed: list[int], width: int) -> list[int]:
    return sorted({value for value in allowed if abs(int(value) - int(base_value)) <= width} | {int(base_value)})


def suggest_params(trial, model: ModelFamily, stage: str, top_configs: list[dict[str, Any]] | None = None) -> dict[str, Any]:
    local = stage in {"stage3", "stage6"} and bool(top_configs)
    base = dict(model.default_params)
    if local and top_configs:
        idx = trial.suggest_int("base_config_index", 0, len(top_configs) - 1)
        base.update(top_configs[idx])
    params = dict(model.default_params)
    if local:
        params.update(base)
    params.update(common_train_params(trial, base=base, local=local))
    if model.key == "neighbormix_scmae":
        if local:
            params.update(
                {
                    "neighbor_k": trial.suggest_categorical("local_neighbor_k_v2", [5, 10, 15, 20, 30]),
                    "mix_neighbors": trial.suggest_categorical("local_mix_neighbors_v2", [2, 4, 6, 8]),
                    "alpha": trial.suggest_float("alpha", max(0.65, float(base.get("alpha", 0.9)) - 0.08), min(0.98, float(base.get("alpha", 0.9)) + 0.08)),
                    "pseudo_weight": trial.suggest_float("pseudo_weight", max(0.03, float(base.get("pseudo_weight", 0.3)) / 2), min(0.8, float(base.get("pseudo_weight", 0.3)) * 2), log=True),
                    "tau": trial.suggest_float("tau", max(0.08, float(base.get("tau", 0.2)) - 0.12), min(0.8, float(base.get("tau", 0.2)) + 0.12)),
                }
            )
        else:
            params.update(
                {
                    "neighbor_k": trial.suggest_categorical("neighbor_k", [5, 10, 15, 20, 30]),
                    "mix_neighbors": trial.suggest_categorical("mix_neighbors", [2, 4, 6, 8]),
                    "alpha": trial.suggest_float("alpha", 0.7, 0.97),
                    "pseudo_weight": trial.suggest_float("pseudo_weight", 0.05, 0.6, log=True),
                    "tau": trial.suggest_float("tau", 0.1, 0.6),
                }
            )
    else:
        if local:
            params.update(
                {
                    "neighbor_k": trial.suggest_categorical("local_neighbor_k_v2", [5, 10, 15, 20, 30]),
                    "mix_neighbors": trial.suggest_categorical("local_mix_neighbors_v2", [2, 4, 6, 8]),
                    "pseudo_weight": trial.suggest_float("pseudo_weight", max(0.03, float(base.get("pseudo_weight", 0.3)) / 2), min(0.8, float(base.get("pseudo_weight", 0.3)) * 2), log=True),
                    "gate_max": trial.suggest_float("gate_max", max(0.03, float(base.get("gate_max", 0.15)) - 0.08), min(0.35, float(base.get("gate_max", 0.15)) + 0.08)),
                    "tau": trial.suggest_float("tau", max(0.08, float(base.get("tau", 0.2)) - 0.12), min(0.8, float(base.get("tau", 0.2)) + 0.12)),
                }
            )
        else:
            params.update(
                {
                    "neighbor_k": trial.suggest_categorical("neighbor_k", [5, 10, 15, 20, 30]),
                    "mix_neighbors": trial.suggest_categorical("mix_neighbors", [2, 4, 6, 8]),
                    "pseudo_weight": trial.suggest_float("pseudo_weight", 0.05, 0.6, log=True),
                    "gate_max": trial.suggest_float("gate_max", 0.05, 0.3),
                    "tau": trial.suggest_float("tau", 0.1, 0.6),
                }
            )
        params["gate_min"] = min(float(params.get("gate_min", 0.0)), float(params["gate_max"]))
        if model.key == "rg_neighbormix_scmae_contrast_safe":
            if local:
                params.update(
                    {
                        "contrast_weight": trial.suggest_float("contrast_weight", max(0.001, float(base.get("contrast_weight", 0.01)) / 2), min(0.1, float(base.get("contrast_weight", 0.01)) * 2), log=True),
                        "contrast_temperature": trial.suggest_float("contrast_temperature", max(0.1, float(base.get("contrast_temperature", 0.5)) - 0.2), min(1.0, float(base.get("contrast_temperature", 0.5)) + 0.2)),
                        "contrast_projection_dim": trial.suggest_categorical("contrast_projection_dim", [32, 64, 128]),
                        "contrast_min_negatives": trial.suggest_categorical("contrast_min_negatives", [8, 16, 32, 64]),
                        "contrast_neighbor_positive_weight": trial.suggest_categorical("contrast_neighbor_positive_weight", [0.0, 0.25, 0.5]),
                        "contrast_hard_negative_weight": trial.suggest_categorical("contrast_hard_negative_weight", [0.0, 0.25, 0.5, 1.0]),
                    }
                )
            else:
                params.update(
                    {
                        "contrast_weight": trial.suggest_float("contrast_weight", 0.003, 0.08, log=True),
                        "contrast_temperature": trial.suggest_float("contrast_temperature", 0.2, 0.8),
                        "contrast_start_epoch": trial.suggest_categorical("contrast_start_epoch", [10, 20, 30]),
                        "contrast_projection_dim": trial.suggest_categorical("contrast_projection_dim", [32, 64, 128]),
                        "contrast_min_negatives": trial.suggest_categorical("contrast_min_negatives", [8, 16, 32, 64]),
                        "contrast_neighbor_positive_weight": trial.suggest_categorical("contrast_neighbor_positive_weight", [0.0, 0.25, 0.5]),
                        "contrast_hard_negative_weight": trial.suggest_categorical("contrast_hard_negative_weight", [0.0, 0.25, 0.5, 1.0]),
                    }
                )
    return params


def optuna_storage(out_root: Path, study_name: str):
    import optuna
    from optuna.storages import JournalStorage
    from optuna.storages.journal import JournalFileBackend

    storage_path = out_root / "optuna" / f"{study_name}.log"
    storage_path.parent.mkdir(parents=True, exist_ok=True)
    return JournalStorage(JournalFileBackend(str(storage_path)))


def load_study(out_root: Path, study_name: str):
    import optuna

    sampler = optuna.samplers.TPESampler(
        seed=20260707,
        multivariate=True,
        group=True,
        constant_liar=True,
    )
    return optuna.create_study(
        study_name=study_name,
        direction="maximize",
        storage=optuna_storage(out_root, study_name),
        sampler=sampler,
        load_if_exists=True,
    )


def local_enqueue_params(model: ModelFamily, params: dict[str, Any], base_index: int) -> dict[str, Any]:
    enqueued = {
        "base_config_index": base_index,
        "hidden_size": int(params.get("hidden_size", model.default_params.get("hidden_size", 128))),
        "dropout": float(params.get("dropout", model.default_params.get("dropout", 0.0))),
        "batch_size": int(params.get("batch_size", model.default_params.get("batch_size", 256))),
        "lr": float(params.get("lr", model.default_params.get("lr", 1e-3))),
        "mask_ratio": float(params.get("mask_ratio", model.default_params.get("mask_ratio", 0.4))),
        "masked_data_weight": float(params.get("masked_data_weight", model.default_params.get("masked_data_weight", 0.75))),
        "mask_loss_weight": float(params.get("mask_loss_weight", model.default_params.get("mask_loss_weight", 0.7))),
        "n_top_genes": int(params.get("n_top_genes", model.default_params.get("n_top_genes", 1000))),
        "local_neighbor_k_v2": int(params.get("neighbor_k", model.default_params.get("neighbor_k", 10))),
        "local_mix_neighbors_v2": int(params.get("mix_neighbors", model.default_params.get("mix_neighbors", 4))),
        "pseudo_weight": float(params.get("pseudo_weight", model.default_params.get("pseudo_weight", 0.3))),
        "tau": float(params.get("tau", model.default_params.get("tau", 0.2))),
    }
    if model.key == "neighbormix_scmae":
        enqueued["alpha"] = float(params.get("alpha", model.default_params.get("alpha", 0.9)))
    else:
        enqueued["gate_max"] = float(params.get("gate_max", model.default_params.get("gate_max", 0.15)))
    if model.key == "rg_neighbormix_scmae_contrast_safe":
        enqueued.update(
            {
                "contrast_weight": float(params.get("contrast_weight", model.default_params.get("contrast_weight", 0.01))),
                "contrast_temperature": float(
                    params.get("contrast_temperature", model.default_params.get("contrast_temperature", 0.5))
                ),
                "contrast_projection_dim": int(
                    params.get("contrast_projection_dim", model.default_params.get("contrast_projection_dim", 64))
                ),
                "contrast_min_negatives": int(
                    params.get("contrast_min_negatives", model.default_params.get("contrast_min_negatives", 16))
                ),
                "contrast_neighbor_positive_weight": float(
                    params.get(
                        "contrast_neighbor_positive_weight",
                        model.default_params.get("contrast_neighbor_positive_weight", 0.0),
                    )
                ),
                "contrast_hard_negative_weight": float(
                    params.get(
                        "contrast_hard_negative_weight",
                        model.default_params.get("contrast_hard_negative_weight", 0.0),
                    )
                ),
            }
        )
    return enqueued


def enqueue_stage6_top_configs(study, model: ModelFamily, top_configs: list[dict[str, Any]]) -> None:
    if study.trials:
        return
    for idx, params in enumerate(top_configs):
        study.enqueue_trial(local_enqueue_params(model, params, idx))


def read_top_configs(out_root: Path, model_key: str, source_stage: str, top_n: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    path = out_root / "global_config_summary.csv"
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("model") == model_key and row.get("stage") == source_stage and row.get("ari_mean") not in {"", None}:
                rows.append(row)
    rows.sort(key=lambda row: float(row.get("ari_mean") or -1), reverse=True)
    configs: list[dict[str, Any]] = []
    for row in rows[:top_n]:
        candidates = sorted((out_root / "runs" / row["stage"] / model_key / row["trial_id"]).glob("*/seed_*/trial_params.json"))
        if candidates:
            payload = read_json_or_none(candidates[0]) or {}
            if isinstance(payload.get("params"), dict):
                configs.append(payload["params"])
    return configs


def read_shortlist_configs(out_root: Path, model_key: str, top_n: int) -> list[dict[str, Any]]:
    path = out_root / "global_config_summary.csv"
    if not path.exists():
        return []
    rows: list[dict[str, Any]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("model") != model_key or row.get("stage") not in {"stage2", "stage3", "stage4"}:
                continue
            if row.get("ari_mean") in {"", None}:
                continue
            try:
                if int(float(row.get("n_success") or 0)) < 24 or int(float(row.get("n_rows") or 0)) < 24:
                    continue
            except Exception:
                continue
            rows.append(row)
    rows.sort(key=lambda row: float(row.get("ari_mean") or -1), reverse=True)
    configs: list[dict[str, Any]] = []
    seen: set[str] = set()
    for row in rows:
        signature = row["parameter_signature"]
        if signature in seen:
            continue
        seen.add(signature)
        params = params_for_signature(out_root, model_key, signature)
        if params:
            configs.append(params)
        if len(configs) >= top_n:
            break
    return configs


def read_selected_configs(out_root: Path, model_key: str, top_n: int) -> list[dict[str, Any]]:
    selected_path = out_root / "selected_params.json"
    selected = read_json_or_none(selected_path)
    if not isinstance(selected, dict):
        return []
    payload = selected.get(model_key)
    if not isinstance(payload, dict):
        return []
    params = payload.get("params")
    if not isinstance(params, dict):
        return []
    return [params][:top_n]


def write_provenance(out_root: Path, args: argparse.Namespace, datasets: list[DatasetSpec], models: list[ModelFamily]) -> None:
    out_root.mkdir(parents=True, exist_ok=True)
    commit_sha, branch = git_info()
    write_json(
        out_root / "objective_definition.json",
        {
            "direction": "maximize",
            "primary_metric": "ari",
            "objective": "mean_ARI over selected datasets and seeds minus failure, catastrophic, runtime, and contrast-validity penalties",
            "seeds": args.seeds,
            "datasets": [dataset.name for dataset in datasets],
            "benchmark_tuned": True,
        },
    )
    write_json(
        out_root / "search_space.json",
        {
            "models": [model.key for model in models],
            "common": ["hidden_size", "dropout", "batch_size", "lr", "mask_ratio", "masked_data_weight", "mask_loss_weight", "n_top_genes"],
            "neighbormix": ["neighbor_k", "mix_neighbors", "alpha", "pseudo_weight", "tau"],
            "rg": ["neighbor_k", "mix_neighbors", "pseudo_weight", "gate_max", "tau"],
            "rg_contrast": ["contrast_weight", "contrast_temperature", "contrast_start_epoch", "contrast_projection_dim", "contrast_min_negatives", "contrast_neighbor_positive_weight", "contrast_hard_negative_weight"],
        },
    )
    budget_rows = []
    for stage, n_trials, epochs in [
        ("stage0", 1, args.epochs),
        ("smoke", 1, args.smoke_epochs),
        ("stage2", args.stage2_trials, args.epochs),
        ("stage3", args.stage3_trials, args.epochs),
        ("stage4", args.stage4_top_n, args.epochs),
        ("stage5", args.stage5_top_n, args.epochs),
        ("stage6", args.stage6_trials, args.epochs),
        ("stage7", args.stage7_top_n, args.epochs),
    ]:
        for model in models:
            n_datasets = 2 if stage == "smoke" else len(datasets)
            n_seeds = 1 if stage == "smoke" else len(args.seeds)
            budget_rows.append(
                {
                    "stage": stage,
                    "model": model.key,
                    "n_trials_or_configs": n_trials,
                    "n_datasets": n_datasets,
                    "n_seeds": n_seeds,
                    "epochs": epochs,
                    "expected_runs": n_trials * n_datasets * n_seeds,
                }
            )
    write_csv(out_root / "budget_estimate.csv", budget_rows)
    write_json(out_root / "dataset_manifest.json", [dataset.__dict__ for dataset in datasets])
    write_json(
        out_root / "model_manifest.json",
        [{"key": model.key, "runner": str(model.runner.relative_to(PROJECT_ROOT)), "variant_name": model.variant_name} for model in models],
    )
    write_json(
        out_root / "suite_config.json",
        {
            "out_root": str(out_root),
            "data_root": str(args.data_root),
            "seeds": args.seeds,
            "gpus": args.gpus,
            "models": [model.key for model in models],
            "commit_sha": commit_sha,
            "branch": branch,
            "python_bin": str(args.python_bin),
        },
    )
    try:
        diff = subprocess.run(["git", "diff"], cwd=PROJECT_ROOT, text=True, capture_output=True, timeout=60)
        (out_root / "git_diff.patch").write_text(diff.stdout, encoding="utf-8")
    except Exception as exc:
        (out_root / "git_diff.patch").write_text(f"git diff unavailable: {exc}\n", encoding="utf-8")
    try:
        freeze = subprocess.run([str(args.python_bin), "-m", "pip", "freeze"], text=True, capture_output=True, timeout=120)
        (out_root / "pip_freeze.txt").write_text(freeze.stdout, encoding="utf-8")
    except Exception as exc:
        (out_root / "pip_freeze.txt").write_text(f"pip freeze unavailable: {exc}\n", encoding="utf-8")
    env_payload: dict[str, Any] = {"python_executable": str(args.python_bin)}
    try:
        probe = subprocess.run(
            [
                str(args.python_bin),
                "-c",
                "import json, optuna, torch, sys; print(json.dumps({'python': sys.version.split()[0], 'torch': torch.__version__, 'optuna': optuna.__version__}))",
            ],
            text=True,
            capture_output=True,
            timeout=30,
        )
        if probe.returncode == 0 and probe.stdout.strip():
            env_payload.update(json.loads(probe.stdout.strip().splitlines()[-1]))
        else:
            env_payload["probe_error"] = probe.stderr[-2000:]
        (out_root / "optuna_version.txt").write_text(str(env_payload.get("optuna", "")) + "\n", encoding="utf-8")
    except Exception as exc:
        env_payload["probe_error"] = repr(exc)
    write_json(out_root / "environment.json", env_payload)
    try:
        conda = subprocess.run(["conda", "list", "--explicit", "-p", str(args.python_bin.parent.parent)], text=True, capture_output=True, timeout=120)
        (out_root / "conda_list_explicit.txt").write_text(conda.stdout + conda.stderr, encoding="utf-8")
    except Exception as exc:
        (out_root / "conda_list_explicit.txt").write_text(f"conda list unavailable: {exc}\n", encoding="utf-8")


def cleanup_smoke(out_root: Path) -> None:
    smoke_root = out_root / "runs" / "smoke"
    artifact_root = out_root / "smoke_artifacts"
    rows = rows_from_status_files(out_root)
    smoke_rows = [row for row in rows if row.get("stage") == "smoke"]
    existing_path = out_root / "smoke_summary.csv"
    if existing_path.exists():
        with existing_path.open(newline="", encoding="utf-8") as handle:
            existing_rows = list(csv.DictReader(handle))
    else:
        existing_rows = []
    merged_smoke: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in existing_rows + smoke_rows:
        key = (str(row.get("model", "")), str(row.get("dataset", "")), str(row.get("seed", "")))
        if all(key):
            merged_smoke[key] = row
    smoke_rows = sorted(merged_smoke.values(), key=lambda row: (row.get("model", ""), row.get("dataset", ""), int(row.get("seed", 0) or 0)))
    write_csv(out_root / "smoke_summary.csv", smoke_rows)
    write_json(out_root / "smoke_summary.json", smoke_rows)
    log_rows = []
    for row in smoke_rows:
        run_dir = Path(row["save_dir"])
        if not run_dir.exists():
            continue
        target = artifact_root / row["model"] / f"{row['dataset']}_seed_{row['seed']}"
        target.mkdir(parents=True, exist_ok=True)
        for name in ["args.json", "metrics.json", "status.json", "command.txt", "contrast_diagnostics.json"]:
            src = run_dir / name
            if src.exists():
                shutil.copy2(src, target / name)
        log_path = run_dir / "run.log"
        if log_path.exists():
            tail = "\n".join(log_path.read_text(encoding="utf-8", errors="ignore").splitlines()[-200:])
            (target / "run.log.tail").write_text(tail + "\n", encoding="utf-8")
        shutil.rmtree(run_dir)
        log_rows.append({"source": str(run_dir), "kept": str(target), "status": row.get("status", "")})
    if smoke_root.exists() and not any(smoke_root.rglob("status.json")):
        shutil.rmtree(smoke_root, ignore_errors=True)
    write_json(out_root / "smoke_cleanup_log.json", log_rows)


def write_selected_params(out_root: Path) -> None:
    path = out_root / "global_config_summary.csv"
    if not path.exists():
        return
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    selected: dict[str, dict[str, Any]] = {}
    for model in MODEL_KEYS:
        candidates = [
            row
            for row in rows
            if row.get("model") == model and row.get("stage") in {"stage4", "stage3", "stage2"} and row.get("ari_mean") not in {"", None}
        ]
        complete_candidates = []
        for row in candidates:
            try:
                n_success = int(float(row.get("n_success") or 0))
                n_rows = int(float(row.get("n_rows") or 0))
            except Exception:
                continue
            expected_rows = 24 if row.get("stage") in {"stage2", "stage3", "stage4"} else n_rows
            if n_rows >= expected_rows and n_success >= expected_rows:
                complete_candidates.append(row)
        candidates = complete_candidates
        candidates.sort(key=lambda row: float(row.get("ari_mean") or -1), reverse=True)
        if not candidates:
            continue
        best = candidates[0]
        param_files = sorted((out_root / "runs" / best["stage"] / model / best["trial_id"]).glob("*/seed_*/trial_params.json"))
        params = read_json_or_none(param_files[0]) if param_files else {}
        selected[model] = {"summary": best, "params": params.get("params", {}) if isinstance(params, dict) else {}}
    write_json(out_root / "selected_params.json", selected)
    lines = ["# Benchmark-Tuned Selection Report", ""]
    for model, payload in selected.items():
        summary = payload["summary"]
        lines.append(
            f"- {model}: trial={summary.get('trial_id')} mean_ARI={summary.get('ari_mean')} "
            f"mean_ACC={summary.get('acc_mean')} mean_NMI={summary.get('nmi_mean')}"
        )
    (out_root / "selection_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_existing_trial(out_root: Path, stage: str, model: ModelFamily, trial_id: str) -> TrialSpec:
    candidates = sorted((out_root / "runs" / stage / model.key / trial_id).glob("*/seed_*/trial_params.json"))
    if not candidates:
        raise SystemExit(f"No trial_params.json found for {stage}/{model.key}/{trial_id}")
    payload = read_json_or_none(candidates[0]) or {}
    params = payload.get("params")
    if not isinstance(params, dict):
        raise SystemExit(f"Invalid params payload in {candidates[0]}")
    optuna_number = None
    if trial_id.startswith("trial_"):
        try:
            optuna_number = int(trial_id.split("_", 1)[1])
        except ValueError:
            optuna_number = None
    return TrialSpec(stage=stage, model=model, trial_id=trial_id, params=params, optuna_number=optuna_number)


def run_optuna_stage(
    stage: str,
    out_root: Path,
    models: list[ModelFamily],
    datasets: list[DatasetSpec],
    seeds: list[int],
    args: argparse.Namespace,
    commit_sha: str,
    branch: str,
) -> None:
    import optuna

    rows = write_master_and_summaries(out_root)
    canonical = canonical_lookup(rows)
    trial_count = args.stage2_trials if stage == "stage2" else args.stage3_trials if stage == "stage3" else args.stage6_trials
    source_stage = "stage2" if stage == "stage3" else ""
    for model in models:
        if stage == "stage3":
            top_configs = read_top_configs(out_root, model.key, source_stage, args.stage3_top_n)
        elif stage == "stage6":
            write_full16_summaries(out_root)
            top_configs = read_top_full16_configs(out_root, model.key, args.stage6_top_n, complete_only=True)
            if not top_configs:
                top_configs = read_shortlist_configs(out_root, model.key, args.stage6_top_n)
        else:
            top_configs = []
        study = load_study(out_root, f"{stage}_{model.key}")
        if stage == "stage6" and top_configs:
            enqueue_stage6_top_configs(study, model, top_configs)
        completed = len([trial for trial in study.trials if trial.state == optuna.trial.TrialState.COMPLETE])
        while completed < trial_count:
            opt_trial = study.ask()
            params = suggest_params(opt_trial, model, stage=stage, top_configs=top_configs)
            trial = TrialSpec(stage=stage, model=model, trial_id=f"trial_{opt_trial.number:04d}", params=params, optuna_number=opt_trial.number)
            write_json(out_root / "trial_state_transitions" / f"{stage}_{model.key}_{opt_trial.number:04d}_asked.json", {"state": "asked", "params": params})
            rows = run_trial(
                out_root=out_root,
                trial=trial,
                datasets=datasets,
                seeds=seeds,
                python_bin=args.python_bin,
                gpus=args.gpus,
                max_workers=args.max_workers,
                epochs=args.epochs,
                commit_sha=commit_sha,
                branch=branch,
                resume=args.resume,
                rerun_failed=args.rerun_failed,
                dry_run=args.dry_run,
            )
            value, details = objective_from_rows(model.key, rows, canonical)
            if args.dry_run:
                write_json(out_root / "trial_state_transitions" / f"{stage}_{model.key}_{opt_trial.number:04d}_dry_run.json", {"state": "dry_run", "objective": details})
                break
            study.tell(opt_trial, value)
            write_json(out_root / "trial_state_transitions" / f"{stage}_{model.key}_{opt_trial.number:04d}_told.json", {"state": "complete", "objective": details})
            for row in rows:
                status_path = Path(row["save_dir"]) / "status.json"
                status = read_json_or_none(status_path) or {}
                status["optuna_value"] = value
                write_json(status_path, status)
            write_master_and_summaries(out_root)
            completed += 1


def build_fixed_trials(stage: str, models: list[ModelFamily], top_n: int, out_root: Path) -> list[TrialSpec]:
    trials: list[TrialSpec] = []
    for model in models:
        if stage in {"stage0", "smoke"}:
            params = dict(model.default_params)
            if stage == "smoke" and model.key == "rg_neighbormix_scmae_contrast_safe":
                params["contrast_start_epoch"] = 1
            trials.append(TrialSpec(stage=stage, model=model, trial_id="trial_canonical", params=params))
        elif stage == "stage4":
            configs = (
                read_selected_configs(out_root, model.key, top_n)
                or read_top_configs(out_root, model.key, "stage3", top_n)
                or read_top_configs(out_root, model.key, "stage2", top_n)
            )
            for idx, params in enumerate(configs[:top_n], start=1):
                trials.append(TrialSpec(stage=stage, model=model, trial_id=f"trial_final_{idx:02d}", params=params))
        elif stage == "stage5":
            configs = read_shortlist_configs(out_root, model.key, top_n)
            for idx, params in enumerate(configs[:top_n], start=1):
                trials.append(TrialSpec(stage=stage, model=model, trial_id=f"trial_full_eval_{idx:02d}", params=params))
        elif stage == "stage7":
            configs = read_top_full16_configs(out_root, model.key, top_n, complete_only=True)
            for idx, params in enumerate(configs[:top_n], start=1):
                trials.append(TrialSpec(stage=stage, model=model, trial_id=f"trial_full_final_{idx:02d}", params=params))
    return trials


def main() -> int:
    parser = argparse.ArgumentParser(description="Run benchmark-tuned Optuna search for NM/RG/RG+CL.")
    parser.add_argument("--data_root", type=Path, default=DATA_ROOT_DEFAULT)
    parser.add_argument("--out_root", type=Path, default=None)
    parser.add_argument("--stages", nargs="+", default=["smoke"], choices=["stage0", "smoke", "stage2", "stage3", "stage4", "stage5", "stage6", "stage7"])
    parser.add_argument("--datasets", nargs="+", default=None)
    parser.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    parser.add_argument("--models", nargs="+", default=MODEL_KEYS)
    parser.add_argument("--gpus", type=int, nargs="+", default=DEFAULT_GPUS)
    parser.add_argument("--max_workers", type=int, default=6)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--smoke_epochs", type=int, default=5)
    parser.add_argument("--stage2_trials", type=int, default=60)
    parser.add_argument("--stage3_trials", type=int, default=80)
    parser.add_argument("--stage3_top_n", type=int, default=5)
    parser.add_argument("--stage4_top_n", type=int, default=3)
    parser.add_argument("--stage5_top_n", type=int, default=12)
    parser.add_argument("--stage6_trials", type=int, default=30)
    parser.add_argument("--stage6_top_n", type=int, default=5)
    parser.add_argument("--stage7_top_n", type=int, default=3)
    parser.add_argument("--python_bin", type=Path, default=PYTHON_BIN_DEFAULT)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--no_resume", dest="resume", action="store_false")
    parser.add_argument("--rerun_failed", action="store_true")
    parser.add_argument("--cleanup_smoke", type=lambda value: value.lower() in {"1", "true", "yes", "y"}, default=True)
    parser.add_argument("--skip_provenance", action="store_true")
    parser.add_argument("--rerun_existing_trial", nargs=3, metavar=("STAGE", "MODEL", "TRIAL_ID"))
    args = parser.parse_args()
    if args.datasets is None:
        args.datasets = FULL16_DATASETS if any(stage in FULL16_STAGES for stage in args.stages) else DEFAULT_DATASETS

    validate_gpus(args.gpus)
    if not args.python_bin.exists():
        raise SystemExit(f"Python executable not found: {args.python_bin}")
    out_root = make_out_root(args.out_root)
    out_root.mkdir(parents=True, exist_ok=True)
    selected_models = [MODEL_FAMILIES[key] for key in args.models]
    for model in selected_models:
        if not model.runner.exists():
            raise SystemExit(f"Runner not found: {model.runner}")

    all_datasets = discover_datasets(args.data_root)
    wanted = set(args.datasets)
    datasets = [dataset for dataset in all_datasets if dataset.name in wanted]
    missing = sorted(wanted - {dataset.name for dataset in datasets})
    if missing:
        raise SystemExit(f"Requested datasets were not discovered: {missing}")
    converted_dir = out_root / "converted_data"
    datasets = [
        planned_prepared_dataset(dataset, converted_dir) if args.dry_run else locked_prepare_dataset(dataset, converted_dir)
        for dataset in datasets
    ]

    commit_sha, branch = git_info()
    if not args.skip_provenance:
        write_provenance(out_root, args, datasets, selected_models)
    print(f"Output root: {out_root}")
    print(f"Datasets: {len(datasets)} Models: {len(selected_models)} Seeds: {args.seeds} Stages: {args.stages}")

    if args.rerun_existing_trial:
        stage, model_key, trial_id = args.rerun_existing_trial
        model = MODEL_FAMILIES[model_key]
        trial = build_existing_trial(out_root, stage, model, trial_id)
        rows = run_trial(
            out_root=out_root,
            trial=trial,
            datasets=datasets,
            seeds=args.seeds,
            python_bin=args.python_bin,
            gpus=args.gpus,
            max_workers=args.max_workers,
            epochs=args.epochs,
            commit_sha=commit_sha,
            branch=branch,
            resume=args.resume,
            rerun_failed=args.rerun_failed,
            dry_run=args.dry_run,
        )
        print(f"rerun {stage} {model.key} {trial.trial_id}: {sum(row['status'] == 'success' for row in rows)}/{len(rows)} success")
        write_master_and_summaries(out_root)
        write_selected_params(out_root)
        write_full16_summaries(out_root)
        return 0

    for stage in args.stages:
        if stage in {"stage0", "smoke", "stage4", "stage5", "stage7"}:
            if stage == "smoke":
                stage_datasets = [dataset for dataset in datasets if dataset.name in {"Pollen", "hrvatin_geo"}]
            elif stage == "stage5":
                stage_datasets = [dataset for dataset in datasets if dataset.name not in set(DEFAULT_DATASETS)]
            else:
                stage_datasets = datasets
            stage_seeds = [42] if stage == "smoke" else args.seeds
            stage_epochs = args.smoke_epochs if stage == "smoke" else args.epochs
            top_n = args.stage4_top_n if stage == "stage4" else args.stage5_top_n if stage == "stage5" else args.stage7_top_n
            for trial in build_fixed_trials(stage, selected_models, top_n, out_root):
                rows = run_trial(
                    out_root=out_root,
                    trial=trial,
                    datasets=stage_datasets,
                    seeds=stage_seeds,
                    python_bin=args.python_bin,
                    gpus=args.gpus,
                    max_workers=args.max_workers,
                    epochs=stage_epochs,
                    commit_sha=commit_sha,
                    branch=branch,
                    resume=args.resume,
                    rerun_failed=args.rerun_failed,
                    dry_run=args.dry_run,
                )
                print(f"{stage} {trial.model.key} {trial.trial_id}: {sum(row['status'] == 'success' for row in rows)}/{len(rows)} success")
                write_master_and_summaries(out_root)
                if stage in FULL16_STAGES:
                    write_full16_summaries(out_root)
            if stage == "smoke" and args.cleanup_smoke and not args.dry_run:
                cleanup_smoke(out_root)
        else:
            if args.dry_run:
                # Import Optuna even in dry-run so missing dependency is caught before long jobs.
                import optuna  # noqa: F401
            run_optuna_stage(stage, out_root, selected_models, datasets, args.seeds, args, commit_sha, branch)
        write_master_and_summaries(out_root)
        write_selected_params(out_root)
        if stage in FULL16_STAGES:
            write_full16_summaries(out_root)

    print(f"Master: {out_root / 'run_master.csv'}")
    print(f"Global summary: {out_root / 'global_config_summary.csv'}")
    print(f"Selected params: {out_root / 'selected_params.json'}")
    if any(stage in FULL16_STAGES for stage in args.stages):
        print(f"Full16 summary: {out_root / 'global_config_summary_full16.csv'}")
        print(f"Full16 selected params: {out_root / 'selected_params_full16.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
