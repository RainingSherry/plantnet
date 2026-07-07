#!/usr/bin/env python3
"""Run selected retired scMAE-family methods on the scMAE benchmark suite."""

from __future__ import annotations

import argparse
import csv
import fcntl
import json
import os
import queue
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
    FORBIDDEN_GPUS,
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
DEFAULT_GPUS = [1, 2]


@dataclass(frozen=True)
class RetiredMethod:
    key: str
    path: Path
    source_path: str
    method_name: str
    variant_name: str
    extra_args: list[str] = field(default_factory=list)
    parameter_signature: str = ""


@dataclass(frozen=True)
class RetiredTask:
    dataset: DatasetSpec
    method: RetiredMethod
    seed: int


RETIRED_METHODS = [
    RetiredMethod(
        key="scmae_dec_stdfloor",
        path=PROJECT_ROOT / "experimental_retired_models" / "scMAE_DEC_StdFloor" / "run.py",
        source_path="experimental_retired_models/scMAE_DEC_StdFloor/run.py",
        method_name="VarFloor-scMAE",
        variant_name="varfloor_scmae",
        parameter_signature=(
            "epochs=80,warmup_epochs=20,n_top_genes=1000,input_mode=auto,"
            "target_sum=10000,scale_input=true,variance_weight=0.02"
        ),
    ),
    RetiredMethod(
        key="rg_neighbormix_scmae",
        path=PROJECT_ROOT / "experimental_retired_models" / "RG_NeighborMix_scMAE" / "run.py",
        source_path="experimental_retired_models/RG_NeighborMix_scMAE/run.py",
        method_name="RG_NeighborMix_scMAE",
        variant_name="rg_nm_v1_reliability",
        extra_args=["--no_save_h5ad"],
        parameter_signature=(
            "epochs=80,n_top_genes=1000,input_mode=auto,target_sum=10000,"
            "scale_input=true,mix_mode=reliability,gate_mode=topology,"
            "edge_reliability_mode=sim_mutual_snn_distance"
        ),
    ),
    RetiredMethod(
        key="neighbormix_scmae",
        path=PROJECT_ROOT / "experimental_retired_models" / "NeighborMix_scMAE" / "run.py",
        source_path="experimental_retired_models/NeighborMix_scMAE/run.py",
        method_name="NeighborMix_scMAE",
        variant_name="nm_scmae_mid",
        extra_args=["--no_save_h5ad"],
        parameter_signature=(
            "epochs=80,n_top_genes=1000,input_mode=auto,target_sum=10000,"
            "scale_input=true,hidden_size=128,dropout=0.0,masked_data_weight=0.75,"
            "mask_loss_weight=0.7,mask_ratio=0.4,use_pseudo=true,pseudo_weight=0.3,"
            "alpha=0.9,neighbor_k=5,mix_neighbors=4,tau=0.2,knn_pca_dim=50"
        ),
    ),
    RetiredMethod(
        key="canm_cut_reweighted_mix",
        path=PROJECT_ROOT / "experimental_retired_models" / "CutAware_NeighborMix_scMAE" / "run.py",
        source_path="experimental_retired_models/CutAware_NeighborMix_scMAE/run.py",
        method_name="CutAware_NeighborMix_scMAE",
        variant_name="canm_cut_reweighted_mix",
        extra_args=["--no_save_h5ad"],
        parameter_signature=(
            "epochs=80,n_top_genes=1000,input_mode=auto,target_sum=10000,"
            "scale_input=true,hidden_size=128,dropout=0.0,masked_data_weight=0.75,"
            "mask_loss_weight=0.7,mask_ratio=0.4,neighbor_k=10,mix_neighbors=4,"
            "mix_alpha=0.8,edge_reliability_mode=sim_mutual_snn_distance,"
            "cut_cross_weight=0.05,pseudo_weight=0.2,cut_weight=0.0,ot_weight=0.0"
        ),
    ),
    RetiredMethod(
        key="canm_cut_reweighted_mix_contrast",
        path=PROJECT_ROOT / "experimental_retired_models" / "CutAware_NeighborMix_scMAE" / "run.py",
        source_path="experimental_retired_models/CutAware_NeighborMix_scMAE/run.py",
        method_name="CutAware_NeighborMix_scMAE",
        variant_name="canm_cut_reweighted_mix_contrast",
        extra_args=[
            "--no_save_h5ad",
            "--contrast_weight",
            "0.05",
            "--contrast_temperature",
            "0.2",
            "--contrast_start_epoch",
            "5",
            "--contrast_neighbor_positive_weight",
            "0.5",
            "--contrast_hard_negative_weight",
            "1.0",
            "--contrast_projection_dim",
            "128",
        ],
        parameter_signature=(
            "epochs=80,n_top_genes=1000,input_mode=auto,target_sum=10000,"
            "scale_input=true,hidden_size=128,dropout=0.0,masked_data_weight=0.75,"
            "mask_loss_weight=0.7,mask_ratio=0.4,neighbor_k=10,mix_neighbors=4,"
            "mix_alpha=0.8,edge_reliability_mode=sim_mutual_snn_distance,"
            "cut_cross_weight=0.05,pseudo_weight=0.2,cut_weight=0.0,ot_weight=0.0,"
            "contrast_weight=0.05,contrast_temperature=0.2,contrast_start_epoch=5,"
            "contrast_neighbor_positive_weight=0.5,contrast_hard_negative_weight=1.0,"
            "contrast_projection_dim=128"
        ),
    ),
]


def make_out_root(out_root: Path | None) -> Path:
    if out_root:
        return out_root
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return PROJECT_ROOT / "result" / f"scmae_retired_methods_{stamp}"


def task_out_dir(out_root: Path, task: RetiredTask) -> Path:
    return out_root / task.dataset.name / task.method.key / f"seed_{task.seed}"


def build_tasks(datasets: list[DatasetSpec], methods: list[RetiredMethod], seeds: list[int]) -> list[RetiredTask]:
    return [RetiredTask(dataset, method, seed) for dataset in datasets for method in methods for seed in seeds]


def locked_prepare_dataset(dataset: DatasetSpec, converted_dir: Path) -> DatasetSpec:
    converted_dir.mkdir(parents=True, exist_ok=True)
    lock_path = converted_dir / f".{dataset.name}.prepare.lock"
    with lock_path.open("w", encoding="utf-8") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        try:
            return prepare_dataset(dataset, converted_dir)
        finally:
            fcntl.flock(lock_handle, fcntl.LOCK_UN)


def command_for(
    python_bin: Path,
    task: RetiredTask,
    out_dir: Path,
    gpu: int | None,
    epochs: int,
) -> list[str]:
    if task.dataset.prepared_path is None or task.dataset.n_clusters is None:
        raise RuntimeError(f"Dataset was not prepared: {task.dataset.name}")
    cmd = [
        str(python_bin),
        str(task.method.path),
        "--data_path",
        str(task.dataset.prepared_path),
        "--save_dir",
        str(out_dir),
        "--dataset_name",
        task.dataset.name,
        "--method_name",
        task.method.method_name,
        "--variant_name",
        task.method.variant_name,
        "--n_clusters",
        str(task.dataset.n_clusters),
        "--seed",
        str(task.seed),
        "--epochs",
        str(epochs),
        "--gpu",
        str(gpu if gpu is not None else 1),
    ]
    cmd.extend(task.method.extra_args)
    return cmd


def load_metrics(out_dir: Path) -> dict[str, Any]:
    path = out_dir / "metrics.json"
    if not path.exists():
        return {}
    return read_json_or_none(path) or {}


def required_artifacts_missing(out_dir: Path) -> list[str]:
    return [
        name
        for name in ["metrics.json", "embedding_final.npy", "labels.npy", "args.json"]
        if not (out_dir / name).exists()
    ]


def write_environment_json(out_dir: Path, python_bin: Path, gpu: int | None, method: RetiredMethod) -> None:
    payload = {
        "python_executable": str(python_bin),
        "runtime_env": "scclubench-sccdcg-h100",
        "framework": "PyTorch",
        "category": "experimental_retired_models",
        "gpu": gpu if gpu is not None else "cpu",
        "retired": True,
        "source_path": method.source_path,
    }
    try:
        probe = subprocess.run(
            [str(python_bin), "-c", "import sys, torch; print(sys.version.split()[0]); print(torch.__version__)"],
            text=True,
            capture_output=True,
            timeout=10,
        )
        if probe.returncode == 0:
            lines = [line.strip() for line in probe.stdout.splitlines() if line.strip()]
            if lines:
                payload["python_version"] = lines[0]
            if len(lines) > 1:
                payload["torch_version"] = lines[1]
        else:
            payload["probe_error"] = probe.stderr[-1000:]
    except Exception as exc:
        payload["probe_error"] = repr(exc)
    write_json(out_dir / "environment.json", payload)


def base_row(task: RetiredTask, out_dir: Path, commit_sha: str, branch: str) -> dict[str, Any]:
    return {
        "dataset": task.dataset.name,
        "model": task.method.key,
        "seed": task.seed,
        "status": "pending",
        **{metric: "" for metric in METRIC_NAMES},
        "runtime_seconds": "",
        "gpu": "",
        "command": "",
        "save_dir": str(out_dir),
        "error": "",
        "fallback": False,
        "substitute_model_used": False,
        "commit_sha": commit_sha,
        "branch": branch,
        "retired": True,
        "source_path": task.method.source_path,
        "variant_name": task.method.variant_name,
        "reused_existing": False,
        "reused_from": "",
        "parameter_signature": task.method.parameter_signature,
        "python_executable": "",
    }


def row_from_status(task: RetiredTask, out_root: Path, commit_sha: str, branch: str) -> dict[str, Any]:
    out_dir = task_out_dir(out_root, task)
    row = base_row(task, out_dir, commit_sha, branch)
    status = read_json_or_none(out_dir / "status.json")
    if status:
        row.update(
            {
                "status": status.get("status", "unknown"),
                "runtime_seconds": status.get("runtime_seconds", ""),
                "gpu": status.get("gpu", ""),
                "command": status.get("command", ""),
                "error": status.get("error", ""),
                "fallback": status.get("fallback", False),
                "substitute_model_used": status.get("substitute_model_used", False),
                "reused_existing": status.get("reused_existing", False),
                "reused_from": status.get("reused_from", ""),
                "python_executable": status.get("python_executable", ""),
            }
        )
    metrics = normalize_metrics(load_metrics(out_dir))
    for metric in METRIC_NAMES:
        row[metric] = metrics.get(metric, "")
    return row


def read_existing_master_rows(out_root: Path) -> list[dict[str, Any]]:
    master_path = out_root / "retired_all_runs_master.csv"
    if not master_path.exists():
        return []
    with master_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def merge_master_rows(existing: list[dict[str, Any]], updates: list[dict[str, Any]]) -> list[dict[str, Any]]:
    merged: dict[tuple[str, str, str], dict[str, Any]] = {}
    for row in existing:
        key = (str(row.get("dataset", "")), str(row.get("model", "")), str(row.get("seed", "")))
        if all(key):
            merged[key] = row
    for row in updates:
        key = (str(row.get("dataset", "")), str(row.get("model", "")), str(row.get("seed", "")))
        merged[key] = row
    return sorted(merged.values(), key=lambda row: (str(row.get("dataset", "")), str(row.get("model", "")), int(row.get("seed", 0) or 0)))


def write_master_tables(
    out_root: Path,
    tasks: list[RetiredTask],
    commit_sha: str,
    branch: str,
    merge_existing_master: bool = False,
) -> None:
    rows = [row_from_status(task, out_root, commit_sha, branch) for task in tasks]
    if merge_existing_master:
        rows = merge_master_rows(read_existing_master_rows(out_root), rows)
    if not rows:
        return
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    with (out_root / "retired_all_runs_master.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    write_json(out_root / "retired_all_runs_master.json", rows)
    write_summary_tables(out_root, rows)


def write_summary_tables(out_root: Path, rows: list[dict[str, Any]]) -> None:
    import statistics
    from collections import defaultdict

    valid = [row for row in rows if row.get("status") == "success"]
    for group_fields, filename in [
        (("model", "dataset"), "summary_by_model_dataset.csv"),
        (("model",), "summary_by_model.csv"),
    ]:
        grouped: dict[tuple[Any, ...], list[dict[str, Any]]] = defaultdict(list)
        for row in valid:
            grouped[tuple(row[field] for field in group_fields)].append(row)
        out_rows: list[dict[str, Any]] = []
        for key, group in sorted(grouped.items()):
            out: dict[str, Any] = {field: value for field, value in zip(group_fields, key)}
            out["n_success"] = len(group)
            for metric in METRIC_NAMES:
                vals = [float(row[metric]) for row in group if row.get(metric) not in {"", None}]
                if vals:
                    out[f"{metric}_mean"] = statistics.mean(vals)
                    out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
            out_rows.append(out)
        if out_rows:
            fields = list(out_rows[0].keys())
            with (out_root / filename).open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields)
                writer.writeheader()
                writer.writerows(out_rows)
        else:
            (out_root / filename).write_text("", encoding="utf-8")


def write_catalogs(out_root: Path, datasets: list[DatasetSpec], methods: list[RetiredMethod], tasks: list[RetiredTask]) -> None:
    write_json(out_root / "dataset_catalog.json", [dataset.__dict__ for dataset in datasets])
    write_json(
        out_root / "method_catalog.json",
        [
            {
                "key": method.key,
                "path": str(method.path.relative_to(PROJECT_ROOT)),
                "source_path": method.source_path,
                "variant_name": method.variant_name,
                "retired": True,
                "parameter_signature": method.parameter_signature,
            }
            for method in methods
        ],
    )


def write_subset_manifest(out_root: Path, datasets: list[DatasetSpec], methods: list[RetiredMethod], tasks: list[RetiredTask]) -> None:
    manifest_dir = out_root / "repair_manifests"
    manifest_dir.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    token = "_".join(method.key for method in methods)
    write_json(
        manifest_dir / f"task_manifest_{token}_{stamp}.json",
        [
            {
                "dataset": task.dataset.name,
                "model": task.method.key,
                "seed": task.seed,
                "save_dir": str(task_out_dir(out_root, task)),
            }
            for task in tasks
        ],
    )
    write_json(
        manifest_dir / f"subset_catalog_{token}_{stamp}.json",
        {
            "datasets": [dataset.__dict__ for dataset in datasets],
            "methods": [
                {
                    "key": method.key,
                    "path": str(method.path.relative_to(PROJECT_ROOT)),
                    "source_path": method.source_path,
                    "variant_name": method.variant_name,
                    "retired": True,
                    "parameter_signature": method.parameter_signature,
                }
                for method in methods
            ],
        },
    )
    write_json(
        out_root / "task_manifest.json",
        [
            {
                "dataset": task.dataset.name,
                "model": task.method.key,
                "seed": task.seed,
                "save_dir": str(task_out_dir(out_root, task)),
            }
            for task in tasks
        ],
    )


def scan_existing_candidates(out_root: Path) -> None:
    roots = [PROJECT_ROOT / "experiment_reports", PROJECT_ROOT / "results"]
    needles = {
        "varfloor": "scmae_dec_stdfloor",
        "stdfloor": "scmae_dec_stdfloor",
        "rg_neighbormix": "rg_neighbormix_scmae",
        "rg_reliability": "rg_neighbormix_scmae",
    }
    rows: list[dict[str, Any]] = []
    for root in roots:
        if not root.exists():
            continue
        try:
            found = subprocess.run(
                ["find", str(root), "-maxdepth", "8", "-type", "f"],
                text=True,
                capture_output=True,
                timeout=60,
            )
            candidates = [Path(line) for line in found.stdout.splitlines() if line.strip()] if found.returncode == 0 else []
        except Exception:
            candidates = []
        for path in candidates:
            lower = str(path).lower()
            matched = [model for token, model in needles.items() if token in lower]
            if not matched:
                continue
            if path.name not in {"metrics.json", "summary.json", "eval_fixed.csv", "args.json", "status.json"}:
                continue
            rows.append({"model_hint": matched[0], "path": str(path), "reused": False, "reason": "candidate_only"})
            if len(rows) >= 5000:
                break
    with (out_root / "existing_result_candidates.csv").open("w", newline="", encoding="utf-8") as handle:
        fieldnames = ["model_hint", "path", "reused", "reason"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def run_task(
    task: RetiredTask,
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
    out_dir = task_out_dir(out_root, task)
    status_path = out_dir / "status.json"
    if resume and status_path.exists():
        previous = read_json_or_none(status_path) or {}
        previous_status = str(previous.get("status", "unknown"))
        if previous_status == "success":
            return row_from_status(task, out_root, commit_sha, branch)
        if previous_status in {"failed", "timeout", "error", "incomplete"} and not rerun_failed:
            return row_from_status(task, out_root, commit_sha, branch)

    out_dir.mkdir(parents=True, exist_ok=True)
    gpu = gpu_queue.get()
    try:
        cmd = command_for(python_bin, task, out_dir, gpu, epochs)
        write_json(
            out_dir / "args.json",
            {
                "dataset": task.dataset.name,
                "model": task.method.key,
                "seed": task.seed,
                "n_clusters": task.dataset.n_clusters,
                "data_path": str(task.dataset.prepared_path),
                "retired": True,
                "variant_name": task.method.variant_name,
                "parameter_signature": task.method.parameter_signature,
                "command": cmd,
            },
        )
        write_environment_json(out_dir, python_bin, gpu, task.method)
        (out_dir / "command.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")
        start = datetime.now().isoformat()
        if dry_run:
            status = {
                "dataset": task.dataset.name,
                "method": task.method.key,
                "seed": task.seed,
                "status": "dry_run",
                "runtime_seconds": 0.0,
                "gpu": gpu,
                "command": " ".join(cmd),
                "error": "",
                "fallback": False,
                "substitute_model_used": False,
                "reused_existing": False,
                "reused_from": "",
                "python_executable": str(python_bin),
                "commit_sha": commit_sha,
                "branch": branch,
                "start_time": start,
                "end_time": datetime.now().isoformat(),
            }
            write_json(status_path, status)
            return row_from_status(task, out_root, commit_sha, branch)

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
        t0 = time.time()
        try:
            proc = subprocess.run(
                cmd,
                cwd=PROJECT_ROOT,
                text=True,
                capture_output=True,
                timeout=TIMEOUT_SECONDS,
                env=env,
            )
            elapsed = round(time.time() - t0, 1)
            log_text = (
                f"Command: {' '.join(cmd)}\n"
                f"Exit code: {proc.returncode}\n\n"
                f"=== STDOUT ===\n{proc.stdout}\n\n"
                f"=== STDERR ===\n{proc.stderr}\n"
            )
            (out_dir / "run.log").write_text(log_text, encoding="utf-8")
            missing = required_artifacts_missing(out_dir)
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
            partial_stdout = exc.stdout.decode() if isinstance(exc.stdout, bytes) else (exc.stdout or "")
            partial_stderr = exc.stderr.decode() if isinstance(exc.stderr, bytes) else (exc.stderr or "")
            (out_dir / "run.log").write_text(
                f"Command timed out after {TIMEOUT_SECONDS}s\nCommand: {' '.join(cmd)}\n"
                f"\n=== STDOUT ===\n{partial_stdout}\n\n=== STDERR ===\n{partial_stderr}\n",
                encoding="utf-8",
            )
            state = "timeout"
            error = f"Execution exceeded {TIMEOUT_SECONDS}s"
            elapsed = round(time.time() - t0, 1)

        status = {
            "dataset": task.dataset.name,
            "method": task.method.key,
            "seed": task.seed,
            "status": state,
            "runtime_seconds": elapsed,
            "gpu": gpu,
            "command": " ".join(cmd),
            "error": error,
            "fallback": False,
            "substitute_model_used": False,
            "reused_existing": False,
            "reused_from": "",
            "python_executable": str(python_bin),
            "commit_sha": commit_sha,
            "branch": branch,
            "start_time": start,
            "end_time": datetime.now().isoformat(),
        }
        write_json(status_path, status)
        return row_from_status(task, out_root, commit_sha, branch)
    finally:
        gpu_queue.put(gpu)


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run retired scMAE-family methods on the scMAE benchmark suite.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data_root", type=Path, default=DATA_ROOT_DEFAULT)
    parser.add_argument("--out_root", type=Path, default=None)
    parser.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    parser.add_argument("--gpus", type=int, nargs="+", default=DEFAULT_GPUS)
    parser.add_argument("--max_workers", type=int, default=2)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--datasets", nargs="+", default=None)
    parser.add_argument("--methods", nargs="+", default=None)
    parser.add_argument("--python_bin", type=Path, default=PYTHON_BIN_DEFAULT)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--no_resume", dest="resume", action="store_false")
    parser.add_argument("--rerun_failed", action="store_true")
    parser.add_argument("--merge_existing_master", action="store_true")
    parser.add_argument("--skip_existing_scan", action="store_true")
    args = parser.parse_args()

    validate_gpus(args.gpus)
    if not args.python_bin.exists():
        raise SystemExit(f"Python executable not found: {args.python_bin}")

    out_root = make_out_root(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    converted_dir = out_root / "converted_data"
    commit_sha, branch = git_info()

    datasets = discover_datasets(args.data_root)
    if args.datasets:
        wanted = set(args.datasets)
        datasets = [dataset for dataset in datasets if dataset.name in wanted]
    if not datasets:
        raise SystemExit(f"No datasets discovered under {args.data_root}")
    datasets = [
        planned_prepared_dataset(dataset, converted_dir) if args.dry_run else locked_prepare_dataset(dataset, converted_dir)
        for dataset in datasets
    ]

    selected = set(args.methods) if args.methods else None
    methods = [method for method in RETIRED_METHODS if selected is None or method.key in selected]
    if not methods:
        raise SystemExit("No retired methods selected.")
    for method in methods:
        if not method.path.exists():
            raise SystemExit(f"Retired method runner not found: {method.path}")

    tasks = build_tasks(datasets, methods, args.seeds)
    if args.merge_existing_master:
        write_subset_manifest(out_root, datasets, methods, tasks)
    else:
        write_catalogs(out_root, datasets, methods, tasks)
    write_json(
        out_root / "suite_config.json",
        {
            "data_root": str(args.data_root),
            "out_root": str(out_root),
            "seeds": args.seeds,
            "gpus": args.gpus,
            "max_workers": args.max_workers,
            "epochs": args.epochs,
            "dry_run": args.dry_run,
            "commit_sha": commit_sha,
            "branch": branch,
            "python_bin": str(args.python_bin),
            "retired": True,
            "merge_existing_master": args.merge_existing_master,
        },
    )
    if not args.skip_existing_scan:
        scan_existing_candidates(out_root)

    print(f"Retired suite output: {out_root}")
    print(f"Datasets: {len(datasets)}")
    print(f"Methods:  {len(methods)}")
    print(f"Seeds:    {args.seeds}")
    print(f"Tasks:    {len(tasks)}")
    print(f"Dry run:  {args.dry_run}")
    write_master_tables(out_root, tasks, commit_sha, branch, args.merge_existing_master)

    gpu_queue: queue.Queue[int] = queue.Queue()
    for gpu in args.gpus:
        gpu_queue.put(gpu)

    worker_count = max(1, min(int(args.max_workers), len(args.gpus), len(tasks)))
    lock = threading.Lock()
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = [
            executor.submit(
                run_task,
                task,
                out_root,
                args.python_bin,
                gpu_queue,
                args.epochs,
                commit_sha,
                branch,
                args.resume,
                args.rerun_failed,
                args.dry_run,
            )
            for task in tasks
        ]
        for index, future in enumerate(as_completed(futures), start=1):
            try:
                row = future.result()
                print(f"[{index}/{len(futures)}] {row['dataset']} / {row['model']} / seed {row['seed']}: {row['status']}")
            except Exception as exc:
                print(f"[{index}/{len(futures)}] task wrapper error: {exc}", file=sys.stderr)
            if index == len(futures) or index % 10 == 0:
                with lock:
                    write_master_tables(out_root, tasks, commit_sha, branch, args.merge_existing_master)

    write_master_tables(out_root, tasks, commit_sha, branch, args.merge_existing_master)
    print(f"Master table: {out_root / 'retired_all_runs_master.csv'}")
    print(f"Summary by model/dataset: {out_root / 'summary_by_model_dataset.csv'}")
    print(f"Summary by model: {out_root / 'summary_by_model.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
