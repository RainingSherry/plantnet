#!/usr/bin/env python3
"""Run every active method runner on the scMAE benchmark data suite.

This is an outer suite runner for the user-requested all-method benchmark.  It
does not replace ``scripts/run_formal_benchmark.py``; it discovers active
``methods/**/run.py`` files, prepares datasets, schedules model runs over GPUs
1-6, and maintains a master table with one row per dataset/model/seed.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import queue
import re
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

import yaml


PROJECT_ROOT = Path(__file__).resolve().parent.parent
METHODS_DIR = PROJECT_ROOT / "methods"
MANIFEST_PATH = METHODS_DIR / "method_manifest.yaml"
DATA_ROOT_DEFAULT = Path("/data/luolie/biopipeline/scCluBench/data/scMAE")
DEFAULT_SEEDS = [42, 3047, 3407]
DEFAULT_GPUS = [1, 2, 3, 4, 5, 6]
FORBIDDEN_GPUS = {0, 7}
TIMEOUT_SECONDS = 6 * 3600
METRIC_NAMES = [
    "acc",
    "nmi",
    "ari",
    "f1_macro",
    "fmi",
    "v_measure",
    "homogeneity",
    "completeness",
]


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    source_path: Path
    kind: str
    prepared_path: Path | None = None
    n_clusters: int | None = None


@dataclass(frozen=True)
class MethodSpec:
    key: str
    path: Path
    category: str
    framework: str = ""
    authenticity: str = "UNREGISTERED"
    runtime_env: str = ""
    extra_args: list[str] = field(default_factory=list)
    source_path: str = ""
    accepted_args: set[str] = field(default_factory=set)
    manifest_entry: dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Task:
    dataset: DatasetSpec
    method: MethodSpec
    seed: int


def _json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    return str(value)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = path.with_name(f".{path.name}.{os.getpid()}.{threading.get_ident()}.tmp")
    with tmp_path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=_json_default)
    tmp_path.replace(path)


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def read_json_or_none(path: Path) -> dict[str, Any] | None:
    try:
        return read_json(path)
    except Exception:
        return None


def load_manifest() -> dict[str, dict[str, Any]]:
    if not MANIFEST_PATH.exists():
        return {}
    with MANIFEST_PATH.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return {item["key"]: item for item in data.get("methods", []) if "key" in item}


def manifest_by_path(manifest: dict[str, dict[str, Any]]) -> dict[str, dict[str, Any]]:
    by_path: dict[str, dict[str, Any]] = {}
    for item in manifest.values():
        path = item.get("path")
        if path:
            by_path[str(Path(path))] = item
    return by_path


def category_from_path(path: Path) -> str:
    rel = path.relative_to(METHODS_DIR)
    return rel.parts[0] if rel.parts else "Unknown"


def key_from_path(path: Path) -> str:
    rel = path.relative_to(METHODS_DIR)
    parts = rel.parts[:-1]
    if len(parts) >= 3 and parts[0] == "GNN" and parts[1] == "scDSC":
        return "_".join(parts[1:]).lower()
    name = parts[-1] if parts else path.parent.name
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def parse_accepted_args(path: Path) -> set[str]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    return set(re.findall(r"add_argument\(\s*[\"'](--[A-Za-z0-9_-]+)[\"']", text))


def get_runtime_env(entry: dict[str, Any]) -> str:
    runtime = entry.get("runtime") or {}
    if isinstance(runtime, dict) and runtime.get("env_name"):
        return str(runtime["env_name"])
    return str(entry.get("runtime_env") or "")


def discover_methods(selected_keys: set[str] | None = None) -> list[MethodSpec]:
    manifest = load_manifest()
    by_path = manifest_by_path(manifest)
    methods: list[MethodSpec] = []
    for run_py in sorted(METHODS_DIR.rglob("run.py")):
        rel = run_py.relative_to(PROJECT_ROOT)
        if "_archive" in rel.parts or "__pycache__" in rel.parts:
            continue
        entry = by_path.get(str(rel), {})
        key = str(entry.get("key") or key_from_path(run_py))
        if selected_keys and key not in selected_keys:
            continue
        category = str(entry.get("category") or category_from_path(run_py))
        methods.append(
            MethodSpec(
                key=key,
                path=run_py,
                category=category,
                framework=str(entry.get("framework") or ""),
                authenticity=str(entry.get("authenticity") or "UNREGISTERED"),
                runtime_env=get_runtime_env(entry),
                extra_args=[str(x) for x in entry.get("extra_args", [])],
                source_path=str(entry.get("source_path") or ""),
                accepted_args=parse_accepted_args(run_py),
                manifest_entry=entry,
            )
        )
    return methods


def discover_datasets(data_root: Path) -> list[DatasetSpec]:
    datasets = [
        DatasetSpec(name=path.stem, source_path=path, kind="h5")
        for path in sorted(data_root.glob("*.h5"))
    ]
    geo_dir = data_root / "hrvatin_geo"
    matrix = geo_dir / "GSE102827_MATRIX.csv.gz"
    labels = geo_dir / "GSE102827_cell_type_assignments.csv.gz"
    if matrix.exists() and labels.exists():
        datasets.append(DatasetSpec(name="hrvatin_geo", source_path=geo_dir, kind="hrvatin_geo"))
    return datasets


def infer_label_col(adata) -> str:
    for candidate in ["resolved_label", "maintype", "celltype", "cell_type", "Celltype", "label", "Y"]:
        if candidate in adata.obs.columns:
            return candidate
    raise ValueError(f"Cannot infer label column; obs columns={list(adata.obs.columns)}")


def infer_h5ad_n_clusters(path: Path) -> int:
    import scanpy as sc

    adata = sc.read_h5ad(path)
    label_col = infer_label_col(adata)
    return int(adata.obs[label_col].astype(str).nunique())


def infer_h5_n_clusters(path: Path) -> int:
    import h5py
    import pandas as pd

    with h5py.File(path, "r") as handle:
        label_key = None
        for candidate in ["Y", "label", "labels", "cell_type", "celltype", "cluster", "type"]:
            if candidate in handle:
                label_key = candidate
                break
        if label_key is None:
            raise ValueError(f"Cannot infer label key from {path}; keys={list(handle.keys())}")
        labels = handle[label_key][...]
    return int(pd.unique(labels.astype(str)).size)


def infer_hrvatin_geo_n_clusters(path: Path) -> int:
    import pandas as pd

    labels = pd.read_csv(path / "GSE102827_cell_type_assignments.csv.gz", index_col=0, usecols=[0, 3])
    return int(labels["maintype"].astype(str).nunique())


def planned_prepared_dataset(dataset: DatasetSpec, converted_dir: Path) -> DatasetSpec:
    if dataset.kind == "h5":
        return DatasetSpec(
            name=dataset.name,
            source_path=dataset.source_path,
            kind=dataset.kind,
            prepared_path=converted_dir / f"{dataset.name}.h5ad",
            n_clusters=infer_h5_n_clusters(dataset.source_path),
        )
    if dataset.kind == "hrvatin_geo":
        return DatasetSpec(
            name=dataset.name,
            source_path=dataset.source_path,
            kind=dataset.kind,
            prepared_path=converted_dir / "hrvatin_geo.h5ad",
            n_clusters=infer_hrvatin_geo_n_clusters(dataset.source_path),
        )
    if dataset.kind == "h5ad":
        return DatasetSpec(
            name=dataset.name,
            source_path=dataset.source_path,
            kind=dataset.kind,
            prepared_path=dataset.source_path,
            n_clusters=infer_h5ad_n_clusters(dataset.source_path),
        )
    raise ValueError(f"Unsupported dataset kind: {dataset.kind}")


def converter_python() -> str:
    candidate = Path("/data/luolie/conda/envs/scclubench-main/bin/python")
    return str(candidate) if candidate.exists() else sys.executable


def prepare_h5_dataset(dataset: DatasetSpec, converted_dir: Path) -> DatasetSpec:
    converted_dir.mkdir(parents=True, exist_ok=True)
    out_path = converted_dir / f"{dataset.name}.h5ad"
    meta_path = converted_dir / f"{dataset.name}.meta.json"
    if out_path.exists() and meta_path.exists():
        meta = read_json(meta_path)
        return DatasetSpec(
            name=dataset.name,
            source_path=dataset.source_path,
            kind=dataset.kind,
            prepared_path=out_path,
            n_clusters=int(meta["n_clusters"]),
        )

    cmd = [
        converter_python(),
        str(PROJECT_ROOT / "scripts" / "prepare_dataset.py"),
        "--input_path",
        str(dataset.source_path),
        "--dataset_name",
        dataset.name,
        "--output_dir",
        str(converted_dir),
        "--force",
    ]
    env = os.environ.copy()
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("OPENBLAS_NUM_THREADS", "4")
    env.setdefault("MKL_NUM_THREADS", "4")
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)
    result = subprocess.run(cmd, cwd=PROJECT_ROOT, text=True, capture_output=True, timeout=1800, env=env)
    if result.returncode != 0:
        raise RuntimeError(
            f"prepare_dataset.py failed for {dataset.name}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    meta = read_json(meta_path)
    return DatasetSpec(
        name=dataset.name,
        source_path=dataset.source_path,
        kind=dataset.kind,
        prepared_path=out_path,
        n_clusters=int(meta["n_clusters"]),
    )


def prepare_hrvatin_geo_dataset(dataset: DatasetSpec, converted_dir: Path) -> DatasetSpec:
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    converted_dir.mkdir(parents=True, exist_ok=True)
    out_path = converted_dir / "hrvatin_geo.h5ad"
    meta_path = converted_dir / "hrvatin_geo.meta.json"
    if out_path.exists() and meta_path.exists():
        meta = read_json(meta_path)
        return DatasetSpec(
            name=dataset.name,
            source_path=dataset.source_path,
            kind=dataset.kind,
            prepared_path=out_path,
            n_clusters=int(meta["n_clusters"]),
        )

    matrix_path = dataset.source_path / "GSE102827_MATRIX.csv.gz"
    label_path = dataset.source_path / "GSE102827_cell_type_assignments.csv.gz"
    labels = pd.read_csv(label_path, index_col=0)
    if "maintype" not in labels.columns:
        raise ValueError(f"hrvatin_geo labels do not contain maintype: {list(labels.columns)}")

    # The source matrix is gene x cell, so transpose it into AnnData's cell x gene shape.
    matrix = pd.read_csv(matrix_path, index_col=0)
    common_cells = labels.index.intersection(matrix.columns)
    if common_cells.empty:
        raise ValueError("hrvatin_geo has no overlapping cells between matrix columns and labels")
    matrix = matrix.loc[:, common_cells]
    labels = labels.loc[common_cells].copy()
    x_csr = sp.csr_matrix(matrix.T.values.astype("float32"))
    adata = ad.AnnData(X=x_csr, obs=labels, var=pd.DataFrame(index=matrix.index.astype(str)))
    adata.obs_names = common_cells.astype(str)
    adata.var_names = matrix.index.astype(str)
    adata.uns["source_format"] = "csv_gz"
    adata.uns["source_file"] = str(dataset.source_path)
    adata.uns["resolved_label_key"] = "maintype"
    adata.uns["matrix_key"] = "GSE102827_MATRIX.csv.gz"
    adata.uns["n_clusters"] = int(labels["maintype"].astype(str).nunique())
    adata.write_h5ad(out_path)

    meta = {
        "dataset_name": dataset.name,
        "source_path": str(dataset.source_path),
        "output_path": str(out_path),
        "label_key": "maintype",
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_clusters": int(adata.uns["n_clusters"]),
    }
    write_json(meta_path, meta)
    return DatasetSpec(
        name=dataset.name,
        source_path=dataset.source_path,
        kind=dataset.kind,
        prepared_path=out_path,
        n_clusters=int(meta["n_clusters"]),
    )


def prepare_dataset(dataset: DatasetSpec, converted_dir: Path) -> DatasetSpec:
    if dataset.kind == "h5":
        return prepare_h5_dataset(dataset, converted_dir)
    if dataset.kind == "hrvatin_geo":
        return prepare_hrvatin_geo_dataset(dataset, converted_dir)
    if dataset.kind == "h5ad":
        return DatasetSpec(
            name=dataset.name,
            source_path=dataset.source_path,
            kind=dataset.kind,
            prepared_path=dataset.source_path,
            n_clusters=infer_h5ad_n_clusters(dataset.source_path),
        )
    raise ValueError(f"Unsupported dataset kind: {dataset.kind}")


def load_runtime_registry(path: Path | None) -> dict[str, dict[str, Any]]:
    registry: dict[str, dict[str, Any]] = {}
    if path and path.exists():
        with path.open("r", encoding="utf-8") as handle:
            data = yaml.safe_load(handle) or {}
        registry.update(data.get("runtimes", {}))

    candidates = {
        "plantnet-core": "/data/luolie/conda/envs/scclubench-main/bin/python",
        "plantnet-tf1": "/data/luolie/conda/envs/scclubench-main/bin/python",
        "plantnet-tf212-cpu": "/data/luolie/conda/envs/scclubench-tf212-cpu/bin/python",
        "plantnet-desc": "/data/luolie/conda/envs/scclubench-main/bin/python",
        "plantnet-scgnn": "/data/luolie/conda/envs/scclubench-main/bin/python",
        "plantnet-attentionae": "/data/luolie/conda/envs/scclubench-main/bin/python",
        "plantnet-sccdcg-h100": "/data/luolie/conda/envs/scclubench-sccdcg-h100/bin/python",
        "plantnet-sccdcg": "/data/luolie/conda/envs/scclubench-sccdcg/bin/python",
        "plantnet-scrcl": "/data/luolie/conda/envs/plantnet-scrcl/bin/python",
        "foundation": "/data/luolie/conda/envs/scclubench-foundation/bin/python",
        "phytocluster": "/data/luolie/conda/envs/PhytoCluster/bin/python",
    }
    for name, python_path in candidates.items():
        registry.setdefault(name, {"python": python_path, "backend": "conda", "auto": True})
    return registry


def resolve_python(method: MethodSpec, registry: dict[str, dict[str, Any]]) -> tuple[str, str]:
    runtime_env = method.runtime_env
    if not runtime_env:
        if method.category == "Foundation":
            runtime_env = "foundation"
        elif method.key == "phytocluster":
            runtime_env = "phytocluster"
        elif method.key == "sccdcg":
            runtime_env = "plantnet-sccdcg"
        else:
            runtime_env = "plantnet-core"
    candidate = registry.get(runtime_env, {}).get("python", "")
    if candidate and Path(str(candidate)).exists():
        return str(candidate), runtime_env
    return sys.executable, runtime_env


def is_cpu_method(method: MethodSpec) -> bool:
    return method.category == "Traditional"


def build_command(
    python_bin: str,
    task: Task,
    prepared_data: DatasetSpec,
    out_dir: Path,
    gpu: int | None,
    epochs: int,
    pretrain_epochs: int,
) -> list[str]:
    accepted = task.method.accepted_args
    cmd = [
        python_bin,
        str(task.method.path),
        "--data_path",
        str(prepared_data.prepared_path),
        "--save_dir",
        str(out_dir),
        "--n_clusters",
        str(prepared_data.n_clusters),
        "--seed",
        str(task.seed),
    ]
    if "--dataset_name" in accepted:
        cmd.extend(["--dataset_name", task.dataset.name])
    if "--epochs" in accepted:
        cmd.extend(["--epochs", str(epochs)])
    if "--pretrain_epochs" in accepted:
        cmd.extend(["--pretrain_epochs", str(pretrain_epochs)])
    if "--cluster_max_iter" in accepted:
        cmd.extend(["--cluster_max_iter", str(epochs)])
    if "--pretrain_max_iter" in accepted:
        cmd.extend(["--pretrain_max_iter", str(pretrain_epochs)])
    if "--maxiter" in accepted:
        cmd.extend(["--maxiter", str(epochs)])
    if gpu is not None and "--gpu" in accepted:
        cmd.extend(["--gpu", str(gpu)])
    if gpu is None and "--no_cuda" in accepted:
        cmd.append("--no_cuda")
    if task.method.extra_args:
        cmd.extend(task.method.extra_args)
    return cmd


def git_info() -> tuple[str, str]:
    try:
        sha = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=PROJECT_ROOT, text=True).strip()
        branch = subprocess.check_output(
            ["git", "rev-parse", "--abbrev-ref", "HEAD"], cwd=PROJECT_ROOT, text=True
        ).strip()
        return sha, branch
    except Exception:
        return "unknown", "unknown"


def normalize_metrics(metrics: dict[str, Any] | None) -> dict[str, Any]:
    if not metrics:
        return {}
    if "kmeans_known_k" in metrics and isinstance(metrics["kmeans_known_k"], dict):
        metrics = metrics["kmeans_known_k"]
    elif "fixed" in metrics and isinstance(metrics["fixed"], dict):
        fixed = metrics["fixed"]
        metrics = fixed.get("kmeans_known_k", fixed) if isinstance(fixed, dict) else metrics
    normalized: dict[str, Any] = {}
    aliases = {
        "ACC": "acc",
        "NMI": "nmi",
        "ARI": "ari",
        "F1_macro": "f1_macro",
        "FMI": "fmi",
        "V_measure": "v_measure",
        "Homogeneity": "homogeneity",
        "Completeness": "completeness",
    }
    for key, value in metrics.items():
        normalized[aliases.get(str(key), str(key).lower())] = value
    return normalized


def load_metrics(out_dir: Path) -> dict[str, Any]:
    path = out_dir / "metrics.json"
    if not path.exists():
        return {}
    try:
        return read_json(path)
    except Exception:
        return {}


def detect_fallback(method: MethodSpec, log_text: str) -> bool:
    if method.category != "Foundation":
        return False
    lowered = log_text.lower()
    return "fallback" in lowered or "pca baseline" in lowered or "using pca" in lowered


def write_environment_json(
    out_dir: Path,
    python_bin: str,
    runtime_env: str,
    gpu: int | None,
    method: MethodSpec,
) -> None:
    payload = {
        "python_executable": python_bin,
        "runtime_env": runtime_env,
        "framework": method.framework,
        "category": method.category,
        "gpu": gpu if gpu is not None else "cpu",
        "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", ""),
    }
    probe = subprocess.run(
        [python_bin, "-c", "import sys; print(sys.version.split()[0])"],
        text=True,
        capture_output=True,
        timeout=5,
    )
    if probe.returncode == 0:
        payload["python_version"] = probe.stdout.strip()
    else:
        payload["python_version_error"] = probe.stderr.strip()[-500:]
    write_json(out_dir / "environment.json", payload)


def task_out_dir(out_root: Path, task: Task) -> Path:
    return out_root / task.dataset.name / task.method.key / f"seed_{task.seed}"


def make_base_row(task: Task, out_dir: Path, commit_sha: str, branch: str) -> dict[str, Any]:
    return {
        "dataset": task.dataset.name,
        "model": task.method.key,
        "seed": task.seed,
        "status": "pending",
        "acc": "",
        "nmi": "",
        "ari": "",
        "f1_macro": "",
        "fmi": "",
        "v_measure": "",
        "homogeneity": "",
        "completeness": "",
        "runtime_seconds": "",
        "gpu": "",
        "command": "",
        "save_dir": str(out_dir),
        "error": "",
        "substitute_model_used": False,
        "fallback": False,
        "authenticity": task.method.authenticity,
        "category": task.method.category,
        "framework": task.method.framework,
        "runtime_env": task.method.runtime_env,
        "python_executable": "",
        "commit_sha": commit_sha,
        "branch": branch,
    }


def row_from_status(task: Task, out_dir: Path, commit_sha: str, branch: str) -> dict[str, Any]:
    row = make_base_row(task, out_dir, commit_sha, branch)
    status_path = out_dir / "status.json"
    if not status_path.exists():
        return row
    status = read_json_or_none(status_path)
    if status is None:
        return row
    row.update(
        {
            "status": status.get("status", "unknown"),
            "runtime_seconds": status.get("runtime_seconds", ""),
            "gpu": status.get("gpu", ""),
            "command": status.get("command", ""),
            "error": status.get("error", ""),
            "substitute_model_used": status.get("substitute_model_used", False),
            "fallback": status.get("fallback", False),
            "runtime_env": status.get("runtime_env", ""),
            "python_executable": status.get("python_executable", ""),
        }
    )
    metrics = normalize_metrics(load_metrics(out_dir))
    for metric in METRIC_NAMES:
        row[metric] = metrics.get(metric, "")
    return row


def read_existing_master_rows(out_root: Path) -> list[dict[str, Any]]:
    master_path = out_root / "all_runs_master.csv"
    if not master_path.exists():
        return []
    with master_path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def merge_master_rows(existing_rows: list[dict[str, Any]], update_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    def row_key(row: dict[str, Any]) -> tuple[str, str, str]:
        return (str(row.get("dataset", "")), str(row.get("model", "")), str(row.get("seed", "")))

    updates = {row_key(row): row for row in update_rows}
    merged: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in existing_rows:
        key = row_key(row)
        merged.append(updates.get(key, row))
        seen.add(key)
    for row in update_rows:
        key = row_key(row)
        if key not in seen:
            merged.append(row)
    return merged


def write_master_tables(
    out_root: Path,
    tasks: list[Task],
    commit_sha: str,
    branch: str,
    merge_existing_master: bool = False,
) -> None:
    rows = [row_from_status(task, task_out_dir(out_root, task), commit_sha, branch) for task in tasks]
    fieldnames = list(make_base_row(tasks[0], task_out_dir(out_root, tasks[0]), commit_sha, branch).keys()) if tasks else []
    if merge_existing_master:
        existing_rows = read_existing_master_rows(out_root)
        if existing_rows:
            rows = merge_master_rows(existing_rows, rows)
            for row in rows:
                for field in row.keys():
                    if field not in fieldnames:
                        fieldnames.append(field)
    with (out_root / "all_runs_master.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    write_json(out_root / "all_runs_master.json", rows)
    write_summary_tables(out_root, rows)


def write_summary_tables(out_root: Path, rows: list[dict[str, Any]]) -> None:
    import statistics
    from collections import defaultdict

    valid = [
        row
        for row in rows
        if row.get("status") == "success" and str(row.get("substitute_model_used")).lower() not in {"true", "1"}
    ]
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
                vals: list[float] = []
                for row in group:
                    value = row.get(metric, "")
                    if value != "":
                        vals.append(float(value))
                if vals:
                    out[f"{metric}_mean"] = statistics.mean(vals)
                    out[f"{metric}_std"] = statistics.stdev(vals) if len(vals) > 1 else 0.0
                else:
                    out[f"{metric}_mean"] = ""
                    out[f"{metric}_std"] = ""
            out_rows.append(out)
        fields = list(out_rows[0].keys()) if out_rows else list(group_fields) + ["n_success"]
        with (out_root / filename).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(out_rows)


def run_task(
    task: Task,
    out_root: Path,
    registry: dict[str, dict[str, Any]],
    gpu_queue: queue.Queue[int],
    epochs: int,
    pretrain_epochs: int,
    commit_sha: str,
    branch: str,
    resume: bool,
    rerun_failed: bool,
    dry_run: bool,
) -> dict[str, Any]:
    out_dir = task_out_dir(out_root, task)
    status_path = out_dir / "status.json"
    if resume and status_path.exists():
        previous = read_json(status_path)
        previous_status = str(previous.get("status", "unknown"))
        if previous_status in {"success", "fallback", "skipped"} or (previous_status not in {"failed", "timeout", "error", "incomplete"}):
            return row_from_status(task, out_dir, commit_sha, branch)
        if previous_status in {"failed", "timeout", "error", "incomplete"} and not rerun_failed:
            return row_from_status(task, out_dir, commit_sha, branch)

    out_dir.mkdir(parents=True, exist_ok=True)
    prepared_data = task.dataset
    if prepared_data.prepared_path is None or prepared_data.n_clusters is None:
        raise RuntimeError(f"Dataset was not prepared before scheduling: {task.dataset.name}")
    python_bin, runtime_env = resolve_python(task.method, registry)

    gpu: int | None = None
    if not is_cpu_method(task.method):
        gpu = gpu_queue.get()
    try:
        cmd = build_command(python_bin, task, prepared_data, out_dir, gpu, epochs, pretrain_epochs)
        write_json(
            out_dir / "args.json",
            {
                "dataset": task.dataset.name,
                "method": task.method.key,
                "seed": task.seed,
                "n_clusters": prepared_data.n_clusters,
                "data_path": str(prepared_data.prepared_path),
                "command": cmd,
            },
        )
        write_json(
            out_dir / "authenticity.json",
            {
                "method": task.method.key,
                "authenticity": task.method.authenticity,
                "source_path": task.method.source_path,
                "framework": task.method.framework,
                "category": task.method.category,
                "substitute_model_used": False,
            },
        )
        write_environment_json(out_dir, python_bin, runtime_env, gpu, task.method)
        (out_dir / "command.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")

        start = datetime.now().isoformat()
        if dry_run:
            status = {
                "dataset": task.dataset.name,
                "method": task.method.key,
                "seed": task.seed,
                "status": "dry_run",
                "runtime_seconds": 0.0,
                "gpu": gpu if gpu is not None else "cpu",
                "command": " ".join(cmd),
                "error": "",
                "fallback": False,
                "substitute_model_used": False,
                "runtime_env": runtime_env,
                "python_executable": python_bin,
                "commit_sha": commit_sha,
                "branch": branch,
                "start_time": start,
                "end_time": datetime.now().isoformat(),
            }
            write_json(status_path, status)
            return row_from_status(task, out_dir, commit_sha, branch)

        env = os.environ.copy()
        env["PYTHONPATH"] = str(METHODS_DIR) + ":" + str(PROJECT_ROOT) + ":" + env.get("PYTHONPATH", "")
        env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
        env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
        env.setdefault("OMP_NUM_THREADS", "4")
        env.setdefault("OPENBLAS_NUM_THREADS", "4")
        env.setdefault("MKL_NUM_THREADS", "4")
        env.setdefault("NUMEXPR_NUM_THREADS", "4")
        if runtime_env in {"plantnet-sccdcg-h100", "plantnet-tf212-cpu", "plantnet-scrcl"}:
            env["PYTHONNOUSERSITE"] = "1"
        if runtime_env == "plantnet-tf212-cpu":
            env["CUDA_VISIBLE_DEVICES"] = ""
            env.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
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
            fallback = detect_fallback(task.method, log_text)
            missing = [
                name
                for name in ["metrics.json", "embedding_final.npy", "labels.npy", "args.json"]
                if not (out_dir / name).exists()
            ]
            if proc.returncode != 0:
                state = "failed"
                error = f"Exit code {proc.returncode}"
            elif fallback:
                state = "fallback"
                error = "Foundation model used fallback path; excluded from success summaries."
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
            fallback = False

        substitute = bool(fallback)
        status = {
            "dataset": task.dataset.name,
            "method": task.method.key,
            "seed": task.seed,
            "status": state,
            "runtime_seconds": elapsed,
            "gpu": gpu if gpu is not None else "cpu",
            "command": " ".join(cmd),
            "error": error,
            "fallback": fallback,
            "substitute_model_used": substitute,
            "runtime_env": runtime_env,
            "python_executable": python_bin,
            "commit_sha": commit_sha,
            "branch": branch,
            "start_time": start,
            "end_time": datetime.now().isoformat(),
        }
        write_json(status_path, status)

        if substitute:
            auth = read_json(out_dir / "authenticity.json")
            auth["substitute_model_used"] = True
            auth["fallback"] = True
            write_json(out_dir / "authenticity.json", auth)
        return row_from_status(task, out_dir, commit_sha, branch)
    finally:
        if gpu is not None:
            gpu_queue.put(gpu)


def validate_gpus(gpus: list[int]) -> None:
    forbidden = sorted(set(gpus) & FORBIDDEN_GPUS)
    if forbidden:
        raise SystemExit(f"Forbidden GPU(s) requested: {forbidden}; use only 1-6.")


def build_tasks(datasets: list[DatasetSpec], methods: list[MethodSpec], seeds: list[int]) -> list[Task]:
    return [Task(dataset, method, seed) for dataset in datasets for method in methods for seed in seeds]


def write_catalogs(out_root: Path, datasets: list[DatasetSpec], methods: list[MethodSpec], tasks: list[Task]) -> None:
    write_json(out_root / "dataset_catalog.json", [dataset.__dict__ for dataset in datasets])
    write_json(
        out_root / "method_catalog.json",
        [
            {
                "key": method.key,
                "path": str(method.path.relative_to(PROJECT_ROOT)),
                "category": method.category,
                "framework": method.framework,
                "authenticity": method.authenticity,
                "runtime_env": method.runtime_env,
                "source_path": method.source_path,
            }
            for method in methods
        ],
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


def make_out_root(base: Path | None) -> Path:
    if base:
        return base
    stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return PROJECT_ROOT / "result" / f"scmae_all_methods_{stamp}"


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run all active /methods runners on all scMAE datasets.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data_root", type=Path, default=DATA_ROOT_DEFAULT)
    parser.add_argument("--out_root", type=Path, default=None)
    parser.add_argument("--seeds", type=int, nargs="+", default=DEFAULT_SEEDS)
    parser.add_argument("--gpus", type=int, nargs="+", default=DEFAULT_GPUS)
    parser.add_argument("--max_workers", type=int, default=6)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--pretrain_epochs", type=int, default=200)
    parser.add_argument("--datasets", nargs="+", default=None, help="Optional subset of discovered dataset names.")
    parser.add_argument("--methods", nargs="+", default=None, help="Optional subset of discovered method keys.")
    parser.add_argument("--runtime_registry", type=Path, default=None)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--resume", action="store_true", default=True)
    parser.add_argument("--no_resume", dest="resume", action="store_false")
    parser.add_argument("--rerun_failed", action="store_true")
    parser.add_argument(
        "--merge_existing_master",
        action="store_true",
        help="Merge subset repair rows into an existing all_runs_master table instead of replacing it.",
    )
    args = parser.parse_args()

    validate_gpus(args.gpus)
    out_root = make_out_root(args.out_root).resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    converted_dir = out_root / "converted_data"
    commit_sha, branch = git_info()

    datasets = discover_datasets(args.data_root)
    if args.datasets:
        wanted_datasets = set(args.datasets)
        datasets = [dataset for dataset in datasets if dataset.name in wanted_datasets]
    if not datasets:
        raise SystemExit(f"No datasets discovered under {args.data_root}")

    prepared_datasets: list[DatasetSpec] = []
    for dataset in datasets:
        prepared = planned_prepared_dataset(dataset, converted_dir) if args.dry_run else prepare_dataset(dataset, converted_dir)
        prepared_datasets.append(prepared)
    datasets = prepared_datasets

    selected_keys = set(args.methods) if args.methods else None
    methods = discover_methods(selected_keys)
    if not methods:
        raise SystemExit("No active methods discovered.")
    tasks = build_tasks(datasets, methods, args.seeds)
    registry = load_runtime_registry(args.runtime_registry)

    suite_config = {
        "data_root": str(args.data_root),
        "out_root": str(out_root),
        "seeds": args.seeds,
        "gpus": args.gpus,
        "max_workers": args.max_workers,
        "epochs": args.epochs,
        "pretrain_epochs": args.pretrain_epochs,
        "dry_run": args.dry_run,
        "rerun_failed": args.rerun_failed,
        "merge_existing_master": args.merge_existing_master,
        "commit_sha": commit_sha,
        "branch": branch,
    }
    if args.merge_existing_master:
        repair_dir = out_root / "repair_manifests"
        repair_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        write_json(
            repair_dir / f"repair_task_manifest_{stamp}.json",
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
        write_json(repair_dir / f"repair_suite_config_{stamp}.json", suite_config)
    else:
        write_catalogs(out_root, datasets, methods, tasks)
        write_json(out_root / "suite_config.json", suite_config)

    print(f"Suite output: {out_root}")
    print(f"Datasets: {len(datasets)}")
    print(f"Methods:  {len(methods)}")
    print(f"Seeds:    {args.seeds}")
    print(f"Tasks:    {len(tasks)}")
    print(f"Dry run:  {args.dry_run}")

    # Prepopulate the master table so dry-runs and interrupted runs still expose
    # the complete dataset x model x seed grid.
    write_master_tables(out_root, tasks, commit_sha, branch, args.merge_existing_master)

    gpu_queue: queue.Queue[int] = queue.Queue()
    for gpu in args.gpus:
        gpu_queue.put(gpu)

    write_lock = threading.Lock()
    worker_count = max(1, min(int(args.max_workers), len(args.gpus)))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = [
            executor.submit(
                run_task,
                task,
                out_root,
                registry,
                gpu_queue,
                args.epochs,
                args.pretrain_epochs,
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
                print(
                    f"[{index}/{len(futures)}] {row['dataset']} / {row['model']} / "
                    f"seed {row['seed']}: {row['status']}"
                )
            except Exception as exc:
                print(f"[{index}/{len(futures)}] task wrapper error: {exc}", file=sys.stderr)
            if index == len(futures) or index % 25 == 0:
                with write_lock:
                    write_master_tables(out_root, tasks, commit_sha, branch, args.merge_existing_master)

    write_master_tables(out_root, tasks, commit_sha, branch, args.merge_existing_master)
    print(f"Master table: {out_root / 'all_runs_master.csv'}")
    print(f"Summary by model/dataset: {out_root / 'summary_by_model_dataset.csv'}")
    print(f"Summary by model: {out_root / 'summary_by_model.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
