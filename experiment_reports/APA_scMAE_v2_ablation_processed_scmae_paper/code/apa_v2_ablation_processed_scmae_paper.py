#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import os
import platform
import queue
import socket
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from statistics import mean, stdev
from typing import Any

import anndata as ad


REPO = Path("/home/luolie/biopipeline/dimension-reduction/plantnet")
DATA_ROOT = Path("/data/luolie/biopipeline/scCluBench/data/processed_scmae")
OUTPUT_ROOT = REPO / "benchmark_outputs" / "APA_scMAE_v2_ablation_processed_scmae_paper"
RUNNER = REPO / "methods" / "DeepLearning" / "APA_scMAE" / "run.py"
SEEDS = [42, 2024, 3407]
GPUS = [1, 2, 3, 4, 6]
PARALLEL_JOBS = 5
CLUSTER_METHODS = ["kmeans_known_k", "leiden_fixed"]
METRIC_NAMES = ["acc", "nmi", "ari", "f1_macro", "fmi", "v_measure", "homogeneity", "completeness"]


DATASETS = {
    "Limb_Muscle": {"file": "Limb_Muscle.h5ad", "n_clusters": 6},
    "Quake_10x_Spleen": {"file": "Quake_10x_Spleen.h5ad", "n_clusters": 5},
    "Guo": {"file": "Guo.h5ad", "n_clusters": 12},
    "Macosko": {"file": "Macosko.h5ad", "n_clusters": 12},
    "Tosches": {"file": "Tosches.h5ad", "n_clusters": 15},
    "Young": {"file": "Young.h5ad", "n_clusters": 11},
}


VARIANTS = {
    "apa_v2_A_random_student": {
        "label": "A_random_student",
        "flags": {
            "use_repr_loss": "false",
            "use_ema_teacher": "false",
            "use_proto_consistency": "false",
            "student_warmup_epochs": "20",
            "enable_generator_after_warmup": "false",
        },
    },
    "apa_v2_B_vicreg": {
        "label": "B_vicreg",
        "flags": {
            "use_repr_loss": "true",
            "use_ema_teacher": "false",
            "use_proto_consistency": "false",
            "student_warmup_epochs": "20",
            "enable_generator_after_warmup": "false",
        },
    },
    "apa_v2_C_ema_teacher": {
        "label": "C_ema_teacher",
        "flags": {
            "use_repr_loss": "true",
            "use_ema_teacher": "true",
            "use_proto_consistency": "false",
            "student_warmup_epochs": "20",
            "enable_generator_after_warmup": "false",
        },
    },
    "apa_v2_D_proto": {
        "label": "D_proto",
        "flags": {
            "use_repr_loss": "true",
            "use_ema_teacher": "true",
            "use_proto_consistency": "true",
            "student_warmup_epochs": "20",
            "enable_generator_after_warmup": "false",
        },
    },
    "apa_v2_E_full_generator": {
        "label": "E_full_generator",
        "flags": {
            "use_repr_loss": "true",
            "use_ema_teacher": "true",
            "use_proto_consistency": "true",
            "student_warmup_epochs": "20",
            "enable_generator_after_warmup": "true",
        },
    },
}


@dataclass
class Task:
    dataset: str
    data_path: Path
    n_clusters: int
    variant: str
    method_name: str
    seed: int
    run_dir: Path


def now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
    tmp.replace(path)


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def run_checked(cmd: list[str], cwd: Path = REPO) -> None:
    print("$ " + " ".join(cmd), flush=True)
    subprocess.run(cmd, cwd=str(cwd), check=True)


def git_text(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], cwd=str(REPO), text=True).strip()


def bool_flag_args(flags: dict[str, str]) -> list[str]:
    args: list[str] = []
    for key, value in flags.items():
        args.extend([f"--{key}", value])
    return args


def build_command(task: Task, gpu: int, batch_size: int, epochs: int) -> list[str]:
    flags = VARIANTS[task.variant]["flags"]
    cmd = [
        sys.executable,
        str(RUNNER),
        "--data_path",
        str(task.data_path),
        "--save_dir",
        str(task.run_dir),
        "--dataset_name",
        task.dataset,
        "--label_key",
        "resolved_label",
        "--input_mode",
        "log1p",
        "--n_top_genes",
        "2000",
        "--epochs",
        str(epochs),
        "--batch_size",
        str(batch_size),
        "--decoder_mode",
        "z_with_stopgrad_h",
        "--n_clusters",
        str(task.n_clusters),
        "--seed",
        str(task.seed),
        "--gpu",
        str(gpu),
        "--method_name",
        task.method_name,
        "--skip_eval",
        "false",
    ]
    cmd.extend(bool_flag_args(flags))
    return cmd


def shell_join(cmd: list[str]) -> str:
    import shlex

    return " ".join(shlex.quote(part) for part in cmd)


def is_complete(run_dir: Path) -> bool:
    manifest = run_dir / "artifact_manifest.json"
    metrics = run_dir / "metrics.json"
    embedding = run_dir / "embedding_final.npy"
    if not (manifest.exists() and metrics.exists() and embedding.exists()):
        return False
    try:
        return read_json(manifest).get("status") == "complete"
    except Exception:
        return False


def preflight_datasets() -> None:
    print("Dataset preflight", flush=True)
    if not DATA_ROOT.exists():
        raise FileNotFoundError(f"DATA_ROOT does not exist: {DATA_ROOT}")
    for name, spec in DATASETS.items():
        path = DATA_ROOT / spec["file"]
        if not path.exists():
            raise FileNotFoundError(f"Missing dataset: {path}")
        adata = ad.read_h5ad(path, backed="r")
        try:
            if "resolved_label" not in adata.obs:
                raise ValueError(f"{name}: adata.obs['resolved_label'] is missing")
            n_unique = int(adata.obs["resolved_label"].nunique())
            planned = int(spec["n_clusters"])
            if n_unique != planned:
                raise ValueError(f"{name}: resolved_label unique labels {n_unique} != planned n_clusters {planned}")
            print(
                f"{name}: n_cells={adata.n_obs} n_genes={adata.n_vars} "
                f"label_key=resolved_label n_unique_labels={n_unique} planned_n_clusters={planned}",
                flush=True,
            )
        finally:
            adata.file.close()


def make_tasks() -> list[Task]:
    tasks: list[Task] = []
    for dataset, spec in DATASETS.items():
        data_path = DATA_ROOT / spec["file"]
        for variant in VARIANTS:
            for seed in SEEDS:
                run_dir = OUTPUT_ROOT / dataset / f"{variant}__seed{seed}"
                tasks.append(
                    Task(
                        dataset=dataset,
                        data_path=data_path,
                        n_clusters=int(spec["n_clusters"]),
                        variant=variant,
                        method_name=variant,
                        seed=seed,
                        run_dir=run_dir,
                    )
                )
    return tasks


def initial_manifest(tasks: list[Task], branch: str, commit_sha: str, batch_size: int) -> list[dict[str, Any]]:
    rows = []
    for task in tasks:
        rows.append(
            {
                "dataset": task.dataset,
                "data_path": str(task.data_path),
                "data_path_resolved": str(task.data_path.resolve()),
                "n_clusters": task.n_clusters,
                "variant": task.variant,
                "method_name": task.method_name,
                "seed": task.seed,
                "gpu": None,
                "command": "",
                "run_dir": str(task.run_dir),
                "status": "pending",
                "start_time": None,
                "end_time": None,
                "elapsed_seconds": None,
                "batch_size_resolved": batch_size,
                "commit_sha": commit_sha,
                "branch": branch,
            }
        )
    return rows


def metric_value(metrics: dict[str, Any], key: str) -> float | None:
    value = metrics.get(key)
    if value is None:
        return None
    try:
        value_f = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(value_f):
        return None
    return value_f


def n_pred_clusters(run_dir: Path, cluster_method: str) -> int | None:
    import numpy as np

    path = run_dir / ("pred_labels.npy" if cluster_method == "kmeans_known_k" else "eval_leiden_fixed.npy")
    if not path.exists():
        return None
    try:
        pred = np.load(path)
        return int(len(set(pred.tolist())))
    except Exception:
        return None


def rows_for_task(manifest_row: dict[str, Any]) -> list[dict[str, Any]]:
    run_dir = Path(manifest_row["run_dir"])
    metrics_path = run_dir / "metrics.json"
    base = {
        "dataset": manifest_row["dataset"],
        "variant": manifest_row["variant"],
        "method_name": manifest_row["method_name"],
        "seed": manifest_row["seed"],
        "run_dir": str(run_dir),
        "metrics_path": str(metrics_path),
    }
    metrics_doc: dict[str, Any] = {}
    if metrics_path.exists():
        try:
            metrics_doc = read_json(metrics_path)
        except Exception:
            metrics_doc = {}
    rows = []
    for cluster_method in CLUSTER_METHODS:
        cluster = metrics_doc.get(cluster_method, {}) if isinstance(metrics_doc, dict) else {}
        status = "success"
        skip_reason = ""
        if manifest_row["status"] != "complete":
            status = manifest_row["status"]
            skip_reason = manifest_row.get("error_tail") or ""
        elif cluster_method not in metrics_doc:
            status = "missing"
            skip_reason = "cluster_method_missing_from_metrics_json"
        elif cluster.get("status") == "skipped":
            status = "skipped"
            skip_reason = str(cluster.get("skip_reason") or cluster.get("reason") or "")
        row = {
            **base,
            "cluster_method": cluster_method,
            "status": status,
            "skip_reason": skip_reason,
            "n_pred_clusters": n_pred_clusters(run_dir, cluster_method),
        }
        for name in METRIC_NAMES:
            row[name] = metric_value(cluster, name) if status == "success" else None
        rows.append(row)
    return rows


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(path)


def fmt_mean_std(avg: float | None, sd: float | None) -> str:
    if avg is None or not math.isfinite(avg):
        return ""
    sd_value = 0.0 if sd is None or not math.isfinite(sd) else sd
    return f"{avg:.3f} \u00b1 {sd_value:.3f}"


def summarize(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    summary_rows: list[dict[str, Any]] = []
    index: dict[tuple[str, str, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        key = (row["dataset"], row["variant"], row["method_name"], row["cluster_method"])
        index.setdefault(key, []).append(row)
    for dataset in DATASETS:
        for variant in VARIANTS:
            for cluster_method in CLUSTER_METHODS:
                key = (dataset, variant, variant, cluster_method)
                group = index.get(key, [])
                success = [r for r in group if r.get("status") == "success"]
                out = {
                    "dataset": dataset,
                    "variant": variant,
                    "method_name": variant,
                    "cluster_method": cluster_method,
                    "count_success": len(success),
                }
                for metric in METRIC_NAMES:
                    vals = [float(r[metric]) for r in success if r.get(metric) is not None]
                    out[f"{metric}_mean"] = mean(vals) if vals else None
                    out[f"{metric}_std"] = stdev(vals) if len(vals) >= 2 else 0.0 if len(vals) == 1 else None
                summary_rows.append(out)

    paper_rows: list[dict[str, Any]] = []
    paper_ari_rows: list[dict[str, Any]] = []
    paper_nmi_rows: list[dict[str, Any]] = []
    for dataset in DATASETS:
        for cluster_method in CLUSTER_METHODS:
            row = {"dataset": dataset, "cluster_method": cluster_method}
            ari_row = {"dataset": dataset, "cluster_method": cluster_method}
            nmi_row = {"dataset": dataset, "cluster_method": cluster_method}
            for variant, spec in VARIANTS.items():
                label = spec["label"]
                summary = next(
                    r for r in summary_rows if r["dataset"] == dataset and r["variant"] == variant and r["cluster_method"] == cluster_method
                )
                ari = fmt_mean_std(summary["ari_mean"], summary["ari_std"])
                nmi = fmt_mean_std(summary["nmi_mean"], summary["nmi_std"])
                row[f"{label}_ARI_mean_std"] = ari
                row[f"{label}_NMI_mean_std"] = nmi
                ari_row[f"{label}_ARI_mean_std"] = ari
                nmi_row[f"{label}_NMI_mean_std"] = nmi
            paper_rows.append(row)
            paper_ari_rows.append(ari_row)
            paper_nmi_rows.append(nmi_row)

    paired_rows: list[dict[str, Any]] = []
    for dataset in DATASETS:
        for cluster_method in CLUSTER_METHODS:
            e_by_seed = {
                int(r["seed"]): r
                for r in rows
                if r["dataset"] == dataset
                and r["variant"] == "apa_v2_E_full_generator"
                and r["cluster_method"] == cluster_method
                and r["status"] == "success"
            }
            d_by_seed = {
                int(r["seed"]): r
                for r in rows
                if r["dataset"] == dataset
                and r["variant"] == "apa_v2_D_proto"
                and r["cluster_method"] == cluster_method
                and r["status"] == "success"
            }
            ari_diffs = []
            nmi_diffs = []
            for seed in SEEDS:
                if seed in e_by_seed and seed in d_by_seed:
                    if e_by_seed[seed].get("ari") is not None and d_by_seed[seed].get("ari") is not None:
                        ari_diffs.append(float(e_by_seed[seed]["ari"]) - float(d_by_seed[seed]["ari"]))
                    if e_by_seed[seed].get("nmi") is not None and d_by_seed[seed].get("nmi") is not None:
                        nmi_diffs.append(float(e_by_seed[seed]["nmi"]) - float(d_by_seed[seed]["nmi"]))
            paired_rows.append(
                {
                    "dataset": dataset,
                    "cluster_method": cluster_method,
                    "E_minus_D_ari_mean": mean(ari_diffs) if ari_diffs else None,
                    "E_minus_D_ari_std": stdev(ari_diffs) if len(ari_diffs) >= 2 else 0.0 if len(ari_diffs) == 1 else None,
                    "E_better_than_D_seed_count": sum(1 for x in ari_diffs if x > 0),
                    "E_minus_D_nmi_mean": mean(nmi_diffs) if nmi_diffs else None,
                    "E_minus_D_nmi_std": stdev(nmi_diffs) if len(nmi_diffs) >= 2 else 0.0 if len(nmi_diffs) == 1 else None,
                    "E_better_than_D_seed_count_nmi": sum(1 for x in nmi_diffs if x > 0),
                }
            )
    return summary_rows, paper_rows, paper_ari_rows, paper_nmi_rows, paired_rows


def refresh_outputs(manifest_rows: list[dict[str, Any]], failures: list[dict[str, Any]]) -> None:
    write_json(OUTPUT_ROOT / "task_manifest.json", manifest_rows)
    by_run_rows = [row for manifest_row in manifest_rows for row in rows_for_task(manifest_row)]
    write_csv(
        OUTPUT_ROOT / "apa_v2_ablation_metrics_by_run.csv",
        by_run_rows,
        [
            "dataset",
            "variant",
            "method_name",
            "seed",
            "cluster_method",
            "status",
            *METRIC_NAMES,
            "n_pred_clusters",
            "skip_reason",
            "run_dir",
            "metrics_path",
        ],
    )
    summary_rows, paper_rows, paper_ari_rows, paper_nmi_rows, paired_rows = summarize(by_run_rows)
    summary_fields = ["dataset", "variant", "method_name", "cluster_method", "count_success"]
    for metric in METRIC_NAMES:
        summary_fields.extend([f"{metric}_mean", f"{metric}_std"])
    write_csv(OUTPUT_ROOT / "apa_v2_ablation_metrics_summary.csv", summary_rows, summary_fields)
    paper_fields = ["dataset", "cluster_method"]
    for spec in VARIANTS.values():
        paper_fields.extend([f"{spec['label']}_ARI_mean_std", f"{spec['label']}_NMI_mean_std"])
    write_csv(OUTPUT_ROOT / "apa_v2_ablation_paper_table.csv", paper_rows, paper_fields)
    ari_fields = ["dataset", "cluster_method", *[f"{spec['label']}_ARI_mean_std" for spec in VARIANTS.values()]]
    nmi_fields = ["dataset", "cluster_method", *[f"{spec['label']}_NMI_mean_std" for spec in VARIANTS.values()]]
    write_csv(OUTPUT_ROOT / "apa_v2_ablation_paper_table_ari.csv", paper_ari_rows, ari_fields)
    write_csv(OUTPUT_ROOT / "apa_v2_ablation_paper_table_nmi.csv", paper_nmi_rows, nmi_fields)
    write_csv(
        OUTPUT_ROOT / "apa_v2_ablation_E_vs_D_paired.csv",
        paired_rows,
        [
            "dataset",
            "cluster_method",
            "E_minus_D_ari_mean",
            "E_minus_D_ari_std",
            "E_better_than_D_seed_count",
            "E_minus_D_nmi_mean",
            "E_minus_D_nmi_std",
            "E_better_than_D_seed_count_nmi",
        ],
    )
    write_json(OUTPUT_ROOT / "failures.json", failures)


def run_one(task: Task, manifest_row: dict[str, Any], gpu: int, batch_size: int, epochs: int, branch: str, commit_sha: str) -> tuple[bool, dict[str, Any] | None]:
    if epochs == 80 and is_complete(task.run_dir):
        manifest_row.update(
            {
                "gpu": gpu,
                "command": shell_join(build_command(task, gpu, batch_size, epochs)),
                "status": "complete",
                "start_time": manifest_row.get("start_time") or now_iso(),
                "end_time": manifest_row.get("end_time") or now_iso(),
                "elapsed_seconds": 0.0,
                "batch_size_resolved": batch_size,
                "commit_sha": commit_sha,
                "branch": branch,
                "skip_reason": "existing_complete",
            }
        )
        return True, None

    task.run_dir.mkdir(parents=True, exist_ok=True)
    cmd = build_command(task, gpu, batch_size, epochs)
    log_path = task.run_dir / "run.log"
    env = os.environ.copy()
    env.update(
        {
            "CUDA_VISIBLE_DEVICES": str(gpu),
            "MPLCONFIGDIR": "/tmp/matplotlib",
            "NUMBA_CACHE_DIR": "/tmp/numba-cache",
            "OMP_NUM_THREADS": "1",
            "MKL_NUM_THREADS": "1",
            "OPENBLAS_NUM_THREADS": "1",
            "NUMBA_NUM_THREADS": "1",
        }
    )
    start = time.time()
    manifest_row.update(
        {
            "gpu": gpu,
            "command": shell_join(cmd),
            "status": "running",
            "start_time": now_iso(),
            "end_time": None,
            "elapsed_seconds": None,
            "batch_size_resolved": batch_size,
            "commit_sha": commit_sha,
            "branch": branch,
        }
    )
    with log_path.open("w", encoding="utf-8") as log:
        log.write("$ " + shell_join(cmd) + "\n")
        log.flush()
        proc = subprocess.run(cmd, cwd=str(REPO), env=env, stdout=log, stderr=subprocess.STDOUT)
    elapsed = time.time() - start
    manifest_row["end_time"] = now_iso()
    manifest_row["elapsed_seconds"] = round(elapsed, 3)
    if proc.returncode == 0 and is_complete(task.run_dir):
        manifest_row["status"] = "complete"
        return True, None
    tail = ""
    try:
        tail = "\n".join(log_path.read_text(encoding="utf-8", errors="replace").splitlines()[-80:])
    except Exception:
        tail = ""
    manifest_row["status"] = "failed"
    manifest_row["error_tail"] = tail
    failure = {
        "dataset": task.dataset,
        "variant": task.variant,
        "seed": task.seed,
        "gpu": gpu,
        "run_dir": str(task.run_dir),
        "command": shell_join(cmd),
        "return_code": proc.returncode,
        "error_tail": tail,
        "log_path": str(log_path),
    }
    return False, failure


def run_vram_preflights(batch_size: int, branch: str, commit_sha: str) -> int:
    preflight_specs = [("Macosko", 42), ("Tosches", 42)]
    for dataset, seed in preflight_specs:
        spec = DATASETS[dataset]
        task = Task(
            dataset=dataset,
            data_path=DATA_ROOT / spec["file"],
            n_clusters=int(spec["n_clusters"]),
            variant="apa_v2_E_full_generator",
            method_name="apa_v2_E_full_generator",
            seed=seed,
            run_dir=OUTPUT_ROOT / "_vram_preflight" / dataset / "apa_v2_E_full_generator__seed42__epochs1",
        )
        row = initial_manifest([task], branch, commit_sha, batch_size)[0]
        ok, failure = run_one(task, row, GPUS[0], batch_size, 1, branch, commit_sha)
        if not ok:
            text = (failure or {}).get("error_tail", "").lower()
            if "out of memory" in text or "cuda error: out of memory" in text or "oom" in text:
                print(f"VRAM preflight OOM at batch_size={batch_size}; formal runs will use batch_size=8.", flush=True)
                return 8
            raise RuntimeError(f"VRAM preflight failed for {dataset} without recognizable OOM. See {task.run_dir / 'run.log'}")
    print(f"VRAM preflight passed at batch_size={batch_size}; formal runs will use batch_size={batch_size}.", flush=True)
    return batch_size


def run_formal(tasks: list[Task], manifest_rows: list[dict[str, Any]], batch_size: int, branch: str, commit_sha: str) -> None:
    failures: list[dict[str, Any]] = []
    lock = threading.Lock()
    gpu_queue: queue.Queue[int] = queue.Queue()
    for gpu in GPUS:
        gpu_queue.put(gpu)

    refresh_outputs(manifest_rows, failures)
    task_to_row = {(r["dataset"], r["variant"], int(r["seed"])): r for r in manifest_rows}

    def worker(task: Task) -> None:
        gpu = gpu_queue.get()
        key = (task.dataset, task.variant, task.seed)
        row = task_to_row[key]
        try:
            with lock:
                row["gpu"] = gpu
                row["command"] = shell_join(build_command(task, gpu, batch_size, 80))
                row["status"] = "running"
                row["start_time"] = now_iso()
                row["batch_size_resolved"] = batch_size
                refresh_outputs(manifest_rows, failures)
            ok, failure = run_one(task, row, gpu, batch_size, 80, branch, commit_sha)
            with lock:
                if not ok and failure is not None:
                    failures.append(failure)
                refresh_outputs(manifest_rows, failures)
                print(f"{row['status']}: {task.dataset} {task.variant} seed={task.seed} gpu={gpu}", flush=True)
        finally:
            gpu_queue.put(gpu)

    with ThreadPoolExecutor(max_workers=PARALLEL_JOBS) as pool:
        futures = [pool.submit(worker, task) for task in tasks]
        for future in as_completed(futures):
            future.result()

    refresh_outputs(manifest_rows, failures)


def final_checks(manifest_rows: list[dict[str, Any]]) -> None:
    failures = read_json(OUTPUT_ROOT / "failures.json")
    by_run_count = sum(1 for _ in csv.DictReader((OUTPUT_ROOT / "apa_v2_ablation_metrics_by_run.csv").open(encoding="utf-8")))
    summary_count = sum(1 for _ in csv.DictReader((OUTPUT_ROOT / "apa_v2_ablation_metrics_summary.csv").open(encoding="utf-8")))
    complete = sum(1 for row in manifest_rows if row.get("status") == "complete")
    if len(manifest_rows) != 90:
        raise RuntimeError(f"Expected 90 tasks, found {len(manifest_rows)}")
    if complete != 90:
        raise RuntimeError(f"Expected 90 complete runs, found {complete}")
    if failures != []:
        raise RuntimeError(f"Expected no failures, found {len(failures)}")
    if by_run_count != 180:
        raise RuntimeError(f"Expected 180 metrics_by_run rows, found {by_run_count}")
    if summary_count != 60:
        raise RuntimeError(f"Expected 60 summary groups, found {summary_count}")
    for row in manifest_rows:
        run_dir = Path(row["run_dir"])
        if not (run_dir / "metrics.json").exists() or not (run_dir / "artifact_manifest.json").exists():
            raise RuntimeError(f"Missing metrics or artifact manifest in {run_dir}")


def print_final_report() -> None:
    summary = list(csv.DictReader((OUTPUT_ROOT / "apa_v2_ablation_metrics_summary.csv").open(encoding="utf-8")))
    paired = list(csv.DictReader((OUTPUT_ROOT / "apa_v2_ablation_E_vs_D_paired.csv").open(encoding="utf-8")))
    failures = read_json(OUTPUT_ROOT / "failures.json")
    tasks = read_json(OUTPUT_ROOT / "task_manifest.json")
    print("\nFINAL COUNTS", flush=True)
    print(f"success={sum(1 for row in tasks if row.get('status') == 'complete')} failed={len(failures)}", flush=True)
    for cluster_method in CLUSTER_METHODS:
        print(f"\n{cluster_method} ARI mean/std", flush=True)
        for dataset in DATASETS:
            parts = []
            for variant, spec in VARIANTS.items():
                row = next(r for r in summary if r["dataset"] == dataset and r["variant"] == variant and r["cluster_method"] == cluster_method)
                parts.append(f"{spec['label']}={fmt_mean_std(float(row['ari_mean']) if row['ari_mean'] else None, float(row['ari_std']) if row['ari_std'] else None)}")
            print(f"{dataset}: " + "; ".join(parts), flush=True)
    print("\nPaired E vs D judgment", flush=True)
    for row in paired:
        ari_mean = row["E_minus_D_ari_mean"]
        nmi_mean = row["E_minus_D_nmi_mean"]
        if ari_mean:
            verdict = "E>D" if float(ari_mean) > 0 and int(row["E_better_than_D_seed_count"]) >= 2 else "E not consistently > D"
        else:
            verdict = "insufficient paired data"
        print(
            f"{row['dataset']} {row['cluster_method']}: "
            f"ARI diff={ari_mean} seeds_better={row['E_better_than_D_seed_count']}; "
            f"NMI diff={nmi_mean} seeds_better={row['E_better_than_D_seed_count_nmi']}; {verdict}",
            flush=True,
        )
    print("\nCSV outputs", flush=True)
    for name in [
        "apa_v2_ablation_metrics_by_run.csv",
        "apa_v2_ablation_metrics_summary.csv",
        "apa_v2_ablation_paper_table.csv",
        "apa_v2_ablation_paper_table_ari.csv",
        "apa_v2_ablation_paper_table_nmi.csv",
        "apa_v2_ablation_E_vs_D_paired.csv",
    ]:
        print(str(OUTPUT_ROOT / name), flush=True)


def main() -> int:
    os.chdir(REPO)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    branch = git_text(["rev-parse", "--abbrev-ref", "HEAD"])
    commit_sha = git_text(["rev-parse", "HEAD"])
    if branch != "apa-v2-objective":
        raise RuntimeError(f"Expected branch apa-v2-objective, got {branch}")
    preflight_datasets()
    write_json(
        OUTPUT_ROOT / "run_environment.json",
        {
            "branch": branch,
            "commit_sha": commit_sha,
            "hostname": socket.gethostname(),
            "python": sys.version,
            "platform": platform.platform(),
            "data_root": str(DATA_ROOT),
            "data_root_resolved": str(DATA_ROOT.resolve()),
            "seeds": SEEDS,
            "datasets": DATASETS,
            "variants": VARIANTS,
            "start_time": now_iso(),
        },
    )
    run_checked([sys.executable, "-m", "compileall", "-q", "methods/DeepLearning/APA_scMAE"])
    run_checked([sys.executable, "-m", "pytest", "-q", "methods/DeepLearning/APA_scMAE/tests"])
    batch_size = run_vram_preflights(16, branch, commit_sha)
    tasks = make_tasks()
    manifest_rows = initial_manifest(tasks, branch, commit_sha, batch_size)
    run_formal(tasks, manifest_rows, batch_size, branch, commit_sha)
    final_manifest = read_json(OUTPUT_ROOT / "task_manifest.json")
    final_checks(final_manifest)
    print_final_report()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
