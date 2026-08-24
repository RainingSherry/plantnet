#!/usr/bin/env python3
"""Run updated DESC/scNAME benchmark after compatibility fixes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PYTHON = Path("/data/luolie/conda/envs/scssl_bench_py310/bin/python")
DEFAULT_DATA_ROOT = Path("/data/luolie/biopipeline/scCluBench/data")
DEFAULT_REPORT_DIR = ROOT / "experiment_reports" / "desc_scname_fix_benchmark_20260628"
DEFAULT_RESULTS_DIR = DEFAULT_REPORT_DIR / "run_artifacts"
METHODS = {
    "desc": ROOT / "methods" / "DeepLearning" / "desc" / "run.py",
    "scname": ROOT / "methods" / "DeepLearning" / "scNAME" / "run.py",
}
METRICS = ["acc", "nmi", "ari", "f1_macro", "fmi", "v_measure", "homogeneity", "completeness"]
VALUE_TAG = "updated_after_desc_scname_fix_20260628"

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
Path(os.environ["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
Path(os.environ["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)


def sanitize_name(path: Path, data_root: Path) -> str:
    rel = path.relative_to(data_root).with_suffix("")
    parts = list(rel.parts)
    if len(parts) == 1:
        name = parts[0]
    else:
        prefix = {
            "processed": "processed",
            "processed_scmae": "",
            "processed_benchmark": "processed_benchmark",
            "smoke": "smoke",
            "其他": "other",
        }.get(parts[0], parts[0])
        name = "__".join([p for p in [prefix, *parts[1:]] if p])
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", name)


def list_datasets(data_root: Path) -> list[tuple[str, Path]]:
    paths = sorted(data_root.rglob("*.h5ad"))
    return [(sanitize_name(path, data_root), path) for path in paths]


def infer_n_clusters(data_path: Path) -> int:
    import scanpy as sc

    adata = sc.read_h5ad(data_path)
    for key in ["resolved_label", "cell_type", "Celltype", "celltype", "cell_label", "label"]:
        if key in adata.obs.columns:
            return int(adata.obs[key].nunique())
    raise ValueError(f"No supported label column found in {data_path}; obs={list(adata.obs.columns)}")


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)


def command_for(method: str, data_path: Path, out_dir: Path, n_clusters: int, seed: int, epochs: int, gpu: int) -> list[str]:
    return [
        str(PYTHON),
        str(METHODS[method]),
        "--data_path",
        str(data_path),
        "--save_dir",
        str(out_dir),
        "--n_clusters",
        str(n_clusters),
        "--seed",
        str(seed),
        "--epochs",
        str(epochs),
        "--pretrain_epochs",
        str(epochs),
        "--gpu",
        str(gpu),
    ]


def run_one(task: dict) -> dict:
    out_dir = task["out_dir"]
    status_path = out_dir / "status.json"
    prior = read_json(status_path)
    if prior.get("status") == "success" and (out_dir / "metrics.json").exists():
        return prior | {"skipped_existing": True}

    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = command_for(
        task["method"],
        task["data_path"],
        out_dir,
        task["n_clusters"],
        task["seed"],
        task["epochs"],
        task["gpu"],
    )
    start = datetime.now().isoformat()
    t0 = time.time()
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "methods") + ":" + env.get("PYTHONPATH", "")
    env["CUDA_VISIBLE_DEVICES"] = str(task["gpu"])
    env["OMP_NUM_THREADS"] = "4"
    env["OPENBLAS_NUM_THREADS"] = "4"
    env["MKL_NUM_THREADS"] = "4"
    env["NUMEXPR_NUM_THREADS"] = "4"
    env["NUMBA_NUM_THREADS"] = "4"
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    Path(env["MPLCONFIGDIR"]).mkdir(parents=True, exist_ok=True)
    Path(env["NUMBA_CACHE_DIR"]).mkdir(parents=True, exist_ok=True)

    status = "unknown"
    error = ""
    return_code = None
    try:
        proc = subprocess.run(
            cmd,
            cwd=str(ROOT),
            capture_output=True,
            text=True,
            timeout=task["timeout_seconds"],
            env=env,
        )
        return_code = proc.returncode
        (out_dir / "run.log").write_text(
            "Command: " + " ".join(cmd) + "\n"
            + f"Exit code: {proc.returncode}\n\n=== STDOUT ===\n{proc.stdout}\n"
            + (f"\n=== STDERR ===\n{proc.stderr}\n" if proc.stderr else ""),
            encoding="utf-8",
        )
        if proc.returncode != 0:
            status = "failed"
            error = f"Exit code {proc.returncode}"
        elif not (out_dir / "metrics.json").exists():
            status = "incomplete"
            error = "metrics.json missing"
        else:
            status = "success"
    except subprocess.TimeoutExpired as exc:
        status = "timeout"
        return_code = -1
        error = f"timeout after {task['timeout_seconds']}s"
        (out_dir / "run.log").write_text(
            "Command: " + " ".join(cmd) + "\n" + str(exc),
            encoding="utf-8",
        )
    except Exception as exc:
        status = "error"
        return_code = -2
        error = str(exc)
        (out_dir / "run.log").write_text(
            "Command: " + " ".join(cmd) + "\n" + repr(exc),
            encoding="utf-8",
        )

    result = {
        "value_tag": VALUE_TAG,
        "dataset": task["dataset"],
        "data_path": str(task["data_path"]),
        "method": task["method"],
        "seed": task["seed"],
        "epochs": task["epochs"],
        "pretrain_epochs": task["epochs"],
        "n_clusters": task["n_clusters"],
        "status": status,
        "return_code": return_code,
        "error": error,
        "gpu": task["gpu"],
        "no_cuda": False,
        "start_time": start,
        "end_time": datetime.now().isoformat(),
        "elapsed_seconds": round(time.time() - t0, 1),
        "command": " ".join(cmd),
        "save_dir": str(out_dir),
    }
    write_json(status_path, result)
    return result


def collect_rows(results_dir: Path) -> list[dict]:
    rows = []
    for status_path in sorted(results_dir.rglob("status.json")):
        row = read_json(status_path)
        metrics = read_json(Path(row.get("save_dir", status_path.parent)) / "metrics.json")
        for metric in METRICS:
            row[metric] = metrics.get(metric, "")
        rows.append(row)
    return rows


def write_tables(results_dir: Path, report_dir: Path) -> None:
    rows = collect_rows(results_dir)
    report_dir.mkdir(parents=True, exist_ok=True)
    fields = [
        "value_tag", "dataset", "data_path", "method", "seed", "epochs", "pretrain_epochs",
        "n_clusters", "status", *METRICS, "elapsed_seconds", "gpu", "no_cuda",
        "return_code", "error", "save_dir", "command",
    ]
    table_path = report_dir / "desc_scname_updated_values_all_runs.csv"
    with table_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    grouped: dict[tuple[str, str], list[dict]] = {}
    for row in rows:
        grouped.setdefault((row.get("dataset", ""), row.get("method", "")), []).append(row)

    summary_fields = [
        "value_tag", "dataset", "method", "n_success", "n_total", "status_summary",
        *[f"{m}_mean" for m in METRICS],
        *[f"{m}_std" for m in METRICS],
    ]
    summary_rows = []
    for (dataset, method), group in sorted(grouped.items()):
        success = [r for r in group if r.get("status") == "success"]
        out = {
            "value_tag": VALUE_TAG,
            "dataset": dataset,
            "method": method,
            "n_success": len(success),
            "n_total": len(group),
            "status_summary": ";".join(f"seed{r.get('seed')}:{r.get('status')}" for r in sorted(group, key=lambda x: x.get("seed", 0))),
        }
        for metric in METRICS:
            vals = []
            for row in success:
                try:
                    val = float(row.get(metric, ""))
                    if math.isfinite(val):
                        vals.append(val)
                except (TypeError, ValueError):
                    pass
            out[f"{metric}_mean"] = sum(vals) / len(vals) if vals else ""
            out[f"{metric}_std"] = (
                (sum((v - out[f"{metric}_mean"]) ** 2 for v in vals) / (len(vals) - 1)) ** 0.5
                if len(vals) > 1
                else (0.0 if vals else "")
            )
        summary_rows.append(out)

    summary_path = report_dir / "desc_scname_updated_values_mean_std.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(summary_rows)


def build_tasks(args: argparse.Namespace) -> list[dict]:
    datasets = list_datasets(args.data_root)
    tasks = []
    for dataset, data_path in datasets:
        try:
            n_clusters = infer_n_clusters(data_path)
        except Exception as exc:
            for method in args.methods:
                for seed in args.seeds:
                    out_dir = args.results_dir / dataset / f"{method}__seed{seed}"
                    write_json(out_dir / "status.json", {
                        "value_tag": VALUE_TAG,
                        "dataset": dataset,
                        "data_path": str(data_path),
                        "method": method,
                        "seed": seed,
                        "epochs": args.epochs,
                        "pretrain_epochs": args.epochs,
                        "n_clusters": "",
                        "status": "skipped",
                        "return_code": None,
                        "error": f"n_clusters inference failed: {exc}",
                        "gpu": "",
                        "no_cuda": False,
                        "elapsed_seconds": 0.0,
                        "save_dir": str(out_dir),
                        "command": "",
                    })
            continue
        for method in args.methods:
            for seed in args.seeds:
                gpu = args.gpus[len(tasks) % len(args.gpus)]
                tasks.append({
                    "dataset": dataset,
                    "data_path": data_path,
                    "method": method,
                    "seed": seed,
                    "epochs": args.epochs,
                    "n_clusters": n_clusters,
                    "gpu": gpu,
                    "out_dir": args.results_dir / dataset / f"{method}__seed{seed}",
                    "timeout_seconds": args.timeout_hours * 3600,
                })
    return tasks


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data_root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--results_dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--report_dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--methods", nargs="+", default=["desc", "scname"], choices=sorted(METHODS))
    parser.add_argument("--seeds", nargs="+", type=int, default=[42, 2024, 3407])
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--gpus", nargs="+", type=int, default=[1, 2, 3, 4, 5, 6])
    parser.add_argument("--max_workers", type=int, default=3)
    parser.add_argument("--timeout_hours", type=int, default=6)
    parser.add_argument("--summary_only", action="store_true")
    args = parser.parse_args()

    forbidden = {0, 7} & set(args.gpus)
    if forbidden:
        raise SystemExit(f"Forbidden GPU(s) requested: {sorted(forbidden)}")
    if not PYTHON.exists():
        raise SystemExit(f"Python executable missing: {PYTHON}")

    if args.summary_only:
        write_tables(args.results_dir, args.report_dir)
        return 0

    tasks = build_tasks(args)
    write_tables(args.results_dir, args.report_dir)
    print(f"Prepared {len(tasks)} runnable tasks under {args.results_dir}")
    print(f"Reports: {args.report_dir}")

    pending = iter(tasks)
    running = {}
    completed = 0
    with ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        while True:
            while len(running) < args.max_workers:
                try:
                    task = next(pending)
                except StopIteration:
                    break
                future = pool.submit(run_one, task)
                running[future] = task
                print(f"[start] {task['dataset']} {task['method']} seed{task['seed']} gpu{task['gpu']}", flush=True)
            if not running:
                break
            done, _ = wait(running.keys(), return_when=FIRST_COMPLETED)
            for future in done:
                task = running.pop(future)
                result = future.result()
                completed += 1
                print(
                    f"[done] {completed}/{len(tasks)} {task['dataset']} {task['method']} "
                    f"seed{task['seed']} status={result.get('status')} "
                    f"elapsed={result.get('elapsed_seconds')}s",
                    flush=True,
                )
                write_tables(args.results_dir, args.report_dir)
    write_tables(args.results_dir, args.report_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
