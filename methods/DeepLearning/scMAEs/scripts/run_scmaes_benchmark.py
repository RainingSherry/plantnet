#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Dict, Iterable, List

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

import pandas as pd

SCMAES_DIR = Path(__file__).resolve().parents[1]
ROOT = SCMAES_DIR.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning.scMAEs.common.runner import detect_free_gpu
from methods.DeepLearning.scMAEs.common.variant_defs import VARIANTS

DATA_ROOT = Path("/data/luolie/biopipeline/scCluBench/data/scMAE")
TARGET_DATASETS = ["Melanoma_5K", "Quake_10x_Spleen", "Macosko"]
FORBIDDEN_GPUS = {0, 7}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run scMAEs smoke/screen/formal benchmarks.")
    parser.add_argument("--stage", choices=["smoke", "screen", "formal", "collect"], default="smoke")
    parser.add_argument("--out_dir", default=str(SCMAES_DIR / "benchmark_runs"))
    parser.add_argument("--datasets", nargs="+", default=TARGET_DATASETS)
    parser.add_argument("--variants", nargs="+", default=list(VARIANTS))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--screen_epochs", type=int, default=25)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--max_parallel", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--append_global", action="store_true")
    parser.add_argument("--rerun_existing", action="store_true")
    parser.add_argument("--per_run_csv", default=None)
    parser.add_argument("--aggregate_csv", default=None)
    return parser.parse_args()


def run_command(cmd: List[str], log_path: Path) -> Dict[str, object]:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    env["OPENBLAS_NUM_THREADS"] = "4"
    env["MKL_NUM_THREADS"] = "4"
    env["NUMEXPR_NUM_THREADS"] = "4"
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, env=env)
    elapsed = time.time() - start
    with open(log_path, "w", encoding="utf-8") as handle:
        handle.write("COMMAND: " + " ".join(cmd) + "\n")
        handle.write(f"EXIT_CODE: {proc.returncode}\n")
        handle.write(f"ELAPSED_SECONDS: {elapsed:.1f}\n\n")
        handle.write("STDOUT\n")
        handle.write(proc.stdout)
        handle.write("\nSTDERR\n")
        handle.write(proc.stderr)
    return {"return_code": proc.returncode, "elapsed_seconds": elapsed, "log": str(log_path)}


def variant_run_py(variant: str) -> Path:
    return SCMAES_DIR / "variants" / variant / "run.py"


def smoke_one(variant: str, out_dir: Path, gpu: int, no_cuda: bool) -> Dict[str, object]:
    cfg = VARIANTS[variant]
    save_dir = out_dir / "smoke" / variant
    cmd = [
        sys.executable,
        str(variant_run_py(variant)),
        "--save_dir", str(save_dir),
        "--data_path", str(DATA_ROOT / "Melanoma_5K.h5"),
        "--n_clusters", "8",
        "--epochs", "1",
        "--batch_size", "16",
        "--seed", "42",
        "--smoke",
    ]
    if no_cuda:
        cmd.append("--no_cuda")
    else:
        cmd.extend(["--gpu", str(gpu)])
    result = run_command(cmd, save_dir / "run.log")
    result.update({"variant": variant, "rank": cfg.rank, "stage": "smoke"})
    return result


def run_one_dataset(args: argparse.Namespace, variant: str, dataset: str, gpu: int, epochs: int) -> Dict[str, object]:
    cfg = VARIANTS[variant]
    save_dir = Path(args.out_dir) / args.stage / dataset / variant
    if not args.rerun_existing and (save_dir / "eval_fixed.csv").exists():
        return {
            "variant": variant,
            "rank": cfg.rank,
            "dataset": dataset,
            "stage": args.stage,
            "return_code": 0,
            "elapsed_seconds": 0.0,
            "log": str(save_dir / "run.log"),
            "status": "skipped_existing",
        }
    cmd = [
        sys.executable,
        str(variant_run_py(variant)),
        "--save_dir", str(save_dir),
        "--data_path", str(DATA_ROOT / f"{dataset}.h5"),
        "--dataset_name", dataset,
        "--n_clusters", "0",
        "--epochs", str(epochs),
        "--batch_size", str(args.batch_size),
        "--seed", str(args.seed),
        "--no_save_h5ad",
    ]
    if args.no_cuda:
        cmd.append("--no_cuda")
    else:
        cmd.extend(["--gpu", str(gpu)])
    result = run_command(cmd, save_dir / "run.log")
    result.update({"variant": variant, "rank": cfg.rank, "dataset": dataset, "stage": args.stage})
    return result


def write_status(rows: Iterable[Dict[str, object]], path: Path) -> None:
    rows = list(rows)
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def collect_rows(out_dir: Path) -> pd.DataFrame:
    rows = []
    for eval_path in out_dir.rglob("eval_fixed.csv"):
        row = pd.read_csv(eval_path).iloc[0].to_dict()
        summary_path = eval_path.parent / "summary.json"
        if summary_path.exists():
            with open(summary_path, encoding="utf-8") as handle:
                summary = json.load(handle)
            row["runtime_seconds"] = summary.get("runtime_seconds", row.get("runtime_seconds", ""))
            row["source_paper"] = summary.get("source_paper", row.get("source_paper", ""))
        row["run_dir"] = str(eval_path.parent)
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def write_outputs(out_dir: Path, append_global: bool, per_run_csv: str | None = None, aggregate_csv: str | None = None) -> None:
    df = collect_rows(out_dir)
    if df.empty:
        print("No eval_fixed.csv files found.")
        return
    per_run = Path(per_run_csv) if per_run_csv else SCMAES_DIR / "新模型单次结果.csv"
    df.to_csv(per_run, index=False)
    agg = (
        df.groupby(["dataset", "method"], dropna=False)
        .agg(
            acc_mean=("acc", "mean"),
            acc_std=("acc", "std"),
            acc_count=("acc", "count"),
            nmi_mean=("nmi", "mean"),
            nmi_std=("nmi", "std"),
            nmi_count=("nmi", "count"),
            ari_mean=("ari", "mean"),
            ari_std=("ari", "std"),
            ari_count=("ari", "count"),
            f1_macro_mean=("f1_macro", "mean"),
            f1_macro_std=("f1_macro", "std"),
            f1_macro_count=("f1_macro", "count"),
            runtime_seconds_mean=("runtime_seconds", "mean"),
            runtime_seconds_std=("runtime_seconds", "std"),
            seed__lambda=("seed", lambda values: ",".join(str(int(v)) for v in sorted(set(values)))),
            n_seeds=("seed", lambda values: len(set(values))),
            method_raw=("variant", "first"),
        )
        .reset_index()
    )
    agg = agg.rename(columns={"seed__lambda": "seed_<lambda>"})
    aggregate_path = Path(aggregate_csv) if aggregate_csv else SCMAES_DIR / "新模型汇总结果.csv"
    agg.to_csv(aggregate_path, index=False)
    if append_global:
        global_path = SCMAES_DIR / "全benchmark结果.csv"
        old = pd.read_csv(global_path) if global_path.exists() else pd.DataFrame()
        combined = pd.concat([old, agg], ignore_index=True, sort=False)
        combined.to_csv(global_path, index=False)
    print(f"Wrote {per_run}")
    print(f"Wrote {aggregate_path}")


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    gpu = detect_free_gpu() if args.gpu is None and not args.no_cuda else (args.gpu or 1)
    if gpu in FORBIDDEN_GPUS and not args.no_cuda:
        raise SystemExit("GPU 0 and GPU 7 are forbidden.")
    variants = [variant for variant in args.variants if variant in VARIANTS]

    if args.stage == "collect":
        write_outputs(out_dir, args.append_global, args.per_run_csv, args.aggregate_csv)
        return

    rows = []
    if args.stage == "smoke":
        with ThreadPoolExecutor(max_workers=max(1, min(args.max_parallel, 4))) as pool:
            futures = [pool.submit(smoke_one, variant, out_dir, gpu, args.no_cuda) for variant in variants]
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                print(f"{row['stage']} {row['variant']} rc={row['return_code']}", flush=True)
        write_status(rows, out_dir / "smoke_status.csv")
        return

    if args.stage == "screen":
        max_workers = max(1, min(args.max_parallel, 2))
        with ThreadPoolExecutor(max_workers=max_workers) as pool:
            futures = [
                pool.submit(run_one_dataset, args, variant, "Melanoma_5K", gpu, args.screen_epochs)
                for variant in variants
            ]
            for future in as_completed(futures):
                row = future.result()
                rows.append(row)
                print(f"{row['dataset']} {row['variant']} rc={row['return_code']}", flush=True)
        write_status(rows, out_dir / "screen_status.csv")
        return

    if args.stage == "formal":
        for variant in variants:
            for dataset in args.datasets:
                row = run_one_dataset(args, variant, dataset, gpu, args.epochs)
                rows.append(row)
                status = row.get("status", "run")
                print(f"{row['dataset']} {row['variant']} rc={row['return_code']} status={status}", flush=True)
                write_status(rows, out_dir / "formal_status.csv")
        write_outputs(out_dir, args.append_global, args.per_run_csv, args.aggregate_csv)


if __name__ == "__main__":
    main()
