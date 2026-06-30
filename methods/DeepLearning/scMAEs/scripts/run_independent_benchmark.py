#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
import time
from pathlib import Path
from typing import Iterable

import pandas as pd

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

SCMAES_DIR = Path(__file__).resolve().parents[1]
ROOT = SCMAES_DIR.parents[2]
DATA_ROOT = Path("/data/luolie/biopipeline/scCluBench/data/scMAE")
TARGET_DATASETS = ["Melanoma_5K", "Quake_10x_Spleen", "Macosko"]
FORBIDDEN_GPUS = {0, 7}
REQUIRED_METHOD_FILES = ["model.py", "loss.py", "run.py", "README.md", "source_manifest.json"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run independent full scMAEs methods only.")
    parser.add_argument("--stage", choices=["screen", "formal", "collect"], default="screen")
    parser.add_argument("--methods", nargs="+", default=None)
    parser.add_argument("--datasets", nargs="+", default=None)
    parser.add_argument("--out_dir", default=str(SCMAES_DIR / "independent_runs"))
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--screen_epochs", type=int, default=80)
    parser.add_argument("--formal_epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--rerun_existing", action="store_true")
    parser.add_argument("--append_global", action="store_true")
    return parser.parse_args()


def discover_methods() -> list[str]:
    methods = []
    for path in sorted(SCMAES_DIR.iterdir()):
        if not path.is_dir():
            continue
        if not path.name.startswith("rank") or not path.name.endswith("_full"):
            continue
        missing = [name for name in REQUIRED_METHOD_FILES if not (path / name).exists()]
        if not missing:
            methods.append(path.name)
    return methods


def validate_method(method: str) -> Path:
    method_dir = SCMAES_DIR / method
    if not method_dir.is_dir():
        raise SystemExit(f"Unknown independent method directory: {method}")
    missing = [name for name in REQUIRED_METHOD_FILES if not (method_dir / name).exists()]
    if missing:
        raise SystemExit(f"{method} is not a complete independent method; missing: {', '.join(missing)}")
    return method_dir


def detect_free_gpu() -> int:
    result = subprocess.run(
        ["nvidia-smi", "--query-gpu=index,memory.used", "--format=csv,noheader,nounits"],
        text=True,
        capture_output=True,
        timeout=5,
    )
    if result.returncode != 0:
        return 1
    candidates = []
    for line in result.stdout.splitlines():
        parts = [part.strip() for part in line.split(",")]
        if len(parts) != 2:
            continue
        idx = int(parts[0])
        if idx in FORBIDDEN_GPUS:
            continue
        candidates.append((int(parts[1]), idx))
    return sorted(candidates)[0][1] if candidates else 1


def build_command(
    method_dir: Path,
    dataset: str,
    save_dir: Path,
    epochs: int,
    batch_size: int,
    n_top_genes: int,
    seed: int,
    gpu: int,
    no_cuda: bool,
) -> list[str]:
    cmd = [
        sys.executable,
        str(method_dir / "run.py"),
        "--data_path",
        str(DATA_ROOT / f"{dataset}.h5"),
        "--dataset_name",
        dataset,
        "--save_dir",
        str(save_dir),
        "--n_clusters",
        "0",
        "--n_top_genes",
        str(n_top_genes),
        "--epochs",
        str(epochs),
        "--batch_size",
        str(batch_size),
        "--seed",
        str(seed),
    ]
    if no_cuda:
        cmd.append("--no_cuda")
    else:
        if gpu in FORBIDDEN_GPUS:
            raise SystemExit("GPU 0 and GPU 7 are forbidden.")
        cmd.extend(["--gpu", str(gpu)])
    return cmd


def run_one(args: argparse.Namespace, method: str, dataset: str, gpu: int) -> dict:
    method_dir = validate_method(method)
    epochs = args.screen_epochs if args.stage == "screen" else args.formal_epochs
    save_dir = Path(args.out_dir) / args.stage / dataset / method
    eval_path = save_dir / "eval_fixed.csv"
    if eval_path.exists() and not args.rerun_existing:
        return {
            "stage": args.stage,
            "dataset": dataset,
            "method": method,
            "status": "skipped_existing",
            "return_code": 0,
            "elapsed_seconds": 0.0,
            "save_dir": str(save_dir),
        }
    save_dir.mkdir(parents=True, exist_ok=True)
    cmd = build_command(method_dir, dataset, save_dir, epochs, args.batch_size, args.n_top_genes, args.seed, gpu, args.no_cuda)
    start = time.time()
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    env["OPENBLAS_NUM_THREADS"] = "4"
    env["MKL_NUM_THREADS"] = "4"
    env["NUMEXPR_NUM_THREADS"] = "4"
    proc = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, env=env)
    elapsed = time.time() - start
    with open(save_dir / "run.log", "w", encoding="utf-8") as handle:
        handle.write("COMMAND: " + " ".join(cmd) + "\n")
        handle.write(f"EXIT_CODE: {proc.returncode}\n")
        handle.write(f"ELAPSED_SECONDS: {elapsed:.1f}\n\nSTDOUT\n")
        handle.write(proc.stdout)
        handle.write("\nSTDERR\n")
        handle.write(proc.stderr)
    return {
        "stage": args.stage,
        "dataset": dataset,
        "method": method,
        "status": "success" if proc.returncode == 0 and eval_path.exists() else "failed",
        "return_code": proc.returncode,
        "elapsed_seconds": elapsed,
        "save_dir": str(save_dir),
    }


def write_status(rows: Iterable[dict], path: Path) -> None:
    rows = list(rows)
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def collect_rows(out_dir: Path) -> pd.DataFrame:
    rows = []
    for eval_path in out_dir.rglob("eval_fixed.csv"):
        row = pd.read_csv(eval_path).iloc[0].to_dict()
        run_dir = eval_path.parent
        row["run_dir"] = str(run_dir)
        manifest_path = run_dir / "source_manifest.json"
        if manifest_path.exists():
            with open(manifest_path, encoding="utf-8") as handle:
                manifest = json.load(handle)
            row["source_manifest_method"] = manifest.get("method", "")
            row["source_manifest_paper"] = manifest.get("paper", "")
        summary_path = run_dir / "summary.json"
        if summary_path.exists():
            with open(summary_path, encoding="utf-8") as handle:
                summary = json.load(handle)
            row["summary_rank"] = summary.get("rank", row.get("rank", ""))
        rows.append(row)
    return pd.DataFrame(rows)


def write_outputs(out_dir: Path, append_global: bool) -> None:
    df = collect_rows(out_dir)
    if df.empty:
        print("No independent eval_fixed.csv files found.", flush=True)
        return
    per_run = SCMAES_DIR / "独立模型单次结果.csv"
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
            seed_values=("seed", lambda values: ",".join(str(int(v)) for v in sorted(set(values)))),
            n_seeds=("seed", lambda values: len(set(values))),
        )
        .reset_index()
    )
    agg = agg.rename(columns={"seed_values": "seed_<lambda>"})
    aggregate = SCMAES_DIR / "独立模型汇总结果.csv"
    agg.to_csv(aggregate, index=False)
    if append_global:
        global_path = SCMAES_DIR / "全benchmark结果.csv"
        old = pd.read_csv(global_path) if global_path.exists() else pd.DataFrame()
        pd.concat([old, agg], ignore_index=True, sort=False).to_csv(global_path, index=False)
    print(f"Wrote {per_run}", flush=True)
    print(f"Wrote {aggregate}", flush=True)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    if args.stage == "collect":
        write_outputs(out_dir, args.append_global)
        return
    methods = args.methods or discover_methods()
    datasets = args.datasets or (["Melanoma_5K"] if args.stage == "screen" else TARGET_DATASETS)
    gpu = args.gpu if args.gpu is not None else (1 if args.no_cuda else detect_free_gpu())
    rows = []
    for method in methods:
        for dataset in datasets:
            row = run_one(args, method, dataset, gpu)
            rows.append(row)
            print(f"{row['stage']} {row['dataset']} {row['method']} {row['status']} rc={row['return_code']}", flush=True)
            write_status(rows, out_dir / f"{args.stage}_status.csv")
    write_outputs(out_dir / args.stage, args.append_global)


if __name__ == "__main__":
    main()
