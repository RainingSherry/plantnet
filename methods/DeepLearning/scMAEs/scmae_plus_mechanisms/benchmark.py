#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import os
import re
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

MECH_DIR = Path(__file__).resolve().parent
ROOT = MECH_DIR.parents[3]
DATA_ROOT = Path("/data/luolie/biopipeline/scCluBench/data/scMAE")
TARGET_DATASETS = ["Melanoma_5K", "Quake_10x_Spleen", "Macosko"]
SCREEN_SEEDS = [42, 2024, 3407]
FORBIDDEN_GPUS = {0, 7}
REQUIRED_VARIANT_FILES = ["model.py", "loss.py", "run.py", "README.md"]
ELAPSED_PATTERN = re.compile(r"^ELAPSED_SECONDS:\s*([0-9.]+)", re.MULTILINE)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Benchmark scMAE-plus mechanism variants.")
    parser.add_argument("--stage", choices=["screen", "formal", "collect"], default="screen")
    parser.add_argument("--variants", nargs="+", default=None)
    parser.add_argument("--datasets", nargs="+", default=None)
    parser.add_argument("--seeds", nargs="+", type=int, default=SCREEN_SEEDS)
    parser.add_argument("--out_dir", default=str(MECH_DIR / "runs"))
    parser.add_argument("--screen_epochs", type=int, default=80)
    parser.add_argument("--formal_epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--rerun_existing", action="store_true")
    parser.add_argument("--skip_eval", action="store_true")
    parser.add_argument("--run_label", default=None, help="Optional suffix used for save_dir and variant_name.")
    parser.add_argument(
        "--extra_args",
        nargs=argparse.REMAINDER,
        default=[],
        help="Arguments passed through to the variant run.py; place this flag last.",
    )
    return parser.parse_args()


def discover_variants() -> list[str]:
    variants_dir = MECH_DIR / "variants"
    names = []
    for path in sorted(variants_dir.iterdir()):
        if not path.is_dir():
            continue
        missing = [name for name in REQUIRED_VARIANT_FILES if not (path / name).exists()]
        if not missing:
            names.append(path.name)
    return names


def validate_variant(name: str) -> Path:
    path = MECH_DIR / "variants" / name
    if not path.is_dir():
        raise SystemExit(f"Unknown scMAE-plus variant: {name}")
    missing = [file_name for file_name in REQUIRED_VARIANT_FILES if not (path / file_name).exists()]
    if missing:
        raise SystemExit(f"{name} is incomplete; missing: {', '.join(missing)}")
    return path


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
    variant_dir: Path,
    dataset: str,
    save_dir: Path,
    variant_name: str,
    epochs: int,
    batch_size: int,
    n_top_genes: int,
    seed: int,
    gpu: int,
    no_cuda: bool,
    skip_eval: bool,
    extra_args: list[str],
) -> list[str]:
    cmd = [
        sys.executable,
        str(variant_dir / "run.py"),
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
    if "--variant_name" not in extra_args:
        cmd.extend(["--variant_name", variant_name])
    if skip_eval:
        cmd.extend(["--skip_eval", "true"])
    if no_cuda:
        cmd.append("--no_cuda")
    else:
        if gpu in FORBIDDEN_GPUS:
            raise SystemExit("GPU 0 and GPU 7 are forbidden.")
        cmd.extend(["--gpu", str(gpu)])
    cmd.extend(extra_args)
    return cmd


def labeled_variant_name(variant: str, run_label: str | None) -> str:
    if not run_label:
        return variant
    safe = "".join(ch if ch.isalnum() or ch in {"_", "-"} else "_" for ch in str(run_label))
    return f"{variant}_{safe}"


def run_one(args: argparse.Namespace, variant: str, dataset: str, seed: int, gpu: int) -> dict:
    variant_dir = validate_variant(variant)
    epochs = args.screen_epochs if args.stage == "screen" else args.formal_epochs
    variant_name = labeled_variant_name(variant, args.run_label)
    save_dir = Path(args.out_dir) / args.stage / dataset / f"seed{seed}" / variant_name
    eval_path = save_dir / "eval_fixed.csv"
    if eval_path.exists() and not args.rerun_existing:
        return {
            "stage": args.stage,
            "dataset": dataset,
            "variant": variant_name,
            "base_variant": variant,
            "seed": seed,
            "status": "skipped_existing",
            "return_code": 0,
            "elapsed_seconds": 0.0,
            "save_dir": str(save_dir),
        }
    save_dir.mkdir(parents=True, exist_ok=True)
    cmd = build_command(
        variant_dir=variant_dir,
        dataset=dataset,
        save_dir=save_dir,
        variant_name=variant_name,
        epochs=epochs,
        batch_size=args.batch_size,
        n_top_genes=args.n_top_genes,
        seed=seed,
        gpu=gpu,
        no_cuda=args.no_cuda,
        skip_eval=args.skip_eval,
        extra_args=args.extra_args,
    )
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = "4"
    env["OPENBLAS_NUM_THREADS"] = "4"
    env["MKL_NUM_THREADS"] = "4"
    env["NUMEXPR_NUM_THREADS"] = "4"
    start = time.time()
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
        "variant": variant_name,
        "base_variant": variant,
        "seed": seed,
        "status": "success" if proc.returncode == 0 and (args.skip_eval or eval_path.exists()) else "failed",
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
        log_path = run_dir / "run.log"
        if log_path.exists():
            match = ELAPSED_PATTERN.search(log_path.read_text(encoding="utf-8", errors="replace"))
            if match:
                row["elapsed_seconds"] = float(match.group(1))
        diag_path = run_dir / "diagnostics.json"
        if diag_path.exists():
            diag = pd.read_json(diag_path, typ="series").to_dict()
            for key in [
                "effective_mask_rate",
                "embedding_variance_mean",
                "embedding_variance_min",
                "embedding_mean_cosine",
                "cluster_mass_max",
                "cluster_mass_min",
                "prototype_confidence_mean",
                "neighbor_reliability_mean",
                "neighbor_reliability_min",
                "neighbor_reliability_all_edge_mean",
                "neighbor_mutual_edge_fraction",
                "neighbor_reliable_edge_fraction",
                "neighbor_reliable_count_mean",
                "neighbor_reliable_count_min",
                "neighbor_score_mean",
                "neighbor_score_p25",
                "neighbor_score_p50",
                "neighbor_score_p75",
                "neighbor_similarity_mean",
                "neighbor_shared_score_mean",
                "neighbor_pseudo_confidence_mean",
                "neighbor_pseudo_same_cluster_edge_fraction",
                "neighbor_pseudo_same_cluster_reliable_fraction",
                "neighbor_consensus_hit_mean",
                "neighbor_adaptive_score_mean",
                "neighbor_adaptive_hit_mean",
                "neighbor_adaptive_core_edge_fraction",
                "neighbor_adaptive_strict_edge_fraction",
            ]:
                row[f"diag_{key}"] = diag.get(key)
        rows.append(row)
    return pd.DataFrame(rows)


def write_outputs(out_dir: Path) -> None:
    df = collect_rows(out_dir)
    if df.empty:
        print("No eval_fixed.csv files found.", flush=True)
        return
    per_run = MECH_DIR / "机制快筛单次结果.csv"
    summary = MECH_DIR / "机制快筛汇总结果.csv"
    attempts = MECH_DIR / "机制尝试记录.csv"
    df.to_csv(per_run, index=False)
    attempt_cols = [col for col in ["dataset", "method", "variant", "seed", "elapsed_seconds", "run_dir"] if col in df.columns]
    attempt_df = df[attempt_cols].copy()
    attempt_df["status"] = "success"
    attempt_df.to_csv(attempts, index=False)
    agg = (
        df.groupby(["dataset", "method", "variant"], dropna=False)
        .agg(
            acc_mean=("acc", "mean"),
            acc_std=("acc", "std"),
            acc_count=("acc", "count"),
            nmi_mean=("nmi", "mean"),
            nmi_std=("nmi", "std"),
            ari_mean=("ari", "mean"),
            ari_std=("ari", "std"),
            ari_count=("ari", "count"),
            f1_macro_mean=("f1_macro", "mean"),
            f1_macro_std=("f1_macro", "std"),
            seeds=("seed", lambda values: ",".join(str(int(v)) for v in sorted(set(values)))),
            n_seeds=("seed", lambda values: len(set(values))),
            effective_mask_rate_mean=("diag_effective_mask_rate", "mean"),
            embedding_variance_mean=("diag_embedding_variance_mean", "mean"),
            prototype_confidence_mean=("diag_prototype_confidence_mean", "mean"),
            neighbor_reliability_mean=("diag_neighbor_reliability_mean", "mean"),
            neighbor_reliable_edge_fraction_mean=("diag_neighbor_reliable_edge_fraction", "mean"),
            neighbor_reliable_count_mean=("diag_neighbor_reliable_count_mean", "mean"),
            neighbor_score_mean=("diag_neighbor_score_mean", "mean"),
            neighbor_pseudo_confidence_mean=("diag_neighbor_pseudo_confidence_mean", "mean"),
            neighbor_pseudo_same_cluster_edge_fraction_mean=("diag_neighbor_pseudo_same_cluster_edge_fraction", "mean"),
            neighbor_pseudo_same_cluster_reliable_fraction_mean=("diag_neighbor_pseudo_same_cluster_reliable_fraction", "mean"),
            neighbor_consensus_hit_mean=("diag_neighbor_consensus_hit_mean", "mean"),
            neighbor_adaptive_score_mean=("diag_neighbor_adaptive_score_mean", "mean"),
            neighbor_adaptive_hit_mean=("diag_neighbor_adaptive_hit_mean", "mean"),
            neighbor_adaptive_core_edge_fraction_mean=("diag_neighbor_adaptive_core_edge_fraction", "mean"),
            neighbor_adaptive_strict_edge_fraction_mean=("diag_neighbor_adaptive_strict_edge_fraction", "mean"),
        )
        .reset_index()
    )
    agg.to_csv(summary, index=False)
    print(f"Wrote {per_run}", flush=True)
    print(f"Wrote {summary}", flush=True)
    print(f"Wrote {attempts}", flush=True)


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    if args.stage == "collect":
        write_outputs(out_dir)
        return
    variants = args.variants or discover_variants()
    datasets = args.datasets or (["Melanoma_5K"] if args.stage == "screen" else TARGET_DATASETS)
    gpu = args.gpu if args.gpu is not None else (1 if args.no_cuda else detect_free_gpu())
    rows = []
    for variant in variants:
        for dataset in datasets:
            for seed in args.seeds:
                row = run_one(args, variant, dataset, seed, gpu)
                rows.append(row)
                print(
                    f"{row['stage']} {row['dataset']} seed={row['seed']} "
                    f"{row['variant']} {row['status']} rc={row['return_code']}",
                    flush=True,
                )
                write_status(rows, MECH_DIR / "机制尝试记录.csv")
                write_status(rows, out_dir / f"{args.stage}_status.csv")
    write_outputs(out_dir / args.stage)


if __name__ == "__main__":
    main()
