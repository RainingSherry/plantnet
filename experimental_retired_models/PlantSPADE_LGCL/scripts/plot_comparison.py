#!/usr/bin/env python
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

from aggregate_results import read_metric_csvs

SCRIPT_DIR = Path(__file__).resolve().parent
PKG_DIR = SCRIPT_DIR.parent
ROOT = next(parent for parent in [SCRIPT_DIR, *SCRIPT_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())

DEFAULT_SELECTION = {
    "traditional_pca": ("pca", "kmeans_known_k"),
    "traditional_leiden": ("leiden", "leiden_fixed"),
    "traditional_louvain": ("louvain", "louvain_fixed"),
    "scmae": ("embedding", "kmeans_known_k"),
    "phytocluster": ("embedding", "kmeans_known_k"),
    "scvi": ("embedding", "kmeans_known_k"),
    "plantspade_lgcl_baseline": ("baseline", "kmeans_known_k"),
    "plantspade_lgcl_support_attention": ("support_attention", "kmeans_known_k"),
    "plantspade_lgcl_attention_no_idf": ("attention_no_idf", "kmeans_known_k"),
}

DISPLAY_NAMES = {
    "traditional_pca": "PCA",
    "traditional_leiden": "Leiden",
    "traditional_louvain": "Louvain",
    "scmae": "scMAE",
    "phytocluster": "PhytoCluster",
    "scvi": "scVI",
    "plantspade_lgcl_baseline": "PlantSPADE",
    "plantspade_lgcl_support_attention": "PlantSPADE+SGA",
    "plantspade_lgcl_attention_no_idf": "PlantSPADE+SGA no-IDF",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Plot a fixed-protocol method comparison from aggregated PlantSPADE results.")
    parser.add_argument("--results_dir", default=None)
    parser.add_argument("--main_config", default=str(PKG_DIR / "configs" / "main_lgcl.yaml"))
    parser.add_argument("--tables_dir", default=None)
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--metrics", default="nmi,ari")
    parser.add_argument("--methods", default=None)
    parser.add_argument("--datasets", default=None)
    parser.add_argument("--seeds", default=None)
    parser.add_argument("--output_name", default="comparison_fixed_protocol.png")
    return parser.parse_args()


def split_csv(value: str | None):
    if value is None:
        return None
    return [item.strip() for item in str(value).split(",") if item.strip()]


def default_results_dir(main_config: str) -> Path:
    with open(main_config, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    return Path(cfg.get("output_dir", ROOT / "results" / "PlantSPADE_LGCL_protocol"))


def load_results(results_dir: Path, tables_dir: Path | None) -> pd.DataFrame:
    table_path = (tables_dir / "all_results_long.csv") if tables_dir else (results_dir / "tables" / "all_results_long.csv")
    if table_path.exists():
        df = pd.read_csv(table_path)
    else:
        df = read_metric_csvs(results_dir)
    if df.empty:
        raise SystemExit(f"No metric CSV rows found under {results_dir}")
    return df


def source_priority(path_text: str) -> tuple[int, float]:
    path = Path(str(path_text))
    name = path.name
    preferred_prefixes = (
        "eval_baseline_fixed",
        "eval_support_attention_fixed",
        "eval_attention_no_idf_fixed",
        "eval_traditional_leiden_fixed",
        "eval_traditional_louvain_fixed",
        "external_eval_fixed",
        "pca_fixed",
    )
    priority = 0 if name.startswith(preferred_prefixes) else 1
    try:
        mtime = path.stat().st_mtime
    except OSError:
        mtime = 0.0
    return priority, -mtime


def infer_suite_method(path_text: str) -> str:
    parts = Path(str(path_text)).parts
    for idx, part in enumerate(parts):
        if part.startswith("seed_") and idx >= 1:
            return parts[idx - 1]
    return ""


def infer_variant(path_text: str) -> str:
    name = Path(str(path_text)).name
    if name.startswith("pca_fixed"):
        return "pca"
    if name.startswith("external_eval_fixed"):
        return "embedding"
    if name.startswith("eval_traditional_leiden_fixed"):
        return "leiden"
    if name.startswith("eval_traditional_louvain_fixed"):
        return "louvain"
    if name.startswith("eval_baseline_fixed"):
        return "baseline"
    if name.startswith("eval_support_attention_fixed"):
        return "support_attention"
    if name.startswith("eval_attention_no_idf_fixed"):
        return "attention_no_idf"
    return ""


def normalize_keys(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    if "source_file" in df.columns:
        source_method = df["source_file"].map(infer_suite_method).astype(str)
        source_variant = df["source_file"].map(infer_variant).astype(str)
    else:
        source_method = pd.Series("", index=df.index)
        source_variant = pd.Series("", index=df.index)

    csv_method = df.get("method", pd.Series("", index=df.index)).astype(str)
    df["_method_key"] = source_method.where(source_method.ne(""), csv_method)

    if "variant" in df.columns:
        csv_variant = df["variant"].astype(str)
        missing_variant = df["variant"].isna() | csv_variant.isin(["", "nan", "None"])
    else:
        csv_variant = pd.Series("", index=df.index)
        missing_variant = pd.Series(True, index=df.index)
    df["_variant_key"] = csv_variant.where(~missing_variant, source_variant)
    df["_variant_key"] = df["_variant_key"].where(df["_variant_key"].ne(""), source_variant)
    return df


def selected_rows(
    df: pd.DataFrame,
    methods: list[str] | None,
    datasets: list[str] | None,
    seeds: list[str] | None,
) -> pd.DataFrame:
    df = df[df.get("protocol", "").astype(str).eq("fixed")].copy()
    df = normalize_keys(df)
    if methods:
        df = df[df["_method_key"].astype(str).isin(methods)].copy()
    if datasets:
        df = df[df["dataset"].astype(str).isin(datasets)].copy()
    if seeds:
        df = df[df["seed"].astype(str).isin(seeds)].copy()
    rows = []
    for method, (variant, cluster_method) in DEFAULT_SELECTION.items():
        block = df[
            df["_method_key"].astype(str).eq(method)
            & df["cluster_method"].astype(str).eq(cluster_method)
            & df["_variant_key"].astype(str).eq(variant)
        ].copy()
        if methods and method not in methods:
            continue
        if block.empty:
            continue
        if "source_file" in block.columns:
            priorities = block["source_file"].map(source_priority)
            block["_priority"] = [item[0] for item in priorities]
            block["_mtime_desc"] = [item[1] for item in priorities]
            block = block.sort_values(["dataset", "method", "seed", "_priority", "_mtime_desc"])
        block = block.drop_duplicates(subset=["dataset", "_method_key", "seed", "cluster_method", "_variant_key"], keep="first")
        block["method"] = method
        block["variant"] = variant
        block["display_method"] = DISPLAY_NAMES.get(method, method)
        rows.append(block)
    if not rows:
        raise SystemExit("No selected fixed-protocol rows matched the requested methods/datasets.")
    return pd.concat(rows, ignore_index=True, sort=False)


def plot_heatmaps(summary: pd.DataFrame, metrics: list[str], output_path: Path) -> None:
    datasets = sorted(summary["dataset"].astype(str).unique())
    methods = list(dict.fromkeys(summary["display_method"].astype(str).tolist()))
    fig, axes = plt.subplots(1, len(metrics), figsize=(max(7.5, 1.1 * len(methods) * len(metrics)), max(4.5, 0.55 * len(datasets) + 1.5)))
    if len(metrics) == 1:
        axes = [axes]

    for ax, metric in zip(axes, metrics):
        pivot = (
            summary.pivot_table(index="dataset", columns="display_method", values=f"{metric}_mean", aggfunc="mean")
            .reindex(index=datasets, columns=methods)
        )
        values = pivot.to_numpy(dtype=float)
        im = ax.imshow(values, vmin=0.0, vmax=1.0, cmap="viridis", aspect="auto")
        ax.set_title(metric.upper())
        ax.set_xticks(np.arange(len(methods)))
        ax.set_xticklabels(methods, rotation=35, ha="right")
        ax.set_yticks(np.arange(len(datasets)))
        ax.set_yticklabels(datasets)
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                value = values[i, j]
                if np.isfinite(value):
                    ax.text(j, i, f"{value:.2f}", ha="center", va="center", color="white" if value < 0.62 else "black", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.suptitle("Fixed-protocol comparison across datasets and seeds", y=1.02)
    fig.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main():
    args = parse_args()
    results_dir = Path(args.results_dir) if args.results_dir else default_results_dir(args.main_config)
    tables_dir = Path(args.tables_dir) if args.tables_dir else results_dir / "tables"
    output_dir = Path(args.output_dir) if args.output_dir else results_dir / "figures"
    metrics = split_csv(args.metrics) or ["nmi", "ari"]
    methods = split_csv(args.methods)
    datasets = split_csv(args.datasets)
    seeds = split_csv(args.seeds)
    output_dir.mkdir(parents=True, exist_ok=True)

    df = load_results(results_dir, tables_dir)
    selected = selected_rows(df, methods=methods, datasets=datasets, seeds=seeds)
    selected.to_csv(output_dir / "comparison_selected_fixed_rows.csv", index=False)

    grouped = selected.groupby(["dataset", "method", "display_method"], dropna=False)
    rows = []
    for keys, block in grouped:
        dataset, method, display_method = keys
        row = {
            "dataset": dataset,
            "method": method,
            "display_method": display_method,
            "n_runs": int(block["seed"].nunique()) if "seed" in block.columns else int(len(block)),
        }
        for metric in metrics:
            values = pd.to_numeric(block[metric], errors="coerce").to_numpy(dtype=float)
            row[f"{metric}_mean"] = float(np.nanmean(values)) if np.any(np.isfinite(values)) else float("nan")
            row[f"{metric}_std"] = float(np.nanstd(values, ddof=1)) if np.sum(np.isfinite(values)) > 1 else 0.0
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(output_dir / "comparison_dataset_mean_std.csv", index=False)

    overall_rows = []
    for method, block in summary.groupby(["method", "display_method"], dropna=False):
        method_name, display_method = method
        row = {"method": method_name, "display_method": display_method, "n_datasets": int(block["dataset"].nunique())}
        for metric in metrics:
            values = pd.to_numeric(block[f"{metric}_mean"], errors="coerce").to_numpy(dtype=float)
            row[f"{metric}_mean"] = float(np.nanmean(values)) if np.any(np.isfinite(values)) else float("nan")
        overall_rows.append(row)
    pd.DataFrame(overall_rows).sort_values([f"{metrics[0]}_mean"], ascending=False).to_csv(
        output_dir / "comparison_overall_mean.csv",
        index=False,
    )

    plot_heatmaps(summary, metrics=metrics, output_path=output_dir / args.output_name)
    print(f"Wrote comparison plot to {output_dir / args.output_name}")


if __name__ == "__main__":
    main()
