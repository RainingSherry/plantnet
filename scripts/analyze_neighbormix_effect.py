#!/usr/bin/env python3
"""
analyze_neighbormix_effect.py
=============================
Step 1: Compute NeighborMix gain (delta metrics vs noMix and vs scMAE).
Step 2: Correlate dataset properties with delta_ari.
Step 3: Label each dataset as positive / negative / neutral.
Step 4: Generate interpretation report.

Usage:
    python scripts/analyze_neighbormix_effect.py \\
        --property_csv results/analysis/dataset_property_summary.csv \\
        --benchmark_csv results/formal/benchmark_summary_mean_std.csv \\
        --out_dir results/analysis

    # Or use per-dataset summary CSVs:
    python scripts/analyze_neighbormix_effect.py \\
        --property_csv results/analysis/dataset_property_summary.csv \\
        --benchmark_dir results/formal \\
        --out_dir results/analysis
"""

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from collections import defaultdict

import numpy as np

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


METRICS = ["acc", "nmi", "ari", "f1_macro"]


def _parse_mean_std(val: str):
    """Parse '0.7123 ± 0.0145' → (0.7123, 0.0145). Returns (None, None) if N/A."""
    if not val or val == "N/A":
        return None, None
    parts = str(val).split(" ± ")
    try:
        mean = float(parts[0])
        std = float(parts[1]) if len(parts) > 1 else 0.0
        return mean, std
    except (ValueError, IndexError):
        return None, None


# ─── Step 1: Load benchmark summaries ──────────────────────────────────────────

def load_benchmark_csv(path: str) -> dict:
    """
    Load a benchmark_summary_mean_std.csv into a dict keyed by (dataset, method).
    Each value is a dict with 'acc_mean', 'nmi_mean', etc.
    """
    results = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            dataset = row.get("dataset", "")
            method = row.get("method", "")
            key = (dataset, method)
            results[key] = {}
            for m in METRICS:
                mean, std = _parse_mean_std(row.get(m, ""))
                if mean is not None:
                    results[key][f"{m}_mean"] = mean
                    results[key][f"{m}_std"] = std
            results[key]["n_success"] = int(row.get("n_success", 0))
            results[key]["n_total"] = int(row.get("n_total", 0))
    return results


def load_all_benchmarks(benchmark_dir: str) -> dict:
    """Scan benchmark_dir for *_mean_std.csv or benchmark_summary_mean_std.csv."""
    results = {}
    benchmark_dir = Path(benchmark_dir)
    for csv_path in sorted(benchmark_dir.rglob("benchmark_summary_mean_std.csv")):
        partial = load_benchmark_csv(str(csv_path))
        results.update(partial)
    # Also try the combined file at the root
    combined = benchmark_dir / "benchmark_summary_combined.csv"
    if combined.exists():
        partial = load_benchmark_csv(str(combined))
        results.update(partial)
    return results


# ─── Step 2: Compute NeighborMix gain ─────────────────────────────────────────

def compute_gain_summary(benchmark_data: dict, out_path: str) -> list:
    """
    For each dataset, compute:
      delta = NeighborMix - noMix
      delta_vs_scmae = NeighborMix - scMAE

    Produces neighbormix_gain_summary.csv.
    """
    datasets = sorted(set(k[0] for k in benchmark_data.keys()))
    rows = []

    for ds in datasets:
        row = {"dataset": ds}

        def get_method_metric(data, dataset, method_key, metric):
            """Return the {metric}_mean value for (dataset, method_key)."""
            key = (dataset, method_key)
            method_data = data.get(key, {})
            return method_data.get(f"{metric}_mean")

        for m in METRICS:
            nm_val = get_method_metric(benchmark_data, ds, "neighbormix_scmae", m)
            nomix_val = get_method_metric(benchmark_data, ds, "nm_scmae_nomix", m)
            scmae_val = get_method_metric(benchmark_data, ds, "scmae", m)

            row[f"{m}_neighbormix"] = round(nm_val, 4) if nm_val is not None else None
            row[f"{m}_nomix"] = round(nomix_val, 4) if nomix_val is not None else None
            row[f"{m}_scmae"] = round(scmae_val, 4) if scmae_val is not None else None

            delta_vs_nomix = round(nm_val - nomix_val, 4) if (nm_val is not None and nomix_val is not None) else None
            delta_vs_scmae = round(nm_val - scmae_val, 4) if (nm_val is not None and scmae_val is not None) else None

            row[f"delta_{m}_vs_nomix"] = delta_vs_nomix
            row[f"delta_{m}_vs_scmae"] = delta_vs_scmae

        # Effect classification (based on delta_ari_vs_nomix)
        delta_ari = row.get("delta_ari_vs_nomix")
        if delta_ari is not None:
            if delta_ari >= 0.02:
                effect = "positive"
            elif delta_ari <= -0.02:
                effect = "negative"
            else:
                effect = "neutral"
        else:
            effect = "insufficient_data"
        row["effect_group"] = effect

        rows.append(row)

    # Write CSV
    fieldnames = ["dataset"] + [f"{m}_neighbormix" for m in METRICS] \
        + [f"{m}_nomix" for m in METRICS] \
        + [f"{m}_scmae" for m in METRICS] \
        + [f"delta_{m}_vs_nomix" for m in METRICS] \
        + [f"delta_{m}_vs_scmae" for m in METRICS] \
        + ["effect_group"]

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {out_path}")
    return rows


# ─── Step 3: Correlation analysis ───────────────────────────────────────────────

def load_properties(path: str) -> dict:
    """Load dataset_property_summary.csv into dict keyed by dataset name."""
    props = {}
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            props[row["dataset"]] = {k: v for k, v in row.items() if k != "dataset"}
    return props


def compute_correlations(gain_rows: list, property_dict: dict, out_path: str) -> list:
    """
    Compute Pearson and Spearman correlations between dataset properties and delta_ari.
    Returns list of correlation records.
    """
    datasets = [r["dataset"] for r in gain_rows]
    delta_ari = [r.get("delta_ari_vs_nomix") for r in gain_rows]

    # Which property columns to correlate with delta_ari
    prop_cols = [
        "knn_purity_k5", "knn_purity_k10", "knn_purity_k20",
        "cross_type_edge_ratio_k5", "cross_type_edge_ratio_k10", "cross_type_edge_ratio_k20",
        "silhouette_pca", "silhouette_hvg_pca",
        "zero_fraction",
        "class_imbalance_ratio", "label_entropy",
        "n_cells", "n_genes", "n_clusters",
        "median_detected_genes",
    ]

    correlation_rows = []
    for prop_col in prop_cols:
        prop_vals = []
        for ds in datasets:
            prop = property_dict.get(ds, {})
            raw_val = prop.get(prop_col, "")
            try:
                val = float(raw_val)
            except (ValueError, TypeError):
                val = None
            prop_vals.append(val)

        # Pearson & Spearman (only for datasets where both values are non-None)
        valid_pairs = [(d, da, pv) for d, da, pv in zip(datasets, delta_ari, prop_vals)
                       if da is not None and pv is not None and not np.isnan(pv)]
        if len(valid_pairs) < 3:
            continue

        da_arr = np.array([p[1] for p in valid_pairs])
        pv_arr = np.array([p[2] for p in valid_pairs])

        # Pearson
        if np.std(da_arr) > 1e-10 and np.std(pv_arr) > 1e-10:
            pearson_r = float(np.corrcoef(pv_arr, da_arr)[0, 1])
        else:
            pearson_r = None

        # Spearman (rank correlation)
        try:
            from scipy.stats import spearmanr, pearsonr as scipy_pearsonr
            spearman_r, spearman_p = spearmanr(pv_arr, da_arr)
            pearson_r2, pearson_p = scipy_pearsonr(pv_arr, da_arr)
            pearson_r = pearson_r2
        except Exception:
            spearman_r = None
            spearman_p = None
            pearson_p = None

        correlation_rows.append({
            "property": prop_col,
            "n_datasets": len(valid_pairs),
            "pearson_r": round(pearson_r, 4) if pearson_r is not None else None,
            "pearson_p": round(pearson_p, 4) if pearson_p is not None else None,
            "spearman_r": round(spearman_r, 4) if spearman_r is not None else None,
            "spearman_p": round(spearman_p, 4) if spearman_p is not None else None,
        })

    # Sort by absolute Spearman r
    correlation_rows.sort(key=lambda x: -abs(x.get("spearman_r") or 0))

    fieldnames = ["property", "n_datasets",
                  "pearson_r", "pearson_p",
                  "spearman_r", "spearman_p"]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(correlation_rows)
    print(f"Saved: {out_path}")
    return correlation_rows


# ─── Step 4: Scatter plots ─────────────────────────────────────────────────────

def plot_correlations(gain_rows: list, property_dict: dict, out_dir: str) -> list:
    """Generate delta_ari vs key property scatter plots."""
    if not HAS_MATPLOTLIB:
        print("matplotlib not available — skipping plots")
        return []

    datasets = [r["dataset"] for r in gain_rows]
    delta_ari = [r.get("delta_ari_vs_nomix") for r in gain_rows]

    plots = []
    plot_pairs = [
        ("knn_purity_k10", "KNN Purity (k=10)"),
        ("cross_type_edge_ratio_k10", "Cross-Type Edge Ratio (k=10)"),
        ("silhouette_pca", "Silhouette Score (PCA)"),
        ("zero_fraction", "Zero Fraction"),
        ("class_imbalance_ratio", "Class Imbalance Ratio"),
        ("label_entropy", "Label Entropy"),
        ("knn_purity_k20", "KNN Purity (k=20)"),
        ("cross_type_edge_ratio_k20", "Cross-Type Edge Ratio (k=20)"),
    ]

    for prop_col, prop_label in plot_pairs:
        prop_vals = []
        for ds in datasets:
            raw = property_dict.get(ds, {}).get(prop_col, "")
            try:
                val = float(raw)
            except (ValueError, TypeError):
                val = None
            prop_vals.append(val)

        valid = [(d, da, pv) for d, da, pv in zip(datasets, delta_ari, prop_vals)
                  if da is not None and pv is not None]
        if len(valid) < 3:
            continue

        fig, ax = plt.subplots(figsize=(6, 5))
        x = [p[2] for p in valid]
        y = [p[1] for p in valid]
        names = [p[0] for p in valid]

        # Color by effect
        name_to_effect = {r["dataset"]: r.get("effect_group", "") for r in gain_rows}
        effect_color = {"positive": "#2ca02c", "negative": "#d62728",
                        "neutral": "#7f7f7f", "insufficient_data": "#7f7f7f"}
        colors = [effect_color.get(name_to_effect.get(n, ""), "#1f77b4") for n in names]

        ax.scatter(x, y, c=colors, s=80, zorder=5)
        for xi, yi, ni in zip(x, y, names):
            ax.annotate(ni, (xi, yi), fontsize=8, xytext=(4, 4),
                       textcoords="offset points")

        ax.axhline(0.02, color="green", linestyle="--", alpha=0.5, label="+0.02 threshold")
        ax.axhline(-0.02, color="red", linestyle="--", alpha=0.5, label="-0.02 threshold")
        ax.axhline(0, color="gray", linestyle="-", alpha=0.3)

        ax.set_xlabel(prop_label, fontsize=10)
        ax.set_ylabel("ΔARI (NeighborMix vs noMix)", fontsize=10)
        ax.set_title(f"delta_ari vs {prop_label}", fontsize=11)
        ax.grid(True, alpha=0.3)

        safe_name = prop_col.replace("/", "_")
        out_path = os.path.join(out_dir, f"delta_ari_vs_{safe_name}.png")
        plt.tight_layout()
        plt.savefig(out_path, dpi=150)
        plt.close()
        plots.append(out_path)
        print(f"Saved: {out_path}")

    return plots


# ─── Step 5: Interpretation ─────────────────────────────────────────────────────

def generate_interpretation(
    gain_rows: list,
    correlation_rows: list,
    property_dict: dict,
    out_path: str,
):
    """
    Generate neighbormix_effect_interpretation.csv with per-dataset
    effect labels and main explanatory properties.
    """
    rows = []
    for r in gain_rows:
        ds = r["dataset"]
        props = property_dict.get(ds, {})
        delta_ari = r.get("delta_ari_vs_nomix")
        effect = r.get("effect_group", "insufficient_data")

        # Find top correlated properties
        sig_corrs = [c for c in correlation_rows if c.get("spearman_r") is not None
                     and abs(c["spearman_r"]) > 0.3]

        # Determine main explanation
        if effect == "positive":
            main_exp = "High neighborhood reliability (KNN purity) — NeighborMix finds same-type neighbors reliably"
        elif effect == "negative":
            main_exp = "High cross-type edge ratio — many nearest neighbors belong to other types, poisoning the mixing"
        elif effect == "neutral":
            main_exp = "Intermediate neighborhood structure — NeighborMix neither helps nor hurts"
        else:
            main_exp = "Insufficient data to determine"

        row = {
            "dataset": ds,
            "delta_ari": round(delta_ari, 4) if delta_ari is not None else None,
            "effect_group": effect,
            "main_explanation": main_exp,
            "knn_purity_k10": props.get("knn_purity_k10", ""),
            "cross_type_edge_ratio_k10": props.get("cross_type_edge_ratio_k10", ""),
            "silhouette_pca": props.get("silhouette_pca", ""),
            "zero_fraction": props.get("zero_fraction", ""),
            "class_imbalance_ratio": props.get("class_imbalance_ratio", ""),
            "label_entropy": props.get("label_entropy", ""),
            "n_cells": props.get("n_cells", ""),
            "n_genes": props.get("n_genes", ""),
            "n_clusters": props.get("n_clusters", ""),
        }
        rows.append(row)

    fieldnames = [
        "dataset", "delta_ari", "effect_group", "main_explanation",
        "knn_purity_k10", "cross_type_edge_ratio_k10", "silhouette_pca",
        "zero_fraction", "class_imbalance_ratio", "label_entropy",
        "n_cells", "n_genes", "n_clusters",
    ]
    with open(out_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {out_path}")


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Analyze NeighborMix effect vs dataset properties",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--property_csv", type=str, required=True,
                        help="Path to dataset_property_summary.csv")
    parser.add_argument("--benchmark_csv", type=str, default=None,
                        help="Path to benchmark_summary_combined.csv")
    parser.add_argument("--benchmark_dir", type=str, default=None,
                        help="Or scan this dir recursively for *_mean_std.csv")
    parser.add_argument("--out_dir", type=str, default="results/analysis",
                        help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # ── Load benchmark data ────────────────────────────────────
    if args.benchmark_csv:
        benchmark_data = load_benchmark_csv(args.benchmark_csv)
    elif args.benchmark_dir:
        benchmark_data = load_all_benchmarks(args.benchmark_dir)
    else:
        # Default: try combined at results/formal/
        benchmark_data = load_all_benchmarks("results/formal")
    print(f"Loaded {len(benchmark_data)} method/dataset entries from benchmark data")

    # ── Step 1: Gain summary ─────────────────────────────────
    gain_path = os.path.join(args.out_dir, "neighbormix_gain_summary.csv")
    gain_rows = compute_gain_summary(benchmark_data, gain_path)

    # ── Load properties ────────────────────────────────────────
    property_dict = load_properties(args.property_csv)
    print(f"Loaded {len(property_dict)} datasets from property CSV")

    # ── Step 2: Correlations ──────────────────────────────────
    corr_path = os.path.join(args.out_dir, "neighbormix_property_correlation.csv")
    correlation_rows = compute_correlations(gain_rows, property_dict, corr_path)

    # ── Step 3: Plots ─────────────────────────────────────────
    plot_paths = plot_correlations(gain_rows, property_dict, args.out_dir)

    # ── Step 4: Interpretation ────────────────────────────────
    interp_path = os.path.join(args.out_dir, "neighbormix_effect_interpretation.csv")
    generate_interpretation(gain_rows, correlation_rows, property_dict, interp_path)

    # ── Summary ───────────────────────────────────────────────
    print(f"\n{'=' * 70}")
    print("NeighborMix Effect Summary")
    print(f"{'=' * 70}")
    print(f"\n{'Dataset':<35} {'delta_ari':>10} {'effect':>12}")
    print("-" * 60)
    for r in gain_rows:
        da = r.get("delta_ari_vs_nomix")
        da_str = f"{da:+.4f}" if da is not None else "N/A"
        print(f"{r['dataset']:<35} {da_str:>10} {r['effect_group']:>12}")

    print(f"\nTop Spearman correlations with delta_ari:")
    for c in correlation_rows[:5]:
        if c["spearman_r"] is not None:
            print(f"  {c['property']:<35} r={c['spearman_r']:+.4f}  p={c['spearman_p']}")

    print(f"\nOutputs in {args.out_dir}:")
    for p in [gain_path, corr_path, interp_path] + plot_paths:
        print(f"  {p}")


if __name__ == "__main__":
    main()
