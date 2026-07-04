#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
GENERATED_ROOT = MANUSCRIPT_ROOT / "generated"
DEFAULT_RUN_ROOT = Path("/tmp/caam_mask_ratio_smoke/dev_20260626_gpu")
MASK_RATIOS = (0.2, 0.4, 0.6)
SEEDS = (42, 2024, 3407)


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True), encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float) -> str:
    return f"{float(value):.6f}"


def run_name(mask_ratio: float, seed: int) -> str:
    suffix = str(mask_ratio).replace(".", "p")
    return f"Quake_Smart-seq2_Lung__mask{suffix}__control__seed{int(seed)}__epochs3"


def flatten_run(run_root: Path, mask_ratio: float, seed: int) -> dict:
    run_dir = run_root / run_name(mask_ratio, seed)
    required = (
        "args.json",
        "metrics.json",
        "corruption_stats.json",
        "mask_stats.json",
        "training_history.json",
        "dataset_profile.json",
        "run_manifest.json",
        "artifact_manifest.json",
    )
    missing = [name for name in required if not (run_dir / name).exists()]
    if missing:
        raise FileNotFoundError(f"Missing required mask-ratio smoke files for {mask_ratio}: {missing} in {run_dir}")

    args = load_json(run_dir / "args.json")
    metrics = load_json(run_dir / "metrics.json")
    corruption = load_json(run_dir / "corruption_stats.json")
    mask = load_json(run_dir / "mask_stats.json")
    history = load_json(run_dir / "training_history.json")
    profile = load_json(run_dir / "dataset_profile.json")
    run_manifest = load_json(run_dir / "run_manifest.json")
    artifact = load_json(run_dir / "artifact_manifest.json")

    selected_mask_ratio = float(mask["selected_mask_ratio"])
    effective_masked = float(corruption["effective_corruption_rate"])
    return {
        "mask_ratio": float(args["mask_ratio"]),
        "dataset": str(artifact["dataset"]),
        "seed": int(artifact["seed"]),
        "epochs": int(args["epochs"]),
        "n_cells": int(profile["n_cells"]),
        "n_genes": int(profile["n_genes"]),
        "student_trainable_params": int(run_manifest["student_trainable_params"]),
        "selected_mask_ratio": selected_mask_ratio,
        "effective_corruption_rate_masked": effective_masked,
        "global_effective_change_rate_estimate": selected_mask_ratio * effective_masked,
        "zero_to_zero_rate_masked": float(corruption["zero_to_zero_rate"]),
        "mean_abs_delta": float(corruption["mean_abs_delta"]),
        "budget_deficit_rate_final": float(history["budget_deficit_rate"][-1]),
        "kmeans_known_k_ari": float(metrics["kmeans_known_k"]["ari"]),
        "kmeans_known_k_nmi": float(metrics["kmeans_known_k"]["nmi"]),
        "kmeans_known_k_f1_macro": float(metrics["kmeans_known_k"]["f1_macro"]),
        "leiden_fixed_ari": float(metrics["leiden_fixed"]["ari"]),
        "leiden_fixed_nmi": float(metrics["leiden_fixed"]["nmi"]),
        "leiden_fixed_f1_macro": float(metrics["leiden_fixed"]["f1_macro"]),
        "run_dir": str(run_dir),
    }


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def sample_std(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    mu = mean(values)
    return (sum((value - mu) ** 2 for value in values) / (len(values) - 1)) ** 0.5


def aggregate_by_ratio(rows: list[dict]) -> list[dict]:
    aggregate_rows = []
    for ratio in sorted({row["mask_ratio"] for row in rows}):
        group = [row for row in rows if row["mask_ratio"] == ratio]
        aggregate_rows.append(
            {
                "mask_ratio": ratio,
                "n_runs": len(group),
                "seeds": ",".join(str(row["seed"]) for row in sorted(group, key=lambda item: item["seed"])),
                "effective_corruption_rate_masked_mean": mean([row["effective_corruption_rate_masked"] for row in group]),
                "effective_corruption_rate_masked_std": sample_std([row["effective_corruption_rate_masked"] for row in group]),
                "global_effective_change_rate_estimate_mean": mean([row["global_effective_change_rate_estimate"] for row in group]),
                "global_effective_change_rate_estimate_std": sample_std(
                    [row["global_effective_change_rate_estimate"] for row in group]
                ),
                "kmeans_known_k_ari_mean": mean([row["kmeans_known_k_ari"] for row in group]),
                "kmeans_known_k_ari_std": sample_std([row["kmeans_known_k_ari"] for row in group]),
                "leiden_fixed_ari_mean": mean([row["leiden_fixed_ari"] for row in group]),
                "leiden_fixed_ari_std": sample_std([row["leiden_fixed_ari"] for row in group]),
                "leiden_fixed_f1_macro_mean": mean([row["leiden_fixed_f1_macro"] for row in group]),
                "leiden_fixed_f1_macro_std": sample_std([row["leiden_fixed_f1_macro"] for row in group]),
            }
        )
    return aggregate_rows


def summarize(rows: list[dict]) -> dict:
    aggregate_rows = aggregate_by_ratio(rows)
    best_known = max(aggregate_rows, key=lambda row: row["kmeans_known_k_ari_mean"])
    best_leiden = max(aggregate_rows, key=lambda row: row["leiden_fixed_ari_mean"])
    baseline = next(row for row in aggregate_rows if abs(row["mask_ratio"] - 0.4) < 1e-9)
    return {
        "best_known_k_mask_ratio": best_known["mask_ratio"],
        "best_known_k_ari_mean": best_known["kmeans_known_k_ari_mean"],
        "best_leiden_mask_ratio": best_leiden["mask_ratio"],
        "best_leiden_ari_mean": best_leiden["leiden_fixed_ari_mean"],
        "baseline_mask_ratio": baseline["mask_ratio"],
        "known_k_mean_range": max(row["kmeans_known_k_ari_mean"] for row in aggregate_rows)
        - min(row["kmeans_known_k_ari_mean"] for row in aggregate_rows),
        "leiden_mean_range": max(row["leiden_fixed_ari_mean"] for row in aggregate_rows)
        - min(row["leiden_fixed_ari_mean"] for row in aggregate_rows),
        "effective_masked_mean_range": max(row["effective_corruption_rate_masked_mean"] for row in aggregate_rows)
        - min(row["effective_corruption_rate_masked_mean"] for row in aggregate_rows),
        "global_effective_change_mean_range": max(row["global_effective_change_rate_estimate_mean"] for row in aggregate_rows)
        - min(row["global_effective_change_rate_estimate_mean"] for row in aggregate_rows),
        "max_known_k_seed_std": max(row["kmeans_known_k_ari_std"] for row in aggregate_rows),
        "max_leiden_seed_std": max(row["leiden_fixed_ari_std"] for row in aggregate_rows),
        "interpretation": (
            "single-dataset development sensitivity; supports reporting mask-ratio protocol sensitivity, "
            "not retuning the frozen mask ratio"
        ),
    }


def write_latex(path: Path, aggregate_rows: list[dict], summary: dict) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Mask-ratio development sensitivity on Quake\_Smart-seq2\_Lung across three seeds. This is a single-dataset diagnostic and must not be used to tune validation.}",
        r"\label{tab:mask-ratio-smoke}",
        r"\scriptsize",
        r"\begin{tabular}{rrrrrr}",
        r"\toprule",
        r"Mask ratio & masked effective & global effective & known-\(K\) ARI & fixed-Leiden ARI & fixed-Leiden F1 \\",
        r"\midrule",
    ]
    for row in aggregate_rows:
        lines.append(
            f"{row['mask_ratio']:.1f} & {fmt(row['effective_corruption_rate_masked_mean'])} & "
            f"{fmt(row['global_effective_change_rate_estimate_mean'])} & {fmt(row['kmeans_known_k_ari_mean'])} & "
            f"{fmt(row['leiden_fixed_ari_mean'])} & {fmt(row['leiden_fixed_f1_macro_mean'])} \\\\"
        )
    lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            "",
            rf"\footnotesize{{Means over seeds 42, 2024, and 3407. Known-\(K\) ARI mean range: {fmt(summary['known_k_mean_range'])}; fixed-Leiden ARI mean range: {fmt(summary['leiden_mean_range'])}.}}",
            r"\end{table}",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def write_figure(path: Path, aggregate_rows: list[dict]) -> None:
    ratios = [row["mask_ratio"] for row in aggregate_rows]
    known = [row["kmeans_known_k_ari_mean"] for row in aggregate_rows]
    known_std = [row["kmeans_known_k_ari_std"] for row in aggregate_rows]
    leiden = [row["leiden_fixed_ari_mean"] for row in aggregate_rows]
    leiden_std = [row["leiden_fixed_ari_std"] for row in aggregate_rows]
    global_effective = [row["global_effective_change_rate_estimate_mean"] for row in aggregate_rows]

    fig, ax_metric = plt.subplots(figsize=(6.4, 4.2))
    ax_effective = ax_metric.twinx()

    ax_metric.errorbar(
        ratios,
        known,
        yerr=known_std,
        marker="o",
        linewidth=2.0,
        capsize=4,
        color="#1f77b4",
        label="known-K ARI",
    )
    ax_metric.errorbar(
        ratios,
        leiden,
        yerr=leiden_std,
        marker="s",
        linewidth=2.0,
        capsize=4,
        color="#2ca02c",
        label="fixed-Leiden ARI",
    )
    ax_effective.plot(
        ratios,
        global_effective,
        marker="^",
        linewidth=2.0,
        color="#d62728",
        label="global effective change",
    )

    ax_metric.set_xlabel("Nominal mask ratio")
    ax_metric.set_ylabel("ARI mean across seeds")
    ax_effective.set_ylabel("Estimated global effective change")
    ax_metric.set_xticks(ratios)
    ax_metric.set_ylim(0.34, 0.54)
    ax_effective.set_ylim(0.0, max(global_effective) * 1.25)
    ax_metric.grid(axis="y", alpha=0.25)
    handles1, labels1 = ax_metric.get_legend_handles_labels()
    handles2, labels2 = ax_effective.get_legend_handles_labels()
    ax_metric.legend(handles1 + handles2, labels1 + labels2, loc="lower right", frameon=False)
    ax_metric.set_title("Mask-ratio sensitivity smoke")
    fig.tight_layout()
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def write_markdown(path: Path, aggregate_rows: list[dict], summary: dict) -> None:
    lines = [
        "# Mask-ratio Smoke Summary",
        "",
        "Status: generated development-only mask-ratio sensitivity across three seeds. This is not validation evidence.",
        "",
        "| mask ratio | runs | masked effective mean | global effective mean | known-K ARI mean | known-K ARI std | fixed-Leiden ARI mean | fixed-Leiden ARI std |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            f"| {row['mask_ratio']:.1f} | {row['n_runs']} | {fmt(row['effective_corruption_rate_masked_mean'])} | "
            f"{fmt(row['global_effective_change_rate_estimate_mean'])} | {fmt(row['kmeans_known_k_ari_mean'])} | "
            f"{fmt(row['kmeans_known_k_ari_std'])} | {fmt(row['leiden_fixed_ari_mean'])} | {fmt(row['leiden_fixed_ari_std'])} |"
        )
    lines.extend(
        [
            "",
            "## Sensitivity Notes",
            "",
            f"- Best mean known-K ARI in this smoke: mask ratio {summary['best_known_k_mask_ratio']:.1f} ({fmt(summary['best_known_k_ari_mean'])}).",
            f"- Best mean fixed-Leiden ARI in this smoke: mask ratio {summary['best_leiden_mask_ratio']:.1f} ({fmt(summary['best_leiden_ari_mean'])}).",
            f"- known-K ARI mean range: {fmt(summary['known_k_mean_range'])}.",
            f"- fixed-Leiden ARI mean range: {fmt(summary['leiden_mean_range'])}.",
            f"- maximum known-K seed standard deviation: {fmt(summary['max_known_k_seed_std'])}.",
            f"- maximum fixed-Leiden seed standard deviation: {fmt(summary['max_leiden_seed_std'])}.",
            f"- Masked-position effective corruption mean range: {fmt(summary['effective_masked_mean_range'])}.",
            f"- Global effective change estimate mean range: {fmt(summary['global_effective_change_mean_range'])}.",
            "",
            "## Claim Boundary",
            "",
            "This supports reporting protocol sensitivity around nominal mask ratio. It must not be used to tune the frozen validation mask ratio or to claim that a different ratio is generally superior.",
            "",
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines), encoding="utf-8")


def build_outputs(args: argparse.Namespace) -> dict:
    rows = [flatten_run(args.run_root, ratio, seed) for ratio in MASK_RATIOS for seed in SEEDS]
    aggregate_rows = aggregate_by_ratio(rows)
    summary = summarize(rows)
    payload = {
        "status": "development_only",
        "run_root": str(args.run_root),
        "n_runs": len(rows),
        "seeds": list(SEEDS),
        "rows": rows,
        "aggregate_rows": aggregate_rows,
        "summary": summary,
        "claim_boundary": "single-dataset mask-ratio sensitivity; not validation; not a tuning decision",
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "mask_ratio_smoke.csv", rows)
    write_json(args.output_dir / "mask_ratio_smoke.json", payload)
    write_markdown(args.output_dir / "mask_ratio_smoke.md", aggregate_rows, summary)
    write_latex(args.output_dir / "mask_ratio_smoke.tex", aggregate_rows, summary)
    write_figure(args.output_dir / "mask_ratio_smoke_sensitivity.png", aggregate_rows)
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize the CAAM/scMAE mask-ratio development smoke.")
    parser.add_argument("--run-root", type=Path, default=DEFAULT_RUN_ROOT)
    parser.add_argument("--output-dir", type=Path, default=GENERATED_ROOT / "mask_ratio_smoke")
    args = parser.parse_args()
    payload = build_outputs(args)
    print(f"mask_ratio_smoke_runs={payload['n_runs']}")
    print(f"mask_ratio_smoke_output={args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
