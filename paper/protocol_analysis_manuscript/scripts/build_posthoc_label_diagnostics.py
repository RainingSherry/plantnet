#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from statistics import mean, stdev

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = MANUSCRIPT_ROOT.parents[1]

PHASE13_ROOT = PROJECT_ROOT / "results/CAAM_scMAE_correction/corruption_triad/formal"
PHASE14_ROOT = PROJECT_ROOT / "results/CAAM_scMAE_correction/advmask_triage/formal"
ATTENTION_ROOT = Path("/tmp/caam_attention_context_smoke/dev_20260626")

EVALUATORS = ("kmeans_known_k", "leiden_fixed")
PREDICTION_FILES = {
    "kmeans_known_k": "eval_kmeans_known_k.npy",
    "leiden_fixed": "eval_leiden_fixed.npy",
}
CORRUPTION_ORDER = ("scmae_shuffle", "matched_donor", "nonzero_aware_donor")
PHASE14_VARIANT_ORDER = ("control", "advmask")
ATTENTION_ROLE_ORDER = ("control", "axial", "mlp_parammatched")
COLORS = {
    "scmae_shuffle": "#1b9e77",
    "matched_donor": "#d95f02",
    "nonzero_aware_donor": "#7570b3",
    "control": "#4c78a8",
    "advmask": "#f58518",
    "axial": "#e45756",
    "mlp_parammatched": "#54a24b",
}
DISPLAY_NAMES = {
    "scmae_shuffle": "scMAE shuffle",
    "matched_donor": "matched donor",
    "nonzero_aware_donor": "nonzero-aware donor",
    "control": "control",
    "advmask": "AdvMask",
    "axial": "Axial",
    "mlp_parammatched": "parameter-matched MLP",
}


def load_json(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_seed(token: str) -> int:
    if not token.startswith("seed"):
        raise ValueError(f"Expected seed token, got {token!r}")
    return int(token.removeprefix("seed"))


def flatten_metrics(path: Path, evaluator: str) -> dict[str, float]:
    metrics = load_json(path)
    group = metrics[evaluator]
    return {
        "ari": float(group["ari"]),
        "nmi": float(group["nmi"]),
        "acc": float(group["acc"]),
        "f1_macro": float(group["f1_macro"]),
    }


def discover_phase13(root: Path) -> list[dict]:
    rows = []
    for run_dir in sorted(path.parent for path in root.glob("*__*__seed*__epochs*/metrics.json")):
        parts = run_dir.name.split("__")
        if len(parts) != 4:
            raise ValueError(f"Unexpected Phase 13 run directory name: {run_dir.name}")
        dataset, corruption, seed_token, epoch_token = parts
        rows.append(
            {
                "phase": "phase13_corruption",
                "dataset": dataset,
                "corruption": corruption,
                "condition": corruption,
                "seed": parse_seed(seed_token),
                "epochs": int(epoch_token.removeprefix("epochs")),
                "run_dir": str(run_dir),
            }
        )
    return rows


def discover_phase14(root: Path) -> list[dict]:
    rows = []
    for run_dir in sorted(path.parent for path in root.glob("*__*__*__seed*__epochs*/metrics.json")):
        parts = run_dir.name.split("__")
        if len(parts) != 5:
            raise ValueError(f"Unexpected Phase 14 run directory name: {run_dir.name}")
        dataset, corruption, variant, seed_token, epoch_token = parts
        rows.append(
            {
                "phase": "phase14_advmask",
                "dataset": dataset,
                "corruption": corruption,
                "condition": variant,
                "variant": variant,
                "seed": parse_seed(seed_token),
                "epochs": int(epoch_token.removeprefix("epochs")),
                "run_dir": str(run_dir),
            }
        )
    return rows


def discover_attention(root: Path) -> list[dict]:
    rows = []
    for run_dir in sorted(path.parent for path in root.glob("*__*__*__seed*__epochs*/metrics.json")):
        parts = run_dir.name.split("__")
        if len(parts) != 5:
            raise ValueError(f"Unexpected attention smoke run directory name: {run_dir.name}")
        dataset, corruption, role, seed_token, epoch_token = parts
        rows.append(
            {
                "phase": "attention_context_smoke",
                "dataset": dataset,
                "corruption": corruption,
                "condition": role,
                "role": role,
                "seed": parse_seed(seed_token),
                "epochs": int(epoch_token.removeprefix("epochs")),
                "run_dir": str(run_dir),
            }
        )
    return rows


def harmonic_mean(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.0
    return 2.0 * a * b / (a + b)


def mean_or_nan(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return mean(finite) if finite else float("nan")


def std_or_zero(values: list[float]) -> float:
    finite = [value for value in values if math.isfinite(value)]
    return stdev(finite) if len(finite) > 1 else 0.0


def label_rows_for_run(run: dict, evaluator: str, rare_fraction: float, rare_min_count: int) -> tuple[list[dict], dict]:
    run_dir = Path(run["run_dir"])
    labels = np.load(run_dir / "labels.npy", allow_pickle=True)
    predictions = np.load(run_dir / PREDICTION_FILES[evaluator], allow_pickle=True)
    if labels.shape != predictions.shape:
        raise ValueError(f"Shape mismatch in {run_dir}: labels {labels.shape}, predictions {predictions.shape}")

    n_cells = int(labels.shape[0])
    rare_threshold = max(int(math.ceil(n_cells * rare_fraction)), int(rare_min_count))
    pred_unique = np.unique(predictions)
    pred_sizes = {cluster: int(np.sum(predictions == cluster)) for cluster in pred_unique}

    rows = []
    for label in np.unique(labels):
        label_mask = labels == label
        label_count = int(np.sum(label_mask))
        label_predictions = predictions[label_mask]
        clusters, counts = np.unique(label_predictions, return_counts=True)
        best_index = int(np.argmax(counts))
        best_cluster = clusters[best_index]
        best_count = int(counts[best_index])
        dominant_recall = best_count / label_count
        cluster_purity = best_count / pred_sizes[best_cluster]
        dominant_f1 = harmonic_mean(dominant_recall, cluster_purity)
        rows.append(
            {
                **run,
                "evaluator": evaluator,
                "n_cells": n_cells,
                "n_clusters": int(len(pred_unique)),
                "label_id": str(label),
                "label_count": label_count,
                "label_fraction": label_count / n_cells,
                "rare_threshold_count": rare_threshold,
                "is_rare_label": label_count <= rare_threshold,
                "best_cluster_id": str(best_cluster),
                "best_cluster_size": pred_sizes[best_cluster],
                "dominant_cluster_recall": dominant_recall,
                "dominant_cluster_purity": cluster_purity,
                "dominant_cluster_f1": dominant_f1,
                "fragmented_cluster_count_10pct": int(np.sum((counts / label_count) >= 0.10)),
            }
        )

    rare_rows = [row for row in rows if row["is_rare_label"]]
    metrics = flatten_metrics(run_dir / "metrics.json", evaluator)
    summary = {
        **run,
        "evaluator": evaluator,
        "n_cells": n_cells,
        "n_labels": int(len(np.unique(labels))),
        "n_clusters": int(len(pred_unique)),
        "rare_threshold_count": rare_threshold,
        "n_rare_labels": len(rare_rows),
        "macro_label_recall": mean(row["dominant_cluster_recall"] for row in rows),
        "macro_label_purity": mean(row["dominant_cluster_purity"] for row in rows),
        "macro_label_f1": mean(row["dominant_cluster_f1"] for row in rows),
        "weighted_label_recall": sum(row["dominant_cluster_recall"] * row["label_count"] for row in rows) / n_cells,
        "rare_label_recall": mean_or_nan([row["dominant_cluster_recall"] for row in rare_rows]),
        "rare_label_purity": mean_or_nan([row["dominant_cluster_purity"] for row in rare_rows]),
        "rare_label_f1": mean_or_nan([row["dominant_cluster_f1"] for row in rare_rows]),
        "worst_label_f1": min(row["dominant_cluster_f1"] for row in rows),
        "ari": metrics["ari"],
        "nmi": metrics["nmi"],
        "acc": metrics["acc"],
        "f1_macro": metrics["f1_macro"],
    }
    return rows, summary


def label_cluster_rows_for_run(run: dict, evaluator: str) -> list[dict]:
    run_dir = Path(run["run_dir"])
    labels = np.load(run_dir / "labels.npy", allow_pickle=True)
    predictions = np.load(run_dir / PREDICTION_FILES[evaluator], allow_pickle=True)
    if labels.shape != predictions.shape:
        raise ValueError(f"Shape mismatch in {run_dir}: labels {labels.shape}, predictions {predictions.shape}")

    rows = []
    label_sizes = {label: int(np.sum(labels == label)) for label in np.unique(labels)}
    cluster_sizes = {cluster: int(np.sum(predictions == cluster)) for cluster in np.unique(predictions)}
    for label in np.unique(labels):
        label_mask = labels == label
        clusters, counts = np.unique(predictions[label_mask], return_counts=True)
        for cluster, count in zip(clusters, counts):
            count = int(count)
            rows.append(
                {
                    **run,
                    "evaluator": evaluator,
                    "label_id": str(label),
                    "cluster_id": str(cluster),
                    "count": count,
                    "label_count": label_sizes[label],
                    "cluster_count": cluster_sizes[cluster],
                    "fraction_of_label": count / label_sizes[label],
                    "fraction_of_cluster": count / cluster_sizes[cluster],
                }
            )
    return rows


def build_label_diagnostics(runs: list[dict], rare_fraction: float, rare_min_count: int) -> tuple[list[dict], list[dict]]:
    label_rows = []
    summaries = []
    for run in runs:
        for evaluator in EVALUATORS:
            run_label_rows, summary = label_rows_for_run(run, evaluator, rare_fraction, rare_min_count)
            label_rows.extend(run_label_rows)
            summaries.append(summary)
    return label_rows, summaries


def build_label_cluster_flows(runs: list[dict]) -> list[dict]:
    rows = []
    for run in runs:
        for evaluator in EVALUATORS:
            rows.extend(label_cluster_rows_for_run(run, evaluator))
    return rows


def aggregate(rows: list[dict], keys: tuple[str, ...], value_keys: tuple[str, ...]) -> list[dict]:
    grouped = defaultdict(list)
    for row in rows:
        grouped[tuple(row[key] for key in keys)].append(row)

    out = []
    for key_values, group_rows in sorted(grouped.items()):
        row = {key: value for key, value in zip(keys, key_values)}
        row["n_runs"] = len(group_rows)
        for value_key in value_keys:
            values = [float(item[value_key]) for item in group_rows]
            row[f"{value_key}.mean"] = mean_or_nan(values)
            row[f"{value_key}.std"] = std_or_zero(values)
        out.append(row)
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"No rows to write for {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def fmt(value: float) -> str:
    value = float(value)
    if not math.isfinite(value):
        return "--"
    return f"{value:.6f}"


def latex_escape(value: object) -> str:
    return str(value).replace("_", r"\_")


def write_summary_table(path: Path, caption: str, label: str, rows: list[dict], condition_key: str) -> None:
    table_rows = [row for row in rows if row["evaluator"] == "kmeans_known_k"]
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        rf"\caption{{{caption}}}",
        rf"\label{{{label}}}",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Dataset & Condition & Label F1 & Rare-label F1 & Worst-label F1 & ARI \\",
        r"\midrule",
    ]
    for row in table_rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row[condition_key]),
                    fmt(row["macro_label_f1.mean"]),
                    fmt(row["rare_label_f1.mean"]),
                    fmt(row["worst_label_f1.mean"]),
                    fmt(row["ari.mean"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def value_map(rows: list[dict], row_key: str, series_key: str, value_key: str) -> dict[str, dict[str, float]]:
    out = defaultdict(dict)
    for row in rows:
        if row["evaluator"] == "kmeans_known_k":
            out[str(row[row_key])][str(row[series_key])] = float(row[value_key])
    return out


def grouped_bar(path: Path, data: dict[str, dict[str, float]], series_order: tuple[str, ...], ylabel: str, title: str) -> None:
    labels = sorted(data)
    x = list(range(len(labels)))
    width = 0.8 / len(series_order)
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    for i, series in enumerate(series_order):
        values = [data[label].get(series, float("nan")) for label in labels]
        offsets = [pos - 0.4 + width / 2 + i * width for pos in x]
        ax.bar(offsets, values, width=width, label=DISPLAY_NAMES.get(series, series), color=COLORS.get(series))
    ax.set_xticks(x)
    ax.set_xticklabels([label.replace("_", " ") for label in labels], rotation=15, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def delta_bar(path: Path, rows: list[dict], positive_condition: str, control_condition: str, ylabel: str, title: str) -> None:
    grouped = value_map(rows, "dataset", "condition", "macro_label_f1.mean")
    labels = sorted(grouped)
    deltas = [grouped[label][positive_condition] - grouped[label][control_condition] for label in labels]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    colors = ["#54a24b" if value > 0 else "#e45756" for value in deltas]
    ax.bar([label.replace("_", " ") for label in labels], deltas, color=colors)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=15)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_run_label_cluster_heatmap(ax: plt.Axes, run: dict, evaluator: str):
    run_dir = Path(run["run_dir"])
    labels = np.load(run_dir / "labels.npy", allow_pickle=True)
    predictions = np.load(run_dir / PREDICTION_FILES[evaluator], allow_pickle=True)
    label_values = np.unique(labels)
    cluster_values = np.unique(predictions)
    matrix = np.zeros((len(label_values), len(cluster_values)), dtype=float)
    for label_index, label in enumerate(label_values):
        label_mask = labels == label
        label_total = float(np.sum(label_mask))
        for cluster_index, cluster in enumerate(cluster_values):
            matrix[label_index, cluster_index] = np.sum(predictions[label_mask] == cluster) / label_total

    image = ax.imshow(matrix, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    ax.set_title(f"{run['dataset'].replace('_', ' ')}\n{DISPLAY_NAMES.get(run['condition'], run['condition'])}", fontsize=8)
    ax.set_xlabel("cluster ID", fontsize=7)
    ax.set_ylabel("label ID", fontsize=7)
    if len(cluster_values) <= 15:
        ax.set_xticks(range(len(cluster_values)))
        ax.set_xticklabels([str(value) for value in cluster_values], fontsize=6, rotation=90)
    else:
        ax.set_xticks([])
    if len(label_values) <= 15:
        ax.set_yticks(range(len(label_values)))
        ax.set_yticklabels([str(value) for value in label_values], fontsize=6)
    else:
        ax.set_yticks([])
    return image


def plot_phase_heatmaps(path: Path, runs: list[dict], condition_order: tuple[str, ...], title: str, evaluator: str = "kmeans_known_k", seed: int = 42) -> None:
    datasets = sorted({run["dataset"] for run in runs})
    selected = {
        (run["dataset"], run["condition"]): run
        for run in runs
        if run["seed"] == seed and run["condition"] in condition_order
    }
    fig, axes = plt.subplots(
        len(datasets),
        len(condition_order),
        figsize=(3.1 * len(condition_order), 2.6 * len(datasets)),
        constrained_layout=True,
    )
    if len(datasets) == 1:
        axes = np.array([axes])
    if len(condition_order) == 1:
        axes = axes.reshape(len(datasets), 1)

    image = None
    for row_index, dataset in enumerate(datasets):
        for col_index, condition in enumerate(condition_order):
            ax = axes[row_index, col_index]
            run = selected.get((dataset, condition))
            if run is None:
                ax.axis("off")
                continue
            image = plot_run_label_cluster_heatmap(ax, run, evaluator)
    fig.suptitle(title, fontsize=11)
    if image is not None:
        fig.colorbar(image, ax=axes.ravel().tolist(), shrink=0.78, label="fraction of label")
    fig.savefig(path, dpi=220)
    plt.close(fig)


def sort_phase_rows(rows: list[dict], condition_order: tuple[str, ...]) -> None:
    rows.sort(key=lambda row: (row["dataset"], row["evaluator"], condition_order.index(row["condition"])))


def build_outputs(args: argparse.Namespace) -> None:
    phase13_runs = discover_phase13(args.phase13_root)
    phase14_runs = discover_phase14(args.phase14_root)
    attention_runs = discover_attention(args.attention_root)

    if len(phase13_runs) != 27:
        raise ValueError(f"Expected 27 Phase 13 runs, found {len(phase13_runs)}")
    if len(phase14_runs) != 18:
        raise ValueError(f"Expected 18 Phase 14 runs, found {len(phase14_runs)}")
    if len(attention_runs) != 9:
        raise ValueError(f"Expected 9 attention smoke runs, found {len(attention_runs)}")

    output_dir = args.output_dir
    data_dir = output_dir / "data"
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    for path in (data_dir, tables_dir, figures_dir):
        path.mkdir(parents=True, exist_ok=True)

    phase13_labels, phase13_summaries = build_label_diagnostics(phase13_runs, args.rare_fraction, args.rare_min_count)
    phase14_labels, phase14_summaries = build_label_diagnostics(phase14_runs, args.rare_fraction, args.rare_min_count)
    attention_labels, attention_summaries = build_label_diagnostics(attention_runs, args.rare_fraction, args.rare_min_count)
    phase13_flows = build_label_cluster_flows(phase13_runs)
    phase14_flows = build_label_cluster_flows(phase14_runs)
    attention_flows = build_label_cluster_flows(attention_runs)

    phase13_summary = aggregate(
        phase13_summaries,
        ("dataset", "corruption", "condition", "evaluator"),
        ("macro_label_f1", "rare_label_f1", "worst_label_f1", "ari", "nmi"),
    )
    phase14_summary = aggregate(
        phase14_summaries,
        ("dataset", "variant", "condition", "evaluator"),
        ("macro_label_f1", "rare_label_f1", "worst_label_f1", "ari", "nmi"),
    )
    attention_summary = aggregate(
        attention_summaries,
        ("dataset", "role", "condition", "evaluator"),
        ("macro_label_f1", "rare_label_f1", "worst_label_f1", "ari", "nmi"),
    )
    sort_phase_rows(phase13_summary, CORRUPTION_ORDER)
    sort_phase_rows(phase14_summary, PHASE14_VARIANT_ORDER)
    sort_phase_rows(attention_summary, ATTENTION_ROLE_ORDER)

    write_csv(data_dir / "phase13_posthoc_label_rows.csv", phase13_labels)
    write_csv(data_dir / "phase13_posthoc_run_summary.csv", phase13_summaries)
    write_csv(data_dir / "phase13_posthoc_condition_summary.csv", phase13_summary)
    write_csv(data_dir / "phase13_label_cluster_flow.csv", phase13_flows)
    write_csv(data_dir / "phase14_posthoc_label_rows.csv", phase14_labels)
    write_csv(data_dir / "phase14_posthoc_run_summary.csv", phase14_summaries)
    write_csv(data_dir / "phase14_posthoc_condition_summary.csv", phase14_summary)
    write_csv(data_dir / "phase14_label_cluster_flow.csv", phase14_flows)
    write_csv(data_dir / "attention_posthoc_label_rows.csv", attention_labels)
    write_csv(data_dir / "attention_posthoc_run_summary.csv", attention_summaries)
    write_csv(data_dir / "attention_posthoc_condition_summary.csv", attention_summary)
    write_csv(data_dir / "attention_label_cluster_flow.csv", attention_flows)

    write_summary_table(
        tables_dir / "phase13_posthoc_label_summary.tex",
        "Post-hoc label diagnostics for Phase 13 corruption variants. Numeric labels are used only after training.",
        "tab:generated-phase13-posthoc-labels",
        phase13_summary,
        "condition",
    )
    write_summary_table(
        tables_dir / "phase14_posthoc_label_summary.tex",
        "Post-hoc label diagnostics for Phase 14 AdvMask triage. Numeric labels are used only after training.",
        "tab:generated-phase14-posthoc-labels",
        phase14_summary,
        "condition",
    )
    write_summary_table(
        tables_dir / "attention_posthoc_label_summary.tex",
        "Post-hoc label diagnostics for the attention/context smoke. Numeric labels are used only after training.",
        "tab:generated-attention-posthoc-labels",
        attention_summary,
        "condition",
    )

    grouped_bar(
        figures_dir / "phase13_posthoc_label_f1.png",
        value_map(phase13_summary, "dataset", "condition", "macro_label_f1.mean"),
        CORRUPTION_ORDER,
        "Post-hoc macro label F1",
        "Phase 13 label recovery by corruption",
    )
    delta_bar(
        figures_dir / "phase14_posthoc_label_f1_delta.png",
        phase14_summary,
        "advmask",
        "control",
        "AdvMask - control label F1",
        "AdvMask does not consistently improve label recovery",
    )
    grouped_bar(
        figures_dir / "attention_posthoc_label_f1.png",
        value_map(attention_summary, "dataset", "condition", "macro_label_f1.mean"),
        ATTENTION_ROLE_ORDER,
        "Post-hoc macro label F1",
        "Attention/context smoke label recovery",
    )
    plot_phase_heatmaps(
        figures_dir / "phase13_label_cluster_heatmaps_seed42.png",
        phase13_runs,
        CORRUPTION_ORDER,
        "Phase 13 label-to-cluster flow, seed 42",
    )
    plot_phase_heatmaps(
        figures_dir / "phase14_label_cluster_heatmaps_seed42.png",
        phase14_runs,
        PHASE14_VARIANT_ORDER,
        "Phase 14 label-to-cluster flow, seed 42",
    )
    plot_phase_heatmaps(
        figures_dir / "attention_label_cluster_heatmaps_seed42.png",
        attention_runs,
        ATTENTION_ROLE_ORDER,
        "Attention/context label-to-cluster flow, seed 42",
    )

    manifest = {
        "phase13_root": str(args.phase13_root),
        "phase14_root": str(args.phase14_root),
        "attention_root": str(args.attention_root),
        "phase13_runs": len(phase13_runs),
        "phase14_runs": len(phase14_runs),
        "attention_runs": len(attention_runs),
        "evaluators": list(EVALUATORS),
        "rare_label_rule": {
            "rare_fraction": args.rare_fraction,
            "rare_min_count": args.rare_min_count,
            "threshold": "label_count <= max(ceil(n_cells * rare_fraction), rare_min_count)",
        },
        "label_cluster_flow": {
            "csv_rows": {
                "phase13": len(phase13_flows),
                "phase14": len(phase14_flows),
                "attention": len(attention_flows),
            },
            "figure_evaluator": "kmeans_known_k",
            "figure_seed": 42,
            "interpretation": "row-normalized label-to-cluster heatmaps; cluster IDs are arbitrary within each run",
        },
        "claim_scope": "post-hoc development diagnostics only; labels are never used for training, model selection, corruption, masks, context, early stopping, or validation",
        "label_semantics": "numeric label IDs from existing artifacts; no marker-gene or cell-type-name claim is made here",
    }
    (output_dir / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build post-hoc label diagnostics from existing CAAM development artifacts.")
    parser.add_argument("--phase13-root", type=Path, default=PHASE13_ROOT)
    parser.add_argument("--phase14-root", type=Path, default=PHASE14_ROOT)
    parser.add_argument("--attention-root", type=Path, default=ATTENTION_ROOT)
    parser.add_argument("--output-dir", type=Path, default=MANUSCRIPT_ROOT / "generated/posthoc_label_diagnostics")
    parser.add_argument("--rare-fraction", type=float, default=0.01)
    parser.add_argument("--rare-min-count", type=int, default=50)
    args = parser.parse_args()
    build_outputs(args)
    print(f"Wrote post-hoc label diagnostics to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
