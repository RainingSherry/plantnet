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


SCRIPT_DIR = Path(__file__).resolve().parent
MANUSCRIPT_ROOT = SCRIPT_DIR.parent
PROJECT_ROOT = MANUSCRIPT_ROOT.parents[1]

PHASE13_ROOT = PROJECT_ROOT / "results/CAAM_scMAE_correction/corruption_triad/formal"
PHASE14_ROOT = PROJECT_ROOT / "results/CAAM_scMAE_correction/advmask_triage/formal"
ATTENTION_ROOT = Path("/tmp/caam_attention_context_smoke/dev_20260626")

METRICS = ("acc", "nmi", "ari", "f1_macro")
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


def flatten_metrics(path: Path) -> dict[str, float]:
    metrics = load_json(path)
    out = {}
    for group in ("kmeans_known_k", "leiden_fixed"):
        for metric in METRICS:
            out[f"{group}.{metric}"] = float(metrics[group][metric])
    return out


def parse_seed(token: str) -> int:
    if not token.startswith("seed"):
        raise ValueError(f"Expected seed token, got {token!r}")
    return int(token.removeprefix("seed"))


def discover_phase13(root: Path) -> list[dict]:
    rows = []
    for metrics_path in sorted(root.glob("*__*__seed*__epochs*/metrics.json")):
        run_dir = metrics_path.parent
        parts = run_dir.name.split("__")
        if len(parts) != 4:
            raise ValueError(f"Unexpected Phase 13 run directory name: {run_dir.name}")
        dataset, corruption, seed_token, epoch_token = parts
        corruption_stats = load_json(run_dir / "corruption_stats.json")
        row = {
            "dataset": dataset,
            "corruption": corruption,
            "seed": parse_seed(seed_token),
            "epochs": int(epoch_token.removeprefix("epochs")),
            "run_dir": str(run_dir),
            "zero_to_zero_rate": float(corruption_stats["zero_to_zero_rate"]),
            "effective_corruption_rate": float(corruption_stats["effective_corruption_rate"]),
            "mean_abs_delta": float(corruption_stats["mean_abs_delta"]),
            "budget_deficit_rate": float(corruption_stats["budget_deficit_rate"]),
        }
        row.update(flatten_metrics(metrics_path))
        rows.append(row)
    return rows


def discover_phase14(root: Path) -> list[dict]:
    rows = []
    for metrics_path in sorted(root.glob("*__*__*__seed*__epochs*/metrics.json")):
        run_dir = metrics_path.parent
        parts = run_dir.name.split("__")
        if len(parts) != 5:
            raise ValueError(f"Unexpected Phase 14 run directory name: {run_dir.name}")
        dataset, corruption, variant, seed_token, epoch_token = parts
        corruption_stats = load_json(run_dir / "corruption_stats.json")
        generator_stats_path = run_dir / "generator_stats.json"
        generator_stats = load_json(generator_stats_path) if generator_stats_path.exists() else {}
        row = {
            "dataset": dataset,
            "corruption": corruption,
            "variant": variant,
            "seed": parse_seed(seed_token),
            "epochs": int(epoch_token.removeprefix("epochs")),
            "run_dir": str(run_dir),
            "effective_corruption_rate": float(corruption_stats["effective_corruption_rate"]),
            "generator_grad_norm": float(generator_stats.get("generator_grad_norm", 0.0)),
            "mask_entropy": float(generator_stats.get("mask_entropy", "nan")),
            "mask_gini": float(generator_stats.get("mask_gini", "nan")),
        }
        row.update(flatten_metrics(metrics_path))
        rows.append(row)
    return rows


def discover_attention(root: Path) -> list[dict]:
    rows = []
    for metrics_path in sorted(root.glob("*__*__*__seed*__epochs*/metrics.json")):
        run_dir = metrics_path.parent
        parts = run_dir.name.split("__")
        if len(parts) != 5:
            raise ValueError(f"Unexpected attention smoke run directory name: {run_dir.name}")
        dataset, corruption, role, seed_token, epoch_token = parts
        run_manifest = load_json(run_dir / "run_manifest.json")
        row = {
            "dataset": dataset,
            "corruption": corruption,
            "role": role,
            "seed": parse_seed(seed_token),
            "epochs": int(epoch_token.removeprefix("epochs")),
            "run_dir": str(run_dir),
            "student_trainable_params": int(run_manifest["student_trainable_params"]),
        }
        row.update(flatten_metrics(metrics_path))
        rows.append(row)
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
            row[f"{value_key}.mean"] = mean(values)
            row[f"{value_key}.std"] = stdev(values) if len(values) > 1 else 0.0
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


def latex_escape(value: object) -> str:
    return str(value).replace("_", r"\_")


def fmt(value: float) -> str:
    if not math.isfinite(float(value)):
        return "--"
    return f"{value:.6f}"


def write_phase13_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Phase 13 corruption triad development summary.}",
        r"\label{tab:generated-phase13}",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Dataset & Corruption & ARI mean & ARI std & Effective corruption & Mean abs. delta \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row["corruption"]),
                    fmt(row["kmeans_known_k.ari.mean"]),
                    fmt(row["kmeans_known_k.ari.std"]),
                    fmt(row["effective_corruption_rate.mean"]),
                    fmt(row["mean_abs_delta.mean"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_phase13_nonoracle_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Phase 13 known-K versus fixed-Leiden summary.}",
        r"\label{tab:generated-phase13-nonoracle}",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Dataset & Corruption & Known-K ARI & Leiden ARI & Known-K NMI & Leiden NMI \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row["corruption"]),
                    fmt(row["kmeans_known_k.ari.mean"]),
                    fmt(row["leiden_fixed.ari.mean"]),
                    fmt(row["kmeans_known_k.nmi.mean"]),
                    fmt(row["leiden_fixed.nmi.mean"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_phase14_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Phase 14 AdvMask triage development summary.}",
        r"\label{tab:generated-phase14}",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Dataset & Variant & ARI mean & ARI std & Generator grad. & Mask entropy \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row["variant"]),
                    fmt(row["kmeans_known_k.ari.mean"]),
                    fmt(row["kmeans_known_k.ari.std"]),
                    fmt(row["generator_grad_norm.mean"]),
                    fmt(row["mask_entropy.mean"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_phase14_nonoracle_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Phase 14 known-K versus fixed-Leiden summary.}",
        r"\label{tab:generated-phase14-nonoracle}",
        r"\begin{tabular}{llrrrr}",
        r"\toprule",
        r"Dataset & Variant & Known-K ARI & Leiden ARI & Known-K NMI & Leiden NMI \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row["variant"]),
                    fmt(row["kmeans_known_k.ari.mean"]),
                    fmt(row["leiden_fixed.ari.mean"]),
                    fmt(row["kmeans_known_k.nmi.mean"]),
                    fmt(row["leiden_fixed.nmi.mean"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_attention_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Development-only attention/context smoke.}",
        r"\label{tab:generated-attention-smoke}",
        r"\begin{tabular}{llrrr}",
        r"\toprule",
        r"Dataset & Role & Student params & ARI & NMI \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row["role"]),
                    str(row["student_trainable_params"]),
                    fmt(row["kmeans_known_k.ari"]),
                    fmt(row["kmeans_known_k.nmi"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def write_attention_nonoracle_table(path: Path, rows: list[dict]) -> None:
    lines = [
        r"\begin{table}[ht]",
        r"\centering",
        r"\caption{Attention/context smoke known-K versus fixed-Leiden summary.}",
        r"\label{tab:generated-attention-smoke-nonoracle}",
        r"\begin{tabular}{llrrr}",
        r"\toprule",
        r"Dataset & Role & Student params & Known-K ARI & Leiden ARI \\",
        r"\midrule",
    ]
    for row in rows:
        lines.append(
            " & ".join(
                [
                    latex_escape(row["dataset"]),
                    latex_escape(row["role"]),
                    str(row["student_trainable_params"]),
                    fmt(row["kmeans_known_k.ari"]),
                    fmt(row["leiden_fixed.ari"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabular}", r"\end{table}", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def value_map(rows: list[dict], row_key: str, series_key: str, value_key: str) -> dict[str, dict[str, float]]:
    out = defaultdict(dict)
    for row in rows:
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
        ax.bar(offsets, values, width=width, label=DISPLAY_NAMES.get(series, series.replace("_", " ")), color=COLORS.get(series))
    ax.set_xticks(x)
    ax.set_xticklabels([label.replace("_", " ") for label in labels], rotation=15, ha="right")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=9)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_phase13_scatter(path: Path, rows: list[dict]) -> None:
    fig, ax = plt.subplots(figsize=(6.5, 4.8))
    for corruption in CORRUPTION_ORDER:
        subset = [row for row in rows if row["corruption"] == corruption]
        ax.scatter(
            [row["effective_corruption_rate.mean"] for row in subset],
            [row["kmeans_known_k.ari.mean"] for row in subset],
            label=corruption.replace("_", " "),
            s=70,
            color=COLORS[corruption],
        )
        for row in subset:
            ax.annotate(row["dataset"].replace("_", " "), (row["effective_corruption_rate.mean"], row["kmeans_known_k.ari.mean"]), fontsize=7, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Effective corruption rate")
    ax.set_ylabel("K-means known-K ARI")
    ax.set_title("Corruption diagnostics do not guarantee clustering gains")
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_phase14_delta(path: Path, rows: list[dict]) -> None:
    grouped = value_map(rows, "dataset", "variant", "kmeans_known_k.ari.mean")
    labels = sorted(grouped)
    deltas = [grouped[label]["advmask"] - grouped[label]["control"] for label in labels]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    colors = ["#54a24b" if value > 0 else "#e45756" for value in deltas]
    ax.bar([label.replace("_", " ") for label in labels], deltas, color=colors)
    ax.axhline(0.0, color="black", linewidth=1.0)
    ax.axhline(0.006518, color="#555555", linewidth=1.0, linestyle="--", label="0.5 x seed-std reference")
    ax.set_ylabel("AdvMask - control ARI")
    ax.set_title("AdvMask gains are inconsistent across development datasets")
    ax.tick_params(axis="x", rotation=15)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def plot_phase14_delta_metric(path: Path, rows: list[dict], metric: str, ylabel: str, title: str, reference: float | None = None) -> None:
    grouped = value_map(rows, "dataset", "variant", metric)
    labels = sorted(grouped)
    deltas = [grouped[label]["advmask"] - grouped[label]["control"] for label in labels]
    fig, ax = plt.subplots(figsize=(7.5, 4.4))
    colors = ["#54a24b" if value > 0 else "#e45756" for value in deltas]
    ax.bar([label.replace("_", " ") for label in labels], deltas, color=colors)
    ax.axhline(0.0, color="black", linewidth=1.0)
    if reference is not None:
        ax.axhline(reference, color="#555555", linewidth=1.0, linestyle="--", label="0.5 x seed-std reference")
        ax.legend(frameon=False, fontsize=8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=15)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def build_outputs(args: argparse.Namespace) -> None:
    phase13_rows = discover_phase13(args.phase13_root)
    phase14_rows = discover_phase14(args.phase14_root)
    attention_rows = discover_attention(args.attention_root)

    if len(phase13_rows) != 27:
        raise ValueError(f"Expected 27 Phase 13 runs, found {len(phase13_rows)}")
    if len(phase14_rows) != 18:
        raise ValueError(f"Expected 18 Phase 14 runs, found {len(phase14_rows)}")
    if len(attention_rows) != 9:
        raise ValueError(f"Expected 9 attention smoke runs, found {len(attention_rows)}")

    output_dir = args.output_dir
    data_dir = output_dir / "data"
    tables_dir = output_dir / "tables"
    figures_dir = output_dir / "figures"
    for path in (data_dir, tables_dir, figures_dir):
        path.mkdir(parents=True, exist_ok=True)

    phase13_summary = aggregate(
        phase13_rows,
        ("dataset", "corruption"),
        (
            "kmeans_known_k.ari",
            "kmeans_known_k.nmi",
            "kmeans_known_k.acc",
            "kmeans_known_k.f1_macro",
            "leiden_fixed.ari",
            "leiden_fixed.nmi",
            "leiden_fixed.acc",
            "leiden_fixed.f1_macro",
            "zero_to_zero_rate",
            "effective_corruption_rate",
            "mean_abs_delta",
        ),
    )
    phase14_summary = aggregate(
        phase14_rows,
        ("dataset", "corruption", "variant"),
        (
            "kmeans_known_k.ari",
            "kmeans_known_k.nmi",
            "kmeans_known_k.acc",
            "kmeans_known_k.f1_macro",
            "leiden_fixed.ari",
            "leiden_fixed.nmi",
            "leiden_fixed.acc",
            "leiden_fixed.f1_macro",
            "generator_grad_norm",
            "mask_entropy",
            "mask_gini",
            "effective_corruption_rate",
        ),
    )
    phase13_summary.sort(key=lambda row: (row["dataset"], CORRUPTION_ORDER.index(row["corruption"])))
    phase14_summary.sort(key=lambda row: (row["dataset"], PHASE14_VARIANT_ORDER.index(row["variant"])))
    attention_rows.sort(key=lambda row: (row["dataset"], ATTENTION_ROLE_ORDER.index(row["role"])))

    write_csv(data_dir / "phase13_corruption_runs.csv", phase13_rows)
    write_csv(data_dir / "phase13_corruption_summary.csv", phase13_summary)
    write_csv(data_dir / "phase14_advmask_runs.csv", phase14_rows)
    write_csv(data_dir / "phase14_advmask_summary.csv", phase14_summary)
    write_csv(data_dir / "attention_context_smoke_runs.csv", attention_rows)

    write_phase13_table(tables_dir / "phase13_corruption_summary.tex", phase13_summary)
    write_phase13_nonoracle_table(tables_dir / "phase13_corruption_nonoracle_summary.tex", phase13_summary)
    write_phase14_table(tables_dir / "phase14_advmask_summary.tex", phase14_summary)
    write_phase14_nonoracle_table(tables_dir / "phase14_advmask_nonoracle_summary.tex", phase14_summary)
    write_attention_table(tables_dir / "attention_context_smoke_summary.tex", attention_rows)
    write_attention_nonoracle_table(tables_dir / "attention_context_smoke_nonoracle_summary.tex", attention_rows)

    grouped_bar(
        figures_dir / "phase13_corruption_ari.png",
        value_map(phase13_summary, "dataset", "corruption", "kmeans_known_k.ari.mean"),
        CORRUPTION_ORDER,
        "K-means known-K ARI",
        "Phase 13 corruption triad",
    )
    grouped_bar(
        figures_dir / "phase13_corruption_leiden_ari.png",
        value_map(phase13_summary, "dataset", "corruption", "leiden_fixed.ari.mean"),
        CORRUPTION_ORDER,
        "Fixed-Leiden ARI",
        "Phase 13 corruption triad, non-oracle view",
    )
    plot_phase13_scatter(figures_dir / "phase13_effective_corruption_vs_ari.png", phase13_summary)
    plot_phase14_delta(figures_dir / "phase14_advmask_delta.png", phase14_summary)
    plot_phase14_delta_metric(
        figures_dir / "phase14_advmask_leiden_delta.png",
        phase14_summary,
        "leiden_fixed.ari.mean",
        "AdvMask - control fixed-Leiden ARI",
        "AdvMask fixed-Leiden gains are not consistently positive",
    )
    grouped_bar(
        figures_dir / "attention_context_smoke_ari.png",
        value_map(attention_rows, "dataset", "role", "kmeans_known_k.ari"),
        ATTENTION_ROLE_ORDER,
        "K-means known-K ARI",
        "Attention/context smoke",
    )
    grouped_bar(
        figures_dir / "attention_context_smoke_leiden_ari.png",
        value_map(attention_rows, "dataset", "role", "leiden_fixed.ari"),
        ATTENTION_ROLE_ORDER,
        "Fixed-Leiden ARI",
        "Attention/context smoke, non-oracle view",
    )

    manifest = {
        "phase13_root": str(args.phase13_root),
        "phase14_root": str(args.phase14_root),
        "attention_root": str(args.attention_root),
        "phase13_runs": len(phase13_rows),
        "phase14_runs": len(phase14_rows),
        "attention_runs": len(attention_rows),
        "claim_scope": "development evidence only; no validation or sealed test data",
        "known_k_metric": "kmeans_known_k",
        "non_oracle_metric": "leiden_fixed",
    }
    (output_dir / "artifact_manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Build protocol-analysis tables and figures from existing CAAM artifacts.")
    parser.add_argument("--phase13-root", type=Path, default=PHASE13_ROOT)
    parser.add_argument("--phase14-root", type=Path, default=PHASE14_ROOT)
    parser.add_argument("--attention-root", type=Path, default=ATTENTION_ROOT)
    parser.add_argument("--output-dir", type=Path, default=MANUSCRIPT_ROOT / "generated")
    args = parser.parse_args()
    build_outputs(args)
    print(f"Wrote generated assets to {args.output_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
