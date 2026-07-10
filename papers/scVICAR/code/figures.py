from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
import numpy as np
import pandas as pd
import seaborn as sns

from .config import DATASETS, PAPER_ROOT, SEEDS, VARIANTS


COLORS = {
    "nomix": "#606060",
    "random_mix": "#D8D8D8",
    "fixed": "#7884B4",
    "topology_edge_only": "#B4C0E4",
    "topology_gate_only": "#E4CCD8",
    "topology_full": "#B64342",
    "anchor": "#0F4D92",
    "neighbor": "#AADCA9",
    "masked": "#F0E0D0",
    "affinity": "#9A4D8E",
}

DISPLAY = {
    "nomix": "NoMix",
    "random_mix": "Random",
    "fixed": "scVICAR-F",
    "topology_edge_only": "T: edge only",
    "topology_gate_only": "T: gate only",
    "topology_full": "scVICAR-T",
}


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans", "Liberation Sans"],
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "font.size": 7,
            "axes.spines.right": False,
            "axes.spines.top": False,
            "axes.linewidth": 0.8,
            "legend.frameon": False,
        }
    )


def export(fig: plt.Figure, target: Path) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    for suffix, kwargs in {
        ".svg": {},
        ".pdf": {},
        ".tiff": {"dpi": 600},
        ".png": {"dpi": 300},
    }.items():
        fig.savefig(target.with_suffix(suffix), bbox_inches="tight", facecolor="white", **kwargs)
    plt.close(fig)


def box(ax, xy, width, height, text, color="#FFFFFF", edge="#4D4D4D", size=7, radius=0.02) -> None:
    patch = FancyBboxPatch(
        xy, width, height,
        boxstyle=f"round,pad=0.012,rounding_size={radius}",
        facecolor=color, edgecolor=edge, linewidth=0.9,
    )
    ax.add_patch(patch)
    ax.text(xy[0] + width / 2, xy[1] + height / 2, text, ha="center", va="center", fontsize=size)


def arrow(ax, start, end, color="#4D4D4D", width=1.0, style="-|>") -> None:
    ax.add_patch(FancyArrowPatch(start, end, arrowstyle=style, mutation_scale=8, linewidth=width, color=color))


def cell(ax, xy, label, color, radius=0.036, edge="#4D4D4D") -> None:
    ax.add_patch(Circle(xy, radius, facecolor=color, edgecolor=edge, linewidth=0.8))
    ax.text(*xy, label, ha="center", va="center", fontsize=6, color="#272727")


def plot_method(target: Path) -> None:
    style()
    fig = plt.figure(figsize=(7.2, 4.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[3.1, 1.0], hspace=0.22, wspace=0.16)
    axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    for ax in axes:
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

    # a: original scMAE corruption and anchor target.
    ax = axes[0]
    ax.text(0, 1.04, "a", fontweight="bold", fontsize=8, transform=ax.transAxes)
    ax.text(0.5, 0.96, "Masked anchor reconstruction", ha="center", fontweight="bold")
    cell(ax, (0.09, 0.64), r"$x_i$", COLORS["anchor"])
    box(ax, (0.22, 0.55), 0.25, 0.18, "gene-value\nreplacement", COLORS["masked"], size=6)
    arrow(ax, (0.13, 0.64), (0.22, 0.64))
    box(ax, (0.57, 0.55), 0.25, 0.18, "masked\nautoencoder", "#E4E4F0", size=6)
    arrow(ax, (0.47, 0.64), (0.57, 0.64))
    cell(ax, (0.93, 0.64), r"$x_i$", COLORS["anchor"])
    arrow(ax, (0.82, 0.64), (0.89, 0.64))
    ax.text(0.5, 0.32, "The original cell is both input anchor\nand reconstruction target.", ha="center", color="#4D4D4D")

    # b: fixed graph-vicinal view.
    ax = axes[1]
    ax.text(0, 1.04, "b", fontweight="bold", fontsize=8, transform=ax.transAxes)
    ax.text(0.5, 0.96, "scVICAR-F: fixed local budget", ha="center", fontweight="bold")
    cell(ax, (0.22, 0.62), r"$x_i$", COLORS["anchor"])
    for pos, lab in [((0.46, 0.78), r"$x_j$"), ((0.52, 0.57), r"$x_k$"), ((0.42, 0.38), r"$x_l$")]:
        cell(ax, pos, lab, COLORS["neighbor"], radius=0.03)
        arrow(ax, pos, (0.275, 0.61), COLORS["neighbor"], 0.8, "-")
    box(ax, (0.64, 0.52), 0.27, 0.20, r"$x_i^v=(1-\beta)x_i$" + "\n" + r"$+\,\beta\sum_j\hat p_{ij}x_j$", "#EEF3FA")
    arrow(ax, (0.28, 0.62), (0.64, 0.62), COLORS["fixed"], 1.2)
    ax.text(0.5, 0.25, r"Fixed $\beta=0.1$; target remains $x_i$.", ha="center", color="#4D4D4D")

    # c: topology-adaptive graph.
    ax = axes[2]
    ax.text(0, 1.04, "c", fontweight="bold", fontsize=8, transform=ax.transAxes)
    ax.text(0.5, 0.96, "scVICAR-T: topology-adaptive budget", ha="center", fontweight="bold")
    cell(ax, (0.17, 0.63), r"$x_i$", COLORS["anchor"])
    neighbor_positions = [(0.42, 0.80), (0.48, 0.60), (0.40, 0.38)]
    widths = [1.8, 1.1, 0.45]
    for pos, lw in zip(neighbor_positions, widths):
        cell(ax, pos, r"$x_j$", COLORS["neighbor"], radius=0.03)
        arrow(ax, pos, (0.215, 0.63), COLORS["affinity"], lw, "-")
    box(ax, (0.61, 0.67), 0.30, 0.15, "topology-informed\naffinity $a_{ij}$", "#F4EAF2")
    box(ax, (0.61, 0.40), 0.30, 0.15, "cell-wise budget\n$0\\leq g_i\\leq g_{max}$", "#F6CFCB")
    arrow(ax, (0.50, 0.65), (0.61, 0.74), COLORS["affinity"], 1.0)
    arrow(ax, (0.50, 0.56), (0.61, 0.48), COLORS["topology_full"], 1.0)
    ax.text(0.5, 0.22, "Affinity selects local evidence; the gate\nlimits its cell-specific displacement.", ha="center", color="#4D4D4D")

    ax_eq = fig.add_subplot(gs[1, :])
    ax_eq.axis("off")
    ax_eq.text(0, 0.92, "d", fontweight="bold", fontsize=8, transform=ax_eq.transAxes)
    box(ax_eq, (0.04, 0.18), 0.26, 0.56, r"$X^v=(I-G)X+G\widehat{P}X$", "#F7F7F7", size=8)
    arrow(ax_eq, (0.31, 0.46), (0.38, 0.46))
    box(ax_eq, (0.39, 0.18), 0.24, 0.56, r"$X^v=X-G(I-\widehat{P})X$" + "\n\none graph-diffusion step", "#EEF3FA", size=8)
    arrow(ax_eq, (0.64, 0.46), (0.71, 0.46))
    box(ax_eq, (0.72, 0.18), 0.24, 0.56, r"$\mathcal{L}=\mathcal{L}_{clean}+\lambda\mathcal{L}_{anchor}$" + "\n" + r"target: $x_i$, not $x_i^v$", "#FBEDEC", size=8)
    export(fig, target)


def require_formal_matrix(frame: pd.DataFrame) -> pd.DataFrame:
    formal = frame[frame["execution_mode"] == "formal"].copy()
    expected = len(DATASETS) * len(VARIANTS) * len(SEEDS)
    keys = formal[["dataset", "variant", "seed"]].drop_duplicates()
    if len(keys) != expected:
        raise RuntimeError(f"Refusing partial confirmatory figure: found {len(keys)}/{expected} formal tasks")
    return formal


def plot_confirmatory(run_master: Path, contrast_csv: Path, target: Path) -> None:
    style()
    formal = require_formal_matrix(pd.read_csv(run_master))
    seed_mean = formal.groupby(["dataset", "variant"], as_index=False)[["ari", "nmi", "f1_macro", "runtime_seconds"]].mean()
    pivot = seed_mean.pivot(index="dataset", columns="variant", values="ari")
    contrasts = pd.read_csv(contrast_csv)
    fig = plt.figure(figsize=(7.2, 5.2))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.45, 1.0, 1.0], hspace=0.38, wspace=0.42)
    ax_a = fig.add_subplot(gs[:, 0])
    y = np.arange(len(pivot.index))
    for offset, variant in [(-0.12, "fixed"), (0.12, "topology_full")]:
        delta = pivot[variant] - pivot["nomix"]
        ax_a.scatter(delta, y + offset, s=24, color=COLORS[variant], label=DISPLAY[variant], zorder=3)
    ax_a.axvline(0, color="#767676", linestyle="--", linewidth=0.8)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels([name.replace("_", " ") for name in pivot.index])
    ax_a.set_xlabel(r"Paired $\Delta$ARI vs NoMix")
    ax_a.legend(loc="lower right")
    ax_a.text(-0.20, 1.02, "a", transform=ax_a.transAxes, fontweight="bold", fontsize=8)

    ax_b = fig.add_subplot(gs[0, 1:])
    order = list(VARIANTS)
    sns.boxplot(data=seed_mean, x="variant", y="ari", order=order, color="white", width=0.55, fliersize=0, ax=ax_b)
    sns.stripplot(data=seed_mean, x="variant", y="ari", hue="variant", order=order, hue_order=order,
                  palette=[COLORS[v] for v in order], size=3.5, jitter=0.16, legend=False, ax=ax_b)
    ax_b.set_xticks(range(len(order)), [DISPLAY[v] for v in order], rotation=25, ha="right")
    ax_b.set_xlabel("")
    ax_b.set_ylabel("ARI")
    ax_b.text(-0.10, 1.04, "b", transform=ax_b.transAxes, fontweight="bold", fontsize=8)

    ax_c = fig.add_subplot(gs[1, 1])
    show = contrasts.set_index("contrast").loc[["F_vs_NoMix", "T_vs_NoMix", "T_vs_F"]]
    ypos = np.arange(len(show))[::-1]
    ax_c.hlines(ypos, show["ci95_low"], show["ci95_high"], color="#606060", linewidth=1.2)
    ax_c.scatter(show["mean_delta"], ypos, color=[COLORS["fixed"], COLORS["topology_full"], COLORS["topology_full"]], s=24)
    ax_c.axvline(0, color="#767676", linestyle="--", linewidth=0.8)
    ax_c.set_yticks(ypos)
    ax_c.set_yticklabels(["F−NoMix", "T−NoMix", "T−F"])
    ax_c.set_xlabel("Mean paired ΔARI (95% CI)")
    ax_c.text(-0.22, 1.04, "c", transform=ax_c.transAxes, fontweight="bold", fontsize=8)

    ax_d = fig.add_subplot(gs[1, 2])
    runtime = formal.groupby("variant")["runtime_seconds"].median().reindex(order)
    ax_d.barh(np.arange(len(order)), runtime, color=[COLORS[v] for v in order])
    ax_d.set_yticks(np.arange(len(order)))
    ax_d.set_yticklabels([DISPLAY[v] for v in order])
    ax_d.set_xlabel("Median runtime (s)")
    ax_d.text(-0.22, 1.04, "d", transform=ax_d.transAxes, fontweight="bold", fontsize=8)
    export(fig, target)


def plot_components(run_master: Path, target: Path) -> None:
    style()
    formal = require_formal_matrix(pd.read_csv(run_master))
    seed_mean = formal.groupby(["dataset", "variant"], as_index=False)[
        ["ari", "geometry_between_within_ratio"]
    ].mean()
    ari = seed_mean.pivot(index="dataset", columns="variant", values="ari")
    variants = ["random_mix", "fixed", "topology_edge_only", "topology_gate_only", "topology_full"]
    delta = ari[variants].subtract(ari["nomix"], axis=0).reset_index().melt(
        id_vars="dataset", var_name="variant", value_name="delta_ari"
    )
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.55), gridspec_kw={"width_ratios": [1.65, 1.0, 1.2]})
    sns.boxplot(data=delta, x="variant", y="delta_ari", order=variants, color="white", fliersize=0, ax=axes[0])
    sns.stripplot(
        data=delta, x="variant", y="delta_ari", hue="variant", order=variants, hue_order=variants,
        palette=[COLORS[v] for v in variants], size=3.5, jitter=0.15, legend=False, ax=axes[0],
    )
    axes[0].axhline(0, color="#777777", linestyle="--", linewidth=0.8)
    axes[0].set_xticks(range(len(variants)), [DISPLAY[v] for v in variants], rotation=30, ha="right")
    axes[0].set(xlabel="", ylabel="Paired ΔARI vs NoMix")
    axes[0].text(-0.12, 1.04, "a", transform=axes[0].transAxes, fontweight="bold")

    interaction = np.array([
        [ari["fixed"].mean(), ari["topology_edge_only"].mean()],
        [ari["topology_gate_only"].mean(), ari["topology_full"].mean()],
    ])
    sns.heatmap(
        interaction, annot=True, fmt=".3f", cmap="vlag", center=float(ari["fixed"].mean()),
        cbar_kws={"label": "Mean ARI"}, xticklabels=["uniform\nweights", "topology\nweights"],
        yticklabels=["fixed gate", "adaptive gate"], ax=axes[1],
    )
    axes[1].set_xlabel("Edge affinity")
    axes[1].set_ylabel("")
    axes[1].tick_params(axis="x", labelrotation=0, labelsize=6)
    axes[1].text(-0.25, 1.04, "b", transform=axes[1].transAxes, fontweight="bold")

    geometry = seed_mean[seed_mean["variant"].isin(["nomix", "fixed", "topology_full"])]
    sns.boxplot(
        data=geometry, x="variant", y="geometry_between_within_ratio",
        order=["nomix", "fixed", "topology_full"], color="white", fliersize=0, ax=axes[2],
    )
    sns.stripplot(
        data=geometry, x="variant", y="geometry_between_within_ratio", hue="variant",
        order=["nomix", "fixed", "topology_full"], hue_order=["nomix", "fixed", "topology_full"],
        palette=[COLORS[v] for v in ["nomix", "fixed", "topology_full"]],
        size=3.5, jitter=0.12, legend=False, ax=axes[2],
    )
    axes[2].set_xticks(range(3), ["NoMix", "F", "T"])
    axes[2].set(xlabel="", ylabel="Between/within-class distance")
    axes[2].text(-0.18, 1.04, "c", transform=axes[2].transAxes, fontweight="bold")
    fig.tight_layout(w_pad=1.1)
    export(fig, target)


def plot_stress(stress_runs: Path, target: Path) -> None:
    style()
    frame = pd.read_csv(stress_runs)
    if len(frame) != 126:
        raise RuntimeError(f"Refusing partial stress figure: {len(frame)}/126 runs")
    current = frame[frame["estimator"] == "current"].copy()
    expected = 3 * 2 * 3 * 5
    if len(current[["dataset", "variant", "seed", "contamination"]].drop_duplicates()) != expected:
        raise RuntimeError("Stress contamination grid is incomplete")
    summary = current.groupby(["dataset", "variant", "contamination"], as_index=False).mean(numeric_only=True)
    fig = plt.figure(figsize=(7.2, 4.3))
    gs = fig.add_gridspec(2, 3, width_ratios=[1.35, 1.0, 1.0], hspace=0.42, wspace=0.42)
    ax_a = fig.add_subplot(gs[:, 0])
    for variant in ["fixed", "topology_full"]:
        subset = summary[summary["variant"] == variant]
        for dataset, group in subset.groupby("dataset"):
            ax_a.plot(group["contamination"] * 100, group["ari"], color=COLORS[variant], alpha=0.28, linewidth=0.8)
        mean = subset.groupby("contamination")["ari"].mean()
        ax_a.plot(mean.index * 100, mean.values, color=COLORS[variant], linewidth=2, marker="o", ms=3, label=DISPLAY[variant])
    ax_a.set(xlabel="Injected cross-class edges (%)", ylabel="ARI")
    ax_a.legend()
    ax_a.text(-0.16, 1.02, "a", transform=ax_a.transAxes, fontweight="bold")

    ax_b = fig.add_subplot(gs[0, 1:])
    pivot = summary.pivot(index=["dataset", "contamination"], columns="variant", values="ari_degradation").reset_index()
    pivot["T_minus_F_degradation"] = pivot["topology_full"] - pivot["fixed"]
    sns.lineplot(
        data=pivot, x="contamination", y="T_minus_F_degradation", hue="dataset",
        marker="o", linewidth=1, ax=ax_b,
    )
    ax_b.axhline(0, color="#777777", linestyle="--", linewidth=0.8)
    ax_b.set(xlabel="Contamination fraction", ylabel="T−F degradation (higher is safer)")
    ax_b.legend(fontsize=5, ncol=3, title=None)
    ax_b.text(-0.10, 1.05, "b", transform=ax_b.transAxes, fontweight="bold")

    ax_c = fig.add_subplot(gs[1, 1])
    diag = current[current["variant"] == "topology_full"]
    sns.boxplot(data=diag, x="contamination", y="affinity_same_edge_auroc", color=COLORS["topology_full"], ax=ax_c)
    ax_c.set(xlabel="Contamination", ylabel="Affinity AUROC")
    ax_c.text(-0.24, 1.05, "c", transform=ax_c.transAxes, fontweight="bold")

    ax_d = fig.add_subplot(gs[1, 2])
    sns.scatterplot(
        data=diag, x="weighted_same_edge_purity", y="gate_purity_spearman",
        hue="contamination", palette="viridis", s=20, ax=ax_d,
    )
    ax_d.set(xlabel="Weighted same-edge purity", ylabel="Gate–purity Spearman")
    ax_d.legend(fontsize=5, title="contam.")
    ax_d.text(-0.22, 1.05, "d", transform=ax_d.transAxes, fontweight="bold")
    export(fig, target)


def plot_downstream(dataset_metrics: Path, target: Path) -> None:
    style()
    frame = pd.read_csv(dataset_metrics)
    expected = len(DATASETS) * len(VARIANTS) * 3
    if len(frame) != expected:
        raise RuntimeError(f"Refusing partial downstream figure: {len(frame)}/{expected} rows")
    variants = ["nomix", "fixed", "topology_full"]
    marker = frame[(frame["task"] == "marker") & frame["variant"].isin(variants)]
    probe = frame[(frame["task"] == "linear_probe") & frame["variant"].isin(variants)]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.6))
    for ax, metric, label, panel in [
        (axes[0], "recovery_recovery_at_100", "Recovery@100", "a"),
        (axes[1], "annotation_f1_macro", "Marker annotation macro-F1", "b"),
    ]:
        sns.boxplot(data=marker, x="variant", y=metric, order=variants, color="white", fliersize=0, ax=ax)
        sns.stripplot(
            data=marker, x="variant", y=metric, hue="variant", order=variants, hue_order=variants,
            palette=[COLORS[v] for v in variants], size=3.5, jitter=0.13, legend=False, ax=ax,
        )
        ax.set_xticks(range(len(variants)), [DISPLAY[v] for v in variants], rotation=20, ha="right")
        ax.set(xlabel="", ylabel=label)
        ax.text(-0.18, 1.04, panel, transform=ax.transAxes, fontweight="bold")
    sns.pointplot(
        data=probe, x="label_fraction", y="probe_f1_macro", hue="variant",
        hue_order=variants, palette=[COLORS[v] for v in variants],
        errorbar=None, dodge=0.16, markers=["o", "s", "D"], ax=axes[2],
    )
    axes[2].set(xlabel="Labeled-cell fraction", ylabel="Frozen-probe macro-F1")
    axes[2].legend(labels=[DISPLAY[v] for v in variants], fontsize=6, title=None)
    axes[2].text(-0.18, 1.04, "c", transform=axes[2].transAxes, fontweight="bold")
    fig.tight_layout(w_pad=1.2)
    export(fig, target)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("figure", choices=["method", "confirmatory", "components", "stress", "downstream"])
    parser.add_argument("--run-master", type=Path, default=PAPER_ROOT / "experiments" / "protocol_v1" / "run_master.csv")
    parser.add_argument("--contrasts", type=Path, default=PAPER_ROOT / "tables" / "confirmatory_contrasts_ari.csv")
    parser.add_argument("--output-dir", type=Path, default=PAPER_ROOT / "figures" / "final")
    parser.add_argument("--stress-runs", type=Path, default=PAPER_ROOT / "experiments/stress_v1/stress_runs.csv")
    parser.add_argument(
        "--downstream-metrics",
        type=Path,
        default=PAPER_ROOT / "experiments/downstream_v1/dataset_variant_metrics.csv",
    )
    args = parser.parse_args()
    if args.figure == "method":
        plot_method(args.output_dir / "fig1_method")
    elif args.figure == "confirmatory":
        plot_confirmatory(args.run_master, args.contrasts, args.output_dir / "fig2_confirmatory")
    elif args.figure == "components":
        plot_components(args.run_master, args.output_dir / "fig3_components")
    elif args.figure == "stress":
        plot_stress(args.stress_runs, args.output_dir / "fig4_graph_stress")
    else:
        plot_downstream(args.downstream_metrics, args.output_dir / "fig5_downstream")


if __name__ == "__main__":
    main()
