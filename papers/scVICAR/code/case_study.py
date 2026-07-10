from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns

from .config import DATASETS, PAPER_ROOT, PROTOCOL_VERSION, REMOTE_RESULT_ROOT, SEEDS, VARIANTS
from .figures import COLORS, DISPLAY, export, style
from .io_utils import sha256_file, verify_checksum_manifest
from .remote_store import RemoteStore


CASE_DATASET = "Human_Pancreas_3"
CASE_VARIANTS = ("nomix", "fixed", "topology_full")
CASE_MODEL_SEED = 42
CASE_SPLIT_SEED = 11


def require_case_runs(run_master: Path) -> dict[str, dict]:
    frame = pd.read_csv(run_master)
    formal = frame[frame["execution_mode"] == "formal"]
    expected = len(DATASETS) * len(VARIANTS) * len(SEEDS)
    if len(formal) != expected:
        raise RuntimeError(f"Case study requires the complete primary matrix: {len(formal)}/{expected}")
    selected = formal[
        (formal["dataset"] == CASE_DATASET)
        & (formal["seed"] == CASE_MODEL_SEED)
        & (formal["variant"].isin(CASE_VARIANTS))
    ]
    if len(selected) != len(CASE_VARIANTS):
        raise RuntimeError("The preregistered Human_Pancreas_3 case runs are incomplete")
    return {row["variant"]: row.to_dict() for _, row in selected.iterrows()}


def fetch_case_sources(run_master: Path, target: Path) -> None:
    runs = require_case_runs(run_master)
    store = RemoteStore()
    for variant, row in runs.items():
        variant_dir = target / variant
        downstream_dir = variant_dir / "downstream"
        remote = f"{REMOTE_RESULT_ROOT}/downstream/{PROTOCOL_VERSION}/{row['run_id']}"
        if not (downstream_dir / "COMPLETED").is_file():
            store.download_directory(remote, downstream_dir)
        verify_checksum_manifest(downstream_dir)
        split = downstream_dir / f"split_{CASE_SPLIT_SEED}"
        for name in ("marker_overlap.csv", "marker_results.json", "split_indices.npz"):
            if not (split / name).is_file():
                raise FileNotFoundError(split / name)
        run_dir = variant_dir / "primary"
        run_dir.mkdir(parents=True, exist_ok=True)
        for name in ("SHA256SUMS", "clusters.npz"):
            store.download_file(f"{row['remote_dir']}/{name}", run_dir / name)
        checksums = {}
        for line in (run_dir / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
            if line.strip():
                digest, name = line.split("  ", 1)
                checksums[name] = digest
        if sha256_file(run_dir / "clusters.npz") != checksums["clusters.npz"]:
            raise ValueError(f"Primary cluster checksum failed for {variant}")
    target.mkdir(parents=True, exist_ok=True)
    (target / "case_protocol.json").write_text(
        json.dumps(
            {
                "dataset": CASE_DATASET,
                "variants": list(CASE_VARIANTS),
                "model_seed": CASE_MODEL_SEED,
                "split_seed": CASE_SPLIT_SEED,
                "selection": "fixed before result inspection; no best-seed selection",
                "run_ids": {variant: row["run_id"] for variant, row in runs.items()},
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
        encoding="utf-8",
    )


def load_marker_panel(path: Path) -> dict[str, list[str]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("verification_status") != "verified":
        raise RuntimeError(
            "The fixed pancreas marker panel is not yet literature-verified; "
            "Figure 6 generation is deliberately blocked"
        )
    return {str(key): [str(gene) for gene in genes] for key, genes in payload["markers"].items()}


def dotplot_data(data_path: Path, panel: dict[str, list[str]]) -> pd.DataFrame:
    adata = ad.read_h5ad(data_path)
    labels = adata.obs["resolved_label"].astype(str).to_numpy()
    matrix = adata.layers["counts"] if "counts" in adata.layers else adata.X
    names = np.asarray(adata.var_names, dtype=str)
    lookup = {name.upper(): index for index, name in enumerate(names)}
    genes = list(dict.fromkeys(gene for values in panel.values() for gene in values))
    present = [gene for gene in genes if gene.upper() in lookup]
    indices = [lookup[gene.upper()] for gene in present]
    totals = np.asarray(matrix.sum(axis=1)).ravel().astype(float)
    selected = matrix[:, indices]
    selected = selected.toarray() if sp.issparse(selected) else np.asarray(selected)
    normalized = np.log1p(selected * np.divide(1e4, totals, out=np.zeros_like(totals), where=totals > 0)[:, None])
    rows = []
    for cell_type in panel:
        mask = labels == cell_type
        for column, gene in enumerate(present):
            values = normalized[mask, column]
            rows.append(
                {
                    "cell_type": cell_type,
                    "gene": gene,
                    "mean_log_expression": float(values.mean()) if values.size else 0.0,
                    "fraction_expressed": float((values > 0).mean()) if values.size else 0.0,
                }
            )
    return pd.DataFrame(rows)


def flow_table(case_root: Path, variant: str, data_path: Path) -> pd.DataFrame:
    split_dir = case_root / variant / "downstream" / f"split_{CASE_SPLIT_SEED}"
    split = np.load(split_dir / "split_indices.npz")
    evaluation = split["evaluation"].astype(int)
    clusters = np.load(case_root / variant / "primary/clusters.npz")["predicted"].astype(str)
    labels = ad.read_h5ad(data_path, backed="r")
    gold = labels.obs["resolved_label"].astype(str).to_numpy()[evaluation]
    labels.file.close()
    result = json.loads((split_dir / "marker_results.json").read_text(encoding="utf-8"))
    mapping = {str(key): str(value) for key, value in result["cluster_annotation"].items()}
    annotation = np.asarray([mapping.get(cluster, "unassigned") for cluster in clusters[evaluation]])
    return pd.DataFrame({"gold": gold, "cluster": clusters[evaluation], "annotation": annotation})


def parallel_flow(ax: plt.Axes, flow: pd.DataFrame, title: str) -> None:
    columns = ["gold", "cluster", "annotation", "gold"]
    x = np.asarray([0.0, 1.0, 2.25, 3.55])
    display = {
        "activated_stellate": "act. stellate",
        "quiescent_stellate": "qui. stellate",
        "unassigned": "unassigned",
    }
    positions = {}
    for index, column in enumerate(columns):
        values = sorted(flow[column].astype(str).unique())
        positions[index] = {value: pos for pos, value in enumerate(values)}
        ax.scatter(np.full(len(values), x[index]), range(len(values)), s=8, color="#555555", zorder=3)
        if index in (0, 2, 3):
            side = "right" if index in (0, 2) else "left"
            offset = -0.055 if index in (0, 2) else 0.055
            for value, pos in positions[index].items():
                ax.text(
                    x[index] + offset, pos, display.get(value, value),
                    fontsize=3.8, ha=side, va="center",
                )
    for index in range(3):
        left, right = columns[index], columns[index + 1]
        counts = flow.groupby([left, right]).size()
        maximum = max(1, int(counts.max()))
        for (source, destination), count in counts.items():
            ax.plot(
                [x[index], x[index + 1]],
                [positions[index][str(source)], positions[index + 1][str(destination)]],
                color=COLORS["topology_full"] if title.endswith("T") else "#6E7FA8",
                alpha=0.10 + 0.35 * count / maximum,
                linewidth=0.25 + 2.0 * count / maximum,
            )
    ax.set_xlim(-0.30, 3.88)
    ax.set_xticks(x)
    ax.set_xticklabels(["gold", "cluster", "marker call", "gold"], fontsize=5)
    ax.set_yticks([])
    ax.set_title(title, fontsize=7)
    for spine in ax.spines.values():
        spine.set_visible(False)


def plot_case(case_root: Path, data_path: Path, marker_panel: Path, target: Path) -> None:
    panel = load_marker_panel(marker_panel)
    dots = dotplot_data(data_path, panel)
    dots.to_csv(case_root / "marker_dotplot_data.csv", index=False)
    style()
    fig = plt.figure(figsize=(7.2, 8.2))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.25, 1.15, 1.7], hspace=0.48, wspace=0.35)
    for index, variant in enumerate(CASE_VARIANTS):
        ax = fig.add_subplot(gs[0, index])
        matrix = pd.read_csv(
            case_root / variant / "downstream" / f"split_{CASE_SPLIT_SEED}/marker_overlap.csv",
            index_col=0,
        )
        sns.heatmap(
            matrix, cmap="mako", cbar=index == 2,
            cbar_kws={"label": "shared top-100 markers"},
            xticklabels=True, yticklabels=index == 0, ax=ax,
        )
        ax.set_title(DISPLAY[variant])
        ax.set(xlabel="predicted cluster", ylabel="reference type" if index == 0 else "")
        ax.tick_params(labelsize=4)
        ax.text(-0.18, 1.06, chr(ord("a") + index), transform=ax.transAxes, fontweight="bold")
    ax_dot = fig.add_subplot(gs[1, :])
    gene_order = list(dict.fromkeys(gene for genes in panel.values() for gene in genes if gene in set(dots["gene"])))
    type_order = list(panel)
    color_max = max(1.0, dots["mean_log_expression"].quantile(0.98))
    for _, row in dots.iterrows():
        ax_dot.scatter(
            gene_order.index(row["gene"]), type_order.index(row["cell_type"]),
            s=3 + 70 * row["fraction_expressed"],
            c=row["mean_log_expression"], cmap="Reds", vmin=0,
            vmax=color_max,
        )
    ax_dot.set_xticks(range(len(gene_order)))
    ax_dot.set_xticklabels(gene_order, rotation=45, ha="right", fontsize=5)
    ax_dot.set_yticks(range(len(type_order)))
    ax_dot.set_yticklabels(type_order, fontsize=5)
    ax_dot.invert_yaxis()
    ax_dot.set_title("Fixed literature marker panel (gold-label expression; interpretation only)", fontsize=7)
    scalar = plt.cm.ScalarMappable(norm=plt.Normalize(0, color_max), cmap="Reds")
    colorbar = fig.colorbar(scalar, ax=ax_dot, fraction=0.012, pad=0.012)
    colorbar.set_label("mean log expression", fontsize=5)
    colorbar.ax.tick_params(labelsize=4)
    size_handles = [
        ax_dot.scatter([], [], s=3 + 70 * fraction, color="#D74B32", alpha=0.75, label=f"{fraction:.0%}")
        for fraction in (0.25, 0.50, 0.75)
    ]
    ax_dot.legend(
        handles=size_handles, title="fraction expressed", ncol=3,
        loc="upper right", bbox_to_anchor=(1.0, 1.13), fontsize=4,
        title_fontsize=4, handletextpad=0.2, columnspacing=0.6,
    )
    ax_dot.text(-0.04, 1.06, "d", transform=ax_dot.transAxes, fontweight="bold")
    for index, variant in enumerate(CASE_VARIANTS):
        ax = fig.add_subplot(gs[2, index])
        flow = flow_table(case_root, variant, data_path)
        flow.to_csv(case_root / f"{variant}_flow.csv", index=False)
        parallel_flow(ax, flow, DISPLAY[variant])
        ax.text(-0.12, 1.04, chr(ord("e") + index), transform=ax.transAxes, fontweight="bold")
    export(fig, target)


def main() -> None:
    parser = argparse.ArgumentParser(description="Prepare and plot the preregistered Human_Pancreas_3 case study")
    parser.add_argument("--run-master", type=Path, default=PAPER_ROOT / "experiments/protocol_v1/run_master.csv")
    parser.add_argument("--data-path", type=Path, default=PAPER_ROOT / ".staging/data/Human_Pancreas_3.h5ad")
    parser.add_argument("--marker-panel", type=Path, default=PAPER_ROOT / "configs/human_pancreas_marker_panel.json")
    parser.add_argument("--case-root", type=Path, default=PAPER_ROOT / "figures/source_data/human_pancreas_3")
    parser.add_argument("--output", type=Path, default=PAPER_ROOT / "figures/final/fig6_pancreas_case")
    parser.add_argument("--fetch", action="store_true")
    args = parser.parse_args()
    if args.fetch:
        fetch_case_sources(args.run_master, args.case_root)
    plot_case(args.case_root, args.data_path, args.marker_panel, args.output)


if __name__ == "__main__":
    main()
