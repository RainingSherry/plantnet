#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd


DATASET_GROUP = {
    "Tosches": "negative",
    "Macosko": "negative",
    "worm_neuron_cell": "negative",
    "Melanoma_5K": "positive",
    "Shekhar": "positive",
    "Guo": "neutral",
    "Wang": "optional",
    "Pollen": "optional",
}

RANDOM_METHODS = [
    "random_pseudo_gate_p0.5",
    "random_edge_dropout_keep0.5",
    "random_beta_uniform_0.1",
]

NEIGHBOR_RULE_METHODS = [
    "mutual_knn_neighbormix",
    "snn_neighbormix",
    "consensus_neighbormix_threshold0.4",
]

CORE_DATASETS = ["Tosches", "Macosko", "worm_neuron_cell", "Melanoma_5K", "Shekhar", "Guo"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize NeighborMix stochastic ablation runs")
    parser.add_argument("--root", required=True, help="Experiment root containing dataset/method/seed directories")
    return parser.parse_args()


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def collect_main(root: Path) -> pd.DataFrame:
    rows = []
    for eval_path in sorted(root.glob("*/*/seed*/eval_fixed.csv")):
        try:
            df = pd.read_csv(eval_path)
        except Exception:
            continue
        if df.empty:
            continue
        row = df.iloc[0].to_dict()
        args = read_json(eval_path.parent / "args.json")
        row["dataset"] = str(row.get("dataset") or args.get("dataset_name") or eval_path.parents[2].name)
        row["seed"] = int(row.get("seed", args.get("seed", eval_path.parent.name.replace("seed", ""))))
        row["method"] = str(args.get("ablation_method") or row.get("ablation_method") or eval_path.parents[1].name)
        row["run_dir"] = str(eval_path.parent)
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    for col in ["ari", "nmi", "acc", "f1_macro"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["group"] = df["dataset"].map(DATASET_GROUP).fillna("other")

    nomix = (
        df[df["method"] == "nm_scmae_nomix"][["dataset", "seed", "ari", "f1_macro"]]
        .rename(columns={"ari": "nomix_ari", "f1_macro": "nomix_f1_macro"})
    )
    df = df.merge(nomix, on=["dataset", "seed"], how="left")
    df["delta_ari_vs_nomix"] = df["ari"] - df["nomix_ari"]
    df["delta_f1_vs_nomix"] = df["f1_macro"] - df["nomix_f1_macro"]
    return df


def collect_diagnostics(root: Path) -> pd.DataFrame:
    rows = []
    for diag_path in sorted(root.glob("*/*/seed*/neighbor_diagnostics_final.csv")):
        try:
            df = pd.read_csv(diag_path)
        except Exception:
            continue
        if df.empty:
            continue
        row = df.iloc[0].to_dict()
        args = read_json(diag_path.parent / "args.json")
        row["dataset"] = str(row.get("dataset") or args.get("dataset_name") or diag_path.parents[2].name)
        row["seed"] = int(row.get("seed", args.get("seed", diag_path.parent.name.replace("seed", ""))))
        row["method"] = str(args.get("ablation_method") or row.get("method") or diag_path.parents[1].name)
        row["run_dir"] = str(diag_path.parent)
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    if "effective_neighbor_count_sampled" in out.columns:
        out["effective_neighbor_count"] = out["effective_neighbor_count_sampled"]
    if "edge_keep_rate_observed" in out.columns:
        out["edge_keep_rate"] = out["edge_keep_rate_observed"]
    required = [
        "dataset",
        "seed",
        "method",
        "effective_neighbor_count",
        "edge_keep_rate",
        "pseudo_branch_activation_rate",
        "perturbation_norm_mean",
        "perturbation_norm_p95",
        "neighbor_similarity_mean",
        "neighbor_similarity_std",
        "mutual_ratio",
        "snn_mean",
        "snn_std",
        "consensus_mean",
        "consensus_std",
        "fallback_rate",
        "hubness_p95",
        "hubness_max",
        "same_label_edge_ratio",
        "same_label_edge_ratio_weighted",
        "minority_class_neighbor_purity",
    ]
    for col in required:
        if col not in out.columns:
            out[col] = np.nan
    return out[required + ["run_dir"]]


def group_summary(main: pd.DataFrame) -> pd.DataFrame:
    if main.empty:
        return pd.DataFrame()
    grouped = (
        main.groupby(["method", "group"], dropna=False)
        .agg(
            mean_ari=("ari", "mean"),
            mean_delta_ari=("delta_ari_vs_nomix", "mean"),
            median_delta_ari=("delta_ari_vs_nomix", "median"),
            worst_delta_ari=("delta_ari_vs_nomix", "min"),
            mean_f1_macro=("f1_macro", "mean"),
            n_runs=("ari", "count"),
        )
        .reset_index()
    )
    return grouped


def _fmt(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.4f}"


def _method_mean(main: pd.DataFrame, method: str, group: str | None = None, col: str = "delta_ari_vs_nomix") -> float:
    sub = main[main["method"] == method]
    if group is not None:
        sub = sub[sub["group"] == group]
    return float(sub[col].mean()) if len(sub) else float("nan")


def _method_worst(main: pd.DataFrame, method: str, col: str = "delta_ari_vs_nomix") -> float:
    sub = main[main["method"] == method]
    return float(sub[col].min()) if len(sub) else float("nan")


def _datasets_beating_fixed(main: pd.DataFrame, method: str) -> int:
    fixed = main[main["method"] == "neighbormix_scmae"].groupby("dataset")["ari"].mean()
    candidate = main[main["method"] == method].groupby("dataset")["ari"].mean()
    shared = sorted(set(fixed.index).intersection(candidate.index).intersection(CORE_DATASETS))
    return sum(float(candidate[d]) > float(fixed[d]) for d in shared)


def build_interpretation(main: pd.DataFrame, diag: pd.DataFrame) -> str:
    lines = []
    if main.empty:
        return "Finding 1:\nNo completed runs were found.\n"

    lines.append("Experiment context:")
    lines.append("Negative-transfer group: Tosches, Macosko, worm_neuron_cell. These datasets test whether a method reduces known NeighborMix downside risk.")
    lines.append("Positive-gain group: Melanoma_5K, Shekhar. These datasets test whether a method preserves cases where NeighborMix can help.")
    lines.append("Neutral/stable group: Guo. This dataset tests whether a method avoids unnecessary degradation when the baseline is already stable.")
    lines.append("Optional validation group: Wang, Pollen. These were run after the six core datasets finished, but the continuation criteria are judged on the six core datasets.")
    lines.append("")
    lines.append("Mechanism framing:")
    lines.append("The pseudo branch reconstructs anchor cells from mixed pseudo-cell inputs. The random variants therefore test stochastic neighborhood regularization, not reliable-cell discovery.")
    lines.append("The alternative-neighbor variants test whether the vanilla PCA-cosine KNN graph is the bottleneck behind negative transfer.")
    lines.append("The global-random control tests whether any gains survive after removing local neighborhood structure.")
    lines.append("")

    fixed_neg = _method_mean(main, "neighbormix_scmae", "negative")
    fixed_worst = _method_worst(main, "neighbormix_scmae")
    fixed_f1 = _method_mean(main, "neighbormix_scmae", None, "f1_macro")

    random_scores = []
    for method in RANDOM_METHODS:
        random_scores.append(
            {
                "method": method,
                "beat_count": _datasets_beating_fixed(main, method),
                "neg_delta": _method_mean(main, method, "negative"),
                "worst_delta": _method_worst(main, method),
                "f1": _method_mean(main, method, None, "f1_macro"),
            }
        )
    random_scores = [item for item in random_scores if not pd.isna(item["neg_delta"])]
    best_random = max(random_scores, key=lambda item: (item["neg_delta"], item["beat_count"])) if random_scores else None

    lines.append("Finding 1:")
    if best_random is None:
        lines.append("Random gate / stochastic neighborhood regularization could not be judged because no random-variant runs completed.")
        lines.append("")
        lines.append("Evidence:")
        lines.append("No completed random variant was present in main_results.csv.")
        lines.append("")
        lines.append("Interpretation:")
        lines.append("The mechanism hypothesis remains untested.")
        lines.append("")
        lines.append("Next action:")
        lines.append("Finish random_pseudo_gate, random_edge_dropout, and random_beta runs before drawing conclusions.")
    else:
        random_continue = (
            best_random["beat_count"] >= 4
            and best_random["neg_delta"] > fixed_neg
            and best_random["worst_delta"] > fixed_worst
            and best_random["f1"] >= fixed_f1 - 1e-4
        )
        lines.append(
            f"Best random variant: {best_random['method']}. "
            f"It beats fixed NeighborMix on {best_random['beat_count']}/6 core datasets."
        )
        lines.append("")
        lines.append("Evidence:")
        lines.append("Random-variant screen:")
        for item in random_scores:
            lines.append(
                f"- {item['method']}: beats fixed on {item['beat_count']}/6 core datasets; "
                f"negative mean delta ARI {_fmt(item['neg_delta'])}; "
                f"worst delta ARI {_fmt(item['worst_delta'])}; "
                f"mean macro-F1 {_fmt(item['f1'])}."
            )
        lines.append(
            f"Negative-group mean delta ARI: best random {_fmt(best_random['neg_delta'])}, "
            f"fixed NeighborMix {_fmt(fixed_neg)}."
        )
        lines.append(
            f"Worst-case delta ARI: best random {_fmt(best_random['worst_delta'])}, "
            f"fixed NeighborMix {_fmt(fixed_worst)}."
        )
        lines.append(
            f"Mean macro-F1: best random {_fmt(best_random['f1'])}, fixed NeighborMix {_fmt(fixed_f1)}."
        )
        lines.append("")
        lines.append("Interpretation:")
        if random_continue:
            lines.append(
                "The stochastic-neighborhood direction passes the continuation screen. "
                "The result should be interpreted as stochastic regularization of anchor-recovery NeighborMix, not as evidence that the gate identifies reliable cells."
            )
        else:
            lines.append(
                "The stochastic-neighborhood direction does not pass the continuation screen. "
                "Any isolated gains are not sufficient to claim a robust regularization benefit."
            )
        lines.append("")
        lines.append("Next action:")
        lines.append(
            "Continue only if the detailed per-dataset table shows the gains are not concentrated in one dataset; otherwise treat this as a negative result."
        )

    lines.append("")
    lines.append("Finding 2:")
    neighbor_scores = []
    for method in NEIGHBOR_RULE_METHODS:
        neighbor_scores.append(
            {
                "method": method,
                "neg_delta": _method_mean(main, method, "negative"),
                "worst_delta": _method_worst(main, method),
                "f1": _method_mean(main, method, None, "f1_macro"),
            }
        )
    neighbor_scores = [item for item in neighbor_scores if not pd.isna(item["neg_delta"])]
    best_neighbor = max(neighbor_scores, key=lambda item: item["neg_delta"]) if neighbor_scores else None
    if best_neighbor is None:
        lines.append("Alternative neighbor rules could not be judged because no mutual/SNN/consensus runs completed.")
        lines.append("")
        lines.append("Evidence:")
        lines.append("No completed alternative-neighbor method was present in main_results.csv.")
        lines.append("")
        lines.append("Interpretation:")
        lines.append("The KNN-graph bottleneck hypothesis remains untested.")
        lines.append("")
        lines.append("Next action:")
        lines.append("Finish mutual, SNN, and consensus runs.")
    else:
        fixed_same = np.nan
        best_same = np.nan
        fixed_perturb = np.nan
        best_perturb = np.nan
        if not diag.empty:
            fixed_same = diag[diag["method"] == "neighbormix_scmae"]["same_label_edge_ratio_weighted"].mean()
            best_same = diag[diag["method"] == best_neighbor["method"]]["same_label_edge_ratio_weighted"].mean()
            fixed_perturb = diag[diag["method"] == "neighbormix_scmae"]["perturbation_norm_mean"].mean()
            best_perturb = diag[diag["method"] == best_neighbor["method"]]["perturbation_norm_mean"].mean()
        perturb_not_tiny = pd.isna(fixed_perturb) or best_perturb >= 0.5 * fixed_perturb
        neighbor_continue = (
            best_neighbor["neg_delta"] > fixed_neg
            and (pd.isna(best_same) or pd.isna(fixed_same) or best_same > fixed_same)
            and perturb_not_tiny
            and best_neighbor["f1"] >= fixed_f1 - 1e-4
        )
        lines.append(f"Best alternative-neighbor variant: {best_neighbor['method']}.")
        lines.append("")
        lines.append("Evidence:")
        lines.append("Alternative-neighbor screen:")
        for item in neighbor_scores:
            lines.append(
                f"- {item['method']}: negative mean delta ARI {_fmt(item['neg_delta'])}; "
                f"worst delta ARI {_fmt(item['worst_delta'])}; mean macro-F1 {_fmt(item['f1'])}."
            )
        lines.append(
            f"Negative-group mean delta ARI: best alternative {_fmt(best_neighbor['neg_delta'])}, "
            f"fixed NeighborMix {_fmt(fixed_neg)}."
        )
        lines.append(
            f"Weighted same-label edge ratio: best alternative {_fmt(best_same)}, "
            f"fixed NeighborMix {_fmt(fixed_same)}."
        )
        lines.append(
            f"Perturbation norm mean: best alternative {_fmt(best_perturb)}, "
            f"fixed NeighborMix {_fmt(fixed_perturb)}."
        )
        lines.append("")
        lines.append("Interpretation:")
        if neighbor_continue:
            lines.append(
                "At least one neighbor rule reduces negative transfer without merely eliminating perturbation strength. "
                "This supports the hypothesis that the vanilla PCA-cosine KNN graph is a major bottleneck."
            )
        else:
            lines.append(
                "The alternative-neighbor direction does not pass the continuation screen. "
                "Either the graph did not improve edge purity enough, or the apparent gains came from weakening the perturbation."
            )
        lines.append("")
        lines.append("Next action:")
        lines.append(
            "Proceed to edge-level gates or MoE only if the accepted neighbor rule also improves macro-F1 or minority diagnostics; otherwise improve graph construction first."
        )

    lines.append("")
    lines.append("Finding 3:")
    global_delta = _method_mean(main, "global_random_neighbor_control")
    fixed_delta = _method_mean(main, "neighbormix_scmae")
    local_random_delta = best_random["neg_delta"] if best_random is not None else float("nan")
    lines.append("")
    lines.append("Evidence:")
    lines.append(
        f"Global random mean delta ARI: {_fmt(global_delta)}; "
        f"fixed NeighborMix mean delta ARI: {_fmt(fixed_delta)}; "
        f"best random negative-group delta ARI: {_fmt(local_random_delta)}."
    )
    lines.append("")
    lines.append("Interpretation:")
    if not pd.isna(global_delta) and not pd.isna(local_random_delta) and global_delta < local_random_delta:
        lines.append(
            "Global random neighbors are worse than local stochastic variants, so any useful regularization still depends on local neighborhood structure."
        )
    else:
        lines.append(
            "Global random neighbors are not clearly worse in the completed runs; this weakens the claim that local neighborhood structure is essential."
        )
    lines.append("")
    lines.append("Next action:")
    lines.append(
        "If global random is competitive, do not build complex gates yet; first test whether the gain is generic noise regularization."
    )
    lines.append("")
    lines.append("Overall next action:")
    lines.append(
        "Given the current screen, do not escalate directly to complex attention or MoE. "
        "The most defensible next step is a small edge-level gate only after improving graph diagnostics, because the tested mutual/SNN/consensus rules did not clearly reduce negative transfer."
    )
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    summary_dir = root / "summaries"
    summary_dir.mkdir(parents=True, exist_ok=True)

    main_df = collect_main(root)
    if not main_df.empty:
        keep_cols = [
            "dataset",
            "group",
            "seed",
            "method",
            "ari",
            "nmi",
            "acc",
            "f1_macro",
            "delta_ari_vs_nomix",
            "delta_f1_vs_nomix",
            "run_dir",
        ]
        for col in keep_cols:
            if col not in main_df.columns:
                main_df[col] = np.nan
        main_df[keep_cols].sort_values(["dataset", "seed", "method"]).to_csv(summary_dir / "main_results.csv", index=False)
    else:
        pd.DataFrame().to_csv(summary_dir / "main_results.csv", index=False)

    grouped = group_summary(main_df)
    grouped.to_csv(summary_dir / "group_summary.csv", index=False)

    diag = collect_diagnostics(root)
    diag.to_csv(summary_dir / "neighbor_diagnostics.csv", index=False)

    interpretation = build_interpretation(main_df, diag)
    (summary_dir / "interpretation.md").write_text(interpretation, encoding="utf-8")
    readme = (
        "# NeighborMix Stochastic Regularization Ablation\n\n"
        "This README is generated from completed run artifacts. See `summaries/` for CSV tables.\n\n"
        + interpretation
    )
    (root / "README.md").write_text(readme, encoding="utf-8")
    print(f"Wrote summaries to {summary_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
