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

CORE_DATASETS = ["Tosches", "Macosko", "worm_neuron_cell", "Melanoma_5K", "Shekhar", "Guo"]
STAGE1_PAIRS = [
    ("fixed_beta_0.025", "random_beta_uniform_0.05", 0.025),
    ("fixed_beta_0.05", "random_beta_uniform_0.1", 0.05),
    ("fixed_beta_0.1", "random_beta_uniform_0.2", 0.1),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize staged NeighborMix beta mechanism experiments")
    parser.add_argument("--root", required=True)
    return parser.parse_args()


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def parse_path(eval_path: Path, root: Path) -> tuple[str, str, str, int]:
    rel = eval_path.parent.relative_to(root)
    parts = rel.parts
    if len(parts) >= 4 and parts[-1].startswith("seed"):
        stage, dataset, variant, seed_text = parts[-4], parts[-3], parts[-2], parts[-1]
        return stage, dataset, variant, int(seed_text.replace("seed", ""))
    dataset = eval_path.parents[2].name
    variant = eval_path.parents[1].name
    seed = int(eval_path.parent.name.replace("seed", ""))
    return "unknown", dataset, variant, seed


def collect_main(root: Path) -> pd.DataFrame:
    rows = []
    for eval_path in sorted(root.glob("*/*/*/seed*/eval_fixed.csv")):
        try:
            df = pd.read_csv(eval_path)
        except Exception:
            continue
        if df.empty:
            continue
        args = read_json(eval_path.parent / "args.json")
        stage, dataset, variant, seed = parse_path(eval_path, root)
        row = df.iloc[0].to_dict()
        row["stage"] = stage
        row["dataset"] = str(row.get("dataset") or args.get("dataset_name") or dataset)
        row["variant"] = str(args.get("variant_name") or row.get("variant") or variant)
        row["method"] = row["variant"]
        row["ablation_method"] = str(args.get("ablation_method") or row.get("ablation_method") or "")
        row["seed"] = int(row.get("seed", args.get("seed", seed)))
        row["run_dir"] = str(eval_path.parent)
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    for col in ["ari", "nmi", "acc", "f1_macro", "beta_mean", "beta_fixed", "beta_max", "bad_edge_ratio"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    df["group"] = df["dataset"].map(DATASET_GROUP).fillna("other")
    nomix = df[df["variant"] == "nm_scmae_nomix"][["stage", "dataset", "seed", "ari", "f1_macro"]].rename(
        columns={"ari": "nomix_ari", "f1_macro": "nomix_f1_macro"}
    )
    df = df.merge(nomix, on=["stage", "dataset", "seed"], how="left")
    missing = df["nomix_ari"].isna()
    if missing.any():
        stage1_nomix = (
            df[(df["stage"] == "stage1") & (df["variant"] == "nm_scmae_nomix")][["dataset", "seed", "ari", "f1_macro"]]
            .drop_duplicates(["dataset", "seed"])
            .rename(columns={"ari": "stage1_nomix_ari", "f1_macro": "stage1_nomix_f1_macro"})
        )
        df = df.merge(stage1_nomix, on=["dataset", "seed"], how="left")
        missing = df["nomix_ari"].isna()
        df.loc[missing, "nomix_ari"] = df.loc[missing, "stage1_nomix_ari"]
        df.loc[missing, "nomix_f1_macro"] = df.loc[missing, "stage1_nomix_f1_macro"]
        df = df.drop(columns=["stage1_nomix_ari", "stage1_nomix_f1_macro"])
    df["delta_ari_vs_nomix"] = df["ari"] - df["nomix_ari"]
    df["delta_f1_vs_nomix"] = df["f1_macro"] - df["nomix_f1_macro"]
    rank_base = df.dropna(subset=["ari"]).copy()
    if not rank_base.empty:
        rank_base["ari_rank"] = rank_base.groupby(["stage", "dataset", "seed"])["ari"].rank(ascending=False, method="average")
        df = df.merge(
            rank_base[["stage", "dataset", "seed", "variant", "ari_rank"]],
            on=["stage", "dataset", "seed", "variant"],
            how="left",
        )
    return df


def collect_diagnostics(root: Path) -> pd.DataFrame:
    rows = []
    for diag_path in sorted(root.glob("*/*/*/seed*/neighbor_diagnostics_final.csv")):
        try:
            df = pd.read_csv(diag_path)
        except Exception:
            continue
        if df.empty:
            continue
        args = read_json(diag_path.parent / "args.json")
        stage, dataset, variant, seed = parse_path(diag_path, root)
        row = df.iloc[0].to_dict()
        row["stage"] = stage
        row["dataset"] = str(row.get("dataset") or args.get("dataset_name") or dataset)
        row["variant"] = str(args.get("variant_name") or row.get("method") or variant)
        row["method"] = row["variant"]
        row["seed"] = int(row.get("seed", args.get("seed", seed)))
        row["run_dir"] = str(diag_path.parent)
        rows.append(row)
    if not rows:
        return pd.DataFrame()
    out = pd.DataFrame(rows)
    for col in out.columns:
        if col not in {"stage", "dataset", "variant", "method", "neighbor_backend", "oracle_neighbor", "target_mode", "noise_mode", "run_dir"}:
            try:
                out[col] = pd.to_numeric(out[col])
            except (TypeError, ValueError):
                pass
    return out


def group_summary(main: pd.DataFrame) -> pd.DataFrame:
    if main.empty:
        return pd.DataFrame()
    return (
        main.groupby(["stage", "variant", "group"], dropna=False)
        .agg(
            mean_ari=("ari", "mean"),
            mean_delta_ari=("delta_ari_vs_nomix", "mean"),
            median_delta_ari=("delta_ari_vs_nomix", "median"),
            worst_delta_ari=("delta_ari_vs_nomix", "min"),
            mean_f1_macro=("f1_macro", "mean"),
            mean_delta_f1=("delta_f1_vs_nomix", "mean"),
            mean_rank=("ari_rank", "mean"),
            n_runs=("ari", "count"),
        )
        .reset_index()
    )


def stage1_pairwise(main: pd.DataFrame) -> pd.DataFrame:
    stage = main[main["stage"] == "stage1"].copy()
    if stage.empty:
        return pd.DataFrame()
    rows = []
    for fixed, random, mean_beta in STAGE1_PAIRS:
        fixed_df = stage[stage["variant"] == fixed][["dataset", "seed", "ari", "f1_macro", "delta_ari_vs_nomix", "delta_f1_vs_nomix"]]
        rand_df = stage[stage["variant"] == random][["dataset", "seed", "ari", "f1_macro", "delta_ari_vs_nomix", "delta_f1_vs_nomix"]]
        merged = fixed_df.merge(rand_df, on=["dataset", "seed"], suffixes=("_fixed", "_random"))
        if merged.empty:
            rows.append(
                {
                    "mean_beta": mean_beta,
                    "fixed_variant": fixed,
                    "random_variant": random,
                    "n_pairs": 0,
                }
            )
            continue
        core = merged[merged["dataset"].isin(CORE_DATASETS)]
        rows.append(
            {
                "mean_beta": mean_beta,
                "fixed_variant": fixed,
                "random_variant": random,
                "n_pairs": int(len(merged)),
                "random_minus_fixed_ari_mean": float((merged["ari_random"] - merged["ari_fixed"]).mean()),
                "random_minus_fixed_f1_mean": float((merged["f1_macro_random"] - merged["f1_macro_fixed"]).mean()),
                "random_wins_run_count": int((merged["ari_random"] > merged["ari_fixed"]).sum()),
                "fixed_wins_run_count": int((merged["ari_random"] < merged["ari_fixed"]).sum()),
                "ties_run_count": int((merged["ari_random"] == merged["ari_fixed"]).sum()),
                "core_random_minus_fixed_ari_mean": float((core["ari_random"] - core["ari_fixed"]).mean()) if len(core) else np.nan,
                "core_random_wins_dataset_count": int(
                    sum(
                        core[core["dataset"] == d]["ari_random"].mean() > core[core["dataset"] == d]["ari_fixed"].mean()
                        for d in sorted(core["dataset"].unique())
                    )
                )
                if len(core)
                else 0,
            }
        )
    return pd.DataFrame(rows)


def stage_metric_table(main: pd.DataFrame, stage_name: str) -> pd.DataFrame:
    stage = main[main["stage"] == stage_name].copy()
    if stage.empty:
        return pd.DataFrame()
    return (
        stage.groupby(["variant"], dropna=False)
        .agg(
            mean_ari=("ari", "mean"),
            mean_delta_ari=("delta_ari_vs_nomix", "mean"),
            worst_delta_ari=("delta_ari_vs_nomix", "min"),
            mean_f1_macro=("f1_macro", "mean"),
            mean_delta_f1=("delta_f1_vs_nomix", "mean"),
            mean_rank=("ari_rank", "mean"),
            n_runs=("ari", "count"),
        )
        .reset_index()
        .sort_values(["mean_ari", "mean_f1_macro"], ascending=False)
    )


def stage4_table(main: pd.DataFrame) -> pd.DataFrame:
    stage = main[main["stage"] == "stage4"].copy()
    if stage.empty:
        return pd.DataFrame()
    stage["base_variant"] = stage["variant"].str.replace(r"_bad[0-9.]+$", "", regex=True)
    return (
        stage.groupby(["base_variant", "bad_edge_ratio"], dropna=False)
        .agg(
            mean_ari=("ari", "mean"),
            mean_delta_ari=("delta_ari_vs_nomix", "mean"),
            worst_delta_ari=("delta_ari_vs_nomix", "min"),
            mean_f1_macro=("f1_macro", "mean"),
            n_runs=("ari", "count"),
        )
        .reset_index()
        .sort_values(["bad_edge_ratio", "mean_ari"], ascending=[True, False])
    )


def fmt(value: float) -> str:
    if pd.isna(value):
        return "NA"
    return f"{value:.4f}"


def build_interpretation(main: pd.DataFrame, pairwise: pd.DataFrame) -> str:
    if main.empty:
        return "Finding 1:\nNo completed runs were found.\n"
    lines = []
    lines.append("Experiment context:")
    lines.append("This package tests whether random_beta_uniform_0.1 works because the average perturbation is lower or because beta is stochastic.")
    lines.append("Stage 1 is the required decision point: fixed_beta_0.05 vs random_beta_uniform_0.1 controls for mean beta = 0.05.")
    lines.append("")
    completed = main.groupby("stage")["run_dir"].count().to_dict()
    lines.append("Completed run counts:")
    for stage, count in sorted(completed.items()):
        lines.append(f"- {stage}: {int(count)} runs")
    lines.append("")
    lines.append("Finding 1:")
    if pairwise.empty or not (pairwise["random_variant"] == "random_beta_uniform_0.1").any():
        lines.append("Stage 1 same-mean stochasticity could not be judged yet because the required pairwise table is incomplete.")
        lines.append("")
        lines.append("Evidence:")
        lines.append("stage1_beta_mean_vs_randomness.csv is empty or missing the fixed_beta_0.05 vs random_beta_uniform_0.1 comparison.")
        lines.append("")
        lines.append("Interpretation:")
        lines.append("Do not claim stochastic beta NeighborMix until this comparison is available.")
        lines.append("")
        lines.append("Next action:")
        lines.append("Finish Stage 1.")
        return "\n".join(lines) + "\n"
    row = pairwise[pairwise["random_variant"] == "random_beta_uniform_0.1"].iloc[0]
    delta = float(row.get("random_minus_fixed_ari_mean", np.nan))
    core_delta = float(row.get("core_random_minus_fixed_ari_mean", np.nan))
    wins = int(row.get("random_wins_run_count", 0))
    fixed_wins = int(row.get("fixed_wins_run_count", 0))
    lines.append(
        "The key same-mean comparison is "
        f"random_beta_uniform_0.1 minus fixed_beta_0.05: mean ARI delta {fmt(delta)}, "
        f"core mean ARI delta {fmt(core_delta)}, run wins {wins} vs {fixed_wins}."
    )
    lines.append("")
    lines.append("Evidence:")
    for _, item in pairwise.iterrows():
        lines.append(
            f"- mean beta {item['mean_beta']}: {item['random_variant']} - {item['fixed_variant']} "
            f"mean ARI {fmt(item.get('random_minus_fixed_ari_mean', np.nan))}, "
            f"mean macro-F1 {fmt(item.get('random_minus_fixed_f1_mean', np.nan))}, "
            f"wins {int(item.get('random_wins_run_count', 0))}/{int(item.get('n_pairs', 0))}."
        )
    lines.append("")
    lines.append("Interpretation:")
    if delta > 0 and wins > fixed_wins:
        lines.append("Stage 1 supports a stochasticity contribution beyond lower mean beta. Continue to Stage 2 variance controls and Stage 3 mechanism controls.")
    else:
        lines.append("Stage 1 does not yet support stochasticity beyond lower mean beta. The likely main line should shift toward low-strength NeighborMix unless later robustness criteria overturn this.")
    lines.append("")
    lines.append("Next action:")
    lines.append("Use Stage 2 only to distinguish random distributions if the key same-mean comparison remains favorable or ambiguous.")
    return "\n".join(lines) + "\n"


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    summary_dir = root / "summaries"
    summary_dir.mkdir(parents=True, exist_ok=True)

    main_df = collect_main(root)
    diag_df = collect_diagnostics(root)
    group_df = group_summary(main_df)
    pairwise_df = stage1_pairwise(main_df)

    main_df.to_csv(summary_dir / "main_results.csv", index=False)
    diag_df.to_csv(summary_dir / "neighbor_diagnostics.csv", index=False)
    group_df.to_csv(summary_dir / "group_summary.csv", index=False)
    pairwise_df.to_csv(summary_dir / "stage1_beta_mean_vs_randomness.csv", index=False)
    stage_metric_table(main_df, "stage2").to_csv(summary_dir / "stage2_beta_variance.csv", index=False)
    stage_metric_table(main_df, "stage3").to_csv(summary_dir / "stage3_local_mix_mechanism.csv", index=False)
    stage4_table(main_df).to_csv(summary_dir / "stage4_bad_edge_robustness.csv", index=False)
    stage_metric_table(main_df, "full").to_csv(summary_dir / "full_benchmark_summary.csv", index=False)
    (summary_dir / "interpretation.md").write_text(build_interpretation(main_df, pairwise_df), encoding="utf-8")

    print(f"main_results rows: {len(main_df)}")
    print(f"neighbor_diagnostics rows: {len(diag_df)}")
    print(f"summaries written to: {summary_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
