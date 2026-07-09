#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr


BASE_DIR = Path(__file__).resolve().parent
PCA_STRONG = {"Bach", "Macosko", "Wang", "Limb_Muscle", "hrvatin"}
PCA_WEAK = {"Pollen", "worm_neuron_cell"}
ORACLE_VARIANTS = ["pca_kmeans", "pca_spectral_kmeans", "rdg_cell_only", "rdg_gene_only", "rdg_concat_kmeans", "rdg_always_on"]
NEGATIVE_CONTROL_VARIANTS = {"neg_random_cell_graph", "neg_degree_shuffle_graph", "neg_shuffled_gene_cell_graph"}


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def collect_stage_a(run_root: Path) -> pd.DataFrame:
    rows = []
    for metrics_path in run_root.glob("*/*/*/metrics.json"):
        variant = metrics_path.parent.name
        seed = int(metrics_path.parent.parent.name.replace("seed_", ""))
        dataset = metrics_path.parent.parent.parent.name
        metrics = read_json(metrics_path)
        gate = read_json(metrics_path.parent / "gate_decision.json")
        rows.append(
            {
                "dataset": dataset,
                "seed": seed,
                "variant": variant,
                "ari": float(metrics.get("ari", np.nan)),
                "nmi": float(metrics.get("nmi", np.nan)),
                "acc": float(metrics.get("acc", np.nan)),
                "f1_macro": float(metrics.get("f1_macro", np.nan)),
                "q_cell": float(gate.get("q_cell", np.nan)),
                "q_gene": float(gate.get("q_gene", np.nan)),
                "q_total": float(gate.get("q_total", np.nan)),
                "graph_enabled": bool(gate.get("graph_enabled", False)),
                "negative_control": gate.get("negative_control", ""),
                "chosen_branch": gate.get("chosen_branch", ""),
                "selector_type": gate.get("selector_type", ""),
                "uses_rdg_features": bool(gate.get("uses_rdg_features", False)),
                "uses_graph_clustering": bool(gate.get("uses_graph_clustering", False)),
            }
        )
    return pd.DataFrame(rows)


def append_calibrated(df: pd.DataFrame, analysis_dir: Path) -> pd.DataFrame:
    parts = [df]
    for path in [analysis_dir / "calibrated_threshold_runs.csv", analysis_dir / "calibrated_logistic_runs.csv"]:
        if path.exists():
            cal = pd.read_csv(path)
            cal["f1_macro"] = np.nan
            cal["negative_control"] = ""
            cal["chosen_branch"] = ""
            cal["selector_type"] = "label_calibrated_meta_selector"
            cal["uses_rdg_features"] = False
            cal["uses_graph_clustering"] = cal["graph_enabled"].fillna(False).astype(bool)
            cal["q_cell"] = np.nan
            cal["q_gene"] = np.nan
            parts.append(cal[["dataset", "seed", "variant", "ari", "nmi", "acc", "f1_macro", "q_cell", "q_gene", "q_total", "graph_enabled", "negative_control", "chosen_branch", "selector_type", "uses_rdg_features", "uses_graph_clustering"]])
    return pd.concat(parts, ignore_index=True)


def add_paired_metrics(df: pd.DataFrame) -> pd.DataFrame:
    pca = df[df["variant"] == "pca_kmeans"][["dataset", "seed", "ari", "nmi", "acc"]].rename(columns={"ari": "pca_ari", "nmi": "pca_nmi", "acc": "pca_acc"})
    out = df.merge(pca, on=["dataset", "seed"], how="left")
    out["delta_ari_vs_pca"] = out["ari"] - out["pca_ari"]
    out["regret_vs_pca"] = np.maximum(0.0, out["pca_ari"] - out["ari"])
    out["negative_transfer"] = out["delta_ari_vs_pca"] < -1e-12
    out["win_vs_pca"] = out["delta_ari_vs_pca"] > 1e-12
    out["tie_vs_pca"] = out["delta_ari_vs_pca"].abs() <= 1e-12
    out["loss_vs_pca"] = out["delta_ari_vs_pca"] < -1e-12
    out["rank_ari_desc"] = out.groupby(["dataset", "seed"])["ari"].rank(method="average", ascending=False)
    return out


def add_oracle_gap(df: pd.DataFrame, oracle_df: pd.DataFrame) -> pd.DataFrame:
    oracle_cols = oracle_df[["dataset", "seed", "oracle_ari", "oracle_source_variant", "oracle_gain_vs_pca"]]
    out = df.merge(oracle_cols, on=["dataset", "seed"], how="left")
    out["oracle_gap"] = out["oracle_ari"] - out["ari"]
    return out


def summarize(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    grouped = df.groupby("variant", dropna=False)
    rows = []
    for variant, g in grouped:
        rows.append(
            {
                "variant": variant,
                "n_runs": len(g),
                "n_datasets": g["dataset"].nunique(),
                "ari_mean": g["ari"].mean(),
                "ari_median": g["ari"].median(),
                "nmi_mean": g["nmi"].mean(),
                "acc_mean": g["acc"].mean(),
                "rank_ari_mean": g["rank_ari_desc"].mean(),
                "rank_ari_median": g["rank_ari_desc"].median(),
                "wins_vs_pca": int(g["win_vs_pca"].sum()),
                "ties_vs_pca": int(g["tie_vs_pca"].sum()),
                "losses_vs_pca": int(g["loss_vs_pca"].sum()),
                "delta_ari_vs_pca_mean": g["delta_ari_vs_pca"].mean(),
                "regret_vs_pca_mean": g["regret_vs_pca"].mean(),
                "oracle_gap_mean": g["oracle_gap"].mean() if "oracle_gap" in g else np.nan,
                "negative_transfer_rate": g["negative_transfer"].mean(),
                "graph_activation_rate": g.get("graph_enabled", pd.Series(False, index=g.index)).fillna(False).astype(bool).mean(),
                "enabled_graph_mean_gain": g.loc[g.get("graph_enabled", pd.Series(False, index=g.index)).fillna(False).astype(bool), "delta_ari_vs_pca"].mean(),
                "chosen_non_pca_rate": (g.get("chosen_branch", pd.Series("", index=g.index)).fillna("") != "").mean() if "chosen_branch" in g else 0.0,
            }
        )
    by_variant = pd.DataFrame(rows).sort_values("ari_mean", ascending=False)
    by_dataset = (
        df.groupby(["variant", "dataset"], dropna=False)
        .agg(
            n_runs=("ari", "size"),
            ari_mean=("ari", "mean"),
            nmi_mean=("nmi", "mean"),
            acc_mean=("acc", "mean"),
            rank_ari_mean=("rank_ari_desc", "mean"),
            wins_vs_pca=("win_vs_pca", "sum"),
            ties_vs_pca=("tie_vs_pca", "sum"),
            losses_vs_pca=("loss_vs_pca", "sum"),
            delta_ari_vs_pca_mean=("delta_ari_vs_pca", "mean"),
            regret_vs_pca_mean=("regret_vs_pca", "mean"),
            oracle_gap_mean=("oracle_gap", "mean"),
            negative_transfer_rate=("negative_transfer", "mean"),
            graph_activation_rate=("graph_enabled", lambda x: x.fillna(False).astype(bool).mean()),
            chosen_non_pca_rate=("chosen_branch", lambda x: (x.fillna("") != "").mean()),
        )
        .reset_index()
    )
    return by_variant, by_dataset


def oracle(df: pd.DataFrame) -> pd.DataFrame:
    base = df[df["variant"].isin(ORACLE_VARIANTS)].copy()
    idx = base.groupby(["dataset", "seed"])["ari"].idxmax()
    best = base.loc[idx, ["dataset", "seed", "variant", "ari", "nmi", "acc"]].rename(columns={"variant": "oracle_source_variant", "ari": "oracle_ari", "nmi": "oracle_nmi", "acc": "oracle_acc"})
    pca = df[df["variant"] == "pca_kmeans"][["dataset", "seed", "ari"]].rename(columns={"ari": "pca_ari"})
    return best.merge(pca, on=["dataset", "seed"], how="left").assign(oracle_gain_vs_pca=lambda x: x["oracle_ari"] - x["pca_ari"])


def q_correlations(df: pd.DataFrame) -> dict:
    rows = df[df["variant"] == "rdg_always_on"].dropna(subset=["q_total", "delta_ari_vs_pca"])
    if len(rows) < 3:
        return {}
    out = {"n": int(len(rows))}
    for q_col in ["q_cell", "q_gene", "q_total"]:
        if q_col not in rows or rows[q_col].isna().all():
            continue
        out[f"spearman_{q_col}_delta_ari"] = float(spearmanr(rows[q_col], rows["delta_ari_vs_pca"]).correlation)
        out[f"pearson_{q_col}_delta_ari"] = float(pearsonr(rows[q_col], rows["delta_ari_vs_pca"])[0])
    return out


def negative_control_report(df: pd.DataFrame) -> dict:
    neg = df[df["variant"].isin(NEGATIVE_CONTROL_VARIANTS)].copy()
    if neg.empty:
        return {"available": False}
    rows = {}
    for variant, g in neg.groupby("variant"):
        rows[variant] = {
            "n_runs": int(len(g)),
            "n_datasets": int(g["dataset"].nunique()),
            "delta_ari_vs_pca_mean": float(g["delta_ari_vs_pca"].mean()),
            "regret_vs_pca_mean": float(g["regret_vs_pca"].mean()),
            "negative_transfer_rate": float(g["negative_transfer"].mean()),
            "positive_gain_rate": float((g["delta_ari_vs_pca"] > 0).mean()),
            "mean_positive_gain_when_any": None if (g["delta_ari_vs_pca"] > 0).sum() == 0 else float(g.loc[g["delta_ari_vs_pca"] > 0, "delta_ari_vs_pca"].mean()),
        }
    return {"available": True, "variants": rows}


def selector_report(df: pd.DataFrame) -> dict:
    rows = df[df["variant"] == "safe_rdg_pca_u"].copy()
    if rows.empty:
        return {"available": False}
    return {
        "available": True,
        "chosen_branch_counts": {str(k): int(v) for k, v in rows["chosen_branch"].fillna("pca_kmeans").replace("", "pca_kmeans").value_counts().items()},
        "selector_type_counts": {str(k): int(v) for k, v in rows["selector_type"].fillna("").value_counts().items()},
        "graph_activation_rate": float(rows["graph_enabled"].fillna(False).astype(bool).mean()),
        "mean_delta_ari_vs_pca": float(rows["delta_ari_vs_pca"].mean()),
        "mean_regret_vs_pca": float(rows["regret_vs_pca"].mean()),
        "negative_transfer_rate": float(rows["negative_transfer"].mean()),
    }


def success(by_variant: pd.DataFrame, by_dataset: pd.DataFrame, oracle_df: pd.DataFrame) -> dict:
    pca_ari = float(by_variant.loc[by_variant["variant"] == "pca_kmeans", "ari_mean"].iloc[0])
    result = {"pca_ari_mean": pca_ari}
    for variant in ["safe_rdg_pca_u", "safe_rdg_heuristic", "safe_rdg_calibrated_threshold"]:
        row = by_variant[by_variant["variant"] == variant]
        if row.empty:
            continue
        row = row.iloc[0]
        ds = by_dataset[by_dataset["variant"] == variant]
        strong_drop = ds[ds["dataset"].isin(PCA_STRONG) & (ds["delta_ari_vs_pca_mean"] < -0.03)]["dataset"].nunique()
        weak_gain = ds[ds["dataset"].isin(PCA_WEAK)]["delta_ari_vs_pca_mean"].max()
        always = by_variant[by_variant["variant"] == "rdg_always_on"].iloc[0]
        result[variant] = {
            "mean_or_median_noninferior_to_pca": bool(row["ari_mean"] >= pca_ari or row["ari_median"] >= pca_ari),
            "strong_dataset_drop_count_gt_0p03": int(strong_drop),
            "weak_dataset_best_gain": float(weak_gain) if pd.notna(weak_gain) else None,
            "graph_activation_count": int(round(row["graph_activation_rate"] * row["n_runs"])),
            "enabled_graph_mean_gain": None if pd.isna(row["enabled_graph_mean_gain"]) else float(row["enabled_graph_mean_gain"]),
            "mean_regret_lt_always_on": bool(row["regret_vs_pca_mean"] < always["regret_vs_pca_mean"]),
            "negative_transfer_lt_always_on": bool(row["negative_transfer_rate"] < always["negative_transfer_rate"]),
        }
    result["oracle_gain_mean_vs_pca"] = float(oracle_df["oracle_gain_vs_pca"].mean()) if not oracle_df.empty else None
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--run_root", type=Path, default=BASE_DIR / "runs")
    parser.add_argument("--out_dir", type=Path, default=BASE_DIR)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    df = collect_stage_a(args.run_root)
    if df.empty:
        raise SystemExit(f"No runs found under {args.run_root}")
    args.out_dir.mkdir(parents=True, exist_ok=True)
    df = append_calibrated(df, args.out_dir)
    df = add_paired_metrics(df)
    oracle_df = oracle(df)
    df = add_oracle_gap(df, oracle_df)
    by_variant, by_dataset = summarize(df)
    df.to_csv(args.out_dir / "all_runs.csv", index=False)
    by_variant.to_csv(args.out_dir / "summary_by_variant.csv", index=False)
    by_dataset.to_csv(args.out_dir / "summary_by_variant_dataset.csv", index=False)
    oracle_df.to_csv(args.out_dir / "oracle_best_runs.csv", index=False)
    report = {"q_correlations": q_correlations(df), "negative_controls": negative_control_report(df), "selector": selector_report(df), "success": success(by_variant, by_dataset, oracle_df)}
    (args.out_dir / "summary_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(by_variant.to_string(index=False))
    print(json.dumps(report, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
