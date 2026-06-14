#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[3]
METRIC_COLUMNS = ["acc", "nmi", "ari", "f1_macro", "fmi"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize CutAware NeighborMix representative ablations.")
    parser.add_argument("--out_dir", default="results/experimental/cutaware_neighbormix_20260615")
    parser.add_argument("--baseline_raw", default="results/formal/rg_phase1_allseeds_e80/rg_phase1_allseeds_raw.csv")
    parser.add_argument("--phase2_raw", default="results/formal/rg_phase2_sensitivity_e80/rg_phase2_all_sweeps_raw.csv")
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def read_json(path: Path) -> dict:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def collect_cutaware(out_root: Path) -> pd.DataFrame:
    rows = []
    for eval_path in sorted(out_root.glob("*/*/seed*/eval_fixed.csv")):
        run_dir = eval_path.parent
        df = pd.read_csv(eval_path)
        if df.empty:
            continue
        row = df.iloc[0].to_dict()
        summary = read_json(run_dir / "summary.json")
        cut = read_json(run_dir / "cut_diagnostics.json")
        sim = read_json(run_dir / "embedding_similarity_diagnostics.json")
        cut_reweight = summary.get("cut_reweight_summary", {})
        row.update(
            {
                "source": "cutaware",
                "dataset": row.get("dataset", eval_path.parts[-4]),
                "method_key": row.get("variant", eval_path.parts[-3]),
                "run_path": str(run_dir.relative_to(ROOT)) if run_dir.is_relative_to(ROOT) else str(run_dir),
                "soft_cut": cut.get("soft_cut", ""),
                "ncut_surrogate": cut.get("ncut_surrogate", ""),
                "cluster_mass_min": cut.get("cluster_mass_min", ""),
                "cluster_mass_max": cut.get("cluster_mass_max", ""),
                "fraction_cosine_gt_0p9": sim.get("fraction_cosine_gt_0p9", ""),
                "max_cluster_fraction": cut.get("max_cluster_fraction", ""),
                "cut_enabled": summary.get("cut_enabled", ""),
                "ot_enabled": summary.get("ot_enabled", ""),
                "pseudo_enabled": summary.get("pseudo_enabled", ""),
                "cut_reweight_enabled": cut_reweight.get("cut_reweight_enabled", ""),
                "cross_edge_mass_before": cut_reweight.get("cross_edge_mass_before", ""),
                "cross_edge_mass_after": cut_reweight.get("cross_edge_mass_after", ""),
                "cross_edge_mass_reduction": cut_reweight.get("cross_edge_mass_reduction", ""),
            }
        )
        rows.append(row)
    return pd.DataFrame(rows)


def collect_rg_baseline(path: Path, seed: int) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "seed" in df.columns:
        df = df[df["seed"] == seed].copy()
    keep_methods = {"rg_none", "rg_fixed", "rg_reliability", "rg_random", "rg_far"}
    if "method" in df.columns:
        df = df[df["method"].isin(keep_methods)].copy()
    if df.empty:
        return df
    df["source"] = "rg_phase1_seed42"
    df["method_key"] = df["method"]
    df["run_path"] = df.get("path", "")
    return df


def collect_rg_phase2_best(path: Path, seed: int) -> pd.DataFrame:
    if not path.exists():
        return pd.DataFrame()
    df = pd.read_csv(path)
    if "seed" in df.columns:
        df = df[df["seed"] == seed].copy()
    if df.empty:
        return df
    idx = df.groupby("dataset")["ari"].idxmax()
    best = df.loc[idx].copy()
    best["source"] = "rg_phase2_best_seed42"
    best["method_key"] = "rg_phase2_best_" + best["sweep"].astype(str) + "_" + best["label"].astype(str)
    best["run_path"] = best.get("path", "")
    return best


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    for col in METRIC_COLUMNS:
        if col not in df.columns:
            df[col] = ""
    wanted = [
        "source",
        "dataset",
        "method_key",
        "seed",
        *METRIC_COLUMNS,
        "soft_cut",
        "ncut_surrogate",
        "cluster_mass_min",
        "cluster_mass_max",
        "max_cluster_fraction",
        "fraction_cosine_gt_0p9",
        "cut_enabled",
        "ot_enabled",
        "pseudo_enabled",
        "cut_reweight_enabled",
        "cross_edge_mass_before",
        "cross_edge_mass_after",
        "cross_edge_mass_reduction",
        "run_path",
    ]
    for col in wanted:
        if col not in df.columns:
            df[col] = ""
    return df[wanted].copy()


def add_deltas(all_rows: pd.DataFrame) -> pd.DataFrame:
    if all_rows.empty:
        return all_rows
    out = all_rows.copy()
    baseline = out[out["method_key"] == "rg_none"][["dataset", "ari"]].rename(columns={"ari": "rg_none_ari"})
    fixed = out[out["method_key"] == "rg_fixed"][["dataset", "ari"]].rename(columns={"ari": "rg_fixed_ari"})
    reliability = out[out["method_key"] == "rg_reliability"][["dataset", "ari"]].rename(columns={"ari": "rg_reliability_ari"})
    out = out.merge(baseline, on="dataset", how="left").merge(fixed, on="dataset", how="left").merge(reliability, on="dataset", how="left")
    for base in ["rg_none", "rg_fixed", "rg_reliability"]:
        out[f"delta_ari_vs_{base}"] = pd.to_numeric(out["ari"], errors="coerce") - pd.to_numeric(out[f"{base}_ari"], errors="coerce")
    return out


def summarize(all_rows: pd.DataFrame) -> pd.DataFrame:
    if all_rows.empty:
        return all_rows
    tmp = all_rows.copy()
    tmp["ari"] = pd.to_numeric(tmp["ari"], errors="coerce")
    tmp["delta_ari_vs_rg_none"] = pd.to_numeric(tmp["delta_ari_vs_rg_none"], errors="coerce")
    tmp["delta_ari_vs_rg_reliability"] = pd.to_numeric(tmp["delta_ari_vs_rg_reliability"], errors="coerce")
    return (
        tmp.groupby(["source", "method_key"], as_index=False)
        .agg(
            n_runs=("ari", "count"),
            mean_ari=("ari", "mean"),
            min_ari=("ari", "min"),
            max_ari=("ari", "max"),
            mean_delta_ari_vs_rg_none=("delta_ari_vs_rg_none", "mean"),
            mean_delta_ari_vs_rg_reliability=("delta_ari_vs_rg_reliability", "mean"),
        )
        .sort_values(["mean_ari", "method_key"], ascending=[False, True])
    )


def main() -> int:
    args = parse_args()
    out_root = Path(args.out_dir)
    if not out_root.is_absolute():
        out_root = ROOT / out_root
    out_root.mkdir(parents=True, exist_ok=True)

    parts = [
        normalize_columns(collect_rg_baseline(ROOT / args.baseline_raw, args.seed)),
        normalize_columns(collect_rg_phase2_best(ROOT / args.phase2_raw, args.seed)),
        normalize_columns(collect_cutaware(out_root)),
    ]
    all_rows = pd.concat([p for p in parts if not p.empty], ignore_index=True) if any(not p.empty for p in parts) else pd.DataFrame()
    all_rows = add_deltas(all_rows)
    summary = summarize(all_rows)
    detail_path = out_root / "cutaware_vs_rg_detail.csv"
    summary_path = out_root / "cutaware_vs_rg_summary.csv"
    all_rows.to_csv(detail_path, index=False)
    summary.to_csv(summary_path, index=False)
    print(f"wrote {len(all_rows)} detail rows: {detail_path}")
    print(f"wrote {len(summary)} summary rows: {summary_path}")
    if not summary.empty:
        print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
