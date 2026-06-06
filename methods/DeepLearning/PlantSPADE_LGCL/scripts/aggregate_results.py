#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
PKG_DIR = SCRIPT_DIR.parent
ROOT = SCRIPT_DIR.parents[3]

METRICS = [
    "acc",
    "nmi",
    "ari",
    "f1_macro",
    "fmi",
    "v_measure",
    "homogeneity",
    "completeness",
    "n_pred_clusters",
    "silhouette",
]


def parse_args():
    parser = argparse.ArgumentParser(description="Aggregate PlantSPADE-LGCL protocol outputs.")
    parser.add_argument("--results_dir", default=None)
    parser.add_argument("--main_config", default=str(PKG_DIR / "configs" / "main_lgcl.yaml"))
    parser.add_argument("--output_dir", default=None)
    return parser.parse_args()


def load_default_results_dir(main_config: str) -> Path:
    with open(main_config, "r", encoding="utf-8") as handle:
        cfg = yaml.safe_load(handle) or {}
    return Path(cfg.get("output_dir", ROOT / "results" / "PlantSPADE_LGCL_protocol"))


def read_metric_csvs(results_dir: Path) -> pd.DataFrame:
    patterns = ["**/*_fixed.csv", "**/*_oracle.csv", "**/*_sweep.csv"]
    frames = []
    legacy_replacements = {
        "plantspade_lgcl_fixed.csv": "eval_baseline_fixed.csv",
        "plantspade_lgcl_sga_fixed.csv": "eval_support_attention_fixed.csv",
    }
    for pattern in patterns:
        for path in sorted(results_dir.glob(pattern)):
            if path.name.startswith("table_") or path.name.startswith("all_results"):
                continue
            if not path.name.startswith(("eval_", "pca_", "external_eval_", "sc3_eval_", "plantspade_lgcl_")):
                continue
            replacement = legacy_replacements.get(path.name)
            if replacement and (path.parent / replacement).exists():
                continue
            try:
                df = pd.read_csv(path)
            except Exception:
                continue
            if df.empty:
                continue
            df["source_file"] = str(path)
            if "protocol" not in df.columns:
                if path.name.endswith("_fixed.csv"):
                    df["protocol"] = "fixed"
                elif path.name.endswith("_oracle.csv"):
                    df["protocol"] = "oracle"
                else:
                    df["protocol"] = "full_sweep"
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def mean_std_table(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    group_cols = [col for col in ["dataset", "method", "variant", "negative_sampler", "cluster_method", "protocol"] if col in df.columns]
    numeric = [metric for metric in METRICS if metric in df.columns]
    grouped = df.groupby(group_cols, dropna=False)
    rows = []
    for keys, block in grouped:
        if not isinstance(keys, tuple):
            keys = (keys,)
        row = dict(zip(group_cols, keys))
        row["n_runs"] = int(block["seed"].nunique()) if "seed" in block.columns else int(len(block))
        for metric in numeric:
            values = pd.to_numeric(block[metric], errors="coerce").to_numpy(dtype=float)
            mean = float(np.nanmean(values)) if np.any(np.isfinite(values)) else float("nan")
            std = float(np.nanstd(values, ddof=1)) if np.sum(np.isfinite(values)) > 1 else 0.0
            row[f"{metric}_mean"] = mean
            row[f"{metric}_std"] = std
            row[metric] = f"{mean:.4f} ± {std:.4f}" if np.isfinite(mean) else "nan"
        rows.append(row)
    return pd.DataFrame(rows)


def text_col(df: pd.DataFrame, column: str) -> pd.Series:
    return df.get(column, pd.Series("", index=df.index)).fillna("").astype(str)


def collect_profiles(results_dir: Path) -> pd.DataFrame:
    rows = []
    for path in sorted(results_dir.glob("**/dataset_profile.json")):
        try:
            with open(path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            payload["source_file"] = str(path)
            rows.append(payload)
        except Exception:
            continue
    if not rows:
        return pd.DataFrame()
    df = pd.json_normalize(rows)
    if "dataset_name" in df.columns:
        df = df.drop_duplicates(subset=["dataset_name"], keep="last")
    return df


def main():
    args = parse_args()
    results_dir = Path(args.results_dir) if args.results_dir else load_default_results_dir(args.main_config)
    output_dir = Path(args.output_dir) if args.output_dir else results_dir / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    long_df = read_metric_csvs(results_dir)
    long_df.to_csv(output_dir / "all_results_long.csv", index=False)

    mean_std = mean_std_table(long_df)
    mean_std.to_csv(output_dir / "all_results_mean_std.csv", index=False)

    fixed_all = long_df[long_df.get("protocol", pd.Series(dtype=str)).eq("fixed")].copy() if not long_df.empty else pd.DataFrame()
    if not fixed_all.empty:
        method_text = text_col(fixed_all, "method")
        variant_text = text_col(fixed_all, "variant")
        is_ablation = (
            method_text.str.contains("attention_no|topk_|neg_|neighbor_conflict|idf_weighted", case=False, na=False)
            | variant_text.str.contains("attention_no|topk_|neg_|neighbor_conflict|idf_weighted", case=False, na=False)
        )
        fixed_long = fixed_all[~is_ablation].copy()
        fixed_long["method"] = fixed_long["method"].replace(
            {
                "plantspade_lgcl": "plantspade_lgcl_baseline",
                "plantspade_lgcl_sga": "plantspade_lgcl_support_attention",
                "plantspade_lgcl_sga_support_attention": "plantspade_lgcl_support_attention",
            }
        )
        if "variant" in fixed_long.columns:
            variant_defaults = {
                "phytocluster": "embedding",
                "scmae": "embedding",
                "scvi": "embedding",
                "traditional_pca": "pca",
                "traditional_leiden": "leiden",
                "traditional_louvain": "louvain",
                "plantspade_lgcl_baseline": "baseline",
                "plantspade_lgcl_support_attention": "support_attention",
                "plantspade_lgcl_gated_fusion": "gated_fusion",
            }
            missing_variant = fixed_long["variant"].isna() | fixed_long["variant"].astype(str).eq("")
            fixed_long.loc[missing_variant, "variant"] = fixed_long.loc[missing_variant, "method"].map(variant_defaults)
        fixed_long = fixed_long[fixed_long["method"].astype(str).ne("plantspade_lgcl_sga_baseline")]
        if "negative_sampler" in fixed_long.columns:
            fixed_long["negative_sampler"] = ""
        dedup_cols = [
            col
            for col in ["dataset", "method", "seed", "variant", "cluster_method", "protocol"]
            if col in fixed_long.columns
        ]
        if dedup_cols:
            fixed_long = fixed_long.drop_duplicates(subset=dedup_cols, keep="last")
        fixed = mean_std_table(fixed_long)
    else:
        fixed = pd.DataFrame()
    fixed.to_csv(output_dir / "table_main_fixed_protocol.csv", index=False)

    oracle = long_df[long_df.get("protocol", pd.Series(dtype=str)).eq("oracle")].copy() if not long_df.empty else pd.DataFrame()
    oracle.to_csv(output_dir / "table_oracle_supplement.csv", index=False)

    if not long_df.empty:
        att_fixed = long_df[long_df.get("protocol", pd.Series(dtype=str)).eq("fixed")].copy()
        att_fixed_text = text_col(att_fixed, "method")
        att_fixed_variant = text_col(att_fixed, "variant")
        attention = mean_std_table(att_fixed[
            att_fixed_text.str.contains("attention|sga|support_attention", case=False, na=False)
            | att_fixed_variant.str.contains("attention|support_attention", case=False, na=False)
        ])
        neg_fixed = long_df[long_df.get("protocol", pd.Series(dtype=str)).eq("fixed")].copy()
        neg_fixed_text = text_col(neg_fixed, "method")
        neg_sampler = text_col(neg_fixed, "negative_sampler")
        negative = mean_std_table(neg_fixed[
            neg_fixed_text.str.contains("neg_|negative", case=False, na=False)
            | neg_sampler.isin(["idf_weighted_zero", "neighbor_conflict_zero"])
        ])
    else:
        attention = pd.DataFrame()
        negative = pd.DataFrame()
    attention.to_csv(output_dir / "table_attention_ablation.csv", index=False)
    negative.to_csv(output_dir / "table_negative_sampling_ablation.csv", index=False)

    profiles = collect_profiles(results_dir)
    profiles.to_csv(output_dir / "dataset_profiles_summary.csv", index=False)
    print(f"Aggregated results under {output_dir}")


if __name__ == "__main__":
    main()
