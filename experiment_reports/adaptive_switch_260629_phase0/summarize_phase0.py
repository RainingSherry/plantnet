#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
BASE = Path(__file__).resolve().parent
REFERENCE = ROOT / "results" / "260629全benchmark结果.csv"
METHOD_RAW = "adaptive_switch_scmae"
METHOD = "scMAE + DEC + StdFloor"


def load_metrics(run_dir: Path) -> dict:
    metrics_path = run_dir / "metrics.json"
    summary_path = run_dir / "summary.json"
    args_path = run_dir / "args.json"
    if not metrics_path.exists():
        return {}

    metrics = json.loads(metrics_path.read_text())
    fixed = metrics.get("kmeans_known_k", metrics)
    summary = json.loads(summary_path.read_text()) if summary_path.exists() else {}
    args = json.loads(args_path.read_text()) if args_path.exists() else {}
    return {
        "dataset": args.get("dataset_name") or summary.get("dataset") or run_dir.parent.name,
        "method": METHOD,
        "method_raw": METHOD_RAW,
        "seed": int(args.get("seed", run_dir.name.replace("seed", "").replace("_gpu", ""))),
        "run_dir": str(run_dir),
        "acc": fixed.get("acc"),
        "nmi": fixed.get("nmi"),
        "ari": fixed.get("ari"),
        "f1_macro": fixed.get("f1_macro"),
        "runtime_seconds": summary.get("runtime_seconds"),
        "final_gate": summary.get("final_gate"),
        "final_kl_ref": summary.get("final_kl_ref"),
    }


def collect_runs() -> pd.DataFrame:
    rows = []
    priority = {"seed42_gpu2": 3, "seed42_gpu": 2, "seed42": 1}
    for metrics_path in sorted(BASE.glob("*/*/metrics.json")):
        run_dir = metrics_path.parent
        seed_name = run_dir.name.split("_")[0]
        siblings = [
            p
            for p in run_dir.parent.iterdir()
            if p.is_dir() and p.name.split("_")[0] == seed_name and (p / "metrics.json").exists()
        ]
        best = max(siblings, key=lambda p: priority.get(p.name, 0))
        if run_dir != best:
            continue
        row = load_metrics(run_dir)
        if row:
            rows.append(row)
    return pd.DataFrame(rows)


def summarize_runs(runs: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for dataset, sub in runs.groupby("dataset", sort=True):
        row = {"dataset": dataset, "method": METHOD, "method_raw": METHOD_RAW}
        seeds = sorted(str(int(s)) for s in sub["seed"].dropna().unique())
        for metric in ["acc", "nmi", "ari", "f1_macro", "runtime_seconds"]:
            vals = pd.to_numeric(sub[metric], errors="coerce").dropna()
            row[f"{metric}_mean"] = vals.mean() if len(vals) else pd.NA
            row[f"{metric}_std"] = vals.std(ddof=1) if len(vals) > 1 else pd.NA
            row[f"{metric}_count"] = int(len(vals))
        row["seed_<lambda>"] = ",".join(seeds)
        row["n_seeds"] = len(seeds)
        rows.append(row)
    return pd.DataFrame(rows)


def build_comparison(summary: pd.DataFrame) -> pd.DataFrame:
    ref = pd.read_csv(REFERENCE)
    baselines = ref[ref["method_raw"].isin(["scmae", "neighbormix_scmae", "nm_scmae_nomix"])]
    merged = summary.merge(
        baselines[
            [
                "dataset",
                "method",
                "method_raw",
                "ari_mean",
                "ari_std",
                "nmi_mean",
                "nmi_std",
                "n_seeds",
            ]
        ],
        on="dataset",
        how="left",
        suffixes=("_adaptive", "_baseline"),
    )
    merged["delta_ari"] = merged["ari_mean_adaptive"] - merged["ari_mean_baseline"]
    merged["delta_nmi"] = merged["nmi_mean_adaptive"] - merged["nmi_mean_baseline"]
    return merged.sort_values(["dataset", "method_raw_baseline"])


def main() -> int:
    runs = collect_runs()
    if runs.empty:
        raise SystemExit(f"No completed metrics.json files found under {BASE}")
    summary = summarize_runs(runs)
    comparison = build_comparison(summary)

    runs.to_csv(BASE / "adaptive_switch_phase0_runs.csv", index=False)
    summary.to_csv(BASE / "adaptive_switch_phase0_summary.csv", index=False)
    comparison.to_csv(BASE / "adaptive_switch_phase0_vs_260629.csv", index=False)

    print("Runs:")
    print(runs[["dataset", "seed", "ari", "nmi", "run_dir"]].to_string(index=False))
    print("\nSummary:")
    print(summary[["dataset", "ari_mean", "ari_std", "nmi_mean", "nmi_std", "n_seeds"]].to_string(index=False))
    print("\nComparison:")
    cols = ["dataset", "method_baseline", "ari_mean_adaptive", "ari_mean_baseline", "delta_ari", "nmi_mean_adaptive", "nmi_mean_baseline", "delta_nmi"]
    print(comparison[cols].to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
