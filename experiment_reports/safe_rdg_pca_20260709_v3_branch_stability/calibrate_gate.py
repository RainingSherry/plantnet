#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression


BASE_DIR = Path(__file__).resolve().parent
STAGE_VARIANTS = ["pca_kmeans", "pca_spectral_kmeans", "rdg_cell_only", "rdg_gene_only", "rdg_concat_kmeans", "rdg_always_on"]


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def collect_runs(run_root: Path) -> pd.DataFrame:
    rows = []
    for metrics_path in run_root.glob("*/*/*/metrics.json"):
        variant = metrics_path.parent.name
        seed_dir = metrics_path.parent.parent.name
        dataset = metrics_path.parent.parent.parent.name
        seed = int(seed_dir.replace("seed_", ""))
        metrics = read_json(metrics_path)
        gate_path = metrics_path.parent / "gate_decision.json"
        diag_path = metrics_path.parent / "diagnostics.json"
        gate = read_json(gate_path) if gate_path.exists() else {}
        diag = read_json(diag_path) if diag_path.exists() else {}
        gene_diag = (diag.get("graphs") or {}).get("A_gene_graph") or {}
        cell_diag = (diag.get("graphs") or {}).get("A_cell_reliable") or {}
        rows.append(
            {
                "dataset": dataset,
                "seed": seed,
                "variant": variant,
                "ari": float(metrics.get("ari", np.nan)),
                "nmi": float(metrics.get("nmi", np.nan)),
                "acc": float(metrics.get("acc", np.nan)),
                "q_cell": float(gate.get("q_cell", np.nan)),
                "q_gene": float(gate.get("q_gene", np.nan)),
                "q_total": float(gate.get("q_total", np.nan)),
                "cell_spectral_gap": float(cell_diag.get("spectral_gap_proxy", np.nan)),
                "hubness": float(cell_diag.get("degree_max", 0.0)) / max(float(cell_diag.get("degree_mean", 0.0)), 1e-8),
                "gene_module_entropy": float(gene_diag.get("module_size_entropy", np.nan)),
            }
        )
    return pd.DataFrame(rows)


def paired_base(df: pd.DataFrame) -> pd.DataFrame:
    wide = df.pivot_table(index=["dataset", "seed"], columns="variant", values=["ari", "nmi", "acc", "q_cell", "q_gene", "q_total", "cell_spectral_gap", "hubness", "gene_module_entropy"], aggfunc="first")
    wide.columns = [f"{a}__{b}" for a, b in wide.columns]
    wide = wide.reset_index()
    for metric in ["ari", "nmi", "acc"]:
        wide[f"delta_{metric}_always_vs_pca"] = wide[f"{metric}__rdg_always_on"] - wide[f"{metric}__pca_kmeans"]
    return wide


def choose_threshold(train: pd.DataFrame, lambda_regret: float) -> float:
    qs = np.unique(np.nan_to_num(train["q_total__rdg_always_on"].to_numpy(dtype=float), nan=0.0))
    candidates = np.unique(np.concatenate(([0.0, 1.0], qs)))
    best_tau = 1.0
    best_utility = -1e9
    for tau in candidates:
        enabled = train["q_total__rdg_always_on"] >= tau
        selected = np.where(enabled, train["ari__rdg_always_on"], train["ari__pca_kmeans"])
        delta = selected - train["ari__pca_kmeans"]
        regret = np.maximum(0.0, -delta)
        utility = float(np.mean(delta) - lambda_regret * np.mean(regret))
        if utility > best_utility:
            best_utility = utility
            best_tau = float(tau)
    return best_tau


def synthesize_threshold(wide: pd.DataFrame, lambda_regret: float) -> pd.DataFrame:
    out = []
    for dataset in sorted(wide["dataset"].unique()):
        train = wide[wide["dataset"] != dataset]
        test = wide[wide["dataset"] == dataset]
        tau = choose_threshold(train, lambda_regret)
        for _, row in test.iterrows():
            enabled = bool(row["q_total__rdg_always_on"] >= tau)
            source = "rdg_always_on" if enabled else "pca_kmeans"
            out.append(
                {
                    "dataset": row["dataset"],
                    "seed": int(row["seed"]),
                    "variant": "safe_rdg_calibrated_threshold",
                    "source_variant": source,
                    "graph_enabled": enabled,
                    "tau": tau,
                    "ari": float(row[f"ari__{source}"]),
                    "nmi": float(row[f"nmi__{source}"]),
                    "acc": float(row[f"acc__{source}"]),
                    "q_total": float(row["q_total__rdg_always_on"]),
                }
            )
    return pd.DataFrame(out)


def synthesize_logistic(wide: pd.DataFrame) -> pd.DataFrame:
    features = ["q_cell__rdg_always_on", "q_gene__rdg_always_on", "q_total__rdg_always_on", "cell_spectral_gap__rdg_always_on", "gene_module_entropy__rdg_always_on"]
    out = []
    for dataset in sorted(wide["dataset"].unique()):
        train = wide[wide["dataset"] != dataset].copy()
        test = wide[wide["dataset"] == dataset].copy()
        x_train = np.nan_to_num(train[features].to_numpy(dtype=float), nan=0.0)
        y_train = (train["delta_ari_always_vs_pca"] > 0).astype(int).to_numpy()
        if len(np.unique(y_train)) < 2:
            probs = np.zeros(len(test))
        else:
            clf = LogisticRegression(C=1.0, max_iter=1000)
            clf.fit(x_train, y_train)
            probs = clf.predict_proba(np.nan_to_num(test[features].to_numpy(dtype=float), nan=0.0))[:, 1]
        for (_, row), prob in zip(test.iterrows(), probs):
            enabled = bool(prob >= 0.5)
            source = "rdg_always_on" if enabled else "pca_kmeans"
            out.append(
                {
                    "dataset": row["dataset"],
                    "seed": int(row["seed"]),
                    "variant": "safe_rdg_calibrated_logistic_explore",
                    "source_variant": source,
                    "graph_enabled": enabled,
                    "prob_graph_gain": float(prob),
                    "ari": float(row[f"ari__{source}"]),
                    "nmi": float(row[f"nmi__{source}"]),
                    "acc": float(row[f"acc__{source}"]),
                    "q_total": float(row["q_total__rdg_always_on"]),
                }
            )
    return pd.DataFrame(out)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--run_root", type=Path, default=BASE_DIR / "runs")
    parser.add_argument("--out_dir", type=Path, default=BASE_DIR)
    parser.add_argument("--lambda_regret", type=float, default=1.0)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    df = collect_runs(args.run_root)
    if df.empty:
        raise SystemExit(f"No Stage A runs found under {args.run_root}")
    wide = paired_base(df)
    threshold = synthesize_threshold(wide, args.lambda_regret)
    logistic = synthesize_logistic(wide)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    wide.to_csv(args.out_dir / "paired_stage_a.csv", index=False)
    threshold.to_csv(args.out_dir / "calibrated_threshold_runs.csv", index=False)
    logistic.to_csv(args.out_dir / "calibrated_logistic_runs.csv", index=False)
    print(f"Wrote {len(threshold)} threshold and {len(logistic)} logistic calibrated rows.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
