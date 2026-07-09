#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


STAGE_A_VARIANTS = [
    "pca_kmeans",
    "pca_spectral_kmeans",
    "rdg_cell_only",
    "rdg_gene_only",
    "rdg_concat_kmeans",
    "rdg_always_on",
    "safe_rdg_heuristic",
]
NEGATIVE_CONTROL_VARIANTS = [
    "neg_random_cell_graph",
    "neg_degree_shuffle_graph",
    "neg_shuffled_gene_cell_graph",
]
RUN_FILES = [
    "metrics.json",
    "eval_fixed.csv",
    "embedding_final.npy",
    "labels.npy",
    "diagnostics.json",
    "gate_decision.json",
    "args.json",
    "preprocess_config.json",
]
METRIC_KEYS = ["ari", "nmi", "acc"]
GATE_KEYS = ["q_cell", "q_gene", "q_total", "graph_enabled"]
SUMMARY_COLUMNS = [
    "variant",
    "ari_mean",
    "ari_median",
    "nmi_mean",
    "acc_mean",
    "rank_ari_mean",
    "wins_vs_pca",
    "ties_vs_pca",
    "losses_vs_pca",
    "delta_ari_vs_pca_mean",
    "regret_vs_pca_mean",
    "oracle_gap_mean",
    "negative_transfer_rate",
    "graph_activation_rate",
    "enabled_graph_mean_gain",
]


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def parse_csv(value: str | None) -> list[str] | None:
    if value is None or value == "":
        return None
    if value in {"auto", "all"}:
        return None
    return [x.strip() for x in value.split(",") if x.strip()]


def fail(errors: list[str], message: str) -> None:
    errors.append(message)


def check_run_dir(run_dir: Path, dataset: str, seed: int, variants: list[str], errors: list[str]) -> None:
    status_path = run_dir / "status.json"
    if not status_path.exists():
        fail(errors, f"missing status.json: {run_dir}")
    else:
        status = read_json(status_path)
        if status.get("status") != "success":
            fail(errors, f"non-success status for {dataset} seed={seed}: {status.get('status')} {status.get('error', '')[:160]}")

    for variant in variants:
        variant_dir = run_dir / variant
        if not variant_dir.exists():
            fail(errors, f"missing variant dir: {variant_dir}")
            continue
        for filename in RUN_FILES:
            if not (variant_dir / filename).exists():
                fail(errors, f"missing {filename}: {variant_dir}")

        metrics_path = variant_dir / "metrics.json"
        if metrics_path.exists():
            metrics = read_json(metrics_path)
            for key in METRIC_KEYS:
                val = metrics.get(key)
                if val is None or not np.isfinite(float(val)):
                    fail(errors, f"bad metric {key} in {metrics_path}: {val}")
            if metrics.get("variant") != variant:
                fail(errors, f"variant mismatch in {metrics_path}: {metrics.get('variant')} != {variant}")

        gate_path = variant_dir / "gate_decision.json"
        if gate_path.exists():
            gate = read_json(gate_path)
            for key in GATE_KEYS:
                if key not in gate:
                    fail(errors, f"missing gate key {key} in {gate_path}")
            for key in ["q_cell", "q_gene", "q_total"]:
                if key in gate and not (0.0 <= float(gate[key]) <= 1.0):
                    fail(errors, f"gate score out of [0,1] for {key} in {gate_path}: {gate[key]}")
            if variant in NEGATIVE_CONTROL_VARIANTS and not gate.get("negative_control"):
                fail(errors, f"negative control marker absent in {gate_path}")

        emb_path = variant_dir / "embedding_final.npy"
        labels_path = variant_dir / "labels.npy"
        if emb_path.exists() and labels_path.exists():
            emb = np.load(emb_path, mmap_mode="r")
            labels = np.load(labels_path, mmap_mode="r")
            if emb.shape[0] != labels.shape[0]:
                fail(errors, f"embedding/label row mismatch in {variant_dir}: {emb.shape[0]} != {labels.shape[0]}")


def discover_dataset_seed_dirs(run_root: Path) -> list[tuple[str, int, Path]]:
    rows = []
    for dataset_dir in sorted(p for p in run_root.iterdir() if p.is_dir()):
        for seed_dir in sorted(p for p in dataset_dir.iterdir() if p.is_dir() and p.name.startswith("seed_")):
            try:
                seed = int(seed_dir.name.replace("seed_", ""))
            except ValueError:
                continue
            rows.append((dataset_dir.name, seed, seed_dir))
    return rows


def check_analysis(analysis_dir: Path, include_negative_controls: bool, errors: list[str]) -> None:
    required = [
        "paired_stage_a.csv",
        "calibrated_threshold_runs.csv",
        "calibrated_logistic_runs.csv",
        "all_runs.csv",
        "summary_by_variant.csv",
        "summary_by_variant_dataset.csv",
        "oracle_best_runs.csv",
        "summary_report.json",
    ]
    for filename in required:
        if not (analysis_dir / filename).exists():
            fail(errors, f"missing analysis file: {analysis_dir / filename}")
    summary_path = analysis_dir / "summary_by_variant.csv"
    if summary_path.exists():
        summary = pd.read_csv(summary_path)
        missing = [col for col in SUMMARY_COLUMNS if col not in summary.columns]
        if missing:
            fail(errors, f"summary missing columns {missing}: {summary_path}")
        expected = set(STAGE_A_VARIANTS + ["safe_rdg_calibrated_threshold", "safe_rdg_calibrated_logistic_explore"])
        if include_negative_controls:
            expected.update(NEGATIVE_CONTROL_VARIANTS)
        absent = sorted(expected - set(summary["variant"].astype(str)))
        if absent:
            fail(errors, f"summary missing variants {absent}: {summary_path}")
    report_path = analysis_dir / "summary_report.json"
    if report_path.exists():
        report = read_json(report_path)
        for key in ["q_correlations", "negative_controls", "success"]:
            if key not in report:
                fail(errors, f"summary_report missing key {key}: {report_path}")
        if include_negative_controls and not report.get("negative_controls", {}).get("available"):
            fail(errors, f"negative_controls report not available: {report_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--run_root", type=Path, required=True)
    parser.add_argument("--analysis_dir", type=Path, default=None)
    parser.add_argument("--datasets", default="auto")
    parser.add_argument("--seeds", default="auto")
    parser.add_argument("--include_negative_controls", action="store_true")
    parser.add_argument("--require_analysis", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    errors: list[str] = []
    variants = STAGE_A_VARIANTS + (NEGATIVE_CONTROL_VARIANTS if args.include_negative_controls else [])

    if not args.run_root.exists():
        raise SystemExit(f"run_root does not exist: {args.run_root}")
    rows = discover_dataset_seed_dirs(args.run_root)
    wanted_datasets = parse_csv(args.datasets)
    wanted_seeds = None if args.seeds == "auto" else [int(x) for x in parse_csv(args.seeds)]
    if wanted_datasets is not None:
        rows = [row for row in rows if row[0] in wanted_datasets]
    if wanted_seeds is not None:
        rows = [row for row in rows if row[1] in wanted_seeds]
    if not rows:
        fail(errors, f"no dataset/seed dirs found under {args.run_root}")

    for dataset, seed, run_dir in rows:
        check_run_dir(run_dir, dataset, seed, variants, errors)

    if args.require_analysis:
        analysis_dir = args.analysis_dir or args.run_root
        check_analysis(analysis_dir, args.include_negative_controls, errors)

    if errors:
        print("VERIFY FAILED")
        for err in errors:
            print(f"- {err}")
        return 1
    print(
        "VERIFY OK "
        f"run_root={args.run_root} runs={len(rows)} variants={len(variants)} "
        f"analysis_checked={bool(args.require_analysis)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
