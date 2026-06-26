#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import mean, stdev
from typing import Any

import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[4]
PHASE13_DATASETS = ("Quake_Smart-seq2_Lung", "Mouse_Pancreas_1", "Limb_Muscle")
PHASE13_SEEDS = (42, 2024, 3407)
CORRUPTION_TYPES = ("scmae_shuffle", "matched_donor", "nonzero_aware_donor")
PRIMARY_METRIC = "kmeans_known_k.ari"
PHASE_GATE_RESULT = "pass"
PHASE_GATE_REASON = "formal Phase 13 grid is complete and a Phase 14 corruption recommendation is recorded"
METRIC_KEYS = (
    "kmeans_known_k.acc",
    "kmeans_known_k.nmi",
    "kmeans_known_k.ari",
    "kmeans_known_k.f1_macro",
    "leiden_fixed.acc",
    "leiden_fixed.nmi",
    "leiden_fixed.ari",
    "leiden_fixed.f1_macro",
)
DIAGNOSTIC_KEYS = (
    "zero_to_zero_rate",
    "effective_corruption_rate",
    "budget_deficit_rate",
    "mean_abs_delta",
    "mean_abs_delta_masked",
    "nonzero_aware_success_rate",
    "fallback_to_matched_rate",
    "fallback_to_scmae_shuffle_rate",
)


def load_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_yaml(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def flatten_metrics(metrics: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for key in METRIC_KEYS:
        group, name = key.split(".", 1)
        value = metrics.get(group, {}).get(name)
        out[key] = float(value) if value is not None else math.nan
    return out


def parse_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def require_phase13_protocol(run_dir: Path, row: dict[str, Any]) -> None:
    expected = {
        "variant": "control",
        "encoder_type": "mlp",
        "mask_selector": "random",
        "benchmark_mode": True,
        "input_mode": "log1p",
        "n_top_genes": 2000,
        "scale_input": False,
        "strict_effective_budget": False,
    }
    bad = {key: (row.get(key), value) for key, value in expected.items() if row.get(key) != value}
    if bad:
        raise ValueError(f"{run_dir}: off-protocol Phase 13 run: {bad}")


def discover_runs(root: Path, run_label: str) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    if not root.exists():
        return rows
    for metrics_path in sorted(root.rglob("metrics.json")):
        run_dir = metrics_path.parent
        artifact = load_json(run_dir / "artifact_manifest.json")
        config = load_yaml(run_dir / "resolved_config.yaml")
        corruption = load_json(run_dir / "corruption_stats.json")
        run_manifest = load_json(run_dir / "run_manifest.json")
        run_status = str(run_manifest.get("status", ""))
        artifact_status = str(artifact.get("status", ""))
        if run_status != "complete" or artifact_status != "complete":
            continue
        row: dict[str, Any] = {
            "run_label": run_label,
            "dataset": str(artifact.get("dataset") or config.get("dataset_name") or run_dir.name),
            "corruption_type": str(corruption.get("corruption_type") or config.get("corruption", {}).get("type")),
            "seed": int(artifact.get("seed", config.get("seed", -1))),
            "epochs": int(config.get("training", {}).get("epochs", -1)),
            "variant": str(config.get("variant", "")),
            "encoder_type": str(config.get("model", {}).get("encoder_type", "")),
            "mask_selector": str(config.get("model", {}).get("mask_selector", "")),
            "benchmark_mode": bool(config.get("benchmark_mode", False)),
            "input_mode": str(config.get("preprocessing", {}).get("input_mode", "")),
            "n_top_genes": int(config.get("preprocessing", {}).get("n_top_genes", -1)),
            "scale_input": bool(config.get("preprocessing", {}).get("scale_input", True)),
            "strict_effective_budget": bool(corruption.get("strict_effective_budget", True)),
            "run_status": run_status,
            "artifact_status": artifact_status,
            "run_dir": str(run_dir),
        }
        require_phase13_protocol(run_dir, row)
        row.update(flatten_metrics(load_json(metrics_path)))
        for key in DIAGNOSTIC_KEYS:
            value = corruption.get(key)
            row[key] = float(value) if value is not None else math.nan
        rows.append(row)
    return rows


def validate_expected_grid(
    rows: list[dict[str, Any]],
    *,
    expected_datasets: list[str],
    expected_corruption_types: tuple[str, ...],
    expected_seeds: list[int],
) -> None:
    expected = {
        (dataset, corruption_type, int(seed))
        for dataset in expected_datasets
        for corruption_type in expected_corruption_types
        for seed in expected_seeds
    }
    seen: dict[tuple[str, str, int], int] = {}
    for row in rows:
        key = (str(row["dataset"]), str(row["corruption_type"]), int(row["seed"]))
        seen[key] = seen.get(key, 0) + 1
    actual = set(seen)
    missing = sorted(expected - actual)
    unexpected = sorted(actual - expected)
    duplicates = sorted(key for key, count in seen.items() if count > 1)
    if missing or unexpected or duplicates:
        details = {
            "missing": missing,
            "unexpected": unexpected,
            "duplicates": duplicates,
            "expected_runs": len(expected),
            "complete_runs": len(rows),
        }
        raise ValueError(f"Formal Phase 13 grid is incomplete or contaminated: {details}")


def mean_std(values: list[float]) -> tuple[float, float]:
    clean = [float(value) for value in values if not math.isnan(float(value))]
    if not clean:
        return math.nan, math.nan
    if len(clean) == 1:
        return clean[0], 0.0
    return mean(clean), stdev(clean)


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str], list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault((row["dataset"], row["corruption_type"]), []).append(row)
    out: list[dict[str, Any]] = []
    for (dataset, corruption_type), group_rows in sorted(groups.items()):
        row: dict[str, Any] = {
            "dataset": dataset,
            "corruption_type": corruption_type,
            "n_runs": len(group_rows),
            "seeds": ",".join(str(row["seed"]) for row in sorted(group_rows, key=lambda item: item["seed"])),
        }
        for key in (*METRIC_KEYS, *DIAGNOSTIC_KEYS):
            avg, spread = mean_std([float(item.get(key, math.nan)) for item in group_rows])
            row[f"{key}.mean"] = avg
            row[f"{key}.std"] = spread
        out.append(row)
    return out


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = sorted({key for row in rows for key in row})
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def by_dataset_corruption(aggregate_rows: list[dict[str, Any]]) -> dict[str, dict[str, dict[str, Any]]]:
    out: dict[str, dict[str, dict[str, Any]]] = {}
    for row in aggregate_rows:
        out.setdefault(row["dataset"], {})[row["corruption_type"]] = row
    return out


def build_differences(aggregate_rows: list[dict[str, Any]]) -> dict[str, Any]:
    grouped = by_dataset_corruption(aggregate_rows)
    out: dict[str, Any] = {}
    for dataset, rows in grouped.items():
        if {"scmae_shuffle", "matched_donor"}.issubset(rows):
            out.setdefault("scmae_shuffle_minus_matched_donor", {})[dataset] = (
                rows["scmae_shuffle"][f"{PRIMARY_METRIC}.mean"] - rows["matched_donor"][f"{PRIMARY_METRIC}.mean"]
            )
        if {"nonzero_aware_donor", "scmae_shuffle"}.issubset(rows):
            out.setdefault("nonzero_aware_minus_scmae_shuffle", {})[dataset] = (
                rows["nonzero_aware_donor"][f"{PRIMARY_METRIC}.mean"] - rows["scmae_shuffle"][f"{PRIMARY_METRIC}.mean"]
            )
    return out


def assess_nonzero_aware(aggregate_rows: list[dict[str, Any]]) -> dict[str, Any]:
    grouped = by_dataset_corruption(aggregate_rows)
    beats_scmae: list[str] = []
    improves_effective_rate: list[str] = []
    for dataset, rows in grouped.items():
        if {"nonzero_aware_donor", "scmae_shuffle"}.issubset(rows):
            if rows["nonzero_aware_donor"][f"{PRIMARY_METRIC}.mean"] > rows["scmae_shuffle"][f"{PRIMARY_METRIC}.mean"]:
                beats_scmae.append(dataset)
            if rows["nonzero_aware_donor"]["effective_corruption_rate.mean"] > rows["scmae_shuffle"]["effective_corruption_rate.mean"]:
                improves_effective_rate.append(dataset)
    return {
        "beats_scmae_on_primary_metric": beats_scmae,
        "improves_effective_corruption_rate": improves_effective_rate,
        "assessment": (
            "improves diagnostics but not consistently clustering"
            if len(improves_effective_rate) >= 2 and len(beats_scmae) < 2
            else "clustering improvement supported on development datasets"
            if len(beats_scmae) >= 2
            else "no consistent diagnostic or clustering advantage"
        ),
    }


def recommend(aggregate_rows: list[dict[str, Any]]) -> dict[str, Any]:
    grouped = by_dataset_corruption(aggregate_rows)
    complete = all(set(rows) >= set(CORRUPTION_TYPES) for rows in grouped.values())
    if not complete:
        return {"recommended_corruption_type": None, "reason": "formal results are incomplete"}

    means_by_type: dict[str, list[float]] = {key: [] for key in CORRUPTION_TYPES}
    stds_by_type: dict[str, list[float]] = {key: [] for key in CORRUPTION_TYPES}
    matched_weaker_count = 0
    nonzero_beats_scmae_count = 0
    for rows in grouped.values():
        for corruption_type in CORRUPTION_TYPES:
            means_by_type[corruption_type].append(float(rows[corruption_type][f"{PRIMARY_METRIC}.mean"]))
            stds_by_type[corruption_type].append(float(rows[corruption_type][f"{PRIMARY_METRIC}.std"]))
        if rows["matched_donor"][f"{PRIMARY_METRIC}.mean"] < rows["scmae_shuffle"][f"{PRIMARY_METRIC}.mean"]:
            matched_weaker_count += 1
        if rows["nonzero_aware_donor"][f"{PRIMARY_METRIC}.mean"] > rows["scmae_shuffle"][f"{PRIMARY_METRIC}.mean"]:
            nonzero_beats_scmae_count += 1

    global_means = {key: mean(values) for key, values in means_by_type.items()}
    avg_seed_std = mean([value for values in stds_by_type.values() for value in values])
    spread = max(global_means.values()) - min(global_means.values())
    if spread <= avg_seed_std:
        return {
            "recommended_corruption_type": "scmae_shuffle",
            "reason": "corruption differences are within average seed variability; choose the simplest baseline",
            "global_primary_metric_means": global_means,
        }
    if nonzero_beats_scmae_count >= 2 and global_means["nonzero_aware_donor"] >= global_means["scmae_shuffle"]:
        return {
            "recommended_corruption_type": "nonzero_aware_donor",
            "reason": "nonzero-aware donor improves primary clustering metric over scMAE shuffle on at least 2/3 datasets",
            "global_primary_metric_means": global_means,
        }
    if matched_weaker_count >= 2:
        return {
            "recommended_corruption_type": "scmae_shuffle",
            "reason": "matched donor is weaker than scMAE shuffle on at least 2/3 datasets",
            "global_primary_metric_means": global_means,
        }
    best = max(global_means, key=global_means.get)
    return {
        "recommended_corruption_type": best,
        "reason": "selected by highest average primary clustering metric across Phase 13 datasets",
        "global_primary_metric_means": global_means,
    }


def write_markdown_report(
    path: Path,
    *,
    smoke_rows: list[dict[str, Any]],
    formal_rows: list[dict[str, Any]],
    aggregate_rows: list[dict[str, Any]],
    differences: dict[str, Any],
    nonzero_aware_assessment: dict[str, Any],
    recommendation: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Phase 13 Corruption Triad Report",
        "",
        "## 1. Smoke validation results",
        "",
        f"- Smoke run count: {len(smoke_rows)}",
        "- Smoke runs are implementation checks only and are not used as scientific evidence.",
        "",
        "## 2. Formal Phase 13 results",
        "",
        f"- Formal run count: {len(formal_rows)}",
        f"- Primary metric: `{PRIMARY_METRIC}`",
        "",
        "| dataset | corruption | ARI mean | ARI std | NMI mean | ACC mean | F1 mean | zero_to_zero mean | effective mean | mean_abs_delta |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            "| {dataset} | {corruption_type} | {ari:.6f} | {ari_std:.6f} | {nmi:.6f} | {acc:.6f} | {f1:.6f} | {zero:.6f} | {eff:.6f} | {delta:.6f} |".format(
                dataset=row["dataset"],
                corruption_type=row["corruption_type"],
                ari=float(row["kmeans_known_k.ari.mean"]),
                ari_std=float(row["kmeans_known_k.ari.std"]),
                nmi=float(row["kmeans_known_k.nmi.mean"]),
                acc=float(row["kmeans_known_k.acc.mean"]),
                f1=float(row["kmeans_known_k.f1_macro.mean"]),
                zero=float(row["zero_to_zero_rate.mean"]),
                eff=float(row["effective_corruption_rate.mean"]),
                delta=float(row["mean_abs_delta.mean"]),
            )
        )
    lines.extend(
        [
            "",
            "## 3. Phase gate",
            "",
            f"- gate_result: `{PHASE_GATE_RESULT}`",
            f"- gate_reason: {PHASE_GATE_REASON}",
            "",
            "## 4. Corruption recommendation for Phase 14",
            "",
            f"- Recommendation: `{recommendation.get('recommended_corruption_type')}`",
            f"- Reason: {recommendation.get('reason')}",
            f"- Differences: `{json.dumps(differences, sort_keys=True)}`",
            f"- Nonzero-aware assessment: `{json.dumps(nonzero_aware_assessment, sort_keys=True)}`",
            "",
            "## 5. Remaining risks",
            "",
            "- Phase 13 is limited to MLP encoder plus random mask; it does not validate AdvMask, Axial, or full CAAM.",
            "- A corruption that improves mask diagnostics without improving clustering must not be claimed as a main method contribution.",
            "- Results must stay development-only; validation and sealed test are not used in Phase 13.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize Phase 13 CAAM corruption triad runs.")
    parser.add_argument("--formal_root", type=Path, required=True)
    parser.add_argument("--smoke_root", type=Path, default=None)
    parser.add_argument("--expected_datasets", type=str, default=",".join(PHASE13_DATASETS))
    parser.add_argument("--expected_seeds", type=str, default=",".join(str(seed) for seed in PHASE13_SEEDS))
    parser.add_argument("--output_dir", type=Path, default=PROJECT_ROOT / "results/CAAM_scMAE_correction")
    parser.add_argument(
        "--report_path",
        type=Path,
        default=PROJECT_ROOT / "methods/DeepLearning/CAAM_scMAE/benchmark/PHASE13_CORRUPTION_TRIAD_REPORT.md",
    )
    args = parser.parse_args()

    formal_rows = discover_runs(args.formal_root, "formal")
    smoke_rows = discover_runs(args.smoke_root, "smoke") if args.smoke_root is not None else []
    validate_expected_grid(
        formal_rows,
        expected_datasets=parse_csv(args.expected_datasets),
        expected_corruption_types=CORRUPTION_TYPES,
        expected_seeds=[int(seed) for seed in parse_csv(args.expected_seeds)],
    )
    aggregate_rows = aggregate(formal_rows)
    differences = build_differences(aggregate_rows)
    nonzero_aware_assessment = assess_nonzero_aware(aggregate_rows)
    recommendation = recommend(aggregate_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "corruption_triad_runs.csv", formal_rows)
    write_csv(args.output_dir / "corruption_triad_summary.csv", aggregate_rows)
    if smoke_rows:
        write_csv(args.output_dir / "corruption_triad_smoke_runs.csv", smoke_rows)
    report = {
        "gate_result": PHASE_GATE_RESULT,
        "gate_reason": PHASE_GATE_REASON,
        "primary_metric": PRIMARY_METRIC,
        "n_smoke_runs": len(smoke_rows),
        "n_formal_runs": len(formal_rows),
        "aggregate": aggregate_rows,
        "differences": differences,
        "nonzero_aware_assessment": nonzero_aware_assessment,
        "recommendation": recommendation,
        "smoke_is_scientific_evidence": False,
    }
    with open(args.output_dir / "corruption_triad_report.json", "w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    write_markdown_report(
        args.report_path,
        smoke_rows=smoke_rows,
        formal_rows=formal_rows,
        aggregate_rows=aggregate_rows,
        differences=differences,
        nonzero_aware_assessment=nonzero_aware_assessment,
        recommendation=recommendation,
    )
    print(f"Wrote {args.output_dir / 'corruption_triad_summary.csv'}")
    print(f"Wrote {args.output_dir / 'corruption_triad_report.json'}")
    print(f"Wrote {args.report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
