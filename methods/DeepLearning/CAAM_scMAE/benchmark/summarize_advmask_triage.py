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
PHASE14_DATASETS = ("Quake_Smart-seq2_Lung", "Mouse_Pancreas_1", "Limb_Muscle")
PHASE14_SEEDS = (42, 2024, 3407)
PHASE14_VARIANTS = ("control", "advmask")
PRIMARY_METRIC = "kmeans_known_k.ari"
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
)
MASK_KEYS = ("mask_entropy", "mask_gini", "top_gene_concentration")
EMBEDDING_KEYS = ("effective_rank", "mean_pairwise_cosine", "min_per_dimension_variance")


def load_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def load_yaml(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def parse_csv(value: str) -> list[str]:
    return [part.strip() for part in value.split(",") if part.strip()]


def flatten_metrics(metrics: dict[str, Any]) -> dict[str, float]:
    out: dict[str, float] = {}
    for key in METRIC_KEYS:
        group, name = key.split(".", 1)
        value = metrics.get(group, {}).get(name)
        out[key] = float(value) if value is not None else math.nan
    return out


def require_phase14_protocol(run_dir: Path, row: dict[str, Any]) -> None:
    expected_common = {
        "encoder_type": "mlp",
        "benchmark_mode": True,
        "input_mode": "log1p",
        "n_top_genes": 2000,
        "scale_input": False,
        "strict_effective_budget": False,
    }
    bad = {key: (row.get(key), value) for key, value in expected_common.items() if row.get(key) != value}
    variant = row.get("variant")
    if variant == "control" and row.get("mask_selector") != "random":
        bad["mask_selector"] = (row.get("mask_selector"), "random")
    elif variant == "advmask" and row.get("mask_selector") != "adversarial":
        bad["mask_selector"] = (row.get("mask_selector"), "adversarial")
    elif variant not in PHASE14_VARIANTS:
        bad["variant"] = (variant, PHASE14_VARIANTS)
    if bad:
        raise ValueError(f"{run_dir}: off-protocol Phase 14 run: {bad}")


def _read_optional(path: Path) -> dict[str, Any]:
    return load_json(path) if path.exists() else {}


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
        if str(run_manifest.get("status", "")) != "complete" or str(artifact.get("status", "")) != "complete":
            continue
        generator_stats = _read_optional(run_dir / "generator_stats.json")
        mask_stats = _read_optional(run_dir / "mask_stats.json")
        embedding_stats = _read_optional(run_dir / "embedding_stats.json")
        row: dict[str, Any] = {
            "run_label": run_label,
            "dataset": str(artifact.get("dataset") or config.get("dataset_name") or run_dir.name),
            "corruption_type": str(corruption.get("corruption_type") or config.get("corruption", {}).get("type")),
            "variant": str(config.get("variant", "")),
            "seed": int(artifact.get("seed", config.get("seed", -1))),
            "epochs": int(config.get("training", {}).get("epochs", -1)),
            "encoder_type": str(config.get("model", {}).get("encoder_type", "")),
            "mask_selector": str(config.get("model", {}).get("mask_selector", "")),
            "benchmark_mode": bool(config.get("benchmark_mode", False)),
            "input_mode": str(config.get("preprocessing", {}).get("input_mode", "")),
            "n_top_genes": int(config.get("preprocessing", {}).get("n_top_genes", -1)),
            "scale_input": bool(config.get("preprocessing", {}).get("scale_input", True)),
            "strict_effective_budget": bool(corruption.get("strict_effective_budget", True)),
            "generator_grad_norm": float(generator_stats.get("generator_grad_norm", 0.0) or 0.0),
            "student_grad_norm_during_generator_step": float(
                generator_stats.get("student_grad_norm_during_generator_step", 0.0) or 0.0
            ),
            "run_dir": str(run_dir),
        }
        row.update(flatten_metrics(load_json(metrics_path)))
        for key in DIAGNOSTIC_KEYS:
            value = corruption.get(key)
            row[key] = float(value) if value is not None else math.nan
        for key in MASK_KEYS:
            value = generator_stats.get(key, mask_stats.get(key))
            row[key] = float(value) if value is not None else math.nan
        for key in EMBEDDING_KEYS:
            value = embedding_stats.get(key)
            row[key] = float(value) if value is not None else math.nan
        require_phase14_protocol(run_dir, row)
        rows.append(row)
    return rows


def validate_expected_grid(
    rows: list[dict[str, Any]],
    *,
    expected_datasets: list[str],
    expected_corruption_types: list[str],
    expected_seeds: list[int],
) -> None:
    expected = {
        (dataset, corruption_type, variant, int(seed))
        for dataset in expected_datasets
        for corruption_type in expected_corruption_types
        for variant in PHASE14_VARIANTS
        for seed in expected_seeds
    }
    seen: dict[tuple[str, str, str, int], int] = {}
    for row in rows:
        key = (str(row["dataset"]), str(row["corruption_type"]), str(row["variant"]), int(row["seed"]))
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
        raise ValueError(f"Formal Phase 14 grid is incomplete or contaminated: {details}")


def mean_std(values: list[float]) -> tuple[float, float]:
    clean = [float(value) for value in values if not math.isnan(float(value))]
    if not clean:
        return math.nan, math.nan
    if len(clean) == 1:
        return clean[0], 0.0
    return mean(clean), stdev(clean)


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    groups: dict[tuple[str, str, str], list[dict[str, Any]]] = {}
    for row in rows:
        groups.setdefault((row["dataset"], row["corruption_type"], row["variant"]), []).append(row)
    out: list[dict[str, Any]] = []
    keys = (*METRIC_KEYS, *DIAGNOSTIC_KEYS, *MASK_KEYS, *EMBEDDING_KEYS, "generator_grad_norm")
    for (dataset, corruption_type, variant), group_rows in sorted(groups.items()):
        row: dict[str, Any] = {
            "dataset": dataset,
            "corruption_type": corruption_type,
            "variant": variant,
            "n_runs": len(group_rows),
            "seeds": ",".join(str(row["seed"]) for row in sorted(group_rows, key=lambda item: item["seed"])),
        }
        for key in keys:
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


def build_differences(aggregate_rows: list[dict[str, Any]]) -> dict[str, Any]:
    grouped: dict[tuple[str, str], dict[str, dict[str, Any]]] = {}
    for row in aggregate_rows:
        grouped.setdefault((row["dataset"], row["corruption_type"]), {})[row["variant"]] = row
    out: dict[str, Any] = {}
    for (dataset, corruption_type), rows in grouped.items():
        if {"control", "advmask"}.issubset(rows):
            key = f"{dataset}|{corruption_type}"
            out[key] = {
                f"advmask_minus_control.{metric}": rows["advmask"][f"{metric}.mean"]
                - rows["control"][f"{metric}.mean"]
                for metric in METRIC_KEYS
            }
    return out


def assess_gate(aggregate_rows: list[dict[str, Any]], differences: dict[str, Any]) -> dict[str, Any]:
    positive_ari = [
        key
        for key, values in differences.items()
        if float(values.get(f"advmask_minus_control.{PRIMARY_METRIC}", math.nan)) > 0.0
    ]
    advmask_rows = [row for row in aggregate_rows if row["variant"] == "advmask"]
    generator_positive = [
        f"{row['dataset']}|{row['corruption_type']}"
        for row in advmask_rows
        if float(row.get("generator_grad_norm.mean", 0.0)) > 0.0
    ]
    concentrated = [
        f"{row['dataset']}|{row['corruption_type']}"
        for row in advmask_rows
        if float(row.get("top_gene_concentration.mean", math.nan)) > 0.8
    ]
    collapsed = [
        f"{row['dataset']}|{row['corruption_type']}"
        for row in advmask_rows
        if float(row.get("mean_pairwise_cosine.mean", 0.0)) > 0.95
        or float(row.get("min_per_dimension_variance.mean", 1.0)) < 1e-4
    ]
    keep = len(positive_ari) >= 2 and len(generator_positive) == len(advmask_rows) and not concentrated and not collapsed
    return {
        "gate_result": "pass" if keep else "fail",
        "recommendation": "keep_advmask" if keep else "drop_or_downgrade_advmask",
        "positive_ari_dataset_corruptions": positive_ari,
        "generator_grad_norm_positive": generator_positive,
        "mask_concentration_flags": concentrated,
        "embedding_collapse_flags": collapsed,
        "note": "Phase 14 development results are not publication claims and do not use validation or sealed test data.",
    }


def write_markdown_report(
    path: Path,
    *,
    rows: list[dict[str, Any]],
    aggregate_rows: list[dict[str, Any]],
    differences: dict[str, Any],
    gate: dict[str, Any],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Phase 14 AdvMask Triage Report",
        "",
        "## 1. Run summary",
        "",
        f"- Complete run count: {len(rows)}",
        f"- Primary metric: `{PRIMARY_METRIC}`",
        "- Scope: MLP encoder only; control random mask vs AdvMask selector.",
        "",
        "## 2. Aggregate results",
        "",
        "| dataset | corruption | variant | ARI mean | ARI std | NMI mean | ACC mean | F1 mean | generator grad | mask entropy | mask gini | effective corruption |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in aggregate_rows:
        lines.append(
            "| {dataset} | {corruption_type} | {variant} | {ari:.6f} | {ari_std:.6f} | {nmi:.6f} | {acc:.6f} | {f1:.6f} | {grad:.6f} | {entropy:.6f} | {gini:.6f} | {eff:.6f} |".format(
                dataset=row["dataset"],
                corruption_type=row["corruption_type"],
                variant=row["variant"],
                ari=float(row["kmeans_known_k.ari.mean"]),
                ari_std=float(row["kmeans_known_k.ari.std"]),
                nmi=float(row["kmeans_known_k.nmi.mean"]),
                acc=float(row["kmeans_known_k.acc.mean"]),
                f1=float(row["kmeans_known_k.f1_macro.mean"]),
                grad=float(row["generator_grad_norm.mean"]),
                entropy=float(row["mask_entropy.mean"]),
                gini=float(row["mask_gini.mean"]),
                eff=float(row["effective_corruption_rate.mean"]),
            )
        )
    lines.extend(
        [
            "",
            "## 3. AdvMask minus control",
            "",
            f"`{json.dumps(differences, sort_keys=True)}`",
            "",
            "## 4. Phase gate",
            "",
            f"- gate_result: `{gate['gate_result']}`",
            f"- recommendation: `{gate['recommendation']}`",
            f"- positive_ari_dataset_corruptions: `{json.dumps(gate['positive_ari_dataset_corruptions'])}`",
            f"- generator_grad_norm_positive: `{json.dumps(gate['generator_grad_norm_positive'])}`",
            f"- mask_concentration_flags: `{json.dumps(gate['mask_concentration_flags'])}`",
            f"- embedding_collapse_flags: `{json.dumps(gate['embedding_collapse_flags'])}`",
            "",
            "## 5. Remaining risks",
            "",
            "- Phase 14 must not run Axial or full model.",
            "- Three-epoch trends must not be written as final publication claims.",
            "- Validation and sealed test datasets are not used for this mechanism decision.",
        ]
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def summarize(
    formal_root: Path,
    *,
    output_dir: Path,
    report_path: Path,
    expected_datasets: list[str],
    expected_corruption_types: list[str],
    expected_seeds: list[int],
) -> tuple[Path, Path, Path]:
    rows = discover_runs(formal_root, "formal")
    validate_expected_grid(
        rows,
        expected_datasets=expected_datasets,
        expected_corruption_types=expected_corruption_types,
        expected_seeds=expected_seeds,
    )
    aggregate_rows = aggregate(rows)
    differences = build_differences(aggregate_rows)
    gate = assess_gate(aggregate_rows, differences)
    output_dir.mkdir(parents=True, exist_ok=True)
    runs_path = output_dir / "advmask_triage_runs.csv"
    summary_path = output_dir / "advmask_triage_summary.csv"
    json_path = output_dir / "advmask_triage_report.json"
    write_csv(runs_path, rows)
    write_csv(summary_path, aggregate_rows)
    with open(json_path, "w", encoding="utf-8") as handle:
        json.dump(
            {
                "primary_metric": PRIMARY_METRIC,
                "n_formal_runs": len(rows),
                "aggregate": aggregate_rows,
                "differences": differences,
                "gate": gate,
                "validation_or_sealed_used": False,
            },
            handle,
            indent=2,
            sort_keys=True,
        )
    write_markdown_report(report_path, rows=rows, aggregate_rows=aggregate_rows, differences=differences, gate=gate)
    return summary_path, json_path, report_path


def main() -> int:
    parser = argparse.ArgumentParser(description="Summarize Phase 14 CAAM AdvMask triage runs.")
    parser.add_argument("--formal_root", type=Path, required=True)
    parser.add_argument("--expected_datasets", type=str, default=",".join(PHASE14_DATASETS))
    parser.add_argument("--expected_corruption_types", type=str, required=True)
    parser.add_argument("--expected_seeds", type=str, default=",".join(str(seed) for seed in PHASE14_SEEDS))
    parser.add_argument("--output_dir", type=Path, default=PROJECT_ROOT / "results/CAAM_scMAE_correction")
    parser.add_argument(
        "--report_path",
        type=Path,
        default=PROJECT_ROOT / "methods/DeepLearning/CAAM_scMAE/benchmark/PHASE14_ADVMASK_TRIAGE_REPORT.md",
    )
    args = parser.parse_args()
    try:
        summary_path, json_path, report_path = summarize(
            args.formal_root,
            output_dir=args.output_dir,
            report_path=args.report_path,
            expected_datasets=parse_csv(args.expected_datasets),
            expected_corruption_types=parse_csv(args.expected_corruption_types),
            expected_seeds=[int(seed) for seed in parse_csv(args.expected_seeds)],
        )
    except Exception as exc:
        print(f"CAAM Phase 14 summary failed: {exc}")
        return 1
    print(f"Wrote {summary_path}")
    print(f"Wrote {json_path}")
    print(f"Wrote {report_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
