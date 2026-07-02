#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import numpy as np
import yaml


def _load_json(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return json.load(handle)


def _load_yaml(path: Path) -> dict[str, Any]:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def _candidate_run_dirs(path: Path) -> list[Path]:
    """Return CAAM run directories under a single run dir, dataset dir, or output root.

    scripts/run_formal_benchmark.py may write either:
      out_dir/run_id/
    or, more commonly:
      out_dir/dataset_name/run_id/

    Therefore this validator searches recursively for artifact_manifest.json while
    still accepting a single run directory directly.
    """
    if (path / "artifact_manifest.json").exists():
        return [path]
    run_dirs = {manifest.parent for manifest in path.rglob("artifact_manifest.json")}
    return sorted(run_dirs)


def _require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def validate_run_dir(run_dir: Path) -> list[str]:
    failures: list[str] = []
    artifact_path = run_dir / "artifact_manifest.json"
    config_path = run_dir / "resolved_config.yaml"
    metrics_path = run_dir / "metrics.json"
    embedding_path = run_dir / "embedding_final.npy"
    labels_path = run_dir / "labels.npy"
    run_manifest_path = run_dir / "run_manifest.json"
    preprocess_path = run_dir / "preprocess_config.json"
    corruption_path = run_dir / "corruption_stats.json"
    selected_gene_indices_path = run_dir / "selected_gene_indices.npy"

    for path in (
        artifact_path,
        config_path,
        metrics_path,
        embedding_path,
        labels_path,
        run_manifest_path,
        preprocess_path,
        corruption_path,
        selected_gene_indices_path,
    ):
        _require(path.exists(), f"{run_dir}: missing {path.name}", failures)
    if failures:
        return failures

    artifact = _load_json(artifact_path)
    config = _load_yaml(config_path)
    metrics = _load_json(metrics_path)
    run_manifest = _load_json(run_manifest_path)
    preprocess_config = _load_json(preprocess_path)
    corruption_stats = _load_json(corruption_path)
    embedding = np.load(embedding_path)
    labels = np.load(labels_path)
    selected_gene_indices = np.load(selected_gene_indices_path)

    preprocessing = config.get("preprocessing", {})
    n_top_genes = preprocessing.get("n_top_genes")
    feature_space_source = preprocessing.get("feature_space_source") or preprocess_config.get("feature_space_source")
    _require(artifact.get("status") == "complete", f"{run_dir}: artifact_manifest.status != complete", failures)
    _require(artifact.get("variant") == "full", f"{run_dir}: artifact_manifest.variant != full", failures)
    _require(config.get("benchmark_mode") is True, f"{run_dir}: resolved_config.benchmark_mode != true", failures)
    _require(config.get("variant") == "full", f"{run_dir}: resolved_config.variant != full", failures)
    _require(preprocessing.get("input_mode") == "log1p", f"{run_dir}: preprocessing.input_mode != log1p", failures)
    _require(n_top_genes in {2000, 3000, 0}, f"{run_dir}: preprocessing.n_top_genes not in {{2000, 3000, 0}}", failures)
    if n_top_genes == 0:
        _require(
            feature_space_source in {"full_gene_stress", "external_hvg"},
            f"{run_dir}: n_top_genes=0 requires feature_space_source full_gene_stress or external_hvg",
            failures,
        )
    else:
        _require(bool(feature_space_source), f"{run_dir}: missing feature_space_source", failures)
    _require(preprocessing.get("scale_input") is False, f"{run_dir}: preprocessing.scale_input != false", failures)
    _require(
        selected_gene_indices.ndim == 1 and selected_gene_indices.shape[0] == int(preprocessing.get("actual_n_genes_after_selection", 0)),
        f"{run_dir}: selected_gene_indices length does not match actual_n_genes_after_selection",
        failures,
    )
    for key in (
        "corruption_type",
        "mask_ratio",
        "n_top_genes",
        "actual_n_genes",
        "zero_to_zero_rate",
        "effective_corruption_rate",
        "budget_deficit_rate",
        "mean_abs_delta",
        "mean_abs_delta_masked",
        "strict_effective_budget",
    ):
        _require(key in corruption_stats, f"{run_dir}: corruption_stats missing {key}", failures)
    _require(
        corruption_stats.get("strict_effective_budget") is False,
        f"{run_dir}: strict_effective_budget should default to false",
        failures,
    )
    _require(
        metrics.get("kmeans_known_k", {}).get("uses_known_k") is True,
        f"{run_dir}: kmeans_known_k.uses_known_k != true",
        failures,
    )
    _require(
        metrics.get("leiden_fixed", {}).get("uses_known_k") is False,
        f"{run_dir}: leiden_fixed.uses_known_k != false",
        failures,
    )
    _require(
        embedding.ndim == 2 and embedding.shape[0] == labels.shape[0],
        f"{run_dir}: embedding_final rows do not match labels",
        failures,
    )
    _require(run_manifest.get("status") == "complete", f"{run_dir}: run_manifest.status != complete", failures)
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate CAAM-scMAE formal smoke outputs.")
    parser.add_argument("output_dir", type=Path, help="Formal smoke root, dataset directory, or a single CAAM run directory.")
    args = parser.parse_args()

    output_dir = args.output_dir
    if not output_dir.exists():
        print(f"ERROR: output directory does not exist: {output_dir}")
        return 2

    run_dirs = _candidate_run_dirs(output_dir)
    if not run_dirs:
        print(f"ERROR: no CAAM formal run directories found under: {output_dir}")
        return 2

    failures: list[str] = []
    for run_dir in run_dirs:
        failures.extend(validate_run_dir(run_dir))

    if failures:
        print("CAAM formal smoke validation failed:")
        for failure in failures:
            print(f"  - {failure}")
        return 1

    print(f"CAAM formal smoke validation passed for {len(run_dirs)} run(s).")
    for run_dir in run_dirs:
        print(f"  - {run_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
