import json

import pytest
import yaml

from methods.DeepLearning.CAAM_scMAE.benchmark.summarize_corruption_triad import (
    CORRUPTION_TYPES,
    discover_runs,
    validate_expected_grid,
    write_markdown_report,
)


def _write_run(
    root,
    *,
    dataset="D1",
    corruption_type="scmae_shuffle",
    seed=42,
    status="complete",
    artifact_status="complete",
    variant="control",
    encoder_type="mlp",
    mask_selector="random",
):
    run_dir = root / f"{dataset}__{corruption_type}__seed{seed}"
    run_dir.mkdir(parents=True)
    (run_dir / "metrics.json").write_text(
        json.dumps(
            {
                "kmeans_known_k": {"acc": 1.0, "nmi": 1.0, "ari": 1.0, "f1_macro": 1.0},
                "leiden_fixed": {"acc": 1.0, "nmi": 1.0, "ari": 1.0, "f1_macro": 1.0},
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "artifact_manifest.json").write_text(
        json.dumps({"status": artifact_status, "dataset": dataset, "seed": seed}),
        encoding="utf-8",
    )
    (run_dir / "run_manifest.json").write_text(json.dumps({"status": status}), encoding="utf-8")
    (run_dir / "resolved_config.yaml").write_text(
        yaml.safe_dump(
            {
                "benchmark_mode": True,
                "dataset_name": dataset,
                "seed": seed,
                "variant": variant,
                "model": {"encoder_type": encoder_type, "mask_selector": mask_selector},
                "training": {"epochs": 10},
                "preprocessing": {"input_mode": "log1p", "n_top_genes": 2000, "scale_input": False},
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "corruption_stats.json").write_text(
        json.dumps(
            {
                "corruption_type": corruption_type,
                "strict_effective_budget": False,
                "zero_to_zero_rate": 0.0,
                "effective_corruption_rate": 1.0,
                "budget_deficit_rate": 0.0,
                "mean_abs_delta": 1.0,
                "mean_abs_delta_masked": 1.0,
            }
        ),
        encoding="utf-8",
    )
    return run_dir


def test_summary_discovery_excludes_failed_runs(tmp_path):
    _write_run(tmp_path, status="failed")
    rows = discover_runs(tmp_path, "formal")
    assert rows == []
    with pytest.raises(ValueError, match="Formal Phase 13 grid"):
        validate_expected_grid(
            rows,
            expected_datasets=["D1"],
            expected_corruption_types=("scmae_shuffle",),
            expected_seeds=[42],
        )


def test_summary_rejects_off_protocol_runs(tmp_path):
    _write_run(tmp_path, variant="axial", encoder_type="axial")
    with pytest.raises(ValueError, match="off-protocol Phase 13 run"):
        discover_runs(tmp_path, "formal")


def test_summary_requires_full_expected_dataset_corruption_seed_grid(tmp_path):
    for corruption_type in CORRUPTION_TYPES:
        _write_run(tmp_path, dataset="D1", corruption_type=corruption_type, seed=42)
    rows = discover_runs(tmp_path, "formal")
    with pytest.raises(ValueError, match="Formal Phase 13 grid"):
        validate_expected_grid(
            rows,
            expected_datasets=["D1", "D2"],
            expected_corruption_types=CORRUPTION_TYPES,
            expected_seeds=[42, 2024],
        )


def test_markdown_report_includes_phase_gate_result(tmp_path):
    report_path = tmp_path / "report.md"
    aggregate_rows = [
        {
            "dataset": "D1",
            "corruption_type": "scmae_shuffle",
            "kmeans_known_k.ari.mean": 1.0,
            "kmeans_known_k.ari.std": 0.0,
            "kmeans_known_k.nmi.mean": 1.0,
            "kmeans_known_k.acc.mean": 1.0,
            "kmeans_known_k.f1_macro.mean": 1.0,
            "zero_to_zero_rate.mean": 0.0,
            "effective_corruption_rate.mean": 1.0,
            "mean_abs_delta.mean": 1.0,
        }
    ]
    write_markdown_report(
        report_path,
        smoke_rows=[],
        formal_rows=[{}],
        aggregate_rows=aggregate_rows,
        differences={},
        nonzero_aware_assessment={},
        recommendation={"recommended_corruption_type": "scmae_shuffle", "reason": "complete formal grid"},
    )

    text = report_path.read_text(encoding="utf-8")
    assert "gate_result: `pass`" in text
    assert "gate_reason:" in text
