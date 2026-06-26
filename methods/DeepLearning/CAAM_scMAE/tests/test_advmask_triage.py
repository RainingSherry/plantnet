import json

import pytest
import yaml

from methods.DeepLearning.CAAM_scMAE.benchmark.run_advmask_triage import build_run_commands, write_phase14_config
from methods.DeepLearning.CAAM_scMAE.benchmark.summarize_advmask_triage import (
    build_differences,
    discover_runs,
    summarize,
)


def _write_run(
    root,
    *,
    dataset="D1",
    corruption_type="scmae_shuffle",
    variant="control",
    seed=42,
    ari=0.5,
    generator_grad_norm=0.0,
):
    run_dir = root / f"{dataset}__{corruption_type}__{variant}__seed{seed}"
    run_dir.mkdir(parents=True)
    metrics = {
        "kmeans_known_k": {"acc": ari + 0.1, "nmi": ari + 0.2, "ari": ari, "f1_macro": ari + 0.3},
        "leiden_fixed": {"acc": ari / 2.0, "nmi": ari / 3.0, "ari": ari / 4.0, "f1_macro": ari / 5.0},
    }
    mask_selector = "adversarial" if variant == "advmask" else "random"
    (run_dir / "metrics.json").write_text(json.dumps(metrics), encoding="utf-8")
    (run_dir / "artifact_manifest.json").write_text(
        json.dumps({"status": "complete", "dataset": dataset, "seed": seed}),
        encoding="utf-8",
    )
    (run_dir / "run_manifest.json").write_text(json.dumps({"status": "complete"}), encoding="utf-8")
    (run_dir / "resolved_config.yaml").write_text(
        yaml.safe_dump(
            {
                "benchmark_mode": True,
                "dataset_name": dataset,
                "seed": seed,
                "variant": variant,
                "model": {"encoder_type": "mlp", "mask_selector": mask_selector},
                "training": {"epochs": 3},
                "preprocessing": {"input_mode": "log1p", "n_top_genes": 2000, "scale_input": False},
                "corruption": {"type": corruption_type},
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "corruption_stats.json").write_text(
        json.dumps(
            {
                "corruption_type": corruption_type,
                "strict_effective_budget": False,
                "zero_to_zero_rate": 0.1,
                "effective_corruption_rate": 0.9,
                "budget_deficit_rate": 0.0,
                "mean_abs_delta": 0.2,
                "mean_abs_delta_masked": 0.3,
            }
        ),
        encoding="utf-8",
    )
    (run_dir / "mask_stats.json").write_text(
        json.dumps({"mask_entropy": 0.7, "mask_gini": 0.2, "top_gene_concentration": 0.3}),
        encoding="utf-8",
    )
    (run_dir / "generator_stats.json").write_text(
        json.dumps({"generator_grad_norm": generator_grad_norm, "mask_entropy": 0.8, "mask_gini": 0.25}),
        encoding="utf-8",
    )
    (run_dir / "embedding_stats.json").write_text(
        json.dumps({"effective_rank": 4.0, "mean_pairwise_cosine": 0.2, "min_per_dimension_variance": 0.01}),
        encoding="utf-8",
    )
    return run_dir


def test_advmask_triage_runner_only_builds_control_and_advmask_commands(tmp_path):
    data_root = tmp_path / "data"
    data_file = data_root / "processed/Quake_Smart-seq2_Lung.h5ad"
    data_file.parent.mkdir(parents=True)
    data_file.write_text("placeholder", encoding="utf-8")

    commands = build_run_commands(
        runner=tmp_path / "run.py",
        config_path=tmp_path / "phase14_config.json",
        run_root=tmp_path / "runs",
        data_root=data_root,
        dataset_names=["Quake_Smart-seq2_Lung"],
        corruption_types=["scmae_shuffle"],
        seeds=[42],
        epochs=3,
        gpu=1,
        no_cuda=True,
    )

    variants = [cmd[3] for cmd in commands]
    assert variants == ["control", "advmask"]
    flattened = "\n".join(" ".join(cmd[-1]) for cmd in commands)
    assert "--variant axial" not in flattened
    assert "--variant full" not in flattened
    assert "method_manifest.yaml" not in flattened
    assert "--config" in flattened


def test_advmask_triage_config_activates_generator_within_three_epochs(tmp_path):
    config_path = write_phase14_config(
        tmp_path / "phase14_advmask_triage_config.json",
        student_warmup_epochs=1,
        generator_update_interval=5,
    )
    config = json.loads(config_path.read_text(encoding="utf-8"))

    assert config["training"]["student_warmup_epochs"] == 1
    assert config["training"]["student_warmup_epochs"] < 3
    assert config["training"]["generator_update_interval"] == 5


def test_advmask_triage_summary_calculates_advmask_minus_control(tmp_path):
    _write_run(tmp_path, variant="control", ari=0.4)
    _write_run(tmp_path, variant="advmask", ari=0.55, generator_grad_norm=1.2)

    rows = discover_runs(tmp_path, "formal")
    summary_path, report_json_path, report_md_path = summarize(
        tmp_path,
        output_dir=tmp_path / "out",
        report_path=tmp_path / "PHASE14_ADVMASK_TRIAGE_REPORT.md",
        expected_datasets=["D1"],
        expected_corruption_types=["scmae_shuffle"],
        expected_seeds=[42],
    )

    assert summary_path.exists()
    assert report_json_path.exists()
    assert report_md_path.exists()
    differences = build_differences([row for row in json.loads(report_json_path.read_text())["aggregate"]])
    value = differences["D1|scmae_shuffle"]["advmask_minus_control.kmeans_known_k.ari"]
    assert round(value, 6) == 0.15
    assert any(row["variant"] == "advmask" and row["generator_grad_norm"] == 1.2 for row in rows)


def test_advmask_triage_gate_requires_mean_delta_above_seed_noise(tmp_path):
    for dataset, delta in (("D1", 0.01), ("D2", 0.01), ("D3", -0.01)):
        _write_run(tmp_path, dataset=dataset, variant="control", seed=1, ari=0.4)
        _write_run(tmp_path, dataset=dataset, variant="control", seed=2, ari=0.8)
        _write_run(tmp_path, dataset=dataset, variant="advmask", seed=1, ari=0.4 + delta, generator_grad_norm=1.0)
        _write_run(tmp_path, dataset=dataset, variant="advmask", seed=2, ari=0.8 + delta, generator_grad_norm=1.0)

    _summary_path, report_json_path, _report_md_path = summarize(
        tmp_path,
        output_dir=tmp_path / "out",
        report_path=tmp_path / "PHASE14_ADVMASK_TRIAGE_REPORT.md",
        expected_datasets=["D1", "D2", "D3"],
        expected_corruption_types=["scmae_shuffle"],
        expected_seeds=[1, 2],
    )
    gate = json.loads(report_json_path.read_text(encoding="utf-8"))["gate"]

    assert len(gate["positive_ari_dataset_corruptions"]) == 2
    assert gate["effect_size_gate_pass"] is False
    assert gate["gate_result"] == "fail"


def test_advmask_triage_rejects_axial_or_full_runs(tmp_path):
    _write_run(tmp_path, variant="full", ari=0.9)
    with pytest.raises(ValueError, match="off-protocol Phase 14 run"):
        discover_runs(tmp_path, "formal")
