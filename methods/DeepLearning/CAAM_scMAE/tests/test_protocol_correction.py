import json
from pathlib import Path

import anndata as ad
import numpy as np
import pytest
import yaml

from methods.DeepLearning.CAAM_scMAE.benchmark.validate_formal_smoke import validate_run_dir
from methods.DeepLearning.CAAM_scMAE.data.preprocessing import load_caam_data
from methods.DeepLearning.CAAM_scMAE.registry import build_arg_parser, resolve_config
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_config
from methods.DeepLearning.CAAM_scMAE.trainers.common import CAAMTrainer


def _args(extra=None):
    parser = build_arg_parser()
    base = [
        "--data_path",
        "toy.h5ad",
        "--save_dir",
        "out",
        "--n_clusters",
        "2",
        "--benchmark_mode",
        "true",
    ]
    return parser.parse_args(base + list(extra or []))


def test_benchmark_mode_defaults_to_hvg_feature_space():
    cfg = resolve_config(_args())
    assert cfg["preprocessing"]["input_mode"] == "log1p"
    assert cfg["preprocessing"]["n_top_genes"] == 2000
    assert cfg["preprocessing"]["scale_input"] is False


def test_cli_n_top_genes_zero_is_not_overwritten():
    cfg = resolve_config(_args(["--n_top_genes", "0"]))
    assert cfg["preprocessing"]["n_top_genes"] == 0


def test_strict_effective_budget_defaults_false():
    cfg = resolve_config(_args())
    assert cfg["corruption"]["strict_effective_budget"] is False


def test_budget_deficit_is_diagnostic_unless_strict():
    trainer = object.__new__(CAAMTrainer)
    trainer.config = toy_config("control")
    trainer.config["runtime"]["fail_fast"] = True
    trainer.config["corruption"]["strict_effective_budget"] = False
    assert trainer._check_budget_deficit({"budget_deficit_rate": 1.0}, where="unit test") == 1.0

    trainer.config["corruption"]["strict_effective_budget"] = True
    with pytest.raises(RuntimeError, match="effective mask budget deficit rate"):
        trainer._check_budget_deficit({"budget_deficit_rate": 1.0}, where="unit test")


def test_selected_gene_indices_reproducible_same_seed(tmp_path):
    rng = np.random.default_rng(7)
    x = np.log1p(rng.poisson(2.0, size=(30, 20)).astype(np.float32))
    adata = ad.AnnData(X=x)
    adata.var_names = [f"gene_{i}" for i in range(20)]
    path = tmp_path / "toy.h5ad"
    adata.write_h5ad(path)

    first = load_caam_data(
        str(path),
        input_mode="log1p",
        target_sum=10000.0,
        n_top_genes=5,
        scale_input=False,
        benchmark_mode=True,
        seed=123,
    )
    second = load_caam_data(
        str(path),
        input_mode="log1p",
        target_sum=10000.0,
        n_top_genes=5,
        scale_input=False,
        benchmark_mode=True,
        seed=123,
    )
    assert np.array_equal(first.selected_gene_indices, second.selected_gene_indices)
    assert first.preprocess_config["feature_space_source"] == "hvg"
    assert first.preprocess_config["actual_n_genes_after_selection"] == 5


def test_validate_formal_smoke_accepts_hvg_protocol(tmp_path):
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    embedding = np.ones((4, 2), dtype=np.float32)
    labels = np.asarray([0, 0, 1, 1], dtype=np.int64)
    np.save(run_dir / "embedding_final.npy", embedding)
    np.save(run_dir / "labels.npy", labels)
    np.save(run_dir / "selected_gene_indices.npy", np.arange(2000, dtype=np.int64))
    (run_dir / "args.json").write_text("{}", encoding="utf-8")
    (run_dir / "artifact_manifest.json").write_text(
        json.dumps({"status": "complete", "variant": "full"}),
        encoding="utf-8",
    )
    (run_dir / "metrics.json").write_text(
        json.dumps({"kmeans_known_k": {"uses_known_k": True}, "leiden_fixed": {"uses_known_k": False}}),
        encoding="utf-8",
    )
    (run_dir / "run_manifest.json").write_text(json.dumps({"status": "complete"}), encoding="utf-8")
    config = {
        "benchmark_mode": True,
        "variant": "full",
        "preprocessing": {
            "input_mode": "log1p",
            "n_top_genes": 2000,
            "scale_input": False,
            "feature_space_source": "hvg",
            "actual_n_genes_after_selection": 2000,
            "selected_gene_indices_path": "selected_gene_indices.npy",
            "selected_gene_names_path": "selected_genes.txt",
        },
    }
    (run_dir / "resolved_config.yaml").write_text(yaml.safe_dump(config), encoding="utf-8")
    (run_dir / "preprocess_config.json").write_text(
        json.dumps({"feature_space_source": "hvg", "actual_n_genes_after_selection": 2000}),
        encoding="utf-8",
    )
    (run_dir / "corruption_stats.json").write_text(
        json.dumps(
            {
                "corruption_type": "matched_donor",
                "mask_ratio": 0.4,
                "n_top_genes": 2000,
                "actual_n_genes": 2000,
                "zero_to_zero_rate": 0.0,
                "effective_corruption_rate": 1.0,
                "budget_deficit_rate": 0.0,
                "mean_abs_delta": 0.1,
                "mean_abs_delta_masked": 0.2,
                "strict_effective_budget": False,
            }
        ),
        encoding="utf-8",
    )
    assert validate_run_dir(run_dir) == []
