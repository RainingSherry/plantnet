from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import torch

from methods.DeepLearning.APA_scMAE.config import build_arg_parser, resolve_config
from methods.DeepLearning.APA_scMAE.data import (
    APAExpressionDataset,
    ScMAEShuffleCorruption,
    assert_no_training_labels,
    compute_gene_stats,
    load_apa_data,
)
from methods.DeepLearning.APA_scMAE.losses import distortion_loss, student_losses
from methods.DeepLearning.APA_scMAE.model import APAModel
from methods.DeepLearning.APA_scMAE.run import (
    artifact_manifest,
    attention_guard,
    hungarian_map,
    mapped_clustering_metrics,
    validate_embedding_inputs,
    validate_required_files,
)
from methods.DeepLearning.APA_scMAE.trainer import APATrainer


def _config(tmp_path):
    return {
        "seed": 7,
        "method_name": "apa_scmae",
        "runtime": {
            "num_workers": 0,
            "fail_fast": True,
            "max_attention_elements": 300_000_000,
            "force_large_attention": False,
        },
        "model": {
            "token_dim": 8,
            "cell_dim": 6,
            "attention_heads": 2,
            "attention_dropout": 0.0,
            "dropout": 0.0,
        },
        "training": {
            "epochs": 1,
            "batch_size": 4,
            "lr_student": 1.0e-3,
            "lr_generator": 1.0e-3,
            "weight_decay": 0.0,
            "student_grad_clip": 5.0,
            "generator_grad_clip": 1.0,
            "generator_update_interval": 1,
            "gamma": 0.5,
        },
        "mask": {"ratio": 0.4, "temperature": 1.0, "masked_data_weight": 0.75, "generator_topk_only_effective": True},
        "corruption": {"type": "scmae_shuffle", "changed_tolerance_abs": 1.0e-6, "changed_tolerance_rel": 1.0e-5},
        "generator_loss": {"lambda_entropy": 0.01, "lambda_balance": 0.01, "lambda_distortion": 0.01, "lambda_coverage": 0.01},
        "save_dir": str(tmp_path),
    }


def _write_h5ad(path, x, obs=None, layers=None, raw_x=None):
    obs = pd.DataFrame(obs or {}, index=[f"c{i}" for i in range(x.shape[0])])
    var = pd.DataFrame(index=[f"g{i}" for i in range(x.shape[1])])
    adata = ad.AnnData(X=np.asarray(x, dtype=np.float32), obs=obs, var=var)
    for key, value in (layers or {}).items():
        adata.layers[key] = np.asarray(value, dtype=np.float32)
    if raw_x is not None:
        adata.raw = ad.AnnData(X=np.asarray(raw_x, dtype=np.float32), obs=obs.copy(), var=var.copy())
    adata.write_h5ad(path)


def test_dataset_returns_no_labels():
    dataset = APAExpressionDataset(np.ones((3, 5), dtype=np.float32))
    assert set(dataset[0]) == {"index", "x"}
    with pytest.raises(RuntimeError):
        assert_no_training_labels({"x": torch.ones(2, 3), "labels": torch.zeros(2)})


@pytest.mark.parametrize("label_key", ["Celltype", "maintype", "resolved_label"])
def test_safe_label_candidates_are_detected(tmp_path, label_key):
    data_path = tmp_path / f"{label_key}.h5ad"
    _write_h5ad(data_path, np.ones((4, 5)), obs={label_key: ["a", "b", "a", "b"]})
    bundle = load_apa_data(str(data_path), input_mode="log1p", target_sum=1.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert bundle.label_key == label_key
    assert bundle.labels is not None


def test_label_errors_and_explicit_key(tmp_path):
    no_label = tmp_path / "unlabeled.h5ad"
    _write_h5ad(no_label, np.ones((4, 5)))
    with pytest.raises(ValueError, match="No label column"):
        load_apa_data(str(no_label), input_mode="auto", target_sum=1.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)

    multi = tmp_path / "multi.h5ad"
    _write_h5ad(multi, np.ones((3, 4)), obs={"cell_type": ["a", "b", "a"], "maintype": ["x", "y", "x"]})
    with pytest.raises(ValueError, match="Multiple candidate"):
        load_apa_data(str(multi), input_mode="log1p", target_sum=1.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    bundle = load_apa_data(
        str(multi),
        input_mode="log1p",
        target_sum=1.0,
        n_top_genes=0,
        scale_input=False,
        n_prototypes=2,
        pca_dim=2,
        seed=1,
        label_key="maintype",
    )
    assert bundle.label_key == "maintype"
    with pytest.raises(ValueError, match="not found"):
        load_apa_data(str(multi), input_mode="log1p", target_sum=1.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1, label_key="missing")

    nan_path = tmp_path / "nan.h5ad"
    _write_h5ad(nan_path, np.ones((3, 4)), obs={"cell_type": ["a", None, "b"]})
    with pytest.raises(ValueError, match="missing values"):
        load_apa_data(str(nan_path), input_mode="log1p", target_sum=1.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)


def test_raw_count_source_modes_and_feature_metadata(tmp_path):
    counts = np.array([[1, 2, 0, 0], [3, 0, 1, 0], [0, 4, 2, 1]], dtype=np.float32)
    log_x = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 10.0)
    obs = {"cell_type": ["a", "b", "a"]}

    counts_path = tmp_path / "counts_priority.h5ad"
    _write_h5ad(counts_path, np.full_like(counts, 99.0), obs=obs, layers={"counts": counts})
    bundle = load_apa_data(str(counts_path), input_mode="auto", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    np.testing.assert_allclose(bundle.x, log_x, rtol=1.0e-6, atol=1.0e-6)
    assert bundle.preprocess_config["raw_count_source"] == 'layers["counts"]'
    assert bundle.preprocess_config["input_mode_resolved"] == "raw"
    assert bundle.preprocess_config["feature_space_source"] == "full_gene"

    raw_path = tmp_path / "raw_priority.h5ad"
    _write_h5ad(raw_path, log_x, obs=obs, raw_x=counts)
    raw_bundle = load_apa_data(str(raw_path), input_mode="auto", target_sum=10.0, n_top_genes=4, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert raw_bundle.preprocess_config["raw_count_source"] == "adata.raw.X"
    assert raw_bundle.preprocess_config["feature_space_source"] == "full_gene_all_available"

    log_path = tmp_path / "log.h5ad"
    _write_h5ad(log_path, log_x, obs=obs)
    with pytest.raises(ValueError, match="no raw count source"):
        load_apa_data(str(log_path), input_mode="raw", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    auto_log = load_apa_data(str(log_path), input_mode="auto", target_sum=10.0, n_top_genes=2, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert auto_log.preprocess_config["input_mode_resolved"] == "log1p"
    assert auto_log.preprocess_config["feature_space_source"] == "hvg_variance"


def test_gene_stats_zero_rate_uses_prescale_values(tmp_path):
    counts = np.array([[0, 1, 0], [0, 2, 3], [4, 0, 0]], dtype=np.float32)
    data_path = tmp_path / "scale.h5ad"
    _write_h5ad(data_path, counts, obs={"cell_type": ["a", "b", "a"]}, layers={"counts": counts})
    bundle = load_apa_data(str(data_path), input_mode="raw", target_sum=10.0, n_top_genes=0, scale_input=True, n_prototypes=2, pca_dim=2, seed=1)
    assert bundle.preprocess_config["gene_stats_source"] == "pre_scale_log1p"
    assert np.any(bundle.gene_stats[:, 2] > 0)


def test_config_defaults_do_not_override_yaml_and_no_cuda_works(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("seed: 123\nmethod_name: yaml_method\nskip_eval: true\nruntime:\n  gpu: 3\n", encoding="utf-8")
    parser = build_arg_parser()
    args = parser.parse_args(["--config", str(config_path), "--data_path", "dummy.h5ad", "--save_dir", str(tmp_path / "out"), "--n_clusters", "2"])
    cfg = resolve_config(args)
    assert cfg["seed"] == 123
    assert cfg["method_name"] == "yaml_method"
    assert cfg["skip_eval"] is True
    assert cfg["runtime"]["gpu"] == 3

    args = parser.parse_args(["--data_path", "dummy.h5ad", "--save_dir", str(tmp_path / "out2"), "--n_clusters", "2", "--no_cuda"])
    cfg = resolve_config(args)
    assert cfg["runtime"]["no_cuda"] is True


def test_cyclic_and_extra_cluster_metrics():
    metrics = mapped_clustering_metrics(np.array([0, 1, 2, 0, 1, 2]), np.array([1, 2, 0, 1, 2, 0]))
    assert metrics["acc"] == pytest.approx(1.0)
    assert metrics["f1_macro"] == pytest.approx(1.0)

    mapped = hungarian_map(np.array([0, 0, 1, 1]), np.array([0, 1, 2, 3]))
    assert np.any(mapped == -1)
    metrics = mapped_clustering_metrics(np.array([0, 0, 1, 1]), np.array([0, 1, 2, 3]))
    assert metrics["acc"] < 1.0
    assert metrics["f1_macro"] < 1.0


def test_embedding_validation_errors():
    with pytest.raises(ValueError, match="NaN"):
        validate_embedding_inputs(np.array([[0.0], [np.nan]]), np.array([0, 1]), 2)
    with pytest.raises(ValueError, match="labels length"):
        validate_embedding_inputs(np.ones((3, 2)), np.array([0, 1]), 2)
    with pytest.raises(ValueError, match="cannot exceed"):
        validate_embedding_inputs(np.ones((3, 2)), np.array([0, 1, 1]), 4)


def test_student_loss_and_distortion_contracts():
    x = torch.zeros(2, 4)
    x_recon = torch.ones(2, 4)
    logits = torch.zeros(2, 4)
    half = torch.tensor([[1, 1, 0, 0], [1, 1, 0, 0]], dtype=torch.float32)
    sparse = torch.tensor([[1, 0, 0, 0], [1, 0, 0, 0]], dtype=torch.float32)
    rec_half = student_losses(x, x_recon, logits, half, masked_data_weight=0.75, gamma=0.0)["loss_rec"]
    rec_sparse = student_losses(x, x_recon, logits, sparse, masked_data_weight=0.75, gamma=0.0)["loss_rec"]
    assert rec_half == pytest.approx(rec_sparse)
    assert torch.isfinite(student_losses(x, x_recon, logits, torch.zeros_like(x), masked_data_weight=0.75, gamma=0.0)["loss_rec"])
    assert torch.isfinite(student_losses(x, x_recon, logits, torch.ones_like(x), masked_data_weight=0.75, gamma=0.0)["loss_rec"])

    delta = torch.tensor([[1.0, 100.0, 1.0, 100.0]])
    assert distortion_loss(torch.tensor([[0.0, 1.0, 0.0, 1.0]]), delta) > distortion_loss(torch.tensor([[1.0, 1.0, 0.0, 0.0]]), delta)


def test_model_shapes_effective_topk_and_budget_deficit():
    x = torch.ones(2, 6)
    replacement = x.clone()
    replacement[0, :3] += 1.0
    replacement[1, :1] += 1.0
    effective = (~torch.isclose(x, replacement)).float()
    gene_stats = torch.as_tensor(compute_gene_stats(x.numpy()), dtype=torch.float32)
    prototypes = torch.randn(2, 3)
    model = APAModel(n_genes=6, token_dim=6, cell_dim=4, proto_dim=3, attention_heads=2, dropout=0.0)
    gene_vec, stat_vec = model.shared_context(gene_stats)
    gen = model.generator(x, replacement, effective, gene_vec, stat_vec, prototypes, mask_ratio=0.5, temperature=1.0)
    assert gen["logits"].shape == (2, 6)
    assert torch.all(gen["mask_hard"] <= effective)
    assert gen["mask_hard"][0].sum() == 3
    assert gen["mask_hard"][1].sum() == 1
    assert gen["budget_deficit"][1] == 2
    out = model.student(x * (1.0 - gen["mask_hard"]) + replacement * gen["mask_hard"], gene_vec, stat_vec, prototypes)
    assert out["z"].shape == (2, 4)
    assert out["mask_logits"].shape == (2, 6)
    assert out["x_recon"].shape == (2, 6)


def test_dynamic_corruption_no_global_matrix_and_reproducible():
    x = torch.arange(20, dtype=torch.float32).view(5, 4)
    corr_a = ScMAEShuffleCorruption(x, seed=3, atol=1e-6, rtol=1e-5)
    corr_b = ScMAEShuffleCorruption(x, seed=3, atol=1e-6, rtol=1e-5)
    assert corr_a.x_full.device.type == "cpu"
    assert not hasattr(corr_a, "permutations")
    sample_a = corr_a.sample(torch.tensor([0, 3]), device=torch.device("cpu"))
    sample_b = corr_b.sample(torch.tensor([0, 3]), device=torch.device("cpu"))
    assert sample_a["replacement"].shape == (2, 4)
    assert sample_a["donor_indices"].shape == (2, 4)
    assert torch.all(sample_a["donor_indices"] != torch.tensor([0, 3]).view(-1, 1))
    assert torch.equal(sample_a["donor_indices"], sample_b["donor_indices"])
    large = ScMAEShuffleCorruption(torch.zeros(100_000, 2), seed=4, atol=1e-6, rtol=1e-5)
    assert not hasattr(large, "permutations")


def test_attention_guard_blocks_large_default_and_allows_force():
    config = _config("/tmp")
    config["training"]["batch_size"] = 256
    config["model"]["attention_heads"] = 4
    with pytest.raises(ValueError, match="attention"):
        attention_guard(config, 2000)
    config["runtime"]["force_large_attention"] = True
    assert attention_guard(config, 2000)["warning"]
    config["training"]["batch_size"] = 16
    config["runtime"]["force_large_attention"] = False
    assert attention_guard(config, 2000)["attention_elements"] <= config["runtime"]["max_attention_elements"]


def test_trainer_forward_backward_and_generator_interval(tmp_path):
    rng = np.random.default_rng(1)
    x = np.log1p(rng.poisson(2.0, size=(8, 12)).astype(np.float32))
    config = _config(tmp_path)
    config["training"]["generator_update_interval"] = 2
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=12, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    assert trainer.full_x.device.type == "cpu"
    history = trainer.train()
    trainer.save_diagnostics()
    assert trainer.extract_embeddings(batch_size=4).shape == (8, 6)
    assert history["student_grad_norm"][-1] > 0
    assert history["generator_grad_norm"][-1] > 0
    assert history["generator_update_count"][-1] == 1.0
    stats = json.loads((tmp_path / "corruption_stats.json").read_text())
    assert "budget_deficit_rate" in stats


def test_artifact_manifest_required_files(tmp_path):
    config = _config(tmp_path)
    config.update({"dataset_name": "toy", "skip_eval": False})
    config["preprocessing"] = {"input_mode_resolved": "raw", "feature_space_source": "hvg_variance", "label_key": "cell_type", "labels_available": True}
    manifest = artifact_manifest(config, (3, 2), "kmeans_known_k")
    for name in manifest["required_files"]:
        if name != "artifact_manifest.json":
            (tmp_path / name).write_text("x", encoding="utf-8")
    validate_required_files(tmp_path, manifest)
    (tmp_path / "metrics.json").unlink()
    with pytest.raises(RuntimeError, match="missing"):
        validate_required_files(tmp_path, manifest)

    config["skip_eval"] = True
    manifest_skip = artifact_manifest(config, (3, 2), "prediction_only_no_labels")
    assert "labels.npy" not in manifest_skip["required_files"]
