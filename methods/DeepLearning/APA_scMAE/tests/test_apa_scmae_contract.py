from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path

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
from methods.DeepLearning.APA_scMAE.losses import balance_loss, distortion_loss, moderate_delta_loss, student_losses
from methods.DeepLearning.APA_scMAE.model import APAModel
from methods.DeepLearning.APA_scMAE.representation_losses import (
    balanced_assignment_loss,
    covariance_loss,
    prototype_kl_loss,
    soft_assignment,
    variance_loss,
    vicreg_losses,
)
from methods.DeepLearning.APA_scMAE.run import (
    artifact_manifest,
    attention_guard,
    embedding_extraction_batch_size,
    hungarian_map,
    mapped_clustering_metrics,
    resolve_device,
    try_leiden_fixed,
    validate_runtime_config,
    validate_embedding_inputs,
    validate_required_files,
)
from methods.DeepLearning.APA_scMAE.trainer import APATrainer
from methods.DeepLearning.scNAME.run import _get_scname_count_input, _shuffle_scname_inputs


_DESC_LABEL_UTILS_PATH = Path(__file__).resolve().parents[2] / "desc" / "desc" / "models" / "label_utils.py"
_DESC_LABEL_SPEC = importlib.util.spec_from_file_location("desc_label_utils_for_test", _DESC_LABEL_UTILS_PATH)
assert _DESC_LABEL_SPEC is not None and _DESC_LABEL_SPEC.loader is not None
_DESC_LABEL_MODULE = importlib.util.module_from_spec(_DESC_LABEL_SPEC)
_DESC_LABEL_SPEC.loader.exec_module(_DESC_LABEL_MODULE)
remap_desc_labels = _DESC_LABEL_MODULE.remap_desc_labels


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
            "decoder_mode": "z_with_stopgrad_h",
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
            "student_warmup_epochs": 0,
            "enable_generator_after_warmup": True,
            "gamma": 0.5,
        },
        "mask": {"ratio": 0.4, "temperature": 1.0, "masked_data_weight": 0.75, "generator_topk_only_effective": True},
        "corruption": {"type": "scmae_shuffle", "changed_tolerance_abs": 1.0e-6, "changed_tolerance_rel": 1.0e-5},
        "representation_loss": {
            "use_repr_loss": True,
            "lambda_invariance": 1.0,
            "lambda_variance": 0.1,
            "lambda_covariance": 0.01,
            "variance_margin": 1.0,
            "eps": 1.0e-4,
        },
        "teacher": {
            "use_ema_teacher": True,
            "teacher_momentum_start": 0.9,
            "teacher_momentum_end": 0.99,
            "teacher_consistency_weight": 1.0,
        },
        "prototype_consistency": {
            "use_proto_consistency": True,
            "num_embedding_prototypes": 3,
            "proto_temperature": 0.2,
            "proto_loss_weight": 0.1,
            "proto_balance_weight": 0.01,
            "proto_ema_momentum": 0.95,
        },
        "generator_loss": {
            "lambda_entropy": 0.01,
            "lambda_balance": 0.01,
            "lambda_distortion": 0.01,
            "lambda_coverage": 0.01,
            "lambda_reconstruction_moderate": 0.1,
            "lambda_representation_moderate": 0.1,
        },
        "n_clusters": 3,
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


def test_skip_eval_without_explicit_label_key_skips_label_discovery(tmp_path):
    multi = tmp_path / "multi_skip_eval.h5ad"
    _write_h5ad(multi, np.ones((3, 4)), obs={"cell_type": ["a", "b", "a"], "maintype": ["x", "y", "x"]})
    bundle = load_apa_data(
        str(multi),
        input_mode="log1p",
        target_sum=1.0,
        n_top_genes=0,
        scale_input=False,
        n_prototypes=2,
        pca_dim=2,
        seed=1,
        require_labels=False,
    )
    assert bundle.labels is None
    assert bundle.label_key is None
    assert bundle.preprocess_config["labels_available"] is False

    explicit = load_apa_data(
        str(multi),
        input_mode="log1p",
        target_sum=1.0,
        n_top_genes=0,
        scale_input=False,
        n_prototypes=2,
        pca_dim=2,
        seed=1,
        require_labels=False,
        label_key="maintype",
    )
    assert explicit.labels is not None
    assert explicit.label_key == "maintype"


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
    with pytest.raises(ValueError, match="no valid count-like raw source"):
        load_apa_data(str(log_path), input_mode="raw", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    auto_log = load_apa_data(str(log_path), input_mode="auto", target_sum=10.0, n_top_genes=2, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert auto_log.preprocess_config["input_mode_resolved"] == "log1p"
    assert auto_log.preprocess_config["feature_space_source"] == "hvg_variance"


def test_raw_source_rejects_log_normalized_raw_in_auto_and_raw_modes(tmp_path):
    counts = np.array([[1, 2], [3, 0], [0, 4]], dtype=np.float32)
    log_x = np.log1p(counts / counts.sum(axis=1, keepdims=True) * 10.0)
    obs = {"cell_type": ["a", "b", "a"]}

    data_path = tmp_path / "raw_log1p.h5ad"
    _write_h5ad(data_path, log_x, obs=obs, raw_x=log_x)
    bundle = load_apa_data(str(data_path), input_mode="auto", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert bundle.preprocess_config["input_mode_resolved"] == "log1p"
    assert bundle.preprocess_config["raw_count_source"] is None

    with pytest.raises(ValueError, match="no valid count-like raw source"):
        load_apa_data(str(data_path), input_mode="raw", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)

    counts_path = tmp_path / "counts_override_bad_raw.h5ad"
    _write_h5ad(counts_path, log_x, obs=obs, layers={"counts": counts}, raw_x=log_x)
    raw_bundle = load_apa_data(str(counts_path), input_mode="raw", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert raw_bundle.preprocess_config["raw_count_source"] == 'layers["counts"]'


def test_log1p_mode_does_not_probe_raw_source(tmp_path, monkeypatch):
    counts = np.array([[1, 2], [3, 0], [0, 4]], dtype=np.float32)
    log_x = np.log1p(counts)
    data_path = tmp_path / "log1p_no_probe.h5ad"
    _write_h5ad(data_path, log_x, obs={"cell_type": ["a", "b", "a"]}, raw_x=counts)

    import methods.DeepLearning.APA_scMAE.data as apa_data

    def fail_if_called(_adata):
        raise AssertionError("_raw_count_source should not be called for input_mode='log1p'")

    monkeypatch.setattr(apa_data, "_raw_count_source", fail_if_called)
    bundle = apa_data.load_apa_data(
        str(data_path),
        input_mode="log1p",
        target_sum=10.0,
        n_top_genes=0,
        scale_input=False,
        n_prototypes=2,
        pca_dim=2,
        seed=1,
    )
    assert bundle.preprocess_config["input_mode_resolved"] == "log1p"


def test_layers_counts_priority_requires_count_like_values(tmp_path):
    counts = np.array([[1, 2], [3, 0], [0, 4]], dtype=np.float32)
    bad_counts_layer = np.log1p(counts)
    obs = {"cell_type": ["a", "b", "a"]}

    valid_path = tmp_path / "valid_counts_layer.h5ad"
    _write_h5ad(valid_path, np.log1p(counts), obs=obs, layers={"counts": counts}, raw_x=counts + 10)
    valid = load_apa_data(str(valid_path), input_mode="auto", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert valid.preprocess_config["raw_count_source"] == 'layers["counts"]'

    invalid_path = tmp_path / "invalid_counts_layer.h5ad"
    _write_h5ad(invalid_path, np.log1p(counts), obs=obs, layers={"counts": bad_counts_layer})
    invalid = load_apa_data(str(invalid_path), input_mode="auto", target_sum=10.0, n_top_genes=0, scale_input=False, n_prototypes=2, pca_dim=2, seed=1)
    assert invalid.preprocess_config["input_mode_resolved"] == "log1p"
    assert invalid.preprocess_config["raw_count_source"] is None


def test_gene_stats_zero_rate_uses_prescale_values(tmp_path):
    counts = np.array([[0, 1, 0], [0, 2, 3], [4, 0, 0]], dtype=np.float32)
    data_path = tmp_path / "scale.h5ad"
    _write_h5ad(data_path, counts, obs={"cell_type": ["a", "b", "a"]}, layers={"counts": counts})
    bundle = load_apa_data(str(data_path), input_mode="raw", target_sum=10.0, n_top_genes=0, scale_input=True, n_prototypes=2, pca_dim=2, seed=1)
    assert bundle.preprocess_config["gene_stats_source"] == "pre_scale_log1p"
    assert np.any(bundle.gene_stats[:, 2] > 0)


def test_config_defaults_do_not_override_yaml_and_no_cuda_works(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text("seed: 123\nmethod_name: yaml_method\nskip_eval: true\nn_clusters: 5\nruntime:\n  gpu: 3\n", encoding="utf-8")
    parser = build_arg_parser()
    args = parser.parse_args(["--config", str(config_path), "--data_path", "dummy.h5ad", "--save_dir", str(tmp_path / "out")])
    cfg = resolve_config(args)
    assert cfg["seed"] == 123
    assert cfg["method_name"] == "yaml_method"
    assert cfg["skip_eval"] is True
    assert cfg["n_clusters"] == 5
    assert cfg["runtime"]["gpu"] == 3

    args = parser.parse_args(["--config", str(config_path), "--data_path", "dummy.h5ad", "--save_dir", str(tmp_path / "out_override"), "--n_clusters", "2"])
    cfg = resolve_config(args)
    assert cfg["n_clusters"] == 2

    args = parser.parse_args(["--data_path", "dummy.h5ad", "--save_dir", str(tmp_path / "out2"), "--n_clusters", "2", "--no_cuda"])
    cfg = resolve_config(args)
    assert cfg["runtime"]["no_cuda"] is True

    args = parser.parse_args(
        [
            "--data_path",
            "dummy.h5ad",
            "--save_dir",
            str(tmp_path / "out3"),
            "--n_clusters",
            "2",
            "--use_repr_loss",
            "false",
            "--use_ema_teacher",
            "false",
            "--use_proto_consistency",
            "false",
            "--student_warmup_epochs",
            "0",
            "--decoder_mode",
            "current",
        ]
    )
    cfg = resolve_config(args)
    assert cfg["representation_loss"]["use_repr_loss"] is False
    assert cfg["teacher"]["use_ema_teacher"] is False
    assert cfg["prototype_consistency"]["use_proto_consistency"] is False
    assert cfg["training"]["student_warmup_epochs"] == 0
    assert cfg["model"]["decoder_mode"] == "current"

    cfg["skip_eval"] = True
    cfg["n_clusters"] = 0
    with pytest.raises(ValueError, match="skip_eval=true requires n_clusters"):
        validate_runtime_config(cfg)

    with pytest.raises(SystemExit):
        parser.parse_args(
            [
                "--data_path",
                "dummy.h5ad",
                "--save_dir",
                str(tmp_path / "out4"),
                "--use_random_mask_during_warmup",
                "false",
            ]
        )


def test_resolve_device_maps_visible_gpu_to_logical_index(monkeypatch):
    import torch

    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "1,2")
    config = resolve_config(
        build_arg_parser().parse_args(["--data_path", "dummy.h5ad", "--save_dir", "/tmp/out", "--n_clusters", "2", "--gpu", "2"])
    )
    device, runtime = resolve_device(config)
    assert str(device) == "cuda:1"
    assert runtime["physical_gpu"] == 2
    assert runtime["logical_device"] == "cuda:1"

    config["runtime"]["gpu"] = 3
    with pytest.raises(ValueError, match="CUDA_VISIBLE_DEVICES"):
        resolve_device(config)


def test_resolve_device_uses_first_visible_gpu_without_explicit_gpu(monkeypatch):
    import torch

    monkeypatch.setattr(torch.cuda, "is_available", lambda: True)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "2")
    config = resolve_config(build_arg_parser().parse_args(["--data_path", "dummy.h5ad", "--save_dir", "/tmp/out", "--n_clusters", "2"]))
    assert config["runtime"]["gpu"] == 1
    assert config["runtime"]["gpu_explicit"] is False
    device, runtime = resolve_device(config)
    assert str(device) == "cuda:0"
    assert runtime["physical_gpu"] == 2
    assert runtime["logical_device"] == "cuda:0"
    assert runtime["gpu_explicit"] is False


def test_cyclic_and_extra_cluster_metrics():
    metrics = mapped_clustering_metrics(np.array([0, 1, 2, 0, 1, 2]), np.array([1, 2, 0, 1, 2, 0]))
    assert metrics["acc"] == pytest.approx(1.0)
    assert metrics["f1_macro"] == pytest.approx(1.0)

    mapped = hungarian_map(np.array([0, 0, 1, 1]), np.array([0, 1, 2, 3]))
    assert np.any(mapped == -1)
    metrics = mapped_clustering_metrics(np.array([0, 0, 1, 1]), np.array([0, 1, 2, 3]))
    assert metrics["acc"] < 1.0
    assert metrics["f1_macro"] < 1.0


def test_leiden_metrics_are_recomputed_with_apa_mapper(monkeypatch):
    labels = np.array([0, 0, 1, 1])
    pred = np.array([0, 1, 2, 3])
    embedding = np.arange(8, dtype=np.float32).reshape(4, 2)
    fake_module = types.ModuleType("methods.DeepLearning.CAAM_scMAE.evaluation.local_metrics")
    fake_module.leiden_fixed = lambda *_args, **_kwargs: ({"acc": 1.0, "f1_macro": 1.0}, pred)
    monkeypatch.setitem(sys.modules, fake_module.__name__, fake_module)

    result = try_leiden_fixed(
        embedding,
        labels,
        {"seed": 7, "evaluation": {"leiden_fixed_resolution": 1.0, "n_neighbors": 3}},
    )

    assert result["status"] == "success"
    assert result["metrics"]["acc"] < 1.0
    assert result["metrics"]["f1_macro"] < 1.0
    assert result["metrics"]["cluster_method"] == "leiden_fixed"


def test_leiden_unexpected_exceptions_surface(monkeypatch):
    fake_module = types.ModuleType("methods.DeepLearning.CAAM_scMAE.evaluation.local_metrics")

    def fail(*_args, **_kwargs):
        raise RuntimeError("evaluation bug")

    fake_module.leiden_fixed = fail
    monkeypatch.setitem(sys.modules, fake_module.__name__, fake_module)
    with pytest.raises(RuntimeError, match="evaluation bug"):
        try_leiden_fixed(
            np.ones((4, 2), dtype=np.float32),
            np.array([0, 0, 1, 1]),
            {"seed": 7, "evaluation": {"leiden_fixed_resolution": 1.0, "n_neighbors": 3}},
        )


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
    assert moderate_delta_loss(torch.tensor([[0.0, 1.0, 0.0, 1.0]]), delta) > moderate_delta_loss(torch.tensor([[1.0, 1.0, 0.0, 0.0]]), delta)
    concentrated = torch.tensor([[1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    balanced = torch.tensor([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    assert balance_loss(concentrated) > balance_loss(balanced)


def test_representation_losses_are_finite_and_penalize_collapse():
    z_clean = torch.randn(8, 5)
    z_masked = z_clean + 0.01 * torch.randn(8, 5)
    losses = vicreg_losses(z_clean, z_masked)
    assert all(torch.isfinite(value) for value in losses.values())
    collapsed = torch.zeros(8, 5)
    spread = torch.randn(8, 5) * 3.0
    assert variance_loss(collapsed) > variance_loss(spread)
    cov = covariance_loss(z_clean)
    assert cov.ndim == 0
    assert torch.isfinite(cov)


def test_prototype_consistency_shapes_gradients_and_balance():
    z_clean = torch.randn(6, 4)
    z_masked = torch.randn(6, 4, requires_grad=True)
    prototypes = torch.randn(3, 4)
    q_clean = soft_assignment(z_clean, prototypes, temperature=0.5)
    q_masked = soft_assignment(z_masked, prototypes, temperature=0.5)
    assert q_clean.shape == (6, 3)
    assert q_masked.shape == (6, 3)
    loss = prototype_kl_loss(q_clean, q_masked)
    loss.backward()
    assert z_masked.grad is not None
    assert torch.isfinite(z_masked.grad).all()
    collapsed = torch.tensor([[1.0, 0.0, 0.0]]).repeat(6, 1)
    balanced = torch.eye(3).repeat(2, 1)
    assert balanced_assignment_loss(collapsed) > balanced_assignment_loss(balanced)


def test_decoder_modes_forward_shapes():
    x = torch.ones(2, 6)
    gene_stats = torch.as_tensor(compute_gene_stats(x.numpy()), dtype=torch.float32)
    prototypes = torch.randn(2, 3)
    for decoder_mode in ["current", "z_only", "z_with_stopgrad_h"]:
        model = APAModel(
            n_genes=6,
            token_dim=6,
            cell_dim=4,
            proto_dim=3,
            attention_heads=2,
            dropout=0.0,
            decoder_mode=decoder_mode,
        )
        gene_vec, stat_vec = model.shared_context(gene_stats)
        out = model.student(x, gene_vec, stat_vec, prototypes)
        assert out["z"].shape == (2, 4)
        assert out["mask_logits"].shape == (2, 6)
        assert out["x_recon"].shape == (2, 6)
        assert model.student.encode_clean(x, gene_vec, stat_vec, prototypes).shape == (2, 4)
        assert model.student.encode_masked(x, gene_vec, stat_vec, prototypes).shape == (2, 4)


def test_decoder_bottleneck_detaches_mask_logits_for_v2_modes():
    x = torch.ones(2, 6)
    gene_stats = torch.as_tensor(compute_gene_stats(x.numpy()), dtype=torch.float32)
    prototypes = torch.randn(2, 3)
    for decoder_mode in ["z_only", "z_with_stopgrad_h"]:
        model = APAModel(n_genes=6, token_dim=6, cell_dim=4, proto_dim=3, attention_heads=2, dropout=0.0, decoder_mode=decoder_mode)
        gene_vec, stat_vec = model.shared_context(gene_stats)
        out = model.student(x, gene_vec, stat_vec, prototypes)
        out["x_recon"].sum().backward()
        assert all(param.grad is None for param in model.student.mask_head.parameters())

    model = APAModel(n_genes=6, token_dim=6, cell_dim=4, proto_dim=3, attention_heads=2, dropout=0.0, decoder_mode="current")
    gene_vec, stat_vec = model.shared_context(gene_stats)
    out = model.student(x, gene_vec, stat_vec, prototypes)
    out["x_recon"].sum().backward()
    assert any(param.grad is not None for param in model.student.mask_head.parameters())


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


def test_random_mask_respects_effective_eligibility_flag(tmp_path):
    rng = np.random.default_rng(6)
    x = np.ones((4, 6), dtype=np.float32)
    config = _config(tmp_path)
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=6, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    effective = torch.tensor([[1, 0, 0, 0, 0, 0]], dtype=torch.float32)
    random_mask = trainer._random_mask(torch.ones_like(effective), torch.zeros_like(effective), effective)
    assert torch.all(random_mask["mask_hard"] <= effective)
    assert random_mask["budget_deficit"][0] > 0

    trainer.config["mask"]["generator_topk_only_effective"] = False
    random_mask = trainer._random_mask(torch.ones_like(effective), torch.zeros_like(effective), effective)
    assert random_mask["mask_hard"].sum() == int(float(config["mask"]["ratio"]) * effective.shape[1])
    assert random_mask["budget_deficit"][0] == 0


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


def test_ema_teacher_is_separate_frozen_and_updates(tmp_path):
    rng = np.random.default_rng(11)
    x = np.log1p(rng.poisson(2.0, size=(8, 10)).astype(np.float32))
    config = _config(tmp_path)
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=10, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    assert trainer.teacher is not None
    student_param = next(trainer.model.student.parameters())
    teacher_param = next(trainer.teacher.student.parameters())
    assert teacher_param is not student_param
    assert all(not param.requires_grad for param in trainer.teacher.parameters())
    before_student = [param.detach().clone() for param in trainer.teacher.student.parameters()]
    before_shared = [param.detach().clone() for param in trainer.teacher.shared.parameters()]
    batch = next(iter(trainer._loader()))
    trainer._student_step(batch, warmup=False, generator_enabled=True)
    assert all(param.grad is None for param in trainer.teacher.parameters())
    after_student = [param.detach().clone() for param in trainer.teacher.student.parameters()]
    after_shared = [param.detach().clone() for param in trainer.teacher.shared.parameters()]
    assert any(not torch.allclose(a, b) for a, b in zip(after_student, before_student))
    assert any(not torch.allclose(a, b) for a, b in zip(after_shared, before_shared))


def test_trainer_forces_generator_update_when_interval_exceeds_batches(tmp_path):
    rng = np.random.default_rng(2)
    x = np.log1p(rng.poisson(2.0, size=(6, 10)).astype(np.float32))
    config = _config(tmp_path)
    config["training"]["generator_update_interval"] = 99
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=10, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    history = trainer.train()
    assert history["generator_update_count"][-1] >= 1.0
    assert history["generator_update_forced"][-1] == 1.0


def test_warmup_skips_generator_and_zero_warmup_enables_it(tmp_path):
    rng = np.random.default_rng(3)
    x = np.log1p(rng.poisson(2.0, size=(8, 10)).astype(np.float32))
    config = _config(tmp_path)
    config["training"]["student_warmup_epochs"] = 1
    config["training"]["epochs"] = 1
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=10, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    assert trainer.embedding_prototypes is None
    history = trainer.train()
    assert history["warmup_epoch"][-1] == 1.0
    assert history["generator_enabled"][-1] == 0.0
    assert history["generator_update_count"][-1] == 0.0
    assert history["loss_generator"][-1] is None
    assert history["generator_grad_norm"][-1] is None
    assert trainer.embedding_prototypes is None

    no_warmup = _config(tmp_path)
    no_warmup["training"]["student_warmup_epochs"] = 0
    no_warmup["training"]["epochs"] = 1
    trainer = APATrainer(
        config=no_warmup,
        model=APAModel(n_genes=10, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    history = trainer.train()
    assert history["warmup_epoch"][-1] == 0.0
    assert history["generator_enabled"][-1] == 1.0
    assert history["generator_update_count"][-1] > 0.0
    assert trainer.embedding_prototypes is not None


def test_warmup_never_calls_generator_forward(tmp_path, monkeypatch):
    rng = np.random.default_rng(8)
    x = np.log1p(rng.poisson(2.0, size=(8, 10)).astype(np.float32))
    config = _config(tmp_path)
    config["training"]["epochs"] = 1
    config["training"]["student_warmup_epochs"] = 1
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=10, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )

    def fail_generator_forward(*_args, **_kwargs):
        raise AssertionError("_generator_forward should not run during warmup")

    monkeypatch.setattr(trainer, "_generator_forward", fail_generator_forward)
    history = trainer.train()
    assert history["warmup_epoch"] == [1.0]
    assert history["generator_update_count"] == [0.0]


def test_generator_disabled_after_warmup_uses_random_mask_without_generator_forward(tmp_path, monkeypatch):
    rng = np.random.default_rng(7)
    x = np.log1p(rng.poisson(2.0, size=(8, 10)).astype(np.float32))
    config = _config(tmp_path)
    config["training"]["epochs"] = 2
    config["training"]["student_warmup_epochs"] = 1
    config["training"]["enable_generator_after_warmup"] = False
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=10, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )

    def fail_generator_forward(*_args, **_kwargs):
        raise AssertionError("_generator_forward should not run when generator is disabled after warmup")

    monkeypatch.setattr(trainer, "_generator_forward", fail_generator_forward)
    history = trainer.train()
    assert history["warmup_epoch"] == [1.0, 0.0]
    assert history["generator_enabled"] == [0.0, 0.0]
    assert history["generator_update_count"] == [0.0, 0.0]
    assert history["loss_generator"] == [None, None]


def test_embedding_prototypes_initialize_after_warmup_and_fallback_to_prototype_count(tmp_path):
    rng = np.random.default_rng(4)
    x = np.log1p(rng.poisson(2.0, size=(10, 8)).astype(np.float32))
    config = _config(tmp_path)
    config["n_clusters"] = 7
    config["prototype"] = {"n_prototypes": 2}
    config["prototype_consistency"]["num_embedding_prototypes"] = 0
    config["training"]["epochs"] = 2
    config["training"]["student_warmup_epochs"] = 1
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=8, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    assert trainer.embedding_prototypes is None
    history = trainer.train()
    assert history["warmup_epoch"] == [1.0, 0.0]
    assert trainer.embedding_prototypes is not None
    assert trainer.embedding_prototypes.shape[0] == 2
    assert history["loss_proto"][0] == pytest.approx(0.0)


def test_all_representation_mechanisms_off_does_not_compute_clean_encode(tmp_path, monkeypatch):
    rng = np.random.default_rng(5)
    x = np.log1p(rng.poisson(2.0, size=(8, 8)).astype(np.float32))
    config = _config(tmp_path)
    config["representation_loss"]["use_repr_loss"] = False
    config["teacher"]["use_ema_teacher"] = False
    config["prototype_consistency"]["use_proto_consistency"] = False
    config["training"]["student_warmup_epochs"] = 0
    config["training"]["epochs"] = 1
    trainer = APATrainer(
        config=config,
        model=APAModel(n_genes=8, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0),
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(compute_gene_stats(x)),
        prototypes=torch.as_tensor(rng.normal(size=(3, 4)).astype(np.float32)),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )

    def fail_encode_clean(*_args, **_kwargs):
        raise AssertionError("encode_clean should not run when repr/teacher/proto are all disabled")

    monkeypatch.setattr(trainer.model.student, "encode_clean", fail_encode_clean)
    history = trainer.train()
    assert history["loss_repr"][-1] == pytest.approx(0.0)
    assert history["loss_teacher"][-1] == pytest.approx(0.0)
    assert history["loss_proto"][-1] == pytest.approx(0.0)


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


def test_embedding_extraction_uses_training_batch_size():
    config = _config("/tmp")
    config["training"]["batch_size"] = 17
    assert embedding_extraction_batch_size(config) == 17


def test_desc_remap_labels_are_dense_integer_without_nan():
    mapped = remap_desc_labels([2, 2, 4])
    np.testing.assert_array_equal(mapped, np.array([0, 0, 1], dtype=np.int64))
    assert np.issubdtype(mapped.dtype, np.integer)
    assert not pd.Series(mapped).isna().any()
    assert pd.Series(mapped, dtype=np.int64).astype(int).tolist() == [0, 0, 1]


def test_scname_rejects_scaled_X_and_norm_log_without_raw_counts():
    counts = np.array([[1, 2], [0, 3], [4, 0]], dtype=np.float32)
    scaled_negative = np.array([[-1.0, 0.2], [0.1, -0.4], [0.3, 0.5]], dtype=np.float32)
    obs = pd.DataFrame(index=[f"c{i}" for i in range(3)])
    var = pd.DataFrame(index=["g0", "g1"])

    adata = ad.AnnData(X=scaled_negative, obs=obs, var=var)
    with pytest.raises(ValueError, match="scaled adata.X"):
        _get_scname_count_input(adata)

    adata = ad.AnnData(X=np.log1p(counts), obs=obs.copy(), var=var.copy())
    adata.layers["norm_log"] = np.log1p(counts)
    with pytest.raises(ValueError, match="norm_log"):
        _get_scname_count_input(adata)


def test_scname_layers_counts_source_is_accepted():
    counts = np.array([[1, 2], [0, 3], [4, 0]], dtype=np.float32)
    scaled_negative = np.array([[-1.0, 0.2], [0.1, -0.4], [0.3, 0.5]], dtype=np.float32)
    obs = pd.DataFrame(index=[f"c{i}" for i in range(3)])
    var = pd.DataFrame(index=["g0", "g1"])

    adata = ad.AnnData(X=scaled_negative, obs=obs, var=var)
    adata.layers["counts"] = counts
    count_input = _get_scname_count_input(adata)
    assert count_input.source == "layers_counts"
    np.testing.assert_array_equal(count_input.matrix, counts)
    np.testing.assert_array_equal(count_input.size_factor_counts, counts)


def test_scname_raw_counts_use_hvg_subset_and_full_counts_for_size_factor():
    counts = np.array([[1, 2], [0, 3], [4, 0]], dtype=np.float32)
    scaled_negative = np.array([[-1.0, 0.2], [0.1, -0.4], [0.3, 0.5]], dtype=np.float32)
    obs = pd.DataFrame(index=[f"c{i}" for i in range(3)])
    hvg_var = pd.DataFrame(index=["g0", "g1"])
    raw_var = pd.DataFrame(index=["g0", "g1", "g2"])
    raw_counts = np.array([[11, 12, 20], [10, 13, 30], [14, 10, 40]], dtype=np.float32)

    adata = ad.AnnData(X=scaled_negative, obs=obs, var=hvg_var.copy())
    adata.raw = ad.AnnData(X=np.log1p(counts), obs=obs.copy(), var=hvg_var.copy())
    with pytest.raises(ValueError, match="not count-like"):
        _get_scname_count_input(adata)

    raw_only = ad.AnnData(X=scaled_negative, obs=obs, var=hvg_var.copy())
    raw_only.raw = ad.AnnData(X=raw_counts, obs=obs.copy(), var=raw_var.copy())
    raw_input = _get_scname_count_input(raw_only)
    assert raw_input.source == "adata_raw_counts"
    np.testing.assert_array_equal(raw_input.matrix, raw_counts[:, :2])
    np.testing.assert_array_equal(raw_input.size_factor_counts, raw_counts)


def test_scname_shuffle_keeps_X_Y_count_like_sf_aligned():
    X = np.array([[10], [20], [30]], dtype=np.float32)
    Y = np.array([1, 2, 3])
    count_like = np.array([[100], [200], [300]], dtype=np.float32)
    sf = np.array([[1.0], [2.0], [3.0]], dtype=np.float32)
    shuffle_ix = np.array([2, 0, 1])
    X_s, Y_s, count_s, sf_s = _shuffle_scname_inputs(X, Y, count_like, sf, shuffle_ix)
    np.testing.assert_array_equal(X_s.ravel(), np.array([30, 10, 20], dtype=np.float32))
    np.testing.assert_array_equal(Y_s, np.array([3, 1, 2]))
    np.testing.assert_array_equal(count_s.ravel(), np.array([300, 100, 200], dtype=np.float32))
    np.testing.assert_array_equal(sf_s.ravel(), np.array([3.0, 1.0, 2.0], dtype=np.float32))
