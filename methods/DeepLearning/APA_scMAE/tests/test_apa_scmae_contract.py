from __future__ import annotations

import json
import anndata as ad
import numpy as np
import pandas as pd
import pytest
import torch

from methods.DeepLearning.APA_scMAE.config import build_arg_parser, resolve_config
from methods.DeepLearning.APA_scMAE.data import APAExpressionDataset, assert_no_training_labels, compute_gene_stats, load_apa_data
from methods.DeepLearning.APA_scMAE.losses import distortion_loss
from methods.DeepLearning.APA_scMAE.model import APAModel
from methods.DeepLearning.APA_scMAE.run import mapped_clustering_metrics
from methods.DeepLearning.APA_scMAE.trainer import APATrainer


def _config(tmp_path):
    return {
        "seed": 7,
        "method_name": "apa_scmae",
        "runtime": {"num_workers": 0, "fail_fast": True},
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
        "mask": {"ratio": 0.4, "temperature": 1.0, "masked_data_weight": 0.75},
        "corruption": {
            "type": "scmae_shuffle",
            "changed_tolerance_abs": 1.0e-6,
            "changed_tolerance_rel": 1.0e-5,
        },
        "generator_loss": {
            "lambda_entropy": 0.01,
            "lambda_balance": 0.01,
            "lambda_distortion": 0.01,
            "lambda_coverage": 0.01,
        },
        "save_dir": str(tmp_path),
    }


def test_dataset_returns_no_labels():
    dataset = APAExpressionDataset(np.ones((3, 5), dtype=np.float32))
    item = dataset[0]
    assert set(item) == {"index", "x"}
    with pytest.raises(RuntimeError):
        assert_no_training_labels({"x": torch.ones(2, 3), "labels": torch.zeros(2)})


def test_unlabeled_h5ad_fails_when_eval_required(tmp_path):
    data_path = tmp_path / "unlabeled.h5ad"
    ad.AnnData(X=np.ones((4, 5), dtype=np.float32), var=pd.DataFrame(index=[f"g{i}" for i in range(5)])).write_h5ad(data_path)
    with pytest.raises(ValueError, match="No label column"):
        load_apa_data(
            str(data_path),
            input_mode="auto",
            target_sum=10000.0,
            n_top_genes=0,
            scale_input=False,
            n_prototypes=2,
            pca_dim=2,
            seed=1,
            require_labels=True,
        )


def test_counts_layer_takes_priority_over_adata_x(tmp_path):
    data_path = tmp_path / "counts_priority.h5ad"
    counts = np.array([[1, 2, 0], [3, 0, 1], [0, 4, 2]], dtype=np.float32)
    misleading_x = np.full_like(counts, 99.0, dtype=np.float32)
    obs = pd.DataFrame({"cell_type": ["a", "b", "a"]})
    adata = ad.AnnData(X=misleading_x, obs=obs, var=pd.DataFrame(index=["g0", "g1", "g2"]))
    adata.layers["counts"] = counts
    adata.write_h5ad(data_path)
    bundle = load_apa_data(
        str(data_path),
        input_mode="auto",
        target_sum=10.0,
        n_top_genes=0,
        scale_input=False,
        n_prototypes=2,
        pca_dim=2,
        seed=1,
        require_labels=True,
    )
    totals = counts.sum(axis=1, keepdims=True)
    expected = np.log1p(counts / totals * 10.0)
    assert bundle.preprocess_config["raw_count_source"] == 'layers["counts"]'
    np.testing.assert_allclose(bundle.x, expected, rtol=1.0e-6, atol=1.0e-6)


def test_cyclic_cluster_labels_use_mapped_acc_and_f1():
    y_true = np.array([0, 1, 2, 0, 1, 2])
    y_pred = np.array([1, 2, 0, 1, 2, 0])
    metrics = mapped_clustering_metrics(y_true, y_pred)
    assert metrics["acc"] == pytest.approx(1.0)
    assert metrics["f1_macro"] == pytest.approx(1.0)


def test_yaml_seed_and_defaults_are_not_overridden_by_parser_defaults(tmp_path):
    config_path = tmp_path / "config.yaml"
    config_path.write_text(
        "\n".join(
            [
                "seed: 123",
                "method_name: yaml_method",
                "skip_eval: true",
                "runtime:",
                "  gpu: 3",
            ]
        ),
        encoding="utf-8",
    )
    parser = build_arg_parser()
    args = parser.parse_args(
        [
            "--config",
            str(config_path),
            "--data_path",
            "dummy.h5ad",
            "--save_dir",
            str(tmp_path / "out"),
            "--n_clusters",
            "2",
        ]
    )
    cfg = resolve_config(args)
    assert cfg["seed"] == 123
    assert cfg["method_name"] == "yaml_method"
    assert cfg["skip_eval"] is True
    assert cfg["runtime"]["gpu"] == 3


def test_distortion_loss_penalizes_extreme_delta_shortcut():
    delta = torch.tensor([[1.0, 100.0, 1.0, 100.0]])
    extreme = torch.tensor([[0.0, 1.0, 0.0, 1.0]])
    balanced = torch.tensor([[1.0, 1.0, 0.0, 0.0]])
    assert distortion_loss(extreme, delta) > distortion_loss(balanced, delta)


def test_model_shapes_and_mask_semantics():
    torch.manual_seed(0)
    x = torch.randn(4, 10).abs()
    replacement = x.roll(1, dims=0)
    effective = (~torch.isclose(x, replacement)).float()
    gene_stats = torch.as_tensor(compute_gene_stats(x.numpy()), dtype=torch.float32)
    prototypes = torch.randn(3, 5)
    model = APAModel(n_genes=10, token_dim=8, cell_dim=6, proto_dim=5, attention_heads=2, dropout=0.0)
    gene_vec, stat_vec = model.shared_context(gene_stats)
    gen = model.generator(
        x,
        replacement,
        effective,
        gene_vec,
        stat_vec,
        prototypes,
        mask_ratio=0.4,
        temperature=1.0,
    )
    assert gen["logits"].shape == (4, 10)
    assert gen["mask_hard"].shape == (4, 10)
    assert torch.all(gen["mask_hard"].sum(dim=1) == 4)
    m_eff = gen["mask_hard"] * effective
    out = model.student(x * (1.0 - gen["mask_hard"]) + replacement * gen["mask_hard"], gene_vec, stat_vec, prototypes)
    assert out["z"].shape == (4, 6)
    assert out["mask_logits"].shape == (4, 10)
    assert out["x_recon"].shape == (4, 10)
    assert torch.all(m_eff <= gen["mask_hard"])


def test_trainer_smoke_forward_backward(tmp_path):
    torch.manual_seed(1)
    rng = np.random.default_rng(1)
    x = rng.poisson(2.0, size=(8, 12)).astype(np.float32)
    x = np.log1p(x)
    gene_stats = compute_gene_stats(x)
    prototypes = rng.normal(size=(3, 4)).astype(np.float32)
    config = _config(tmp_path)
    model = APAModel(n_genes=12, token_dim=8, cell_dim=6, proto_dim=4, attention_heads=2, dropout=0.0)
    trainer = APATrainer(
        config=config,
        model=model,
        train_dataset=APAExpressionDataset(x),
        full_x=torch.as_tensor(x),
        gene_stats=torch.as_tensor(gene_stats),
        prototypes=torch.as_tensor(prototypes),
        device=torch.device("cpu"),
        save_dir=tmp_path,
    )
    history = trainer.train()
    trainer.save_diagnostics()
    embedding = trainer.extract_embeddings(batch_size=4)
    assert embedding.shape == (8, 6)
    assert history["student_grad_norm"][-1] > 0
    assert history["generator_grad_norm"][-1] > 0
    stats = json.loads((tmp_path / "corruption_stats.json").read_text())
    assert stats["corruption_type"] == "scmae_shuffle"
    assert 0.0 <= stats["effective_corruption_rate"] <= 1.0
