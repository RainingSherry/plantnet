import inspect

import numpy as np
import torch

from methods.DeepLearning.CAAM_scMAE.data.donor_candidates import DonorCandidateProvider
from methods.DeepLearning.CAAM_scMAE.registry import build_arg_parser, resolve_config
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_arrays, toy_config, toy_provider
from methods.DeepLearning.CAAM_scMAE.trainers.common import CAAMTrainer, build_corruption


def test_cli_accepts_phase13_corruption_types():
    parser = build_arg_parser()
    for corruption_type in ("scmae_shuffle", "matched_donor", "nonzero_aware_donor"):
        args = parser.parse_args(
            [
                "--data_path",
                "toy.h5ad",
                "--save_dir",
                "out",
                "--n_clusters",
                "2",
                "--benchmark_mode",
                "true",
                "--variant",
                "control",
                "--corruption_type",
                corruption_type,
            ]
        )
        cfg = resolve_config(args)
        assert cfg["variant"] == "control"
        assert cfg["model"]["encoder_type"] == "mlp"
        assert cfg["model"]["mask_selector"] == "random"
        assert cfg["corruption"]["type"] == corruption_type


def test_scmae_shuffle_preserves_each_gene_marginal_distribution():
    x, batch, lib, zero = toy_arrays(n=10, g=6)
    provider = toy_provider(x, batch, lib, zero)
    full = torch.as_tensor(x)
    idx = torch.arange(x.shape[0])
    out = provider.sample_scmae_shuffle_batch(idx, full, torch.device("cpu"))
    replacement = out["replacement"].numpy()
    for gene_id in range(x.shape[1]):
        assert np.allclose(np.sort(replacement[:, gene_id]), np.sort(x[:, gene_id]))
    assert out["replacement_info"]["per_gene_permutation_seed"] == provider.scmae_permutation_seed


def test_matched_donor_never_selects_self_donor():
    x, batch, lib, zero = toy_arrays()
    provider = toy_provider(x, batch, lib, zero)
    idx = torch.tensor([0, 1, 2, 3])
    out = provider.sample_batch(idx, torch.as_tensor(x), torch.device("cpu"))
    assert torch.all(out["donor_indices"] != idx.view(-1, 1))
    assert torch.equal(out["mask_eligibility"], out["eligibility"])


def test_nonzero_aware_prioritizes_changed_donor():
    x = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 2.0],
            [0.0, 3.0, 0.0],
            [4.0, 5.0, 6.0],
        ],
        dtype=np.float32,
    )
    provider = DonorCandidateProvider(
        x,
        np.zeros(x.shape[0], dtype=np.int64),
        x.sum(axis=1).astype(np.float32),
        (x == 0).mean(axis=1).astype(np.float32),
        candidate_pool_size=2,
        library_size_bins=1,
        zero_ratio_bins=1,
        atol=1.0e-6,
        rtol=1.0e-5,
        seed=0,
    )
    provider.candidates[0] = np.asarray([1, 2], dtype=np.int64)
    out = provider.sample_nonzero_aware_batch(torch.tensor([0]), torch.as_tensor(x), torch.device("cpu"))
    assert torch.all(out["replacement"] != 0.0)
    assert torch.all(out["replacement_info"]["nonzero_aware_success"])
    assert not torch.any(out["replacement_info"]["fallback_to_matched"])
    assert not torch.any(out["replacement_info"]["fallback_to_scmae_shuffle"])


def test_corruption_fallback_apis_do_not_accept_labels():
    for method_name in ("sample_batch", "sample_scmae_shuffle_batch", "sample_nonzero_aware_batch"):
        signature = inspect.signature(getattr(DonorCandidateProvider, method_name))
        forbidden = {"label", "labels", "cell_type", "n_clusters", "true_cluster"}
        assert forbidden.isdisjoint(signature.parameters)


def test_all_phase13_corruptions_emit_required_corruption_stats_fields():
    generic = {
        "corruption_type",
        "zero_to_zero_rate",
        "effective_corruption_rate",
        "budget_deficit_rate",
        "mean_abs_delta",
        "mean_abs_delta_masked",
        "strict_effective_budget",
    }
    for corruption_type in ("scmae_shuffle", "matched_donor", "nonzero_aware_donor"):
        x, batch, lib, zero = toy_arrays()
        trainer = object.__new__(CAAMTrainer)
        trainer.config = toy_config("control")
        trainer.config["corruption"]["type"] = corruption_type
        trainer.full_x = torch.as_tensor(x)
        trainer.donor_provider = toy_provider(x, batch, lib, zero)
        trainer.corruption = build_corruption(corruption_type)
        trainer.corruption_totals = {
            "selected": 4.0,
            "changed_selected": 3.0,
            "zero_to_zero_selected": 1.0,
            "abs_delta_all": 8.0,
            "abs_delta_masked": 4.0,
            "positions": 16.0,
            "deficit_cells": 0.0,
            "cells": 4.0,
        }
        trainer.corruption_source_totals = {
            "nonzero_aware_success": 2.0,
            "fallback_to_matched": 1.0,
            "fallback_to_scmae_shuffle": 1.0,
        }
        stats = trainer.corruption_stats()
        assert generic.issubset(stats)
        if corruption_type == "scmae_shuffle":
            assert "per_gene_permutation_seed" in stats
        if corruption_type == "matched_donor":
            assert {"donor_fallback_matched", "donor_fallback_batch", "donor_fallback_global"}.issubset(stats)
        if corruption_type == "nonzero_aware_donor":
            assert {"nonzero_aware_success_rate", "fallback_to_matched_rate", "fallback_to_scmae_shuffle_rate"}.issubset(stats)
