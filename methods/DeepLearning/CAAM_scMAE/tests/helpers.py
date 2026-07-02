from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import torch

from methods.DeepLearning.CAAM_scMAE.data.donor_candidates import DonorCandidateProvider
from methods.DeepLearning.CAAM_scMAE.data.gene_modules import normalized_assignment_dense
from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student


def toy_arrays(n: int = 12, g: int = 8):
    rng = np.random.default_rng(0)
    x = rng.poisson(2.0, size=(n, g)).astype(np.float32)
    x = np.log1p(x)
    batch = np.arange(n) % 2
    lib = x.sum(axis=1).astype(np.float32)
    zero = (x == 0).mean(axis=1).astype(np.float32)
    return x, batch.astype(np.int64), lib, zero


def toy_provider(x, batch, lib, zero):
    return DonorCandidateProvider(
        x,
        batch,
        lib,
        zero,
        candidate_pool_size=4,
        library_size_bins=2,
        zero_ratio_bins=2,
        atol=1.0e-6,
        rtol=1.0e-5,
        seed=0,
    )


def toy_config(variant: str = "control"):
    from methods.DeepLearning.CAAM_scMAE.registry import DEFAULT_CONFIG, apply_variant
    import copy

    cfg = copy.deepcopy(DEFAULT_CONFIG)
    cfg["seed"] = 0
    cfg["variant"] = variant
    cfg["training"]["epochs"] = 1
    cfg["training"]["batch_size"] = 6
    cfg["training"]["student_warmup_epochs"] = 0
    cfg["training"]["generator_update_interval"] = 1
    cfg["model"]["latent_dim"] = 4
    cfg["model"]["mlp_hidden_dim"] = 8
    cfg["axial"]["n_gene_modules"] = 4
    cfg["axial"]["token_dim"] = 8
    cfg["axial"]["gene_attention_heads"] = 2
    cfg["axial"]["gene_attention_layers"] = 1
    cfg["axial"]["context_size"] = 4
    cfg["generator"]["hidden_dim"] = 8
    cfg["mask"]["ratio"] = 0.25
    cfg["runtime"]["fail_fast"] = True
    return apply_variant(cfg, variant)


def toy_assignment(g: int = 8, m: int = 4):
    rows = np.arange(g)
    cols = rows % m
    mat = sp.csr_matrix((np.ones(g), (rows, cols)), shape=(g, m), dtype=np.float32)
    return normalized_assignment_dense(mat)

