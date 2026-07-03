import torch

from methods.DeepLearning.CAAM_scMAE.mask_generator.random_mask import RandomFixedBudgetMask
from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_assignment, toy_config


def test_mlp_model_shapes():
    cfg = toy_config("control")
    model = build_student(n_genes=8, config=cfg)
    x = torch.randn(3, 8)
    mask = torch.zeros_like(x)
    out = model(x, mask=mask)
    assert out["z"].shape == (3, 4)
    assert out["x_hat"].shape == (3, 8)
    assert out["mask_logits"].shape == (3, 8)


def test_axial_model_shapes():
    cfg = toy_config("axial")
    model = build_student(n_genes=8, config=cfg, assignment=toy_assignment())
    context_idx = torch.tensor([0, 1, 2, 3])
    model.refresh_context_cache(torch.randn(4, 8), context_idx)
    x = torch.randn(3, 8)
    out = model(x, mask=torch.zeros_like(x), indices=torch.tensor([0, 5, 6]))
    assert out["z"].shape == (3, 4)
    assert out["module_tokens"].shape == (3, 4, 8)
    assert out["gene_attn"].shape[-2:] == (4, 4)
    assert out["cell_attn"].shape[:2] == (3, 4)


def test_random_mask_budget_shape():
    x = torch.zeros(4, 10)
    eligibility = torch.ones_like(x, dtype=torch.bool)
    _, mask, _ = RandomFixedBudgetMask(0.3)(x, eligibility)
    assert mask.shape == x.shape
    assert torch.all(mask.sum(dim=1) == 3)

