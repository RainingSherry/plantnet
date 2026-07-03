import torch

from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student
from methods.DeepLearning.CAAM_scMAE.tests.helpers import toy_assignment, toy_config


def test_context_self_exclusion_attention_zero():
    cfg = toy_config("axial")
    model = build_student(n_genes=8, config=cfg, assignment=toy_assignment())
    context_idx = torch.tensor([0, 1, 2, 3])
    model.refresh_context_cache(torch.randn(4, 8), context_idx)
    out = model(torch.randn(2, 8), mask=torch.zeros(2, 8), indices=torch.tensor([1, 5]))
    attn = out["cell_attn"]
    assert torch.allclose(attn[0, :, :, 1], torch.zeros_like(attn[0, :, :, 1]), atol=1.0e-6)

