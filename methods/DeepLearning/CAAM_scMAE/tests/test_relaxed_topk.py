import torch

from methods.DeepLearning.CAAM_scMAE.mask_generator.relaxed_topk import relaxed_topk_straight_through


def test_relaxed_topk_hard_forward_soft_backward():
    logits = torch.randn(2, 6, requires_grad=True)
    eligibility = torch.ones_like(logits, dtype=torch.bool)
    k = torch.tensor([2, 3])
    hard, soft, st = relaxed_topk_straight_through(logits, k, 0.7, eligibility, add_gumbel=False)
    assert torch.all(hard.sum(dim=1) == k)
    assert torch.all((soft >= 0) & (soft <= 1.0001))
    assert torch.allclose(st.detach(), hard)
    assert st.requires_grad
    loss = (st * logits).sum()
    loss.backward()
    assert logits.grad is not None

