import torch

from methods.DeepLearning.CAAM_scMAE.mask_generator.relaxed_topk import relaxed_topk_straight_through


def test_relaxed_topk_hard_forward_soft_backward():
    logits = torch.randn(2, 6, requires_grad=True)
    eligibility = torch.ones_like(logits, dtype=torch.bool)
    k = torch.tensor([2, 3])
    hard, soft, st = relaxed_topk_straight_through(logits, k, 0.7, eligibility, add_gumbel=False)
    assert torch.all(hard.sum(dim=1) == k)
    assert torch.all((soft >= 0) & (soft <= 1.0001))
    assert torch.allclose(soft.sum(dim=1), k.float(), atol=1.0e-3)
    assert torch.allclose(st.detach(), hard)
    assert st.requires_grad
    loss = (st * logits).sum()
    loss.backward()
    assert logits.grad is not None


def test_relaxed_topk_respects_partial_eligibility_and_deficit():
    logits = torch.tensor(
        [[0.1, 0.2, 9.0, 8.0, 0.3, -1.0], [0.6, 0.4, 10.0, 9.0, 0.2, 8.0]],
        requires_grad=True,
    )
    eligibility = torch.tensor(
        [[True, True, False, False, True, False], [True, True, False, False, True, False]]
    )
    k = torch.tensor([2, 4])
    hard, soft, st = relaxed_topk_straight_through(logits, k, 0.7, eligibility, add_gumbel=False)
    expected_budget = torch.tensor([2.0, 3.0])
    assert torch.all(hard[:, ~eligibility[0]] == 0)
    assert torch.all(soft[:, ~eligibility[0]] == 0)
    assert torch.allclose(hard.sum(dim=1), expected_budget)
    assert torch.allclose(soft.sum(dim=1), expected_budget, atol=1.0e-3)
    loss = (st * logits).sum()
    loss.backward()
    assert logits.grad is not None
