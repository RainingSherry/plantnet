import torch

from methods.DeepLearning.CAAM_scMAE.mask_generator.random_mask import RandomFixedBudgetMask


def test_ineligible_positions_not_masked_and_deficit_recorded():
    x = torch.zeros(2, 10)
    eligibility = torch.zeros_like(x, dtype=torch.bool)
    eligibility[0, :3] = True
    eligibility[1, :10] = True
    _, mask, info = RandomFixedBudgetMask(0.5)(x, eligibility)
    assert torch.all(mask[~eligibility] == 0)
    assert mask[0].sum() == 3
    assert mask[1].sum() == 5
    assert info["budget_deficit"][0] == 2

