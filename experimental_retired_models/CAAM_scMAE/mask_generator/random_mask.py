from __future__ import annotations

import torch


class RandomFixedBudgetMask:
    def __init__(self, mask_ratio: float) -> None:
        self.mask_ratio = float(mask_ratio)

    def __call__(self, x: torch.Tensor, eligibility: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, dict]:
        b, g = x.shape
        base_k = int(round(self.mask_ratio * g))
        eligible_count = eligibility.long().sum(dim=1)
        k_i = torch.minimum(eligible_count, torch.full_like(eligible_count, base_k))
        scores = torch.rand((b, g), device=x.device).masked_fill(~eligibility.bool(), -1.0)
        mask = torch.zeros((b, g), dtype=torch.float32, device=x.device)
        for row in range(b):
            k = int(k_i[row].item())
            if k > 0:
                top = torch.topk(scores[row], k=k).indices
                mask[row, top] = 1.0
        logits = torch.zeros_like(x)
        deficit = (base_k - k_i).clamp_min(0)
        info = {
            "mask_type": "random_fixed_budget",
            "budget_per_cell": int(base_k),
            "budget_deficit": deficit.detach(),
            "budget_deficit_rate": float((deficit > 0).float().mean().detach().cpu()),
        }
        return logits, mask, info

