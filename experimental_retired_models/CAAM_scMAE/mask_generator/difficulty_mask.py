from __future__ import annotations

import torch


class DifficultyMaskGenerator:
    def __init__(self, mask_ratio: float, eta: float = 0.5) -> None:
        self.mask_ratio = float(mask_ratio)
        self.eta = float(eta)

    def __call__(self, x: torch.Tensor, x_hat: torch.Tensor, eligibility: torch.Tensor):
        difficulty = (x - x_hat).abs().detach()
        difficulty = difficulty.masked_fill(~eligibility.bool(), 0.0)
        probs = difficulty / (difficulty.sum(dim=1, keepdim=True) + 1.0e-8)
        uniform = eligibility.float() / (eligibility.float().sum(dim=1, keepdim=True) + 1.0e-8)
        probs = (1.0 - self.eta) * uniform + self.eta * probs
        g = x.shape[1]
        k = int(round(self.mask_ratio * g))
        scores = probs + 1.0e-6 * torch.rand_like(probs)
        mask = torch.zeros_like(x)
        for i in range(x.shape[0]):
            ki = min(k, int(eligibility[i].sum().item()))
            if ki > 0:
                mask[i, torch.topk(scores[i], k=ki).indices] = 1.0
        return probs.log().clamp_min(-30.0), mask, {"mask_type": "difficulty", "difficulty_detached": True}

