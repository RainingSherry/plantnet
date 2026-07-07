from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class JOAOMaskPolicyLoss(nn.Module):
    def __init__(self, masked_data_weight: float = 0.75, mask_weight: float = 0.7, consistency_weight: float = 0.05):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.consistency_weight = float(consistency_weight)

    def forward(
        self,
        out1: dict[str, torch.Tensor],
        out2: dict[str, torch.Tensor],
        target: torch.Tensor,
        mask1: torch.Tensor,
        mask2: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec1 = self._scmae_loss(out1, target, mask1)
        rec2 = self._scmae_loss(out2, target, mask2)
        z1 = F.normalize(out1["projection"], dim=1)
        z2 = F.normalize(out2["projection"], dim=1)
        consistency = F.mse_loss(z1, z2.detach()) + F.mse_loss(z2, z1.detach())
        loss = 0.5 * (rec1 + rec2) + self.consistency_weight * consistency
        parts = {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float((0.5 * (rec1 + rec2)).detach().cpu()),
            "consistency_loss": float(consistency.detach().cpu()),
        }
        return loss, parts

    def _scmae_loss(self, outputs: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = torch.mul(weights, F.mse_loss(outputs["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask.float())
        return (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

