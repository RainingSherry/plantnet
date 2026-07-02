from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class AsymmetricMAELoss(nn.Module):
    """Original scMAE-style weighted MSE + mask BCE for asymmetric decoder."""

    def __init__(self, masked_data_weight: float = 0.75, mask_loss_weight: float = 0.7, variance_weight: float = 0.0):
        super().__init__()
        self.masked_data_weight = masked_data_weight
        self.mask_loss_weight = mask_loss_weight
        self.variance_weight = variance_weight

    def forward(self, outputs: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.mse_loss(outputs["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        var_loss = F.relu(0.5 - std).mean()
        total = (1.0 - self.mask_loss_weight) * rec + self.mask_loss_weight * mask_loss + self.variance_weight * var_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "recon_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
        }
