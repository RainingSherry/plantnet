from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class CMTLoss(nn.Module):
    """scMAE loss plus collaborative teacher/student feature target."""

    def __init__(self, masked_data_weight: float = 0.75, mask_loss_weight: float = 0.7, target_weight: float = 0.08, variance_weight: float = 0.0):
        super().__init__()
        self.masked_data_weight = masked_data_weight
        self.mask_loss_weight = mask_loss_weight
        self.target_weight = target_weight
        self.variance_weight = variance_weight

    def forward(self, outputs: dict[str, torch.Tensor], target_expr: torch.Tensor, mask: torch.Tensor, target_feature: torch.Tensor) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.mse_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        pred = F.normalize(outputs["target_pred"], dim=-1)
        tgt = F.normalize(target_feature.detach(), dim=-1)
        target_loss = F.smooth_l1_loss(pred, tgt)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        var_loss = F.relu(0.5 - std).mean()
        scmae = (1.0 - self.mask_loss_weight) * rec + self.mask_loss_weight * mask_loss
        total = scmae + self.target_weight * target_loss + self.variance_weight * var_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "recon_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "target_loss": float(target_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
