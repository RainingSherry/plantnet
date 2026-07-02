from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class ConditionalMultiMaskLoss(nn.Module):
    """scMAE loss with conditional latent-mask prediction."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        latent_mask_weight: float = 0.05,
        condition_weight: float = 0.02,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.latent_mask_weight = float(latent_mask_weight)
        self.condition_weight = float(condition_weight)
        self.variance_weight = float(variance_weight)

    def forward(self, outputs: dict[str, torch.Tensor], target_expr: torch.Tensor, input_mask: torch.Tensor, clean_latent_target: torch.Tensor) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = input_mask * self.masked_data_weight + (1.0 - input_mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], input_mask)
        latent_loss = F.binary_cross_entropy_with_logits(outputs["latent_mask_logits"], outputs["latent_mask"])
        condition_loss = F.smooth_l1_loss(outputs["condition_pred"], clean_latent_target.detach())
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.latent_mask_weight * latent_loss + self.condition_weight * condition_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "latent_mask_loss": float(latent_loss.detach().cpu()),
            "condition_loss": float(condition_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
