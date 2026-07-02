from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class ScDiffusionLoss(nn.Module):
    """scMAE loss with conditional DDPM epsilon prediction."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        diffusion_weight: float = 0.06,
        x0_weight: float = 0.01,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.diffusion_weight = float(diffusion_weight)
        self.x0_weight = float(x0_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        clean_z: torch.Tensor,
        pred_noise: torch.Tensor,
        true_noise: torch.Tensor,
        pred_x0: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        diffusion_loss = F.mse_loss(pred_noise, true_noise)
        x0_loss = F.smooth_l1_loss(pred_x0, clean_z.detach())
        z = outputs["embedding"]
        std = torch.sqrt(z.var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.diffusion_weight * diffusion_loss + self.x0_weight * x0_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "diffusion_loss": float(diffusion_loss.detach().cpu()),
            "x0_loss": float(x0_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }
