from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class QDiffusionLoss(nn.Module):
    """scMAE loss with q-kernel latent denoising and context consistency."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        denoise_weight: float = 0.08,
        context_weight: float = 0.015,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.denoise_weight = float(denoise_weight)
        self.context_weight = float(context_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        clean_z: torch.Tensor,
        denoised_z: torch.Tensor,
        q_context: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        denoise_loss = F.smooth_l1_loss(denoised_z, clean_z.detach())
        context_loss = 1.0 - F.cosine_similarity(denoised_z, q_context.detach(), dim=1).mean()
        z = outputs["embedding"]
        std = torch.sqrt(z.var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.denoise_weight * denoise_loss + self.context_weight * context_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "denoise_loss": float(denoise_loss.detach().cpu()),
            "context_loss": float(context_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "latent_noise_mse": float(F.mse_loss(denoised_z.detach(), clean_z.detach()).cpu()),
        }
