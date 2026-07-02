from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class IJEPAScMAELoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.55,
        jepa_weight: float = 0.18,
        variance_weight: float = 0.02,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.jepa_weight = float(jepa_weight)
        self.variance_weight = float(variance_weight)

    def forward(self, out: dict[str, torch.Tensor], target_expr: torch.Tensor) -> tuple[torch.Tensor, dict[str, float]]:
        mask = out["gene_mask"].float()
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        reconstruction_loss = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask)
        pred = F.layer_norm(out["prediction"], out["prediction"].shape[-1:])
        target = F.layer_norm(out["target"], out["target"].shape[-1:])
        jepa_loss = F.smooth_l1_loss(pred, target)
        z = out["latent"] - out["latent"].mean(dim=0, keepdim=True)
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        variance_loss = torch.mean(F.relu(0.2 - std))
        scmae_loss = (1.0 - self.mask_weight) * reconstruction_loss + self.mask_weight * mask_loss
        total = scmae_loss + self.jepa_weight * jepa_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "reconstruction_loss": float(reconstruction_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "jepa_loss": float(jepa_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }
