from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class FlowMatchingLatentLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        flow_weight: float = 0.08,
        denoise_weight: float = 0.04,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.flow_weight = float(flow_weight)
        self.denoise_weight = float(denoise_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        flow_pred: torch.Tensor,
        flow_target: torch.Tensor,
        z1_hat: torch.Tensor,
        z_clean: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        flow_loss = F.smooth_l1_loss(flow_pred, flow_target)
        denoise_loss = F.smooth_l1_loss(z1_hat, z_clean.detach())
        var_loss = self.variance_loss(out["latent"])
        total = scmae + self.flow_weight * flow_loss + self.denoise_weight * denoise_loss + self.variance_weight * var_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "flow_loss": float(flow_loss.detach().cpu()),
            "denoise_loss": float(denoise_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "flow_target_norm": float(flow_target.norm(dim=1).mean().detach().cpu()),
            "flow_pred_norm": float(flow_pred.norm(dim=1).mean().detach().cpu()),
        }
