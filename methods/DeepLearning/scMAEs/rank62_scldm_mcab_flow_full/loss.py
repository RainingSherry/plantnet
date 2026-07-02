from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class ScLdmMcabFlowLoss(nn.Module):
    """scMAE loss with linear-interpolant latent flow matching."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        flow_weight: float = 0.035,
        norm_weight: float = 0.005,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.flow_weight = float(flow_weight)
        self.norm_weight = float(norm_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        flow_pred: torch.Tensor,
        flow_target: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        flow_loss = F.mse_loss(flow_pred, flow_target)
        z = outputs["embedding"]
        norm_loss = torch.square(z.norm(dim=1).mean() - 1.0)
        std = torch.sqrt(z.var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.flow_weight * flow_loss + self.norm_weight * norm_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "flow_loss": float(flow_loss.detach().cpu()),
            "norm_loss": float(norm_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }
