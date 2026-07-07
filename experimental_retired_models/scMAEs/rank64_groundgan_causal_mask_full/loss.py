from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class GroundGANCausalMaskLoss(nn.Module):
    """scMAE loss with masked causal dependency prediction."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        causal_weight: float = 0.05,
        regulator_dropout_weight: float = 0.005,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.causal_weight = float(causal_weight)
        self.regulator_dropout_weight = float(regulator_dropout_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        causal_pred: torch.Tensor,
        causal_target: torch.Tensor,
        causal_pred_dropout: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        causal_loss = F.smooth_l1_loss(causal_pred, causal_target)
        dropout_loss = torch.tensor(0.0, device=causal_pred.device)
        if causal_pred_dropout is not None:
            dropout_loss = F.smooth_l1_loss(causal_pred_dropout, causal_pred.detach())
        z = outputs["embedding"]
        std = torch.sqrt(z.var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.causal_weight * causal_loss + self.regulator_dropout_weight * dropout_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "causal_loss": float(causal_loss.detach().cpu()),
            "regulator_dropout_loss": float(dropout_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }
