from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class AnatomyPromptLoss(nn.Module):
    """scMAE loss with replaced-expression and gene rank-token objectives."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.6,
        replaced_weight: float = 0.15,
        token_weight: float = 0.05,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.replaced_weight = float(replaced_weight)
        self.token_weight = float(token_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        replaced_label: torch.Tensor,
        token_target: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        replaced_loss = F.binary_cross_entropy_with_logits(outputs["replaced_logits"], replaced_label)
        token_loss_raw = F.cross_entropy(
            outputs["token_logits"].reshape(-1, outputs["token_logits"].shape[-1]),
            token_target.reshape(-1).long(),
            reduction="none",
        ).view_as(mask)
        token_loss = ((0.5 + mask) * token_loss_raw).mean()
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.replaced_weight * replaced_loss + self.token_weight * token_loss + self.variance_weight * variance_loss
        parts = {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "replaced_loss": float(replaced_loss.detach().cpu()),
            "token_loss": float(token_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }
        return total, parts, rec_raw.detach().mean(dim=0)
