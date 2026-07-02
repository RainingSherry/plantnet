from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


def _masked_mean(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    return (values * mask).sum() / mask.sum().clamp_min(1.0)


class MaskFeatModuleLoss(nn.Module):
    """scMAE reconstruction + mask prediction + module semantic feature loss."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.3,
        semantic_weight: float = 0.5,
        huber_beta: float = 1.0,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.semantic_weight = float(semantic_weight)
        self.huber_beta = float(huber_beta)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        expression_target: torch.Tensor,
        semantic_target: torch.Tensor,
        gene_mask: torch.Tensor,
        module_mask: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec_raw = F.smooth_l1_loss(
            outputs["reconstruction"],
            expression_target,
            reduction="none",
            beta=self.huber_beta,
        )
        reconstruction_loss = _masked_mean(rec_raw, gene_mask.float())
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], gene_mask.float())
        sem_raw = F.smooth_l1_loss(
            outputs["semantic"],
            semantic_target,
            reduction="none",
            beta=self.huber_beta,
        )
        semantic_loss = _masked_mean(sem_raw, module_mask.float().unsqueeze(-1))
        loss = (
            self.reconstruction_weight * reconstruction_loss
            + self.mask_weight * mask_loss
            + self.semantic_weight * semantic_loss
        )
        parts = {
            "loss": float(loss.detach().cpu()),
            "reconstruction_loss": float(reconstruction_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "semantic_loss": float(semantic_loss.detach().cpu()),
        }
        return loss, parts

