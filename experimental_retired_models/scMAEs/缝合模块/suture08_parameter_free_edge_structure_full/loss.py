from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ParameterFreeEdgeStructureLoss(nn.Module):
    """scMAE objective with reliability-weighted structure consistency."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        structure_align_weight: float = 0.025,
        delta_penalty_weight: float = 0.005,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.structure_align_weight = float(structure_align_weight)
        self.delta_penalty_weight = float(delta_penalty_weight)
        self.rec = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        reliability = outputs["reliability"].clamp(0.0, 1.0)
        boundary_weight = 0.75 + 0.50 * (1.0 - reliability)
        rec = self.rec(outputs["reconstruction"], target)
        rec_loss = (rec * (1.0 + 3.0 * mask) * boundary_weight).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        align_each = F.smooth_l1_loss(outputs["base_latent"], outputs["context_z"].detach(), reduction="none").mean(dim=1, keepdim=True)
        structure_align = (align_each * reliability).mean()
        delta_penalty = (outputs["adapter_delta"].pow(2).mean(dim=1, keepdim=True) * (1.0 - reliability)).mean()
        total = (
            self.reconstruction_weight * rec_loss
            + self.mask_weight * mask_loss
            + self.structure_align_weight * structure_align
            + self.delta_penalty_weight * delta_penalty
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "structure_align": float(structure_align.detach().cpu()),
            "delta_penalty": float(delta_penalty.detach().cpu()),
            "reliability_mean": float(reliability.mean().detach().cpu()),
        }
        return total, parts
