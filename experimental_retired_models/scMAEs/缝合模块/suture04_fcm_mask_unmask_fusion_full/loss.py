from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class FCMMaskUnmaskLoss(nn.Module):
    """scMAE loss with conservative masked/unmasked latent correction."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        align_weight: float = 0.05,
        gate_budget_weight: float = 0.04,
        target_gate: float = 0.12,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.align_weight = float(align_weight)
        self.gate_budget_weight = float(gate_budget_weight)
        self.target_gate = float(target_gate)
        self.rec = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        reconstruction = (self.rec(outputs["reconstruction"], target) * (1.0 + 3.0 * mask)).mean()
        mask_bce = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        align = F.smooth_l1_loss(outputs["masked_z"], outputs["clean_z"].detach())
        gate = outputs["correction_gate"]
        gate_budget = (gate.mean() - self.target_gate).pow(2)
        total = (
            self.reconstruction_weight * reconstruction
            + self.mask_weight * mask_bce
            + self.align_weight * align
            + self.gate_budget_weight * gate_budget
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(reconstruction.detach().cpu()),
            "mask_bce": float(mask_bce.detach().cpu()),
            "align": float(align.detach().cpu()),
            "gate_mean": float(gate.mean().detach().cpu()),
            "gate_budget": float(gate_budget.detach().cpu()),
        }
        return total, parts
