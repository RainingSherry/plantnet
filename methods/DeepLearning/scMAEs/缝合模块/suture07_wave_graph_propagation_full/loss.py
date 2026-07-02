from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class GraphWaveLoss(nn.Module):
    """scMAE objective with weak graph-wave latent alignment."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        wave_align_weight: float = 0.03,
        gate_budget_weight: float = 0.02,
        target_gate: float = 0.10,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.wave_align_weight = float(wave_align_weight)
        self.gate_budget_weight = float(gate_budget_weight)
        self.target_gate = float(target_gate)
        self.rec = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        rec_loss = (self.rec(outputs["reconstruction"], target) * (1.0 + 3.0 * mask)).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        wave_align = F.smooth_l1_loss(outputs["base_latent"], outputs["wave_z"].detach())
        gate = outputs["wave_gate"]
        gate_budget = (gate.mean() - self.target_gate).pow(2)
        total = (
            self.reconstruction_weight * rec_loss
            + self.mask_weight * mask_loss
            + self.wave_align_weight * wave_align
            + self.gate_budget_weight * gate_budget
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "wave_align": float(wave_align.detach().cpu()),
            "wave_gate_mean": float(gate.mean().detach().cpu()),
            "gate_budget": float(gate_budget.detach().cpu()),
        }
        return total, parts
