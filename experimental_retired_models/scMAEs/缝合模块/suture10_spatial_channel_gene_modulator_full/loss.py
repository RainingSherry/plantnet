from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class SpatialChannelGeneModulatorLoss(nn.Module):
    """scMAE objective with weak gate-budgeted latent modulation regularization."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        delta_weight: float = 0.004,
        gate_budget_weight: float = 0.02,
        target_cell_gate: float = 0.12,
        target_channel_gate: float = 0.50,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.delta_weight = float(delta_weight)
        self.gate_budget_weight = float(gate_budget_weight)
        self.target_cell_gate = float(target_cell_gate)
        self.target_channel_gate = float(target_channel_gate)
        self.rec = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        rec_loss = (self.rec(outputs["reconstruction"], target) * (1.0 + 3.0 * mask)).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        delta_penalty = outputs["channel_delta"].pow(2).mean() + outputs["cell_delta"].pow(2).mean()
        cell_budget = (outputs["cell_gate"].mean() - self.target_cell_gate).pow(2)
        channel_budget = (outputs["channel_gate"].mean() - self.target_channel_gate).pow(2)
        total = (
            self.reconstruction_weight * rec_loss
            + self.mask_weight * mask_loss
            + self.delta_weight * delta_penalty
            + self.gate_budget_weight * (cell_budget + channel_budget)
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "delta_penalty": float(delta_penalty.detach().cpu()),
            "cell_gate_mean": float(outputs["cell_gate"].mean().detach().cpu()),
            "channel_gate_mean": float(outputs["channel_gate"].mean().detach().cpu()),
            "gate_budget": float((cell_budget + channel_budget).detach().cpu()),
        }
        return total, parts
