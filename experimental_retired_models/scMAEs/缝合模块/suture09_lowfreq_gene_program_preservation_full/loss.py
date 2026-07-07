from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class LowFreqGeneProgramLoss(nn.Module):
    """scMAE reconstruction objective plus weak low-frequency program preservation."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        program_weight: float = 0.05,
        latent_align_weight: float = 0.01,
        gate_budget_weight: float = 0.02,
        target_gate: float = 0.08,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.program_weight = float(program_weight)
        self.latent_align_weight = float(latent_align_weight)
        self.gate_budget_weight = float(gate_budget_weight)
        self.target_gate = float(target_gate)
        self.rec = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor, lowfreq_program: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        rec_loss = (self.rec(outputs["reconstruction"], target) * (1.0 + 3.0 * mask)).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        program_loss = F.smooth_l1_loss(outputs["program_pred"], lowfreq_program)
        latent_align = F.smooth_l1_loss(outputs["base_latent"], outputs["program_z"].detach())
        gate_budget = (outputs["program_gate"].mean() - self.target_gate).pow(2)
        total = (
            self.reconstruction_weight * rec_loss
            + self.mask_weight * mask_loss
            + self.program_weight * program_loss
            + self.latent_align_weight * latent_align
            + self.gate_budget_weight * gate_budget
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "program_loss": float(program_loss.detach().cpu()),
            "latent_align": float(latent_align.detach().cpu()),
            "program_gate_mean": float(outputs["program_gate"].mean().detach().cpu()),
            "gate_budget": float(gate_budget.detach().cpu()),
        }
        return total, parts
