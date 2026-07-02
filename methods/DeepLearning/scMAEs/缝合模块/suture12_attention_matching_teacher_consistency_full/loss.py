from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class AttentionMatchingTeacherLoss(nn.Module):
    """scMAE objective with weak clean/masked attention matching consistency."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        match_weight: float = 0.03,
        delta_weight: float = 0.004,
        gate_budget_weight: float = 0.02,
        target_gate: float = 0.08,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.match_weight = float(match_weight)
        self.delta_weight = float(delta_weight)
        self.gate_budget_weight = float(gate_budget_weight)
        self.target_gate = float(target_gate)
        self.rec = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        rec_loss = (self.rec(outputs["reconstruction"], target) * (1.0 + 3.0 * mask)).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        match_loss = 1.0 - F.cosine_similarity(outputs["masked_proj"], outputs["clean_proj"].detach(), dim=1, eps=1e-6).mean()
        delta_penalty = outputs["matched_delta"].pow(2).mean()
        gate_budget = (outputs["match_gate"].mean() - self.target_gate).pow(2)
        total = (
            self.reconstruction_weight * rec_loss
            + self.mask_weight * mask_loss
            + self.match_weight * match_loss
            + self.delta_weight * delta_penalty
            + self.gate_budget_weight * gate_budget
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "match_loss": float(match_loss.detach().cpu()),
            "delta_penalty": float(delta_penalty.detach().cpu()),
            "match_gate_mean": float(outputs["match_gate"].mean().detach().cpu()),
            "match_attention_mean": float(outputs["match_attention"].mean().detach().cpu()),
            "match_similarity_mean": float(outputs["match_similarity"].mean().detach().cpu()),
            "gate_budget": float(gate_budget.detach().cpu()),
        }
        return total, parts
