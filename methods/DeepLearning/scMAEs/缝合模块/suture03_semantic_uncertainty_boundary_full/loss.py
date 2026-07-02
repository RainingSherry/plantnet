from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class SemanticUncertaintyLoss(nn.Module):
    """scMAE loss with SID-style core/boundary/rare-risk gate supervision."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        gate_weight: float = 0.08,
        boundary_guard_weight: float = 0.03,
        entropy_weight: float = 0.005,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.gate_weight = float(gate_weight)
        self.boundary_guard_weight = float(boundary_guard_weight)
        self.entropy_weight = float(entropy_weight)
        self.reconstruction = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor, gate_target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        recon = outputs["reconstruction"]
        rec = (self.reconstruction(recon, target) * (1.0 + 3.0 * mask)).mean()
        mask_bce = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        gates = outputs["gates"].clamp(1e-6, 1.0)
        gate_loss = -(gate_target * torch.log(gates)).sum(dim=1).mean()
        boundary_strength = gates[:, 1]
        guard_loss = (boundary_strength * outputs["delta"].pow(2).mean(dim=1)).mean()
        entropy = -(gates * torch.log(gates)).sum(dim=1).mean()
        total = (
            self.reconstruction_weight * rec
            + self.mask_weight * mask_bce
            + self.gate_weight * gate_loss
            + self.boundary_guard_weight * guard_loss
            - self.entropy_weight * entropy
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec.detach().cpu()),
            "mask_bce": float(mask_bce.detach().cpu()),
            "gate": float(gate_loss.detach().cpu()),
            "guard": float(guard_loss.detach().cpu()),
            "gate_core": float(gates[:, 0].mean().detach().cpu()),
            "gate_boundary": float(gates[:, 1].mean().detach().cpu()),
            "gate_rare": float(gates[:, 2].mean().detach().cpu()),
        }
        return total, parts
