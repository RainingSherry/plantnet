from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class PrototypeGraphLoss(nn.Module):
    """scMAE loss plus prototype assignment and separation objectives."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        consistency_weight: float = 0.04,
        entropy_weight: float = 0.02,
        separation_weight: float = 0.03,
        target_entropy: float = 0.45,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.consistency_weight = float(consistency_weight)
        self.entropy_weight = float(entropy_weight)
        self.separation_weight = float(separation_weight)
        self.target_entropy = float(target_entropy)
        self.rec = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        rec_loss = (self.rec(outputs["reconstruction"], target) * (1.0 + 3.0 * mask)).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        q = outputs["assignment"].clamp(1e-8, 1.0)
        q_clean = outputs["clean_assignment"].detach().clamp(1e-8, 1.0)
        consistency = F.kl_div(torch.log(q), q_clean, reduction="batchmean")
        entropy = -(q * torch.log(q)).sum(dim=1).mean() / torch.log(torch.tensor(float(q.shape[1]), device=q.device))
        entropy_budget = (entropy - self.target_entropy).pow(2)
        proto = F.normalize(outputs["prototypes"], dim=1)
        sim = torch.matmul(proto, proto.t())
        offdiag = sim - torch.eye(sim.shape[0], device=sim.device, dtype=sim.dtype)
        separation = F.relu(offdiag - 0.15).pow(2).mean()
        total = (
            self.reconstruction_weight * rec_loss
            + self.mask_weight * mask_loss
            + self.consistency_weight * consistency
            + self.entropy_weight * entropy_budget
            + self.separation_weight * separation
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "proto_consistency": float(consistency.detach().cpu()),
            "proto_entropy": float(entropy.detach().cpu()),
            "proto_entropy_budget": float(entropy_budget.detach().cpu()),
            "proto_separation": float(separation.detach().cpu()),
            "proto_gate_mean": float(outputs["proto_gate"].mean().detach().cpu()),
        }
        return total, parts
