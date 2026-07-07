from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class StructuralMaskLoss(nn.Module):
    """scMAE reconstruction plus mask prediction with structural policy regularization."""

    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.2,
        protection_weight: float = 0.02,
        budget_weight: float = 0.75,
        entropy_weight: float = 0.005,
        target_mask_prob: float = 0.4,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.protection_weight = float(protection_weight)
        self.budget_weight = float(budget_weight)
        self.entropy_weight = float(entropy_weight)
        self.target_mask_prob = float(target_mask_prob)
        self.reconstruction = nn.SmoothL1Loss(reduction="none")

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        recon = outputs["reconstruction"]
        mask_logits = outputs["mask_logits"]
        mask_prob = outputs["mask_prob"]
        marker_risk = outputs["marker_risk"].detach()

        rec_element = self.reconstruction(recon, target)
        weighted_rec = rec_element * (1.0 + 3.0 * mask)
        rec_loss = weighted_rec.mean()
        mask_loss = F.binary_cross_entropy_with_logits(mask_logits, mask)

        # High-risk marker genes should be sampled less often for destructive zero masking.
        protection_loss = (mask_prob * marker_risk).mean()
        budget_loss = (mask_prob.mean() - self.target_mask_prob).pow(2)
        entropy = -(mask_prob.clamp(1e-5, 1.0 - 1e-5) * torch.log(mask_prob.clamp(1e-5, 1.0 - 1e-5))
                    + (1.0 - mask_prob).clamp(1e-5, 1.0) * torch.log((1.0 - mask_prob).clamp(1e-5, 1.0))).mean()
        entropy_loss = -entropy

        total = (
            self.reconstruction_weight * rec_loss
            + self.mask_weight * mask_loss
            + self.protection_weight * protection_loss
            + self.budget_weight * budget_loss
            + self.entropy_weight * entropy_loss
        )
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "protection": float(protection_loss.detach().cpu()),
            "budget": float(budget_loss.detach().cpu()),
            "mask_prob_mean": float(mask_prob.mean().detach().cpu()),
            "marker_risk_mean": float(marker_risk.mean().detach().cpu()),
        }
        return total, parts
