from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class SAINTRowColLoss(nn.Module):
    """scMAE reconstruction plus SAINT-style row consistency and gate sparsity."""

    def __init__(self, recon_weight: float = 1.0, mask_weight: float = 0.05, contrast_weight: float = 0.05, gate_weight: float = 0.005, variance_weight: float = 0.01):
        super().__init__()
        self.recon_weight = recon_weight
        self.mask_weight = mask_weight
        self.contrast_weight = contrast_weight
        self.gate_weight = gate_weight
        self.variance_weight = variance_weight

    def forward(self, outputs: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor, view_projection: torch.Tensor | None = None) -> tuple[torch.Tensor, dict[str, float]]:
        recon_raw = F.smooth_l1_loss(outputs["reconstruction"], target, reduction="none")
        recon_loss = (recon_raw * (1.0 + 2.0 * mask)).sum() / (target.numel() + 2.0 * mask.sum().clamp_min(1.0))
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        if view_projection is None:
            contrast_loss = outputs["embedding"].new_tensor(0.0)
        else:
            p = F.normalize(outputs["projection"], dim=-1)
            q = F.normalize(view_projection.detach(), dim=-1)
            contrast_loss = 2.0 - 2.0 * (p * q).sum(dim=-1).mean()
        gate_loss = outputs["feature_gate"].mean()
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        total = self.recon_weight * recon_loss + self.mask_weight * mask_loss + self.contrast_weight * contrast_loss + self.gate_weight * gate_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "recon_loss": float(recon_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "contrast_loss": float(contrast_loss.detach().cpu()),
            "gate_loss": float(gate_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }
