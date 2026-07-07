from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class PlanetGrnAttentionLoss(nn.Module):
    """scMAE loss with Planet-style regulatory edge denoising."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        edge_weight: float = 0.025,
        sparsity_weight: float = 0.002,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.edge_weight = float(edge_weight)
        self.sparsity_weight = float(sparsity_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        edge_logits: torch.Tensor,
        clean_edge_target: torch.Tensor,
        edge_keep: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        edge_loss_raw = F.binary_cross_entropy_with_logits(edge_logits, clean_edge_target[None, :].expand_as(edge_logits), reduction="none")
        edge_loss = (edge_loss_raw * edge_keep).sum() / edge_keep.sum().clamp_min(1.0)
        edge_prob = torch.sigmoid(edge_logits)
        sparsity_loss = edge_prob.mean()
        z = outputs["embedding"]
        std = torch.sqrt(z.var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.edge_weight * edge_loss + self.sparsity_weight * sparsity_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "sparsity_loss": float(sparsity_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "edge_confidence_mean": float(edge_prob.detach().mean().cpu()),
            "edge_survival": float(edge_keep.detach().mean().cpu()),
        }
