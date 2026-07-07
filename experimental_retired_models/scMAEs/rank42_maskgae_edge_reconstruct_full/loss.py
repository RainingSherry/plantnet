from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class MaskGAEEdgeLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        edge_weight: float = 0.08,
        degree_weight: float = 0.01,
        neighbor_weight: float = 0.02,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.edge_weight = float(edge_weight)
        self.degree_weight = float(degree_weight)
        self.neighbor_weight = float(neighbor_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        edge_logits: tuple[torch.Tensor, torch.Tensor],
        degree_target: torch.Tensor,
        neighbor_latent: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        pos_logit, neg_logit = edge_logits
        edge_loss = 0.5 * (
            F.binary_cross_entropy_with_logits(pos_logit, torch.ones_like(pos_logit))
            + F.binary_cross_entropy_with_logits(neg_logit, torch.zeros_like(neg_logit))
        )
        degree_loss = F.smooth_l1_loss(out["degree_pred"], degree_target)
        neighbor_loss = F.smooth_l1_loss(out["latent"], neighbor_latent.detach())
        var_loss = self.variance_loss(out["latent"])
        total = scmae + self.edge_weight * edge_loss + self.degree_weight * degree_loss + self.neighbor_weight * neighbor_loss + self.variance_weight * var_loss
        edge_acc = 0.5 * ((pos_logit.sigmoid() > 0.5).float().mean() + (neg_logit.sigmoid() < 0.5).float().mean())
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "degree_loss": float(degree_loss.detach().cpu()),
            "neighbor_loss": float(neighbor_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "edge_confidence": float(pos_logit.sigmoid().mean().detach().cpu()),
            "edge_negative_confidence": float(neg_logit.sigmoid().mean().detach().cpu()),
            "edge_proxy_accuracy": float(edge_acc.detach().cpu()),
        }
