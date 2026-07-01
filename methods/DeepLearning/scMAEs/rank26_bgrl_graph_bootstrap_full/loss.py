from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class BGRLGraphBootstrapLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        bootstrap_weight: float = 0.08,
        edge_weight: float = 0.04,
        neighbor_weight: float = 0.03,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.bootstrap_weight = float(bootstrap_weight)
        self.edge_weight = float(edge_weight)
        self.neighbor_weight = float(neighbor_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def negative_cosine(pred: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
        pred = F.normalize(pred, dim=1)
        target = F.normalize(target.detach(), dim=1)
        return 2.0 - 2.0 * (pred * target).sum(dim=1).mean()

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        bootstrap_pairs: list[tuple[torch.Tensor, torch.Tensor]],
        edge_logits: tuple[torch.Tensor, torch.Tensor],
        neighbor_latent: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        bootstrap_loss = torch.stack([self.negative_cosine(pred, target) for pred, target in bootstrap_pairs]).mean()
        pos_logit, neg_logit = edge_logits
        edge_loss = 0.5 * (
            F.binary_cross_entropy_with_logits(pos_logit, torch.ones_like(pos_logit))
            + F.binary_cross_entropy_with_logits(neg_logit, torch.zeros_like(neg_logit))
        )
        neighbor_loss = F.smooth_l1_loss(out["latent"], neighbor_latent.detach())
        var_loss = self.variance_loss(out["latent"])
        total = scmae + self.bootstrap_weight * bootstrap_loss + self.edge_weight * edge_loss + self.neighbor_weight * neighbor_loss + self.variance_weight * var_loss
        edge_acc = 0.5 * ((pos_logit.sigmoid() > 0.5).float().mean() + (neg_logit.sigmoid() < 0.5).float().mean())
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "bootstrap_loss": float(bootstrap_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "neighbor_loss": float(neighbor_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "edge_confidence": float(pos_logit.sigmoid().mean().detach().cpu()),
            "edge_negative_confidence": float(neg_logit.sigmoid().mean().detach().cpu()),
            "edge_proxy_accuracy": float(edge_acc.detach().cpu()),
        }
