from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class GraphormerStructuralLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        structure_weight: float = 0.04,
        edge_weight: float = 0.04,
        anchor_weight: float = 0.05,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.structure_weight = float(structure_weight)
        self.edge_weight = float(edge_weight)
        self.anchor_weight = float(anchor_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    @staticmethod
    def cosine_loss(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
        a = F.normalize(a, dim=1)
        b = F.normalize(b.detach(), dim=1)
        return 1.0 - (a * b).sum(dim=1).mean()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        clean_out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        anchor: torch.Tensor,
        edge_logits: tuple[torch.Tensor, torch.Tensor],
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        struct_loss = self.cosine_loss(out["latent"], clean_out["latent"])
        attn_loss = F.smooth_l1_loss(out["attention"], clean_out["attention"].detach())
        anchor_loss = F.smooth_l1_loss(out["anchor_pred"], anchor)
        pos_logit, neg_logit = edge_logits
        edge_loss = 0.5 * (
            F.binary_cross_entropy_with_logits(pos_logit, torch.ones_like(pos_logit))
            + F.binary_cross_entropy_with_logits(neg_logit, torch.zeros_like(neg_logit))
        )
        var_loss = self.variance_loss(out["latent"])
        total = scmae + self.structure_weight * (struct_loss + 0.25 * attn_loss) + self.edge_weight * edge_loss + self.anchor_weight * anchor_loss + self.variance_weight * var_loss
        edge_acc = 0.5 * ((pos_logit.sigmoid() > 0.5).float().mean() + (neg_logit.sigmoid() < 0.5).float().mean())
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "structure_loss": float(struct_loss.detach().cpu()),
            "attention_stability_loss": float(attn_loss.detach().cpu()),
            "anchor_loss": float(anchor_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "edge_confidence": float(pos_logit.sigmoid().mean().detach().cpu()),
            "edge_negative_confidence": float(neg_logit.sigmoid().mean().detach().cpu()),
            "edge_proxy_accuracy": float(edge_acc.detach().cpu()),
        }
