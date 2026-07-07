from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class GraphGPSLoss(nn.Module):
    """Loss for the GraphGPS local-global scMAE candidate."""

    def __init__(
        self,
        recon_weight: float = 1.0,
        mask_weight: float = 0.05,
        pe_weight: float = 0.05,
        edge_weight: float = 0.02,
        consistency_weight: float = 0.05,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.recon_weight = recon_weight
        self.mask_weight = mask_weight
        self.pe_weight = pe_weight
        self.edge_weight = edge_weight
        self.consistency_weight = consistency_weight
        self.variance_weight = variance_weight

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target: torch.Tensor,
        mask: torch.Tensor,
        pe_target: torch.Tensor,
        pos_edge_logits: torch.Tensor,
        neg_edge_logits: torch.Tensor,
        consistency_z: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        recon_raw = F.smooth_l1_loss(outputs["reconstruction"], target, reduction="none")
        mask_denom = mask.sum().clamp_min(1.0)
        recon_loss = (recon_raw * (1.0 + 2.0 * mask)).sum() / (target.numel() + 2.0 * mask_denom)
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        pe_loss = F.smooth_l1_loss(outputs["pe_reconstruction"], pe_target)

        pos_loss = F.binary_cross_entropy_with_logits(pos_edge_logits, torch.ones_like(pos_edge_logits))
        neg_loss = F.binary_cross_entropy_with_logits(neg_edge_logits, torch.zeros_like(neg_edge_logits))
        edge_loss = 0.5 * (pos_loss + neg_loss)

        if consistency_z is None:
            consistency_loss = outputs["embedding"].new_tensor(0.0)
        else:
            consistency_loss = F.smooth_l1_loss(outputs["embedding"], consistency_z.detach())

        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()

        total = (
            self.recon_weight * recon_loss
            + self.mask_weight * mask_loss
            + self.pe_weight * pe_loss
            + self.edge_weight * edge_loss
            + self.consistency_weight * consistency_loss
            + self.variance_weight * variance_loss
        )
        with torch.no_grad():
            edge_conf = torch.sigmoid(pos_edge_logits).mean()
            neg_conf = torch.sigmoid(neg_edge_logits).mean()
            edge_acc = 0.5 * ((torch.sigmoid(pos_edge_logits) > 0.5).float().mean() + (torch.sigmoid(neg_edge_logits) < 0.5).float().mean())
        return total, {
            "loss": float(total.detach().cpu()),
            "recon_loss": float(recon_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "pe_loss": float(pe_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "consistency_loss": float(consistency_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "edge_confidence": float(edge_conf.detach().cpu()),
            "edge_negative_confidence": float(neg_conf.detach().cpu()),
            "edge_proxy_accuracy": float(edge_acc.detach().cpu()),
        }
