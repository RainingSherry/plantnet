from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class SoftGraphLoss(nn.Module):
    """scMAE loss with soft graph edge confidence and residual consistency."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        graph_weight: float = 0.025,
        edge_weight: float = 0.02,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.graph_weight = float(graph_weight)
        self.edge_weight = float(edge_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        model,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        neighbor_z: torch.Tensor,
        edge_target: torch.Tensor,
        edge_keep: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.mse_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        z = outputs["embedding"]
        edge_logits = model.edge_logits(z, neighbor_z.detach())
        edge_loss_raw = F.binary_cross_entropy_with_logits(edge_logits, edge_target, reduction="none")
        edge_loss = (edge_loss_raw * edge_keep).sum() / edge_keep.sum().clamp_min(1.0)
        z_norm = F.normalize(z, dim=1)
        neigh_norm = F.normalize(neighbor_z.detach(), dim=2)
        cosine_dist = 1.0 - (z_norm[:, None, :] * neigh_norm).sum(dim=2)
        graph_loss = (cosine_dist * edge_target * edge_keep).sum() / edge_keep.sum().clamp_min(1.0)
        std = torch.sqrt(z.var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.graph_weight * graph_loss + self.edge_weight * edge_loss + self.variance_weight * variance_loss
        edge_conf = torch.sigmoid(edge_logits).detach()
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "graph_loss": float(graph_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "edge_confidence_mean": float(edge_conf.mean().cpu()),
            "edge_survival": float(edge_keep.mean().detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
