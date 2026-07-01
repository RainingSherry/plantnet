from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ScAGCGraphLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        edge_weight: float = 0.05,
        graph_weight: float = 0.05,
        dropedge_weight: float = 0.02,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.edge_weight = float(edge_weight)
        self.graph_weight = float(graph_weight)
        self.dropedge_weight = float(dropedge_weight)

    def forward(
        self,
        out: dict[str, torch.Tensor],
        neighbor_out: dict[str, torch.Tensor],
        dropped_out: dict[str, torch.Tensor],
        target: torch.Tensor,
        mask: torch.Tensor,
        pos_edge_logits: torch.Tensor,
        neg_edge_logits: torch.Tensor,
        edge_conf: torch.Tensor,
        mix_gate: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        pos_target = torch.ones_like(pos_edge_logits)
        neg_target = torch.zeros_like(neg_edge_logits)
        pos_loss = F.binary_cross_entropy_with_logits(pos_edge_logits, pos_target, reduction="none")
        neg_loss = F.binary_cross_entropy_with_logits(neg_edge_logits, neg_target, reduction="none")
        edge_loss = ((pos_loss * edge_conf.detach()).mean() + neg_loss.mean()) * 0.5

        z = F.normalize(out["latent"], dim=1)
        zn = F.normalize(neighbor_out["latent"], dim=1)
        graph_cons = ((z - zn).pow(2).sum(dim=1) * edge_conf.detach() * (1.0 - 0.5 * mix_gate.detach())).mean()

        zd = F.normalize(dropped_out["latent"], dim=1)
        dropedge = (z - zd).pow(2).sum(dim=1).mean()
        total = scmae + self.edge_weight * edge_loss + self.graph_weight * graph_cons + self.dropedge_weight * dropedge
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "graph_loss": float(graph_cons.detach().cpu()),
            "dropedge_loss": float(dropedge.detach().cpu()),
            "mixed_cell_fraction": float(mix_gate.detach().mean().cpu()),
        }

