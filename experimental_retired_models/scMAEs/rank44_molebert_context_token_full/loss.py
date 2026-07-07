from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class MoleBERTContextLoss(nn.Module):
    """scMAE loss plus Mole-BERT-style masked context-token modeling."""

    def __init__(self, recon_weight: float = 1.0, mask_weight: float = 0.05, token_weight: float = 0.12, edge_weight: float = 0.02, variance_weight: float = 0.01):
        super().__init__()
        self.recon_weight = recon_weight
        self.mask_weight = mask_weight
        self.token_weight = token_weight
        self.edge_weight = edge_weight
        self.variance_weight = variance_weight

    def forward(self, outputs: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor, token_target: torch.Tensor, pos_logits: torch.Tensor, neg_logits: torch.Tensor) -> tuple[torch.Tensor, dict[str, float]]:
        recon_raw = F.smooth_l1_loss(outputs["reconstruction"], target, reduction="none")
        recon_loss = (recon_raw * (1.0 + 2.0 * mask)).sum() / (target.numel() + 2.0 * mask.sum().clamp_min(1.0))
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        token_loss = F.cross_entropy(outputs["token_logits"], token_target)
        pos_loss = F.binary_cross_entropy_with_logits(pos_logits, torch.ones_like(pos_logits))
        neg_loss = F.binary_cross_entropy_with_logits(neg_logits, torch.zeros_like(neg_logits))
        edge_loss = 0.5 * (pos_loss + neg_loss)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        total = self.recon_weight * recon_loss + self.mask_weight * mask_loss + self.token_weight * token_loss + self.edge_weight * edge_loss + self.variance_weight * variance_loss
        with torch.no_grad():
            token_acc = (outputs["token_logits"].argmax(dim=-1) == token_target).float().mean()
            edge_acc = 0.5 * ((torch.sigmoid(pos_logits) > 0.5).float().mean() + (torch.sigmoid(neg_logits) < 0.5).float().mean())
        return total, {
            "loss": float(total.detach().cpu()),
            "recon_loss": float(recon_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "token_loss": float(token_loss.detach().cpu()),
            "edge_loss": float(edge_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "token_accuracy_proxy": float(token_acc.detach().cpu()),
            "edge_confidence": float(torch.sigmoid(pos_logits).mean().detach().cpu()),
            "edge_negative_confidence": float(torch.sigmoid(neg_logits).mean().detach().cpu()),
            "edge_proxy_accuracy": float(edge_acc.detach().cpu()),
        }
