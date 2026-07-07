from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class SelfGuidedMaskLoss(nn.Module):
    """scMAE reconstruction plus mask detection and gene-rank token prediction."""

    def __init__(self, recon_weight: float = 1.0, mask_weight: float = 0.08, token_weight: float = 0.08, difficulty_weight: float = 0.02, variance_weight: float = 0.01):
        super().__init__()
        self.recon_weight = recon_weight
        self.mask_weight = mask_weight
        self.token_weight = token_weight
        self.difficulty_weight = difficulty_weight
        self.variance_weight = variance_weight

    def forward(self, outputs: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor, rank_target: torch.Tensor) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        recon_raw = F.smooth_l1_loss(outputs["reconstruction"], target, reduction="none")
        recon_loss = (recon_raw * (1.0 + 2.0 * mask)).sum() / (target.numel() + 2.0 * mask.sum().clamp_min(1.0))
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        token_loss_all = F.cross_entropy(outputs["rank_logits"].reshape(-1, outputs["rank_logits"].shape[-1]), rank_target.reshape(-1), reduction="none").view_as(mask)
        token_loss = (token_loss_all * (1.0 + mask)).mean()
        gene_error = recon_raw.detach().mean(dim=0)
        diff_target = (gene_error / gene_error.mean().clamp_min(1e-6)).clamp(0.0, 3.0) / 3.0
        difficulty_loss = F.smooth_l1_loss(outputs["difficulty"].mean(dim=0), diff_target)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        total = self.recon_weight * recon_loss + self.mask_weight * mask_loss + self.token_weight * token_loss + self.difficulty_weight * difficulty_loss + self.variance_weight * variance_loss
        with torch.no_grad():
            token_acc = (outputs["rank_logits"].argmax(dim=-1) == rank_target).float()
            token_acc_masked = (token_acc * mask).sum() / mask.sum().clamp_min(1.0)
        return total, {
            "loss": float(total.detach().cpu()),
            "recon_loss": float(recon_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "token_loss": float(token_loss.detach().cpu()),
            "difficulty_loss": float(difficulty_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "masked_token_accuracy": float(token_acc_masked.detach().cpu()),
        }, gene_error
