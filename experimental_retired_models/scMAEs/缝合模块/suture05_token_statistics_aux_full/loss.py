from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class TokenStatisticsAuxLoss(nn.Module):
    """scMAE loss plus a lightweight cell-statistics auxiliary target."""

    def __init__(self, reconstruction_weight: float = 1.0, mask_weight: float = 0.2, stat_weight: float = 0.04, entropy_weight: float = 0.002):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.stat_weight = float(stat_weight)
        self.entropy_weight = float(entropy_weight)
        self.rec = nn.SmoothL1Loss(reduction="none")

    @staticmethod
    def stat_target(log_expr: torch.Tensor) -> torch.Tensor:
        mean = log_expr.mean(dim=1)
        std = log_expr.std(dim=1, unbiased=False)
        dropout = (log_expr <= 1e-6).float().mean(dim=1)
        return torch.stack([mean, std, dropout], dim=1)

    def forward(self, outputs: dict, target: torch.Tensor) -> tuple[torch.Tensor, dict]:
        mask = outputs["mask"]
        rec_loss = (self.rec(outputs["reconstruction"], target) * (1.0 + 3.0 * mask)).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        stat_loss = F.smooth_l1_loss(outputs["stat_pred"], self.stat_target(target))
        imp = outputs["token_importance"].clamp(1e-8, 1.0)
        entropy = -(imp * torch.log(imp)).sum(dim=1).mean() / torch.log(torch.tensor(float(imp.shape[1]), device=imp.device))
        total = self.reconstruction_weight * rec_loss + self.mask_weight * mask_loss + self.stat_weight * stat_loss - self.entropy_weight * entropy
        parts = {
            "total": float(total.detach().cpu()),
            "reconstruction": float(rec_loss.detach().cpu()),
            "mask_bce": float(mask_loss.detach().cpu()),
            "stat": float(stat_loss.detach().cpu()),
            "token_entropy": float(entropy.detach().cpu()),
        }
        return total, parts
