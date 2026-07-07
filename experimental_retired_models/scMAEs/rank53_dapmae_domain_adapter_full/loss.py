from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class DAPMAELoss(nn.Module):
    """scMAE loss plus DAP-MAE domain feature generator guidance."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        domain_ce_weight: float = 0.08,
        domain_feature_weight: float = 0.08,
        contrast_weight: float = 0.02,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.domain_ce_weight = float(domain_ce_weight)
        self.domain_feature_weight = float(domain_feature_weight)
        self.contrast_weight = float(contrast_weight)
        self.variance_weight = float(variance_weight)

    def forward(self, outputs: dict[str, torch.Tensor], target_expr: torch.Tensor, mask: torch.Tensor, domain_id: torch.Tensor, domain_target: torch.Tensor) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        domain_ce = F.cross_entropy(outputs["domain_logits"], domain_id.long())
        pred = F.normalize(outputs["domain_feature"], dim=-1)
        target = F.normalize(domain_target.detach(), dim=-1)
        domain_feature = F.smooth_l1_loss(pred, target)
        sim = pred @ target.T
        same = (domain_id.view(-1, 1) == domain_id.view(1, -1)).float()
        same = same / same.sum(dim=1, keepdim=True).clamp_min(1.0)
        contrast = -(same * F.log_softmax(sim / 0.2, dim=1)).sum(dim=1).mean()
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.domain_ce_weight * domain_ce + self.domain_feature_weight * domain_feature + self.contrast_weight * contrast + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "domain_ce_loss": float(domain_ce.detach().cpu()),
            "domain_feature_loss": float(domain_feature.detach().cpu()),
            "contrast_loss": float(contrast.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
