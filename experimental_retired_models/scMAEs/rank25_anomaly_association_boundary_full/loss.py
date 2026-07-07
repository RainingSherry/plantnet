from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class AnomalyAssociationBoundaryLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        assoc_weight: float = 0.04,
        latent_consistency_weight: float = 0.03,
        boundary_weight: float = 0.02,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.assoc_weight = float(assoc_weight)
        self.latent_consistency_weight = float(latent_consistency_weight)
        self.boundary_weight = float(boundary_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    @staticmethod
    def _risk_target(clean_assoc: torch.Tensor) -> torch.Tensor:
        score = clean_assoc.mean(dim=1)
        center = score.median()
        spread = torch.quantile((score - center).abs(), 0.75).clamp_min(1e-4)
        return torch.sigmoid((score - center) / spread).detach()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        clean_out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        assoc_loss = F.smooth_l1_loss(out["assoc_profile"], clean_out["assoc_profile"].detach())
        latent_loss = F.smooth_l1_loss(out["latent"], clean_out["latent"].detach())
        boundary_target = self._risk_target(clean_out["assoc_profile"])
        boundary_loss = F.binary_cross_entropy_with_logits(out["boundary_logit"], boundary_target)
        var_loss = self.variance_loss(out["latent"])
        total = (
            scmae
            + self.assoc_weight * assoc_loss
            + self.latent_consistency_weight * latent_loss
            + self.boundary_weight * boundary_loss
            + self.variance_weight * var_loss
        )
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "association_stability_loss": float(assoc_loss.detach().cpu()),
            "latent_consistency_loss": float(latent_loss.detach().cpu()),
            "boundary_loss": float(boundary_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "association_risk_mean": float(out["association_risk"].mean().detach().cpu()),
            "association_risk_std": float(out["association_risk"].std(unbiased=False).detach().cpu()),
        }
