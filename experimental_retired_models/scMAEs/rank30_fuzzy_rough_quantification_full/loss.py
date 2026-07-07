from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class FuzzyRoughQuantificationLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        lower_weight: float = 0.07,
        width_weight: float = 0.03,
        upper_entropy_weight: float = 0.02,
        balance_weight: float = 0.02,
        separation_weight: float = 0.01,
        anchor_weight: float = 0.04,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.lower_weight = float(lower_weight)
        self.width_weight = float(width_weight)
        self.upper_entropy_weight = float(upper_entropy_weight)
        self.balance_weight = float(balance_weight)
        self.separation_weight = float(separation_weight)
        self.anchor_weight = float(anchor_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    @staticmethod
    def center_separation_loss(centers: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(centers, centers)
        eye = torch.eye(centers.shape[0], dtype=torch.bool, device=centers.device)
        dist = dist.masked_fill(eye, 10.0)
        return F.relu(1.0 - dist).pow(2).mean()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        anchor: torch.Tensor,
        centers: torch.Tensor,
        rough_active: bool,
        confidence_threshold: float,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        q = out["membership"]
        p = out["target_distribution"]
        lower = out["rough_lower"]
        upper = out["rough_upper"]
        confidence = out["membership_confidence"].detach()
        core = (confidence >= float(confidence_threshold)).float()
        if rough_active and float(core.sum().detach().cpu()) > 0:
            lower_q = lower / lower.sum(dim=1, keepdim=True).clamp_min(1e-8)
            lower_kl = (p * (p.clamp_min(1e-8).log() - lower_q.clamp_min(1e-8).log())).sum(dim=1)
            lower_loss = (lower_kl * core).sum() / core.sum().clamp_min(1.0)
        else:
            lower_loss = q.new_tensor(0.0)
        width = out["rough_width"]
        width_loss = (width * core).sum() / core.sum().clamp_min(1.0) if float(core.sum().detach().cpu()) > 0 else q.new_tensor(0.0)
        upper_q = upper / upper.sum(dim=1, keepdim=True).clamp_min(1e-8)
        upper_entropy = -(upper_q * upper_q.clamp_min(1e-8).log()).sum(dim=1)
        boundary = 1.0 - core
        upper_entropy_loss = -((upper_entropy * boundary).sum() / boundary.sum().clamp_min(1.0)) if float(boundary.sum().detach().cpu()) > 0 else q.new_tensor(0.0)
        avg_q = q.mean(dim=0)
        uniform = torch.full_like(avg_q, 1.0 / avg_q.numel())
        balance_loss = F.kl_div(avg_q.clamp_min(1e-8).log(), uniform, reduction="batchmean")
        sep_loss = self.center_separation_loss(centers)
        anchor_loss = F.smooth_l1_loss(out["anchor_pred"], anchor)
        var_loss = self.variance_loss(out["latent"])
        total = (
            scmae
            + self.lower_weight * lower_loss
            + self.width_weight * width_loss
            + self.upper_entropy_weight * upper_entropy_loss
            + self.balance_weight * balance_loss
            + self.separation_weight * sep_loss
            + self.anchor_weight * anchor_loss
            + self.variance_weight * var_loss
        )
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "lower_kl_loss": float(lower_loss.detach().cpu()),
            "rough_width_loss": float(width_loss.detach().cpu()),
            "upper_entropy_loss": float(upper_entropy_loss.detach().cpu()),
            "balance_loss": float(balance_loss.detach().cpu()),
            "center_separation_loss": float(sep_loss.detach().cpu()),
            "anchor_loss": float(anchor_loss.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "core_fraction": float(core.mean().detach().cpu()),
            "rough_width": float(width.mean().detach().cpu()),
            "rough_certainty": float(out["rough_certainty"].mean().detach().cpu()),
            "membership_confidence": float(confidence.mean().detach().cpu()),
            "membership_entropy": float(out["membership_entropy"].mean().detach().cpu()),
        }
