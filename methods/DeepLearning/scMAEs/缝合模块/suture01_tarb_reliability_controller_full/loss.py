from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class SutureTARBLoss(nn.Module):
    """scMAE loss + confidence-gated DEC + TARB controller regularizers."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        cluster_weight: float = 0.35,
        consistency_weight: float = 0.05,
        confidence_threshold: float = 0.35,
        balance_weight: float = 0.01,
        variance_weight: float = 0.02,
        conservative_weight: float = 0.01,
        variance_floor: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.cluster_weight = float(cluster_weight)
        self.consistency_weight = float(consistency_weight)
        self.confidence_threshold = float(confidence_threshold)
        self.balance_weight = float(balance_weight)
        self.variance_weight = float(variance_weight)
        self.conservative_weight = float(conservative_weight)
        self.variance_floor = float(variance_floor)

    def _scmae_terms(self, out: dict, target: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_bce = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        return (1.0 - self.mask_weight) * rec + self.mask_weight * mask_bce, rec, mask_bce

    def forward(
        self,
        out: dict,
        weak_out: dict,
        target: torch.Tensor,
        mask: torch.Tensor,
        p_target: torch.Tensor | None,
        cluster_scale: float,
        reliability: torch.Tensor,
    ) -> tuple[torch.Tensor, dict]:
        scmae, rec, mask_bce = self._scmae_terms(out, target, mask)

        cluster = target.new_tensor(0.0)
        conf_frac = target.new_tensor(0.0)
        if p_target is not None and cluster_scale > 0:
            q = out["cluster_q"].clamp_min(1e-8)
            conf = p_target.max(dim=1).values
            gate = (conf >= self.confidence_threshold).float()
            conf_frac = gate.mean()
            kl = (p_target * (torch.log(p_target.clamp_min(1e-8)) - torch.log(q))).sum(dim=1)
            cluster = (kl * gate).sum() / gate.sum().clamp_min(1.0)

        z1 = F.normalize(out["latent"], dim=1)
        z2 = F.normalize(weak_out["latent"], dim=1)
        consistency = (z1 - z2).pow(2).sum(dim=1).mean()

        opw = out["operation_weights"].clamp_min(1e-8)
        mean_w = opw.mean(dim=0)
        target_w = torch.full_like(mean_w, 1.0 / mean_w.numel())
        balance = F.kl_div(mean_w.log(), target_w, reduction="sum")
        variance = out["latent"].std(dim=0).mean()
        variance_loss = F.relu(out["latent"].new_tensor(self.variance_floor) - variance)
        risky = opw[:, 3:].sum(dim=1)
        conservative = (risky * (1.0 - reliability.clamp(0.0, 1.0))).mean()

        loss = (
            scmae
            + float(cluster_scale) * self.cluster_weight * cluster
            + self.consistency_weight * consistency
            + self.balance_weight * balance
            + self.variance_weight * variance_loss
            + self.conservative_weight * conservative
        )
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_bce.detach().cpu()),
            "cluster_loss": float(cluster.detach().cpu()),
            "consistency_loss": float(consistency.detach().cpu()),
            "operation_balance_loss": float(balance.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "conservative_loss": float(conservative.detach().cpu()),
            "confidence_fraction": float(conf_frac.detach().cpu()),
            "latent_std": float(variance.detach().cpu()),
        }
