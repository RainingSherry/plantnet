from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class AdaptiveGranularityLoss(nn.Module):
    """scMAE reconstruction + mask BCE + strong/weak consistency + adaptive-granularity KL.

    The KL target p is the per-cell ADAPTIVE target (built outside as
    c_i*sharpen(q_i) + (1-c_i)*q_i). Boundary/rare cells therefore have
    p_i ≈ q_i and contribute ~0 KL automatically — no explicit loss-weight gate
    needed. A light confidence floor still skips near-uniform rows.
    """

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        cluster_weight: float = 0.35,
        consistency_weight: float = 0.05,
        anchor_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.cluster_weight = float(cluster_weight)
        self.consistency_weight = float(consistency_weight)
        self.anchor_weight = float(anchor_weight)

    def forward(
        self,
        out: dict,
        weak_out: dict,
        target: torch.Tensor,
        mask: torch.Tensor,
        p_target: torch.Tensor | None,
        cluster_scale: float,
    ) -> tuple[torch.Tensor, dict]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss

        cluster = target.new_tensor(0.0)
        eff_frac = target.new_tensor(0.0)
        if p_target is not None and cluster_scale > 0:
            q = out["cluster_q"].clamp_min(1e-8)
            kl = (p_target * (torch.log(p_target.clamp_min(1e-8)) - torch.log(q))).sum(dim=1)
            cluster = kl.mean()
            eff_frac = (kl > 1e-4).float().mean()

        z1 = F.normalize(out["latent"], dim=1)
        z2 = F.normalize(weak_out["latent"], dim=1)
        consistency = (z1 - z2).pow(2).sum(dim=1).mean()

        loss = scmae + float(cluster_scale) * self.cluster_weight * cluster + self.consistency_weight * consistency
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "cluster_loss": float(cluster.detach().cpu()) if torch.is_tensor(cluster) else 0.0,
            "consistency_loss": float(consistency.detach().cpu()),
            "active_kl_fraction": float(eff_frac.detach().cpu()) if torch.is_tensor(eff_frac) else 0.0,
        }
