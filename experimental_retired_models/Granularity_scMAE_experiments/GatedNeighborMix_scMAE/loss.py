from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class GatedFusionLoss(nn.Module):
    """rank13's MaskedScClusterLoss (reproduced) + a reliability-gated NeighborMix term.

    rank13 core (unchanged):
        scmae   = (1-mask_weight)*weighted_smooth_l1(recon, unscaled_target) + mask_weight*BCE(mask)
        cluster = confidence-gated KL(p||q)                 # gated by CONFIDENCE only
        consistency = ||norm(z_strong) - norm(z_weak)||^2
        loss = scmae + cluster_scale*cluster_weight*cluster + consistency_weight*consistency

    Added NeighborMix branch:
        pseudo  = scmae_term(neighbor-mixed view, target=REAL cell), weighted per-cell by r_i
        loss   += pseudo_weight * pseudo
    r_i gates ONLY the NeighborMix smoothing (the risky part). DEC is NOT gated
    by r_i (throttling it on imbalanced data kills the DEC gain).
    """

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        cluster_weight: float = 0.35,
        consistency_weight: float = 0.05,
        confidence_threshold: float = 0.35,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.cluster_weight = float(cluster_weight)
        self.consistency_weight = float(consistency_weight)
        self.confidence_threshold = float(confidence_threshold)

    def _scmae_percell(self, out: dict, target: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean(dim=1)
        mask_bce = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float(), reduction="none").mean(dim=1)
        return (1.0 - self.mask_weight) * rec + self.mask_weight * mask_bce

    def forward(
        self,
        out: dict,
        weak_out: dict,
        pseudo_out: dict | None,
        target: torch.Tensor,
        mask: torch.Tensor,
        pseudo_mask: torch.Tensor | None,
        reliability: torch.Tensor,
        pseudo_weight: float,
        p_target: torch.Tensor | None,
        cluster_scale: float,
    ) -> tuple[torch.Tensor, dict]:
        # rank13 scmae (mean over cells) — kept exactly
        scmae = self._scmae_percell(out, target, mask).mean()

        # rank13 confidence-gated DEC KL
        cluster = target.new_tensor(0.0)
        conf_frac = target.new_tensor(0.0)
        if p_target is not None and cluster_scale > 0:
            q = out["cluster_q"].clamp_min(1e-8)
            conf = p_target.max(dim=1).values
            gate = (conf >= self.confidence_threshold).float()
            conf_frac = gate.mean()
            kl = (p_target * (torch.log(p_target.clamp_min(1e-8)) - torch.log(q))).sum(dim=1)
            cluster = (kl * gate).sum() / gate.sum().clamp_min(1.0)

        # rank13 strong/weak consistency
        z1 = F.normalize(out["latent"], dim=1)
        z2 = F.normalize(weak_out["latent"], dim=1)
        consistency = (z1 - z2).pow(2).sum(dim=1).mean()

        # NEW: reliability-gated NeighborMix pseudo term (target = real cell)
        pseudo = target.new_tensor(0.0)
        if pseudo_out is not None and pseudo_mask is not None and pseudo_weight > 0:
            per_cell = self._scmae_percell(pseudo_out, target, pseudo_mask)
            denom = reliability.sum().clamp_min(1e-6)
            pseudo = (per_cell * reliability).sum() / denom

        loss = (
            scmae
            + float(cluster_scale) * self.cluster_weight * cluster
            + self.consistency_weight * consistency
            + float(pseudo_weight) * pseudo
        )
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "cluster_loss": float(cluster.detach().cpu()) if torch.is_tensor(cluster) else 0.0,
            "consistency_loss": float(consistency.detach().cpu()),
            "pseudo_loss": float(pseudo.detach().cpu()) if torch.is_tensor(pseudo) else 0.0,
            "confidence_fraction": float(conf_frac.detach().cpu()) if torch.is_tensor(conf_frac) else 0.0,
        }
