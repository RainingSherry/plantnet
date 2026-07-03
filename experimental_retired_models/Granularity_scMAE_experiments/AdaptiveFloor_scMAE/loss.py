from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class AdaptiveFloorLoss(nn.Module):
    """scMAE recon+mask + consistency + confidence-gated DEC KL + adaptive variance floor.

    L = scmae + consistency_w·consistency + cluster_scale·cluster_w·KL_gated
        + floor_w · L_floor

    L_floor (per-dimension std hinge toward a target):
      fixed : mean_d relu(1 - std_d)                 # legacy VICReg (baseline)
      ref   : mean_d relu(floor_scale·ref_std_d - std_d)   # THEORY-DRIVEN adaptive
              ref_std_d = pre-DEC (post-warmup) per-dim std, set once and frozen.
    """

    def __init__(self, masked_data_weight=0.75, mask_weight=0.65, cluster_weight=0.35,
                 consistency_weight=0.05, floor_weight=0.02, confidence_threshold=0.35,
                 floor_mode="ref", floor_scale=1.0):
        super().__init__()
        self.mdw = float(masked_data_weight)
        self.mw = float(mask_weight)
        self.cw = float(cluster_weight)
        self.consw = float(consistency_weight)
        self.floorw = float(floor_weight)
        self.conf_th = float(confidence_threshold)
        self.floor_mode = str(floor_mode)
        self.floor_scale = float(floor_scale)
        self.ref_std = None  # (hidden,) set by run.py after warmup for floor_mode='ref'

    def set_reference(self, ref_std: torch.Tensor) -> None:
        self.ref_std = ref_std.detach().clone()

    def _floor_loss(self, z):
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        if self.floor_mode == "ref" and self.ref_std is not None:
            target = (self.floor_scale * self.ref_std).to(std.device)
            return F.relu(target - std).mean()
        # fixed VICReg std>=1 (scaled)
        return F.relu(self.floor_scale - std).mean()

    def forward(self, out, weak_out, target, mask, p_target, cluster_scale):
        w = mask * self.mdw + (1.0 - mask) * (1.0 - self.mdw)
        rec = (w * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mw) * rec + self.mw * mask_loss

        z1 = F.normalize(out["latent"], dim=1)
        z2 = F.normalize(weak_out["latent"], dim=1)
        consistency = (z1 - z2).pow(2).sum(dim=1).mean()

        cluster = target.new_tensor(0.0)
        if p_target is not None and cluster_scale > 0:
            q = out["cluster_q"].clamp_min(1e-8)
            kl = (p_target * (torch.log(p_target.clamp_min(1e-8)) - torch.log(q))).sum(dim=1)
            gate = (p_target.max(dim=1).values >= self.conf_th).float()
            cluster = (kl * gate).sum() / gate.sum().clamp_min(1.0)

        floor = self._floor_loss(out["latent"])
        loss = scmae + self.consw * consistency + float(cluster_scale) * self.cw * cluster + self.floorw * floor
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "cluster_loss": float(cluster.detach().cpu()) if torch.is_tensor(cluster) else 0.0,
            "consistency_loss": float(consistency.detach().cpu()),
            "floor_loss": float(floor.detach().cpu()),
        }
