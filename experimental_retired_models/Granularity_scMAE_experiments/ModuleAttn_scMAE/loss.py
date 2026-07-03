from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ModuleAttnLoss(nn.Module):
    """Identical objective to AdaptiveSwitch's winning config (force_gate=1):
    scmae(recon + mask) + confidence-gated DEC KL + per-dim variance floor.
    Only the encoder architecture differs (gene-module attention), so gains are
    attributable to the module structure.
    """

    def __init__(self, masked_data_weight=0.75, mask_weight=0.65, cluster_weight=0.35,
                 variance_weight=0.02, confidence_threshold=0.35):
        super().__init__()
        self.mdw = float(masked_data_weight)
        self.mw = float(mask_weight)
        self.cw = float(cluster_weight)
        self.varw = float(variance_weight)
        self.conf_th = float(confidence_threshold)

    @staticmethod
    def _variance_floor(z):
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(1.0 - std).mean()

    def forward(self, out, target, mask, p_target, cluster_scale):
        w = mask * self.mdw + (1.0 - mask) * (1.0 - self.mdw)
        rec = (w * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mw) * rec + self.mw * mask_loss

        cluster = target.new_tensor(0.0)
        if p_target is not None and cluster_scale > 0:
            q = out["cluster_q"].clamp_min(1e-8)
            kl = (p_target * (torch.log(p_target.clamp_min(1e-8)) - torch.log(q))).sum(dim=1)
            conf = p_target.max(dim=1).values
            gate = (conf >= self.conf_th).float()
            cluster = (kl * gate).sum() / gate.sum().clamp_min(1.0)

        var_loss = self._variance_floor(out["latent"])
        loss = scmae + float(cluster_scale) * self.cw * cluster + self.varw * var_loss
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "cluster_loss": float(cluster.detach().cpu()) if torch.is_tensor(cluster) else 0.0,
            "variance_loss": float(var_loss.detach().cpu()),
        }
