from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


def compute_gate(sharp_p: torch.Tensor, q: torch.Tensor, kappa: float) -> float:
    """Dataset-level sharp gate from DEC-KL reference.

    kl_ref = mean KL(sharpen(q) || q): how far the sharpened target is from the
    current soft assignment == how confidently the data clusters.
      Quake  -> kl_ref ~0.01 -> gate ~1 (sharp)
      Macosko-> kl_ref high  -> gate ~0 (soft)
      gate = 1 / (1 + (kl_ref/kappa)^2)
    """
    p = torch.as_tensor(sharp_p).clamp_min(1e-8)
    qq = torch.as_tensor(q).clamp_min(1e-8)
    kl_ref = float((p * (p.log() - qq.log())).sum(dim=1).mean())
    gate = 1.0 / (1.0 + (kl_ref / max(kappa, 1e-6)) ** 2)
    return float(gate), kl_ref


class AdaptiveSwitchLoss(nn.Module):
    """scMAE recon+mask + consistency + gate-blended (sharp DEC | soft fuzzy/entropy) + variance.

    L = scmae + consistency
        + cluster_scale * cluster_weight * [ gate*L_sharp + (1-gate)*L_soft ]
        + variance_weight * L_variance          (always on, cheap anti-collapse)

    L_sharp : confidence-gated KL(p_sharp || q)          (rank13)
    L_soft  : fuzzy-core KL(p_sharp || q) on core cells only
              + entropy maximization on boundary cells   (rank29 soft levers)
    Both use the SAME sharpened target p_sharp; the difference is WHERE KL is
    applied (all-confident vs core-only) and the boundary-entropy push.
    NO balance loss (dissection: it hurts imbalanced scRNA).
    """

    def __init__(self, masked_data_weight=0.75, mask_weight=0.65, cluster_weight=0.35,
                 consistency_weight=0.05, variance_weight=0.02, entropy_weight=0.10,
                 confidence_threshold=0.35):
        super().__init__()
        self.mdw = float(masked_data_weight)
        self.mw = float(mask_weight)
        self.cw = float(cluster_weight)
        self.consw = float(consistency_weight)
        self.varw = float(variance_weight)
        self.entw = float(entropy_weight)
        self.conf_th = float(confidence_threshold)

    def __init__(self, masked_data_weight=0.75, mask_weight=0.65, cluster_weight=0.35,
                 consistency_weight=0.05, variance_weight=0.02, entropy_weight=0.10,
                 confidence_threshold=0.35, var_mode="hinge"):
        super().__init__()
        self.mdw = float(masked_data_weight)
        self.mw = float(mask_weight)
        self.cw = float(cluster_weight)
        self.consw = float(consistency_weight)
        self.varw = float(variance_weight)
        self.entw = float(entropy_weight)
        self.conf_th = float(confidence_threshold)
        self.var_mode = str(var_mode)

    def _anticollapse_loss(self, z):
        """Anti-collapse regularizers (mechanism-comparison):
          hinge : VICReg std hinge  relu(1 - std_d).mean()      (spread each dim)
          cov   : VICReg covariance off-diagonal decorrelation  (whiten dims)
          koleo : mean -log(nn distance) uniformity (DINOv2-style)
          both  : hinge + cov
        """
        if self.var_mode in ("hinge", "both"):
            std = torch.sqrt(z.var(dim=0) + 1e-4)
            hinge = F.relu(1.0 - std).mean()
            if self.var_mode == "hinge":
                return hinge
        if self.var_mode in ("cov", "both"):
            zc = z - z.mean(dim=0, keepdim=True)
            n, d = zc.shape
            cov = (zc.T @ zc) / max(1, n - 1)
            off = cov - torch.diag(torch.diag(cov))
            cov_loss = off.pow(2).sum() / d
            return cov_loss if self.var_mode == "cov" else hinge + cov_loss
        if self.var_mode == "koleo":
            zn = F.normalize(z, dim=1)
            d = torch.cdist(zn, zn)
            d = d + torch.eye(zn.shape[0], device=z.device) * 1e6
            nn_d = d.min(dim=1).values.clamp_min(1e-6)
            return -torch.log(nn_d).mean()
        return z.new_tensor(0.0)

    def forward(self, out, weak_out, target, mask, p_target, clusterability, cluster_scale, gate):
        w = mask * self.mdw + (1.0 - mask) * (1.0 - self.mdw)
        rec = (w * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mw) * rec + self.mw * mask_loss

        z1 = F.normalize(out["latent"], dim=1)
        z2 = F.normalize(weak_out["latent"], dim=1)
        consistency = (z1 - z2).pow(2).sum(dim=1).mean()
        var_loss = self._anticollapse_loss(out["latent"])

        sharp = target.new_tensor(0.0)
        soft = target.new_tensor(0.0)
        if p_target is not None and cluster_scale > 0:
            q = out["cluster_q"].clamp_min(1e-8)
            kl = (p_target * (torch.log(p_target.clamp_min(1e-8)) - torch.log(q))).sum(dim=1)
            conf = p_target.max(dim=1).values
            # SHARP: KL on all confident cells (rank13)
            g_sharp = (conf >= self.conf_th).float()
            sharp = (kl * g_sharp).sum() / g_sharp.sum().clamp_min(1.0)
            # SOFT: KL only on high-clusterability CORE cells + entropy push on boundary
            if clusterability is not None:
                core = (clusterability >= 0.5).float()
                soft_kl = (kl * core).sum() / core.sum().clamp_min(1.0)
                ent = -(q * q.log()).sum(dim=1)                 # per-cell entropy
                boundary = (1.0 - clusterability).clamp(0.0, 1.0)
                ent_term = -(ent * boundary).sum() / boundary.sum().clamp_min(1.0)  # maximize
                soft = soft_kl + self.entw * ent_term
            else:
                soft = sharp

        cluster = gate * sharp + (1.0 - gate) * soft
        loss = scmae + self.consw * consistency + float(cluster_scale) * self.cw * cluster + self.varw * var_loss
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "sharp_loss": float(sharp.detach().cpu()) if torch.is_tensor(sharp) else 0.0,
            "soft_loss": float(soft.detach().cpu()) if torch.is_tensor(soft) else 0.0,
            "variance_loss": float(var_loss.detach().cpu()),
            "consistency_loss": float(consistency.detach().cpu()),
        }
