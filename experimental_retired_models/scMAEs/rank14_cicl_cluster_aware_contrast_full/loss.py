from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class CICLClusterAwareLoss(nn.Module):
    """scMAE plus filtered instance and cluster-aware contrastive losses."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        instance_weight: float = 0.1,
        cluster_weight: float = 0.08,
        temperature: float = 0.2,
        confidence_threshold: float = 0.45,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.instance_weight = float(instance_weight)
        self.cluster_weight = float(cluster_weight)
        self.temperature = float(temperature)
        self.confidence_threshold = float(confidence_threshold)

    def scmae_loss(self, out: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        return scmae, rec, mask_loss

    def filtered_instance_loss(
        self,
        z1: torch.Tensor,
        z2: torch.Tensor,
        pseudo: torch.Tensor | None,
        confident: torch.Tensor | None,
    ) -> torch.Tensor:
        b = z1.shape[0]
        z = torch.cat([z1, z2], dim=0)
        sim = torch.matmul(z, z.T) / self.temperature
        eye = torch.eye(2 * b, device=z.device, dtype=torch.bool)
        pos = torch.arange(2 * b, device=z.device)
        pos = (pos + b) % (2 * b)
        denom_mask = ~eye
        if pseudo is not None and confident is not None:
            lab = torch.cat([pseudo, pseudo], dim=0)
            conf = torch.cat([confident, confident], dim=0)
            same = lab[:, None].eq(lab[None, :]) & conf[:, None] & conf[None, :]
            pos_mask = torch.zeros_like(same)
            pos_mask[torch.arange(2 * b, device=z.device), pos] = True
            denom_mask = denom_mask & (~same | pos_mask)
        logits = sim.masked_fill(~denom_mask, -1e9)
        return F.cross_entropy(logits, pos)

    def cluster_aware_loss(
        self,
        z1: torch.Tensor,
        z2: torch.Tensor,
        pseudo: torch.Tensor | None,
        confident: torch.Tensor | None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if pseudo is None or confident is None:
            return z1.new_tensor(0.0), z1.new_tensor(0.0)
        b = z1.shape[0]
        z = torch.cat([z1, z2], dim=0)
        labels = torch.cat([pseudo, pseudo], dim=0)
        conf = torch.cat([confident, confident], dim=0)
        sim = torch.matmul(z, z.T) / self.temperature
        eye = torch.eye(2 * b, device=z.device, dtype=torch.bool)
        same = labels[:, None].eq(labels[None, :])
        pos_mask = same & (~eye) & conf[:, None] & conf[None, :]
        denom_mask = (~eye) & conf[:, None]
        log_prob = sim - torch.logsumexp(sim.masked_fill(~denom_mask, -1e9), dim=1, keepdim=True)
        pos_count = pos_mask.sum(dim=1)
        valid = pos_count > 0
        if not bool(valid.any()):
            return z1.new_tensor(0.0), z1.new_tensor(0.0)
        loss = -(log_prob * pos_mask.float()).sum(dim=1)[valid] / pos_count[valid].float()
        return loss.mean(), valid.float().mean()

    def forward(
        self,
        out1: dict[str, torch.Tensor],
        out2: dict[str, torch.Tensor],
        target: torch.Tensor,
        mask: torch.Tensor,
        pseudo: torch.Tensor | None,
        confidence: torch.Tensor | None,
        contrast_scale: float,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        scmae, rec, mask_loss = self.scmae_loss(out1, target, mask)
        confident = None if confidence is None else confidence >= self.confidence_threshold
        instance = self.filtered_instance_loss(out1["projection"], out2["projection"], pseudo, confident)
        cluster, cluster_valid = self.cluster_aware_loss(out1["projection"], out2["projection"], pseudo, confident)
        loss = scmae + float(contrast_scale) * (self.instance_weight * instance + self.cluster_weight * cluster)
        conf_frac = target.new_tensor(0.0) if confident is None else confident.float().mean()
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "instance_loss": float(instance.detach().cpu()),
            "cluster_loss": float(cluster.detach().cpu()),
            "confidence_fraction": float(conf_frac.detach().cpu()),
            "cluster_valid_fraction": float(cluster_valid.detach().cpu()),
        }

