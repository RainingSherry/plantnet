from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class AudioMAEModulePatchLoss(nn.Module):
    """scMAE reconstruction/mask loss plus AudioMAE-style module patch target and view consistency."""

    def __init__(
        self,
        masked_data_weight: float = 0.8,
        mask_weight: float = 0.6,
        module_weight: float = 0.12,
        consistency_weight: float = 0.04,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.module_weight = float(module_weight)
        self.consistency_weight = float(consistency_weight)

    def forward(
        self,
        out_a: dict[str, torch.Tensor],
        out_b: dict[str, torch.Tensor],
        target: torch.Tensor,
        module_target: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        mask = out_a["gene_mask"].float()
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        reconstruction_loss = (weights * F.smooth_l1_loss(out_a["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out_a["mask_logits"], mask)
        module_loss_raw = F.smooth_l1_loss(out_a["module_reconstruction"], module_target, reduction="none").mean(dim=2)
        module_loss = (module_loss_raw * out_a["patch_mask"]).sum() / out_a["patch_mask"].sum().clamp_min(1.0)
        consistency_loss = 1.0 - F.cosine_similarity(out_a["latent"], out_b["latent"].detach(), dim=1).mean()
        scmae_loss = (1.0 - self.mask_weight) * reconstruction_loss + self.mask_weight * mask_loss
        total = scmae_loss + self.module_weight * module_loss + self.consistency_weight * consistency_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "reconstruction_loss": float(reconstruction_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "module_loss": float(module_loss.detach().cpu()),
            "consistency_loss": float(consistency_loss.detach().cpu()),
        }
