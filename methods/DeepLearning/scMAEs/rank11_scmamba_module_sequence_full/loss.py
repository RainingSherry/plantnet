from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ScMambaModuleLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        module_weight: float = 0.15,
        consistency_weight: float = 0.05,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.module_weight = float(module_weight)
        self.consistency_weight = float(consistency_weight)

    def forward(
        self,
        out: dict[str, torch.Tensor],
        weak_out: dict[str, torch.Tensor],
        target: torch.Tensor,
        module_target: torch.Tensor,
        mask: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        module_loss = F.smooth_l1_loss(out["module_reconstruction"], module_target)
        z1 = F.normalize(out["latent"], dim=1)
        z2 = F.normalize(weak_out["latent"], dim=1)
        consistency = (z1 - z2).pow(2).sum(dim=1).mean()
        loss = scmae + self.module_weight * module_loss + self.consistency_weight * consistency
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "module_loss": float(module_loss.detach().cpu()),
            "consistency_loss": float(consistency.detach().cpu()),
        }

