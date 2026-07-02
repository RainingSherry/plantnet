from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class VideoMAEGeneModuleLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.80,
        mask_weight: float = 0.55,
        module_weight: float = 0.12,
        consistency_weight: float = 0.03,
        variance_weight: float = 0.01,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.module_weight = float(module_weight)
        self.consistency_weight = float(consistency_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        gene_mask: torch.Tensor,
        second_latent: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = gene_mask * self.masked_data_weight + (1.0 - gene_mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], gene_mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        patch_mask = out["patch_mask"]
        module_loss = F.smooth_l1_loss(out["module_pred"], out["module_target"], reduction="none")
        module_loss = (module_loss * patch_mask).sum() / patch_mask.sum().clamp_min(1.0)
        if second_latent is None:
            consistency = out["latent"].new_tensor(0.0)
        else:
            consistency = F.smooth_l1_loss(F.normalize(out["latent"], dim=-1), F.normalize(second_latent.detach(), dim=-1))
        var_loss = self.variance_loss(out["latent"])
        total = scmae + self.module_weight * module_loss + self.consistency_weight * consistency + self.variance_weight * var_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "module_loss": float(module_loss.detach().cpu()),
            "consistency_loss": float(consistency.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "masked_patch_fraction": float(patch_mask.mean().detach().cpu()),
        }
