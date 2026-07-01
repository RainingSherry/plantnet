from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class Data2VecLatentLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.65,
        distill_weight: float = 0.08,
        variance_weight: float = 0.02,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.distill_weight = float(distill_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def normalize_target(target: torch.Tensor) -> torch.Tensor:
        target = target - target.mean(dim=0, keepdim=True)
        return target / target.std(dim=0, keepdim=True).clamp_min(1e-4)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(1.0 - std).mean()

    def forward(
        self,
        student_out: dict[str, torch.Tensor],
        teacher_target: torch.Tensor,
        target: torch.Tensor,
        mask: torch.Tensor,
        distill_scale: float,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(student_out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(student_out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        teacher_target = self.normalize_target(teacher_target.detach())
        pred = self.normalize_target(student_out["prediction"])
        distill = F.smooth_l1_loss(pred, teacher_target)
        var_loss = self.variance_loss(student_out["latent"])
        total = scmae + float(distill_scale) * self.distill_weight * distill + self.variance_weight * var_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "distill_loss": float(distill.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
        }

