from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class DinoBloomDistillLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.7,
        distill_weight: float = 0.05,
        student_temp: float = 0.1,
        teacher_temp: float = 0.04,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.distill_weight = float(distill_weight)
        self.student_temp = float(student_temp)
        self.teacher_temp = float(teacher_temp)

    def forward(
        self,
        student_out: dict[str, torch.Tensor],
        teacher_logits: torch.Tensor,
        center: torch.Tensor,
        target: torch.Tensor,
        mask: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = torch.mul(weights, F.mse_loss(student_out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(student_out["mask_logits"], mask.float())
        scmae_loss = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        teacher_probs = F.softmax((teacher_logits.detach() - center) / self.teacher_temp, dim=1)
        student_log_probs = F.log_softmax(student_out["proto_logits"] / self.student_temp, dim=1)
        distill = -(teacher_probs * student_log_probs).sum(dim=1).mean()
        loss = scmae_loss + self.distill_weight * distill
        return loss, {
            "loss": float(loss.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "distill_loss": float(distill.detach().cpu()),
        }

