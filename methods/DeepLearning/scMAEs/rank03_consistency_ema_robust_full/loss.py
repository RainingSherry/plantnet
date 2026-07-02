from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


def masked_mean(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    return (values * mask).sum() / mask.sum().clamp_min(1.0)


def pseudo_huber(diff: torch.Tensor, c: float) -> torch.Tensor:
    c_t = torch.as_tensor(float(c), dtype=diff.dtype, device=diff.device)
    return c_t * c_t * (torch.sqrt(1.0 + (diff / c_t).pow(2)) - 1.0)


class ConsistencyEMALoss(nn.Module):
    def __init__(
        self,
        reconstruction_weight: float = 1.0,
        mask_weight: float = 0.3,
        consistency_weight: float = 0.5,
        huber_beta: float = 1.0,
        pseudo_huber_c: float = 0.03,
    ):
        super().__init__()
        self.reconstruction_weight = float(reconstruction_weight)
        self.mask_weight = float(mask_weight)
        self.consistency_weight = float(consistency_weight)
        self.huber_beta = float(huber_beta)
        self.pseudo_huber_c = float(pseudo_huber_c)

    def forward(
        self,
        student: dict[str, torch.Tensor],
        teacher_projection: torch.Tensor,
        expression_target: torch.Tensor,
        mask: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        rec = F.smooth_l1_loss(
            student["reconstruction"],
            expression_target,
            reduction="none",
            beta=self.huber_beta,
        )
        reconstruction_loss = masked_mean(rec, mask.float())
        mask_loss = F.binary_cross_entropy_with_logits(student["mask_logits"], mask.float())
        s = F.normalize(student["projection"], dim=1)
        t = F.normalize(teacher_projection.detach(), dim=1)
        consistency_loss = pseudo_huber(s - t, self.pseudo_huber_c).mean()
        loss = (
            self.reconstruction_weight * reconstruction_loss
            + self.mask_weight * mask_loss
            + self.consistency_weight * consistency_loss
        )
        parts = {
            "loss": float(loss.detach().cpu()),
            "reconstruction_loss": float(reconstruction_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "consistency_loss": float(consistency_loss.detach().cpu()),
        }
        return loss, parts

