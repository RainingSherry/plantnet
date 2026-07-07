from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


def _masked_mean(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    denom = mask.sum().clamp_min(1.0)
    return (values * mask).sum() / denom


class ScDiVaLoss(nn.Module):
    """Dual denoising loss: expression regression + absorbing-mask + token CE."""

    def __init__(
        self,
        expression_weight: float = 1.0,
        mask_weight: float = 0.3,
        token_weight: float = 0.5,
        huber_beta: float = 1.0,
    ):
        super().__init__()
        self.expression_weight = float(expression_weight)
        self.mask_weight = float(mask_weight)
        self.token_weight = float(token_weight)
        self.huber_beta = float(huber_beta)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        expression_target: torch.Tensor,
        token_target: torch.Tensor,
        mask: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        mask = mask.float()
        expr_raw = F.smooth_l1_loss(
            outputs["reconstruction"],
            expression_target,
            reduction="none",
            beta=self.huber_beta,
        )
        expression_loss = _masked_mean(expr_raw, mask)
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)

        masked_positions = mask.bool()
        if bool(masked_positions.any()):
            token_loss = F.cross_entropy(
                outputs["token_logits"][masked_positions],
                token_target[masked_positions].long(),
            )
        else:
            token_loss = outputs["token_logits"].sum() * 0.0

        loss = (
            self.expression_weight * expression_loss
            + self.mask_weight * mask_loss
            + self.token_weight * token_loss
        )
        parts = {
            "loss": float(loss.detach().cpu()),
            "expression_loss": float(expression_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "token_loss": float(token_loss.detach().cpu()),
        }
        return loss, parts

