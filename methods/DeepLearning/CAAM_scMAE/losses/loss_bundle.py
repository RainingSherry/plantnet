from __future__ import annotations

import torch

from .mask_prediction import mask_prediction_loss
from .reconstruction import reconstruction_losses


def student_loss_bundle(
    x: torch.Tensor,
    x_hat: torch.Tensor,
    mask_logits: torch.Tensor,
    mask: torch.Tensor,
    *,
    lambda_visible: float,
    lambda_mask: float,
) -> dict[str, torch.Tensor]:
    rec = reconstruction_losses(x, x_hat, mask, lambda_visible)
    loss_mask = mask_prediction_loss(mask_logits, mask)
    loss = rec["loss_rec"] + float(lambda_mask) * loss_mask
    return {**rec, "loss_mask": loss_mask, "loss_student": loss}

