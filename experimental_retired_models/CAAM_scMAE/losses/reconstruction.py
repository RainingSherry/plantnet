from __future__ import annotations

import torch


def reconstruction_losses(x: torch.Tensor, x_hat: torch.Tensor, mask: torch.Tensor, lambda_visible: float, eps: float = 1.0e-8):
    sq = (x - x_hat).pow(2)
    mask = mask.float()
    visible = 1.0 - mask
    loss_masked = (mask * sq).sum() / (mask.sum() + eps)
    loss_visible = (visible * sq).sum() / (visible.sum() + eps)
    loss_rec = loss_masked + float(lambda_visible) * loss_visible
    return {"loss_rec_masked": loss_masked, "loss_rec_visible": loss_visible, "loss_rec": loss_rec}

