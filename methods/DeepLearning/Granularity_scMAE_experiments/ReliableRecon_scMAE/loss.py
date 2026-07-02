from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ReliableReconLoss(nn.Module):
    """scMAE loss with per-cell-per-gene reliability precision weighting.

    w = w_pos * ((1-lambda) + lambda * r_hat),  r_hat = r/mean(r) (mean-1 normalized
    so total loss scale is preserved). lambda=0 -> exact vanilla scMAE.
    """

    def __init__(self, masked_data_weight=0.75, mask_weight=0.65, reliability_lambda=1.0):
        super().__init__()
        self.mdw = float(masked_data_weight)
        self.mw = float(mask_weight)
        self.lam = float(reliability_lambda)

    def forward(self, out, target, mask, reliability):
        w_pos = mask * self.mdw + (1.0 - mask) * (1.0 - self.mdw)
        if reliability is not None and self.lam > 0:
            r_hat = reliability / reliability.mean().clamp_min(1e-6)
            w = w_pos * ((1.0 - self.lam) + self.lam * r_hat)
        else:
            w = w_pos
        rec = (w * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        loss = (1.0 - self.mw) * rec + self.mw * mask_loss
        return loss, {"loss": float(loss.detach().cpu()), "reconstruction_loss": float(rec.detach().cpu()),
                      "mask_loss": float(mask_loss.detach().cpu())}
