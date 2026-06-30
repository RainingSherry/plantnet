from __future__ import annotations

import torch
import torch.nn.functional as F


def poisson_nll_like(pred_nonnegative: torch.Tensor, target_nonnegative: torch.Tensor) -> torch.Tensor:
    """Small count-aware auxiliary loss for future variants."""

    pred = F.softplus(pred_nonnegative).clamp_min(1e-6)
    target = target_nonnegative.clamp_min(0.0)
    return (pred - target * pred.log()).mean()

