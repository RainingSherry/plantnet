from __future__ import annotations

import torch
import torch.nn.functional as F


def mask_prediction_loss(mask_logits: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    return F.binary_cross_entropy_with_logits(mask_logits, mask.float())

