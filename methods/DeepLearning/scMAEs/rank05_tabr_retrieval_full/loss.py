from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class TabRLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    entropy: torch.Tensor
    mask_rate: torch.Tensor


def random_mask_expression(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    if x.ndim != 2:
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    mask = (torch.rand_like(x) < float(mask_ratio)).float()
    return x * (1.0 - mask), mask


def tabr_retrieval_loss(
    reconstruction: torch.Tensor,
    clean: torch.Tensor,
    mask_logits: torch.Tensor,
    mask: torch.Tensor,
    context_probs: torch.Tensor,
    *,
    mask_weight: float,
    entropy_weight: float,
) -> TabRLossParts:
    if reconstruction.shape != clean.shape or mask_logits.shape != mask.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, mask_logits, and mask must share [batch, genes] dimensions")
    if context_probs.ndim != 2 or context_probs.shape[0] != clean.shape[0]:
        raise ValueError("context_probs must be [batch, context_size]")
    denom = mask.sum().clamp_min(1.0)
    reconstruction_loss = (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum() / denom
    mask_loss = F.binary_cross_entropy_with_logits(mask_logits, mask)
    entropy = -(context_probs.clamp_min(1e-12) * context_probs.clamp_min(1e-12).log()).sum(dim=1).mean()
    total = reconstruction_loss + float(mask_weight) * mask_loss - float(entropy_weight) * entropy
    return TabRLossParts(total, reconstruction_loss, mask_loss, entropy, mask.mean())
