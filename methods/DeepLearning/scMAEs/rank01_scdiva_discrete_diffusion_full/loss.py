from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F

from model import cosine_mask_schedule


@dataclass
class ScDiVaLossParts:
    total: torch.Tensor
    token: torch.Tensor
    value: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor


def q_sample_discrete_diffusion(
    tokens: torch.Tensor,
    t: torch.Tensor,
    steps: int,
    max_mask_ratio: float,
    mask_token_id: int,
) -> tuple[torch.Tensor, torch.Tensor]:
    """Forward masked discrete diffusion.

    Mask semantics: returned `mask` is 1.0 exactly where the clean expression
    token was replaced by the absorbing mask token. The denominator for masked
    losses is the number of masked gene tokens, clamped only to avoid division
    by zero in degenerate smoke tests.
    """
    if tokens.ndim != 2:
        raise ValueError(f"tokens must be [batch, genes], got {tuple(tokens.shape)}")
    if t.ndim != 1 or t.shape[0] != tokens.shape[0]:
        raise ValueError(f"t must be [batch], got {tuple(t.shape)} for tokens {tuple(tokens.shape)}")
    prob = cosine_mask_schedule(t, steps, max_mask_ratio).view(-1, 1)
    mask_bool = torch.bernoulli(prob.expand_as(tokens).float()).bool()
    corrupted = torch.where(mask_bool, torch.full_like(tokens, int(mask_token_id)), tokens)
    return corrupted, mask_bool.float()


def scdiva_loss(
    token_logits: torch.Tensor,
    clean_tokens: torch.Tensor,
    value_pred: torch.Tensor,
    clean_values: torch.Tensor,
    mask_logits: torch.Tensor,
    mask: torch.Tensor,
    n_bins: int,
    value_weight: float = 0.35,
    mask_weight: float = 0.15,
) -> ScDiVaLossParts:
    if token_logits.ndim != 3:
        raise ValueError(f"token_logits must be [batch, genes, bins], got {tuple(token_logits.shape)}")
    if token_logits.shape[:2] != clean_tokens.shape:
        raise ValueError("token_logits and clean_tokens batch/gene dimensions differ")
    if token_logits.shape[2] != int(n_bins):
        raise ValueError(f"token_logits bins must be {n_bins}, got {token_logits.shape[2]}")
    if value_pred.shape != clean_values.shape or mask_logits.shape != mask.shape:
        raise ValueError(
            "value_pred, clean_values, mask_logits, and mask must share [batch, genes] dimensions"
        )
    mask = mask.to(dtype=value_pred.dtype, device=value_pred.device)
    denom = mask.sum().clamp_min(1.0)

    ce = F.cross_entropy(
        token_logits.reshape(-1, int(n_bins)),
        clean_tokens.reshape(-1),
        reduction="none",
    ).view_as(mask)
    token_loss = (ce * mask).sum() / denom
    value_loss = (F.smooth_l1_loss(value_pred, clean_values, reduction="none") * mask).sum() / denom
    mask_loss = F.binary_cross_entropy_with_logits(mask_logits, mask)
    total = token_loss + float(value_weight) * value_loss + float(mask_weight) * mask_loss
    return ScDiVaLossParts(total, token_loss, value_loss, mask_loss, mask.mean())

