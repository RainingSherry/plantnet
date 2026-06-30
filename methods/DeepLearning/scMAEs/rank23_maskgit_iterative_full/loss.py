from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import torch
import torch.nn.functional as F


@dataclass
class MaskGITLossParts:
    total: torch.Tensor
    token: torch.Tensor
    reconstruction: torch.Tensor
    confidence: torch.Tensor
    mask_rate: torch.Tensor


def patchify_expression(x: torch.Tensor, patch_size: int, num_patches: int) -> torch.Tensor:
    if x.ndim != 2:
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    target_genes = int(patch_size) * int(num_patches)
    pad = target_genes - x.shape[1]
    if pad < 0:
        raise ValueError("patch_size * num_patches cannot be smaller than gene count")
    padded = F.pad(x, (0, pad)) if pad else x
    return padded.view(x.shape[0], int(num_patches), int(patch_size))


def fit_patch_quantiles(x: np.ndarray, patch_size: int, vocab_size: int) -> np.ndarray:
    if x.ndim != 2:
        raise ValueError("x must be [cells, genes]")
    num_patches = int(math.ceil(x.shape[1] / int(patch_size)))
    pad = num_patches * int(patch_size) - x.shape[1]
    padded = np.pad(x, ((0, 0), (0, pad))) if pad else x
    patch_means = padded.reshape(x.shape[0], num_patches, int(patch_size)).mean(axis=-1).reshape(-1)
    qs = np.linspace(0.0, 1.0, int(vocab_size) + 1)[1:-1]
    edges = np.quantile(patch_means, qs).astype(np.float32)
    return np.unique(edges).astype(np.float32)


def quantize_patches(x: torch.Tensor, patch_size: int, num_patches: int, edges: torch.Tensor) -> torch.Tensor:
    patches = patchify_expression(x, patch_size, num_patches)
    patch_means = patches.mean(dim=-1)
    return torch.bucketize(patch_means.contiguous(), edges.to(device=x.device, dtype=x.dtype)).long()


def cosine_mask_ratio(ratio: torch.Tensor) -> torch.Tensor:
    return torch.cos(torch.clamp(ratio, 0.0, 1.0) * (math.pi / 2.0)).clamp(1e-6, 1.0)


def make_maskgit_inputs(
    target_ids: torch.Tensor,
    mask_token_id: int,
    min_mask_ratio: float,
    max_mask_ratio: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    if target_ids.ndim != 2:
        raise ValueError("target_ids must be [batch, patches]")
    if not 0.0 < min_mask_ratio <= max_mask_ratio < 1.0:
        raise ValueError("mask ratios must satisfy 0 < min <= max < 1")
    batch, num_patches = target_ids.shape
    ratio = torch.rand(batch, 1, device=target_ids.device)
    scheduled = cosine_mask_ratio(ratio)
    mask_ratio = float(min_mask_ratio) + (float(max_mask_ratio) - float(min_mask_ratio)) * scheduled
    mask = torch.zeros(batch, num_patches, device=target_ids.device)
    for row in range(batch):
        n_mask = max(1, min(num_patches, int(round(float(mask_ratio[row].item()) * num_patches))))
        chosen = torch.randperm(num_patches, device=target_ids.device)[:n_mask]
        mask[row, chosen] = 1.0
    input_ids = target_ids.clone()
    input_ids[mask.bool()] = int(mask_token_id)
    return input_ids, mask, mask_ratio.squeeze(1)


def masked_token_cross_entropy(logits: torch.Tensor, target_ids: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if logits.ndim != 3 or logits.shape[:2] != target_ids.shape or target_ids.shape != mask.shape:
        raise ValueError("logits [batch, patches, vocab], target_ids, and mask shapes are inconsistent")
    loss = F.cross_entropy(logits.transpose(1, 2), target_ids, reduction="none")
    denom = mask.sum().clamp_min(1.0)
    return (loss * mask).sum() / denom


def masked_patch_reconstruction_loss(
    patch_reconstruction: torch.Tensor,
    clean: torch.Tensor,
    mask: torch.Tensor,
    patch_size: int,
) -> torch.Tensor:
    patches = patchify_expression(clean, patch_size, patch_reconstruction.shape[1])
    if patch_reconstruction.shape != patches.shape or mask.shape != patch_reconstruction.shape[:2]:
        raise ValueError("patch reconstruction, clean patches, and mask shapes are inconsistent")
    weight = mask.to(dtype=patch_reconstruction.dtype).unsqueeze(-1)
    denom = weight.sum().mul(patch_reconstruction.shape[-1]).clamp_min(1.0)
    return (F.smooth_l1_loss(patch_reconstruction, patches, reduction="none") * weight).sum() / denom


def confidence_regularizer(logits: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    probs = torch.softmax(logits, dim=-1)
    entropy = -(probs * probs.clamp_min(1e-8).log()).sum(dim=-1)
    denom = mask.sum().clamp_min(1.0)
    return (entropy * mask).sum() / denom


def maskgit_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    target_ids: torch.Tensor,
    mask: torch.Tensor,
    *,
    patch_size: int,
    reconstruction_weight: float,
    confidence_weight: float,
) -> MaskGITLossParts:
    token = masked_token_cross_entropy(outputs["token_logits"], target_ids, mask)
    reconstruction = masked_patch_reconstruction_loss(outputs["patch_reconstruction"], clean, mask, patch_size)
    confidence = confidence_regularizer(outputs["token_logits"], mask)
    total = token + float(reconstruction_weight) * reconstruction + float(confidence_weight) * confidence
    return MaskGITLossParts(total, token, reconstruction, confidence, mask.mean())
