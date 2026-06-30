from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class MaskFeatLossParts:
    total: torch.Tensor
    feature: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor


def gene_patch_feature_target(patches: torch.Tensor) -> torch.Tensor:
    """Deterministic MaskFeat target for gene patches.

    The visual MaskFeat paper predicts HOG features. For 1D gene-expression
    patches, the analogous fixed target concatenates normalized expression,
    local first-difference gradients, and within-patch ranks. These targets are
    computed from clean patches and detached by the caller/loss.
    """
    if patches.ndim != 3:
        raise ValueError(f"patches must be [batch, patches, patch_size], got {tuple(patches.shape)}")
    mean = patches.mean(dim=-1, keepdim=True)
    std = patches.std(dim=-1, unbiased=False, keepdim=True).clamp_min(1e-4)
    normalized = (patches - mean) / std
    gradient = torch.zeros_like(patches)
    gradient[..., 1:] = patches[..., 1:] - patches[..., :-1]
    ranks = torch.argsort(torch.argsort(patches, dim=-1), dim=-1).to(dtype=patches.dtype)
    if patches.shape[-1] > 1:
        ranks = ranks / float(patches.shape[-1] - 1)
    return torch.cat([normalized, gradient, ranks], dim=-1).detach()


def maskfeat_loss(
    feature_pred: torch.Tensor,
    target_features: torch.Tensor,
    reconstruction_pred: torch.Tensor,
    clean_patches: torch.Tensor,
    mask_logits: torch.Tensor,
    patch_mask: torch.Tensor,
    feature_weight: float = 1.0,
    reconstruction_weight: float = 0.2,
    mask_weight: float = 0.1,
) -> MaskFeatLossParts:
    if feature_pred.shape != target_features.shape:
        raise ValueError(f"feature_pred shape {tuple(feature_pred.shape)} != target {tuple(target_features.shape)}")
    if reconstruction_pred.shape != clean_patches.shape:
        raise ValueError(
            f"reconstruction_pred shape {tuple(reconstruction_pred.shape)} != clean patches {tuple(clean_patches.shape)}"
        )
    if mask_logits.shape != patch_mask.shape:
        raise ValueError(f"mask_logits shape {tuple(mask_logits.shape)} != patch_mask {tuple(patch_mask.shape)}")
    patch_mask = patch_mask.to(dtype=feature_pred.dtype, device=feature_pred.device)
    denom = patch_mask.sum().clamp_min(1.0)

    feature_per_patch = F.mse_loss(feature_pred, target_features, reduction="none").mean(dim=-1)
    feature_loss = (feature_per_patch * patch_mask).sum() / denom
    recon_per_patch = F.smooth_l1_loss(reconstruction_pred, clean_patches, reduction="none").mean(dim=-1)
    recon_loss = (recon_per_patch * patch_mask).sum() / denom
    mask_loss = F.binary_cross_entropy_with_logits(mask_logits, patch_mask)
    total = (
        float(feature_weight) * feature_loss
        + float(reconstruction_weight) * recon_loss
        + float(mask_weight) * mask_loss
    )
    return MaskFeatLossParts(total, feature_loss, recon_loss, mask_loss, patch_mask.mean())

