from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class AudioMAELossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    mask_rate: torch.Tensor


def audiomae_reconstruction_loss(
    prediction: torch.Tensor,
    target_patches: torch.Tensor,
    mask: torch.Tensor,
    valid_patch_weights: torch.Tensor,
    norm_patch_loss: bool = False,
) -> torch.Tensor:
    if prediction.shape != target_patches.shape or prediction.ndim != 3:
        raise ValueError("prediction and target_patches must share [batch, patches, patch_dim] shape")
    if mask.shape != prediction.shape[:2]:
        raise ValueError(f"mask must be [batch, patches], got {tuple(mask.shape)}")
    if valid_patch_weights.shape != (1, prediction.shape[1], prediction.shape[2]):
        raise ValueError(
            "valid_patch_weights must be [1, patches, patch_dim], "
            f"got {tuple(valid_patch_weights.shape)} for prediction {tuple(prediction.shape)}"
        )
    target = target_patches
    valid = valid_patch_weights.to(dtype=prediction.dtype, device=prediction.device)
    if norm_patch_loss:
        denom = valid.sum(dim=2, keepdim=True).clamp_min(1.0)
        mean = (target * valid).sum(dim=2, keepdim=True) / denom
        var = (((target - mean) * valid) ** 2).sum(dim=2, keepdim=True) / denom
        target = (target - mean) / torch.sqrt(var + 1e-6)
    weighted_mask = mask.to(dtype=prediction.dtype).unsqueeze(-1) * valid
    denom = weighted_mask.sum().clamp_min(1.0)
    return (F.mse_loss(prediction, target, reduction="none") * weighted_mask).sum() / denom


def audiomae_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    model,
    norm_patch_loss: bool,
) -> AudioMAELossParts:
    target = model.patchify(clean)
    valid = model.valid_patch_weights(clean.device)
    recon = audiomae_reconstruction_loss(outputs["prediction"], target, outputs["mask"], valid, norm_patch_loss)
    return AudioMAELossParts(recon, recon, outputs["mask"].to(dtype=recon.dtype).mean())
