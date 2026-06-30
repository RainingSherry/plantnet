from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class IjepaLossParts:
    total: torch.Tensor
    latent: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    target_patch_rate: torch.Tensor


def make_context_target_masks(
    batch_size: int,
    num_patches: int,
    target_ratio: float,
    context_ratio: float,
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor]:
    if batch_size <= 0 or num_patches <= 1:
        raise ValueError("batch_size must be positive and num_patches must be > 1")
    if not 0.0 < float(target_ratio) < 1.0:
        raise ValueError("target_ratio must be in (0, 1)")
    if not 0.0 < float(context_ratio) <= 1.0:
        raise ValueError("context_ratio must be in (0, 1]")
    target_len = max(1, min(num_patches - 1, int(round(num_patches * float(target_ratio)))))
    context_len = max(1, min(num_patches - target_len, int(round(num_patches * float(context_ratio)))))
    target_mask = torch.zeros(batch_size, num_patches, device=device)
    context_mask = torch.zeros(batch_size, num_patches, device=device)
    for row in range(batch_size):
        start = int(torch.randint(0, num_patches - target_len + 1, (1,), device=device).item())
        target_mask[row, start : start + target_len] = 1.0
        available = torch.nonzero(target_mask[row] == 0.0, as_tuple=False).flatten()
        perm = available[torch.randperm(available.numel(), device=device)[:context_len]]
        context_mask[row, perm] = 1.0
    return context_mask, target_mask


def patch_mask_to_gene_mask(patch_mask: torch.Tensor, patch_size: int, num_genes: int) -> torch.Tensor:
    if patch_mask.ndim != 2:
        raise ValueError(f"patch_mask must be [batch, patches], got {tuple(patch_mask.shape)}")
    return patch_mask.repeat_interleave(int(patch_size), dim=1)[:, : int(num_genes)]


def target_block_smooth_l1(prediction: torch.Tensor, target: torch.Tensor, target_mask: torch.Tensor, beta: float) -> torch.Tensor:
    if prediction.shape != target.shape or prediction.ndim != 3:
        raise ValueError("prediction and target must share [batch, patches, hidden] shape")
    if target_mask.shape != prediction.shape[:2]:
        raise ValueError("target_mask must be [batch, patches]")
    mask = target_mask.to(dtype=prediction.dtype).unsqueeze(-1)
    denom = mask.sum().mul(prediction.shape[-1]).clamp_min(1.0)
    return (F.smooth_l1_loss(prediction, target.detach(), beta=float(beta), reduction="none") * mask).sum() / denom


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, gene_mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or gene_mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and gene_mask must share [batch, genes] shape")
    denom = gene_mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * gene_mask).sum() / denom


def ijepa_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    target_features: torch.Tensor,
    target_mask: torch.Tensor,
    gene_target_mask: torch.Tensor,
    *,
    beta: float,
    reconstruction_weight: float,
    mask_weight: float,
) -> IjepaLossParts:
    latent = target_block_smooth_l1(outputs["predicted_targets"], target_features, target_mask, beta)
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean, gene_target_mask)
    mask = F.binary_cross_entropy_with_logits(outputs["mask_logits"], gene_target_mask)
    total = latent + float(reconstruction_weight) * reconstruction + float(mask_weight) * mask
    return IjepaLossParts(total, latent, reconstruction, mask, target_mask.mean())
