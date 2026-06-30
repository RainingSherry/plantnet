from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class Data2VecLossParts:
    total: torch.Tensor
    latent: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor


def make_patch_mask(batch_size: int, num_patches: int, mask_ratio: float, device: torch.device) -> torch.Tensor:
    if batch_size <= 0 or num_patches <= 0:
        raise ValueError("batch_size and num_patches must be positive")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    mask = (torch.rand(batch_size, num_patches, device=device) < float(mask_ratio)).float()
    empty = mask.sum(dim=1) == 0
    if empty.any():
        chosen = torch.randint(0, num_patches, (int(empty.sum()),), device=device)
        mask[empty, chosen] = 1.0
    return mask


def patch_mask_to_gene_mask(patch_mask: torch.Tensor, patch_size: int, num_genes: int) -> torch.Tensor:
    if patch_mask.ndim != 2:
        raise ValueError(f"patch_mask must be [batch, patches], got {tuple(patch_mask.shape)}")
    return patch_mask.repeat_interleave(int(patch_size), dim=1)[:, : int(num_genes)]


def masked_latent_smooth_l1(prediction: torch.Tensor, target: torch.Tensor, patch_mask: torch.Tensor, beta: float) -> torch.Tensor:
    if prediction.shape != target.shape or prediction.ndim != 3:
        raise ValueError("prediction and target must share [batch, patches, hidden] shape")
    if patch_mask.shape != prediction.shape[:2]:
        raise ValueError("patch_mask must be [batch, patches]")
    mask = patch_mask.to(dtype=prediction.dtype).unsqueeze(-1)
    denom = mask.sum().mul(prediction.shape[-1]).clamp_min(1.0)
    return (F.smooth_l1_loss(prediction, target.detach(), beta=float(beta), reduction="none") * mask).sum() / denom


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, gene_mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or gene_mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and gene_mask must share [batch, genes] shape")
    denom = gene_mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * gene_mask).sum() / denom


def data2vec_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    teacher_target: torch.Tensor,
    patch_mask: torch.Tensor,
    gene_mask: torch.Tensor,
    *,
    beta: float,
    reconstruction_weight: float,
    mask_weight: float,
) -> Data2VecLossParts:
    latent = masked_latent_smooth_l1(outputs["patch_prediction"], teacher_target, patch_mask, beta)
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean, gene_mask)
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], gene_mask)
    total = latent + float(reconstruction_weight) * reconstruction + float(mask_weight) * mask_loss
    return Data2VecLossParts(total, latent, reconstruction, mask_loss, gene_mask.mean())
