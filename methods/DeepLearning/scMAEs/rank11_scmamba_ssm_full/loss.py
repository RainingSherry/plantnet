from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class ScMambaLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    token_smoothness: torch.Tensor
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
    if patch_size <= 0 or num_genes <= 0:
        raise ValueError("patch_size and num_genes must be positive")
    return patch_mask.repeat_interleave(int(patch_size), dim=1)[:, : int(num_genes)]


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, gene_mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or gene_mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and gene_mask must share [batch, genes] shape")
    denom = gene_mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * gene_mask).sum() / denom


def token_smoothness_loss(tokens: torch.Tensor) -> torch.Tensor:
    if tokens.ndim != 3:
        raise ValueError(f"tokens must be [batch, patches, hidden], got {tuple(tokens.shape)}")
    if tokens.shape[1] < 2:
        return tokens.new_zeros(())
    return F.smooth_l1_loss(tokens[:, 1:], tokens[:, :-1].detach())


def scmamba_masked_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    gene_mask: torch.Tensor,
    *,
    mask_weight: float,
    smoothness_weight: float,
) -> ScMambaLossParts:
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean, gene_mask)
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], gene_mask)
    smoothness = token_smoothness_loss(outputs["tokens"])
    total = reconstruction + float(mask_weight) * mask_loss + float(smoothness_weight) * smoothness
    return ScMambaLossParts(total, reconstruction, mask_loss, smoothness, gene_mask.mean())
