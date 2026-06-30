from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class BEiTLossParts:
    total: torch.Tensor
    token: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor
    token_accuracy: torch.Tensor


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


def masked_token_cross_entropy(logits: torch.Tensor, targets: torch.Tensor, patch_mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    if logits.ndim != 3 or targets.ndim != 2 or patch_mask.ndim != 2:
        raise ValueError("logits must be [batch, patches, vocab], targets/mask [batch, patches]")
    if logits.shape[:2] != targets.shape or targets.shape != patch_mask.shape:
        raise ValueError("logits, targets, and patch_mask patch dimensions must match")
    active = patch_mask.bool()
    if active.sum() == 0:
        return logits.new_zeros(()), logits.new_zeros(())
    active_logits = logits[active]
    active_targets = targets[active].long()
    loss = F.cross_entropy(active_logits, active_targets)
    acc = (active_logits.argmax(dim=1) == active_targets).float().mean()
    return loss, acc


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, gene_mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or gene_mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and gene_mask must share [batch, genes] shape")
    denom = gene_mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * gene_mask).sum() / denom


def beit_expression_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    patch_targets: torch.Tensor,
    patch_mask: torch.Tensor,
    gene_mask: torch.Tensor,
    *,
    reconstruction_weight: float,
    mask_weight: float,
) -> BEiTLossParts:
    token_loss, token_accuracy = masked_token_cross_entropy(outputs["token_logits"], patch_targets, patch_mask)
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean, gene_mask)
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], gene_mask)
    total = token_loss + float(reconstruction_weight) * reconstruction + float(mask_weight) * mask_loss
    return BEiTLossParts(total, token_loss, reconstruction, mask_loss, gene_mask.mean(), token_accuracy)
