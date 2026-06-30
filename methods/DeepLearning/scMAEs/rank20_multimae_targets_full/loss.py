from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class MultiMAELossParts:
    total: torch.Tensor
    expr: torch.Tensor
    rank: torch.Tensor
    stat: torch.Tensor
    masked_fraction: torch.Tensor


def patchify_expression(x: torch.Tensor, num_genes: int, patch_size: int) -> torch.Tensor:
    if x.ndim != 2 or x.shape[1] != int(num_genes):
        raise ValueError(f"x must be [batch, {num_genes}], got {tuple(x.shape)}")
    num_patches = (int(num_genes) + int(patch_size) - 1) // int(patch_size)
    pad = num_patches * int(patch_size) - int(num_genes)
    padded = F.pad(x, (0, pad)) if pad else x
    return padded.view(x.shape[0], num_patches, int(patch_size))


def per_cell_rank_targets(x: torch.Tensor, num_genes: int, patch_size: int) -> torch.Tensor:
    if x.ndim != 2 or x.shape[1] != int(num_genes):
        raise ValueError(f"x must be [batch, {num_genes}], got {tuple(x.shape)}")
    order = torch.argsort(x, dim=1, stable=True)
    ranks = torch.argsort(order, dim=1, stable=True).to(dtype=x.dtype)
    denom = max(1, int(num_genes) - 1)
    ranks = ranks / float(denom)
    return patchify_expression(ranks, num_genes, patch_size)


def patch_stat_targets(x_nonnegative: torch.Tensor, num_genes: int, patch_size: int) -> torch.Tensor:
    patches = patchify_expression(x_nonnegative.clamp_min(0.0), num_genes, patch_size)
    mean = patches.mean(dim=2)
    std = patches.std(dim=2, unbiased=False)
    zero_fraction = (patches <= 1e-8).to(dtype=patches.dtype).mean(dim=2)
    return torch.stack([mean, std, zero_fraction], dim=2)


def build_multimae_tasks(
    x_scaled: torch.Tensor,
    x_nonnegative: torch.Tensor,
    num_genes: int,
    patch_size: int,
) -> dict[str, torch.Tensor]:
    if x_scaled.shape != x_nonnegative.shape:
        raise ValueError("x_scaled and x_nonnegative must have identical [batch, genes] shape")
    return {
        "expr": patchify_expression(x_scaled, num_genes, patch_size),
        "rank": per_cell_rank_targets(x_scaled, num_genes, patch_size),
        "stat": patch_stat_targets(x_nonnegative, num_genes, patch_size),
    }


def masked_smooth_l1(pred: torch.Tensor, target: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if pred.shape != target.shape or pred.ndim != 3:
        raise ValueError(f"pred and target must share [batch, patches, features], got {tuple(pred.shape)} and {tuple(target.shape)}")
    if mask.shape != pred.shape[:2]:
        raise ValueError(f"mask must be [batch, patches], got {tuple(mask.shape)} for prediction {tuple(pred.shape)}")
    mask_f = mask.to(dtype=pred.dtype).unsqueeze(-1)
    denom = mask_f.sum().mul(pred.shape[-1]).clamp_min(1.0)
    return (F.smooth_l1_loss(pred, target, reduction="none") * mask_f).sum() / denom


def multimae_targets_loss(
    outputs: dict[str, torch.Tensor | dict[str, torch.Tensor]],
    targets: dict[str, torch.Tensor],
    *,
    rank_weight: float,
    stat_weight: float,
) -> MultiMAELossParts:
    masks = outputs.get("task_masks")
    if not isinstance(masks, dict):
        raise ValueError("outputs must contain task_masks dict")
    required = ("expr", "rank", "stat")
    for task in required:
        if task not in outputs or task not in targets or task not in masks:
            raise ValueError(f"missing MultiMAE task {task!r} in outputs, targets, or masks")
    expr = masked_smooth_l1(outputs["expr"], targets["expr"], masks["expr"])
    rank = masked_smooth_l1(outputs["rank"], targets["rank"], masks["rank"])
    stat = masked_smooth_l1(outputs["stat"], targets["stat"], masks["stat"])
    total = expr + float(rank_weight) * rank + float(stat_weight) * stat
    mask_values = [masks[task].to(dtype=expr.dtype).mean() for task in required]
    return MultiMAELossParts(total, expr, rank, stat, torch.stack(mask_values).mean())
