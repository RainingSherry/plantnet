from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class GraphBarlowLossParts:
    total: torch.Tensor
    barlow: torch.Tensor
    on_diag: torch.Tensor
    off_diag: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor
    edge_density1: torch.Tensor
    edge_density2: torch.Tensor


def build_knn_adjacency(x: torch.Tensor, k: int) -> torch.Tensor:
    if x.ndim != 2:
        raise ValueError(f"x must be [cells, genes], got {tuple(x.shape)}")
    cells = x.shape[0]
    if cells <= 1:
        return x.new_zeros(cells, cells)
    k = min(int(k), cells - 1)
    if k <= 0:
        raise ValueError("k must be positive when cells > 1")
    dist = torch.cdist(x, x)
    dist = dist.masked_fill(torch.eye(cells, dtype=torch.bool, device=x.device), float("inf"))
    nn_idx = torch.topk(dist, k=k, largest=False, dim=1).indices
    directed = torch.zeros(cells, cells, dtype=x.dtype, device=x.device).scatter_(1, nn_idx, 1.0)
    symmetric = torch.maximum(directed, directed.T)
    return symmetric.masked_fill(torch.eye(cells, dtype=torch.bool, device=x.device), 0.0)


def drop_edges(adjacency: torch.Tensor, drop_prob: float) -> torch.Tensor:
    if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
        raise ValueError("adjacency must be [cells, cells]")
    if float(drop_prob) <= 0.0:
        return adjacency
    if not 0.0 < float(drop_prob) < 1.0:
        raise ValueError("drop_prob must be in [0, 1)")
    cells = adjacency.shape[0]
    upper = torch.rand(cells, cells, device=adjacency.device) >= float(drop_prob)
    upper = torch.triu(upper, diagonal=1)
    keep = (upper | upper.T).to(dtype=adjacency.dtype)
    return adjacency * keep


def apply_feature_mask(x: torch.Tensor, feature_drop_prob: float, element_mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    if x.ndim != 2:
        raise ValueError(f"x must be [cells, genes], got {tuple(x.shape)}")
    if not 0.0 <= float(feature_drop_prob) < 1.0:
        raise ValueError("feature_drop_prob must be in [0, 1)")
    if not 0.0 <= float(element_mask_ratio) < 1.0:
        raise ValueError("element_mask_ratio must be in [0, 1)")
    feature_mask = (torch.rand(x.shape[1], device=x.device) < float(feature_drop_prob)).float()
    element_mask = (torch.rand_like(x) < float(element_mask_ratio)).float() if element_mask_ratio > 0.0 else torch.zeros_like(x)
    mask = torch.maximum(element_mask, feature_mask.unsqueeze(0).expand_as(x))
    return x.masked_fill(mask.bool(), 0.0), mask


def off_diagonal(x: torch.Tensor) -> torch.Tensor:
    if x.ndim != 2 or x.shape[0] != x.shape[1]:
        raise ValueError("x must be a square matrix")
    n = x.shape[0]
    return x.flatten()[:-1].view(n - 1, n + 1)[:, 1:].flatten()


def barlow_twins_loss(z1: torch.Tensor, z2: torch.Tensor, lambda_offdiag: float | None, eps: float) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    if z1.shape != z2.shape or z1.ndim != 2:
        raise ValueError("z1 and z2 must share [cells, features] shape")
    batch_size, feature_dim = z1.shape
    if batch_size < 2:
        zero = z1.new_zeros(())
        return zero, zero, zero
    lam = (1.0 / float(feature_dim)) if lambda_offdiag is None else float(lambda_offdiag)
    z1_norm = (z1 - z1.mean(dim=0)) / (z1.std(dim=0, unbiased=False) + float(eps))
    z2_norm = (z2 - z2.mean(dim=0)) / (z2.std(dim=0, unbiased=False) + float(eps))
    c = (z1_norm.T @ z2_norm) / float(batch_size)
    on_diag = (c.diagonal() - 1.0).pow(2).sum()
    off_diag_loss = off_diagonal(c).pow(2).sum()
    return on_diag + lam * off_diag_loss, on_diag, off_diag_loss


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and mask must share [cells, genes] shape")
    denom = mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum() / denom


def graph_barlow_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    mask1: torch.Tensor,
    mask2: torch.Tensor,
    adjacency1: torch.Tensor,
    adjacency2: torch.Tensor,
    *,
    lambda_offdiag: float | None,
    eps: float,
    barlow_weight: float,
    reconstruction_weight: float,
    mask_weight: float,
) -> GraphBarlowLossParts:
    barlow, on_diag, off_diag_loss = barlow_twins_loss(outputs["projection1"], outputs["projection2"], lambda_offdiag, eps)
    reconstruction = 0.5 * (
        masked_reconstruction_loss(outputs["reconstruction1"], clean, mask1)
        + masked_reconstruction_loss(outputs["reconstruction2"], clean, mask2)
    )
    mask_loss = 0.5 * (
        F.binary_cross_entropy_with_logits(outputs["mask_logits1"], mask1)
        + F.binary_cross_entropy_with_logits(outputs["mask_logits2"], mask2)
    )
    total = float(barlow_weight) * barlow + float(reconstruction_weight) * reconstruction + float(mask_weight) * mask_loss
    return GraphBarlowLossParts(
        total,
        barlow,
        on_diag,
        off_diag_loss,
        reconstruction,
        mask_loss,
        0.5 * (mask1.mean() + mask2.mean()),
        adjacency1.detach().mean(),
        adjacency2.detach().mean(),
    )
