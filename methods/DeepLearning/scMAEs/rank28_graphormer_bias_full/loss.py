from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class GraphormerLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    graph: torch.Tensor
    mask_rate: torch.Tensor
    edge_density: torch.Tensor


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


def apply_mask_corruption(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    if x.ndim != 2:
        raise ValueError(f"x must be [cells, genes], got {tuple(x.shape)}")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    mask = (torch.rand_like(x) < float(mask_ratio)).float()
    replacement = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else torch.zeros_like(x)
    corrupted = torch.where(mask.bool(), replacement, x)
    effective_mask = (corrupted != x).float()
    return corrupted, effective_mask


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and mask must share [cells, genes] shape")
    denom = mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum() / denom


def graph_reconstruction_loss(predicted: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    if predicted.shape != target.shape or predicted.ndim != 2:
        raise ValueError("predicted and target adjacency must share [cells, cells] shape")
    return F.binary_cross_entropy(predicted.clamp(1e-6, 1.0 - 1e-6), target)


def graphormer_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    mask: torch.Tensor,
    *,
    reconstruction_weight: float,
    mask_weight: float,
    graph_weight: float,
) -> GraphormerLossParts:
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean, mask)
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
    graph = graph_reconstruction_loss(outputs["adjacency_reconstruction"], outputs["adjacency"].detach())
    total = (
        float(reconstruction_weight) * reconstruction
        + float(mask_weight) * mask_loss
        + float(graph_weight) * graph
    )
    return GraphormerLossParts(total, reconstruction, mask_loss, graph, mask.mean(), outputs["adjacency"].detach().mean())
