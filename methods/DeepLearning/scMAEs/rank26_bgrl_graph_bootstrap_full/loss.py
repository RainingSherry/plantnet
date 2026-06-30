from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class BGRLLossParts:
    total: torch.Tensor
    bootstrap: torch.Tensor
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
    upper_mask = torch.rand(cells, cells, device=adjacency.device) >= float(drop_prob)
    upper_mask = torch.triu(upper_mask, diagonal=1)
    keep = (upper_mask | upper_mask.T).to(dtype=adjacency.dtype)
    return adjacency * keep


def apply_feature_mask(x: torch.Tensor, mask_ratio: float, column_drop_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
    if x.ndim != 2:
        raise ValueError(f"x must be [cells, genes], got {tuple(x.shape)}")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    if not 0.0 <= float(column_drop_prob) < 1.0:
        raise ValueError("column_drop_prob must be in [0, 1)")
    element_mask = (torch.rand_like(x) < float(mask_ratio)).float()
    if column_drop_prob > 0.0:
        column_mask = (torch.rand(x.shape[1], device=x.device) < float(column_drop_prob)).float()
        element_mask = torch.maximum(element_mask, column_mask.unsqueeze(0).expand_as(x))
    corrupted = x.masked_fill(element_mask.bool(), 0.0)
    return corrupted, element_mask


def bgrl_bootstrap_loss(prediction: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    if prediction.shape != target.shape or prediction.ndim != 2:
        raise ValueError("prediction and target must share [cells, latent] shape")
    prediction = F.normalize(prediction, dim=1)
    target = F.normalize(target.detach(), dim=1)
    return (2.0 - 2.0 * (prediction * target).sum(dim=1)).mean()


def symmetric_bgrl_loss(outputs: dict[str, torch.Tensor]) -> torch.Tensor:
    return 0.5 * (
        bgrl_bootstrap_loss(outputs["prediction1"], outputs["target_z2"])
        + bgrl_bootstrap_loss(outputs["prediction2"], outputs["target_z1"])
    )


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and mask must share [cells, genes] shape")
    denom = mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum() / denom


def bgrl_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    mask1: torch.Tensor,
    mask2: torch.Tensor,
    adjacency1: torch.Tensor,
    adjacency2: torch.Tensor,
    *,
    bootstrap_weight: float,
    reconstruction_weight: float,
    mask_weight: float,
) -> BGRLLossParts:
    bootstrap = symmetric_bgrl_loss(outputs)
    reconstruction = 0.5 * (
        masked_reconstruction_loss(outputs["reconstruction1"], clean, mask1)
        + masked_reconstruction_loss(outputs["reconstruction2"], clean, mask2)
    )
    mask_loss = 0.5 * (
        F.binary_cross_entropy_with_logits(outputs["mask_logits1"], mask1)
        + F.binary_cross_entropy_with_logits(outputs["mask_logits2"], mask2)
    )
    total = float(bootstrap_weight) * bootstrap + float(reconstruction_weight) * reconstruction + float(mask_weight) * mask_loss
    return BGRLLossParts(
        total,
        bootstrap,
        reconstruction,
        mask_loss,
        0.5 * (mask1.mean() + mask2.mean()),
        adjacency1.detach().mean(),
        adjacency2.detach().mean(),
    )
