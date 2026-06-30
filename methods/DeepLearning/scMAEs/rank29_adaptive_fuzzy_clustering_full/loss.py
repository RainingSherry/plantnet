from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class AdaptiveFuzzyLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    fuzzy_reconstruction: torch.Tensor
    fuzzy_compactness: torch.Tensor
    adaptive_entropy: torch.Tensor
    partition_balance: torch.Tensor
    center_separation: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor
    membership_entropy: torch.Tensor
    fuzzifier: torch.Tensor


def apply_mask_corruption(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    """Return corrupted expression and mask where 1 = replaced/masked target."""
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


def fuzzy_compactness_loss(weighted_membership: torch.Tensor, distance_sq: torch.Tensor) -> torch.Tensor:
    if weighted_membership.shape != distance_sq.shape or weighted_membership.ndim != 2:
        raise ValueError("weighted_membership and distance_sq must share [cells, clusters] shape")
    denom = weighted_membership.sum().clamp_min(1e-8)
    return (weighted_membership * distance_sq).sum() / denom


def weighted_adaptive_entropy(membership: torch.Tensor, distance_sq: torch.Tensor, lambda1: float, lambda2: float) -> torch.Tensor:
    if membership.shape != distance_sq.shape or membership.ndim != 2:
        raise ValueError("membership and distance_sq must share [cells, clusters] shape")
    adaptive_weight = torch.softmax(-distance_sq.detach(), dim=1).clamp_min(1e-8)
    membership_sum_penalty = (1.0 - membership.sum(dim=1)).square().mean()
    entropy = -(adaptive_weight * adaptive_weight.log()).sum(dim=1).mean()
    return float(lambda1) * membership_sum_penalty + float(lambda2) * entropy


def partition_balance_loss(membership: torch.Tensor) -> torch.Tensor:
    if membership.ndim != 2:
        raise ValueError(f"membership must be [cells, clusters], got {tuple(membership.shape)}")
    cluster_mass = membership.mean(dim=0)
    uniform = torch.full_like(cluster_mass, 1.0 / membership.shape[1])
    return F.mse_loss(cluster_mass, uniform)


def center_separation_loss(centers: torch.Tensor) -> torch.Tensor:
    if centers.ndim != 2:
        raise ValueError(f"centers must be [clusters, latent], got {tuple(centers.shape)}")
    clusters = centers.shape[0]
    if clusters <= 1:
        return centers.new_zeros(())
    dist = torch.cdist(centers, centers).square()
    off_diag = ~torch.eye(clusters, dtype=torch.bool, device=centers.device)
    return torch.exp(-dist[off_diag]).mean()


def membership_entropy(membership: torch.Tensor) -> torch.Tensor:
    if membership.ndim != 2:
        raise ValueError(f"membership must be [cells, clusters], got {tuple(membership.shape)}")
    return -(membership.clamp_min(1e-8) * membership.clamp_min(1e-8).log()).sum(dim=1).mean()


def adaptive_fuzzy_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    mask: torch.Tensor,
    *,
    reconstruction_weight: float,
    fuzzy_reconstruction_weight: float,
    fuzzy_weight: float,
    entropy_weight: float,
    balance_weight: float,
    separation_weight: float,
    mask_weight: float,
    entropy_lambda1: float,
    entropy_lambda2: float,
) -> AdaptiveFuzzyLossParts:
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean, mask)
    fuzzy_reconstruction = masked_reconstruction_loss(outputs["fuzzy_reconstruction"], clean, mask)
    compactness = fuzzy_compactness_loss(outputs["weighted_membership"], outputs["distance_sq"])
    adaptive_entropy = weighted_adaptive_entropy(
        outputs["membership"],
        outputs["distance_sq"],
        entropy_lambda1,
        entropy_lambda2,
    )
    balance = partition_balance_loss(outputs["membership"])
    separation = center_separation_loss(outputs["centers"])
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
    total = (
        float(reconstruction_weight) * reconstruction
        + float(fuzzy_reconstruction_weight) * fuzzy_reconstruction
        + float(fuzzy_weight) * compactness
        + float(entropy_weight) * adaptive_entropy
        + float(balance_weight) * balance
        + float(separation_weight) * separation
        + float(mask_weight) * mask_loss
    )
    return AdaptiveFuzzyLossParts(
        total,
        reconstruction,
        fuzzy_reconstruction,
        compactness,
        adaptive_entropy,
        balance,
        separation,
        mask_loss,
        mask.mean(),
        membership_entropy(outputs["membership"]),
        outputs["fuzzifier"].detach(),
    )

