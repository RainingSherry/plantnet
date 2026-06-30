from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class FuzzyRoughLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    prototype: torch.Tensor
    lower_consistency: torch.Tensor
    boundary: torch.Tensor
    balance: torch.Tensor
    separation: torch.Tensor
    mask_rate: torch.Tensor
    relation_density: torch.Tensor
    mean_boundary_width: torch.Tensor
    mean_core_strength: torch.Tensor


def apply_mask_corruption(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    """Return corrupted expression and mask where 1 = corrupted/replaced/masked target."""
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


def rim_s_quantifier(p: torch.Tensor, alpha: float, beta: float) -> torch.Tensor:
    if not 0.0 <= float(alpha) < float(beta) <= 1.0:
        raise ValueError("RIM quantifier requires 0 <= alpha < beta <= 1")
    p = p.clamp(0.0, 1.0)
    alpha_t = torch.as_tensor(float(alpha), device=p.device, dtype=p.dtype)
    beta_t = torch.as_tensor(float(beta), device=p.device, dtype=p.dtype)
    midpoint = (alpha_t + beta_t) * 0.5
    denom = (beta_t - alpha_t).square().clamp_min(1e-8)
    rising = 2.0 * (p - alpha_t).square() / denom
    falling = 1.0 - 2.0 * (p - beta_t).square() / denom
    return torch.where(p <= alpha_t, torch.zeros_like(p), torch.where(p <= midpoint, rising, torch.where(p < beta_t, falling, torch.ones_like(p))))


def fuzzy_relation(embedding: torch.Tensor, relation_sigma: float) -> torch.Tensor:
    if embedding.ndim != 2:
        raise ValueError(f"embedding must be [cells, latent], got {tuple(embedding.shape)}")
    z = F.normalize(embedding, dim=1)
    dist_sq = torch.cdist(z, z).square()
    sigma = max(float(relation_sigma), 1e-4)
    relation = torch.exp(-dist_sq / (2.0 * sigma * sigma))
    eye = torch.eye(relation.shape[0], dtype=torch.bool, device=relation.device)
    return torch.where(eye, torch.ones_like(relation), relation)


def kleene_dienes_implication(a: torch.Tensor, b: torch.Tensor) -> torch.Tensor:
    return torch.maximum(1.0 - a, b)


def ywi_lower_approximation(
    relation: torch.Tensor,
    concept: torch.Tensor,
    quantifier_alpha: float,
    quantifier_beta: float,
) -> torch.Tensor:
    if relation.ndim != 2 or relation.shape[0] != relation.shape[1]:
        raise ValueError("relation must be [cells, cells]")
    if concept.ndim != 2 or concept.shape[0] != relation.shape[0]:
        raise ValueError("concept must be [cells, clusters] and align with relation")
    cells, clusters = concept.shape
    implication = kleene_dienes_implication(relation.unsqueeze(2), concept.unsqueeze(0).expand(cells, cells, clusters))
    implication_sorted = torch.sort(implication, dim=1, descending=True).values
    relation_sorted = torch.sort(relation, dim=1, descending=False).values
    cumulative = relation_sorted.cumsum(dim=1) / relation_sorted.sum(dim=1, keepdim=True).clamp_min(1e-8)
    previous = F.pad(cumulative[:, :-1], (1, 0), value=0.0)
    weights = rim_s_quantifier(cumulative, quantifier_alpha, quantifier_beta) - rim_s_quantifier(previous, quantifier_alpha, quantifier_beta)
    weights = weights / weights.sum(dim=1, keepdim=True).clamp_min(1e-8)
    lower = (implication_sorted * weights.unsqueeze(2)).sum(dim=1)
    if lower.shape != (cells, clusters):
        raise RuntimeError(f"lower approximation shape mismatch: {tuple(lower.shape)}")
    return lower.clamp(0.0, 1.0)


def unary_upper_approximation(
    relation: torch.Tensor,
    concept: torch.Tensor,
    quantifier_alpha: float,
    quantifier_beta: float,
) -> torch.Tensor:
    if relation.ndim != 2 or concept.ndim != 2 or relation.shape[0] != concept.shape[0]:
        raise ValueError("relation [cells, cells] and concept [cells, clusters] are incompatible")
    overlap = relation.unsqueeze(2) * concept.unsqueeze(0)
    proportion = overlap.mean(dim=1)
    upper = rim_s_quantifier(proportion, quantifier_alpha, quantifier_beta)
    return upper.clamp(0.0, 1.0)


def fuzzy_rough_approximations(
    embedding: torch.Tensor,
    concept: torch.Tensor,
    *,
    relation_sigma: float,
    lower_alpha: float,
    lower_beta: float,
    upper_alpha: float,
    upper_beta: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    relation = fuzzy_relation(embedding, relation_sigma)
    lower = ywi_lower_approximation(relation, concept, lower_alpha, lower_beta)
    upper = unary_upper_approximation(relation, concept, upper_alpha, upper_beta)
    upper = torch.maximum(upper, lower)
    return lower, upper, relation


def prototype_compactness(embedding: torch.Tensor, centers: torch.Tensor, membership: torch.Tensor, core_strength: torch.Tensor) -> torch.Tensor:
    if embedding.ndim != 2 or centers.ndim != 2 or membership.ndim != 2:
        raise ValueError("embedding, centers, and membership must be rank-2 tensors")
    if membership.shape != (embedding.shape[0], centers.shape[0]) or embedding.shape[1] != centers.shape[1]:
        raise ValueError("membership [cells, clusters] must align with embedding and centers")
    dist_sq = torch.sum((embedding.unsqueeze(1) - centers.unsqueeze(0)) ** 2, dim=2)
    weights = membership.detach().square() * core_strength.detach().clamp_min(0.05)
    return (weights * dist_sq).sum() / weights.sum().clamp_min(1e-8)


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


def fuzzy_rough_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    mask: torch.Tensor,
    *,
    reconstruction_weight: float,
    mask_weight: float,
    prototype_weight: float,
    lower_weight: float,
    boundary_weight: float,
    balance_weight: float,
    separation_weight: float,
    relation_sigma: float,
    lower_alpha: float,
    lower_beta: float,
    upper_alpha: float,
    upper_beta: float,
) -> FuzzyRoughLossParts:
    membership = outputs["membership"]
    boundary_membership = outputs["boundary_membership"]
    lower, upper, relation = fuzzy_rough_approximations(
        outputs["embedding"],
        membership.detach(),
        relation_sigma=relation_sigma,
        lower_alpha=lower_alpha,
        lower_beta=lower_beta,
        upper_alpha=upper_alpha,
        upper_beta=upper_beta,
    )
    boundary_width = (upper - lower).clamp_min(0.0)
    core_strength = lower.max(dim=1, keepdim=True).values
    target = (0.7 * lower + 0.3 * upper).detach()
    lower_consistency = (F.kl_div((boundary_membership + 1e-8).log(), target / target.sum(dim=1, keepdim=True).clamp_min(1e-8), reduction="batchmean"))
    boundary = (boundary_width * membership.detach()).sum() / membership.sum().clamp_min(1e-8)
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean, mask)
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
    prototype = prototype_compactness(outputs["embedding"], outputs["centers"], membership, core_strength)
    balance = partition_balance_loss(membership)
    separation = center_separation_loss(outputs["centers"])
    total = (
        float(reconstruction_weight) * reconstruction
        + float(mask_weight) * mask_loss
        + float(prototype_weight) * prototype
        + float(lower_weight) * lower_consistency
        + float(boundary_weight) * boundary
        + float(balance_weight) * balance
        + float(separation_weight) * separation
    )
    return FuzzyRoughLossParts(
        total,
        reconstruction,
        mask_loss,
        prototype,
        lower_consistency,
        boundary,
        balance,
        separation,
        mask.mean(),
        relation.detach().mean(),
        boundary_width.detach().mean(),
        core_strength.detach().mean(),
    )
