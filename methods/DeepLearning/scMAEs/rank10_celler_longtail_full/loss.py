from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class CellerLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    ginf: torch.Tensor
    compactness: torch.Tensor
    hard_weight_mean: torch.Tensor
    mask_rate: torch.Tensor


def apply_mask_corruption(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    if x.ndim != 2:
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    mask = (torch.rand_like(x) < float(mask_ratio)).float()
    replacement = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else torch.zeros_like(x)
    corrupted = torch.where(mask.bool(), replacement, x)
    effective_mask = (corrupted != x).float()
    return corrupted, effective_mask


def gaussian_inflation_loss(
    logits: torch.Tensor,
    pseudo_labels: torch.Tensor,
    class_counts: torch.Tensor,
    sigma: float,
    mu: float = 0.0,
    eps: float = 1e-6,
) -> torch.Tensor:
    if logits.ndim != 2 or pseudo_labels.ndim != 1 or logits.shape[0] != pseudo_labels.shape[0]:
        raise ValueError("logits must be [batch, classes] and pseudo_labels must be [batch]")
    if class_counts.ndim != 1 or class_counts.shape[0] != logits.shape[1]:
        raise ValueError("class_counts must be [classes]")
    counts = class_counts.to(device=logits.device, dtype=logits.dtype).clamp_min(1.0)
    delta_base = torch.log(counts.max() + eps) - torch.log(counts + eps)
    noise = torch.randn_like(logits) * float(sigma) + float(mu)
    inflated_logits = logits + delta_base.view(1, -1) * noise
    return F.cross_entropy(inflated_logits, pseudo_labels.long())


def longtail_hard_weights(
    reconstruction: torch.Tensor,
    clean: torch.Tensor,
    mask: torch.Tensor,
    logits: torch.Tensor,
    pseudo_labels: torch.Tensor,
    class_counts: torch.Tensor,
    rare_strength: float,
    hard_strength: float,
) -> torch.Tensor:
    if reconstruction.shape != clean.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and mask must share [batch, genes] dimensions")
    denom = mask.sum(dim=1).clamp_min(1.0)
    per_sample_error = (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum(dim=1) / denom
    probs = F.softmax(logits.detach(), dim=1)
    uncertainty = 1.0 - probs.max(dim=1).values
    counts = class_counts.to(device=clean.device, dtype=clean.dtype).clamp_min(1.0)
    rare = torch.sqrt(counts.max() / counts[pseudo_labels.long()]).clamp(1.0, 4.0)
    hard = (per_sample_error.detach() / per_sample_error.detach().mean().clamp_min(1e-6)).clamp(0.5, 3.0)
    weights = 1.0 + float(rare_strength) * (rare - 1.0) + float(hard_strength) * (hard + uncertainty - 1.0)
    return weights.clamp(0.25, 6.0)


def weighted_masked_reconstruction(
    reconstruction: torch.Tensor,
    clean: torch.Tensor,
    mask: torch.Tensor,
    sample_weights: torch.Tensor,
) -> torch.Tensor:
    if sample_weights.ndim != 1 or sample_weights.shape[0] != clean.shape[0]:
        raise ValueError("sample_weights must be [batch]")
    weighted_mask = mask * sample_weights.view(-1, 1)
    denom = weighted_mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * weighted_mask).sum() / denom


def prototype_compactness(
    embedding: torch.Tensor,
    pseudo_labels: torch.Tensor,
    prototype_centers: torch.Tensor,
) -> torch.Tensor:
    if prototype_centers.ndim != 2 or prototype_centers.shape[1] != embedding.shape[1]:
        raise ValueError("prototype_centers must be [classes, hidden]")
    centers = prototype_centers.to(device=embedding.device, dtype=embedding.dtype)[pseudo_labels.long()]
    return F.mse_loss(F.normalize(embedding, dim=1), F.normalize(centers, dim=1))


def celler_longtail_loss(
    outputs: dict[str, torch.Tensor],
    clean: torch.Tensor,
    mask: torch.Tensor,
    pseudo_labels: torch.Tensor,
    class_counts: torch.Tensor,
    prototype_centers: torch.Tensor,
    *,
    ginf_sigma: float,
    rare_strength: float,
    hard_strength: float,
    mask_weight: float,
    ginf_weight: float,
    compactness_weight: float,
) -> CellerLossParts:
    weights = longtail_hard_weights(
        outputs["reconstruction"],
        clean,
        mask,
        outputs["prototype_logits"],
        pseudo_labels,
        class_counts,
        rare_strength,
        hard_strength,
    )
    reconstruction = weighted_masked_reconstruction(outputs["reconstruction"], clean, mask, weights)
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
    ginf = gaussian_inflation_loss(outputs["prototype_logits"], pseudo_labels, class_counts, ginf_sigma)
    compactness = prototype_compactness(outputs["embedding"], pseudo_labels, prototype_centers)
    total = (
        reconstruction
        + float(mask_weight) * mask_loss
        + float(ginf_weight) * ginf
        + float(compactness_weight) * compactness
    )
    return CellerLossParts(total, reconstruction, mask_loss, ginf, compactness, weights.mean(), mask.mean())
