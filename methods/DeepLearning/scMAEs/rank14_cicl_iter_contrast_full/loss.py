from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class CICLLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    instance_contrastive: torch.Tensor
    cluster_contrastive: torch.Tensor
    cluster_kl: torch.Tensor
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


def gaussian_augment(x: torch.Tensor, noise_std: float, feature_drop: float) -> torch.Tensor:
    if x.ndim != 2:
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    keep = (torch.rand_like(x) > float(feature_drop)).to(dtype=x.dtype)
    return x * keep + torch.randn_like(x) * float(noise_std)


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and mask must share [batch, genes] shape")
    denom = mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum() / denom


def student_t_distribution(embedding: torch.Tensor, centers: torch.Tensor, alpha: float) -> torch.Tensor:
    if embedding.ndim != 2 or centers.ndim != 2 or embedding.shape[1] != centers.shape[1]:
        raise ValueError("embedding must be [batch, hidden] and centers [clusters, hidden]")
    norm_squared = torch.sum((embedding.unsqueeze(1) - centers.unsqueeze(0)) ** 2, dim=2)
    numerator = (1.0 + norm_squared / float(alpha)).pow(-0.5 * (float(alpha) + 1.0))
    return numerator / numerator.sum(dim=1, keepdim=True).clamp_min(1e-8)


def target_distribution(q: torch.Tensor) -> torch.Tensor:
    if q.ndim != 2:
        raise ValueError(f"q must be [batch, clusters], got {tuple(q.shape)}")
    weight = (q ** 2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
    return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)


def instance_contrastive_loss(z1: torch.Tensor, z2: torch.Tensor, temperature: float) -> torch.Tensor:
    if z1.ndim != 2 or z2.ndim != 2 or z1.shape != z2.shape:
        raise ValueError("z1 and z2 must share [batch, dim] shape")
    features = torch.cat([F.normalize(z1, dim=1), F.normalize(z2, dim=1)], dim=0)
    batch = z1.shape[0]
    logits = torch.matmul(features, features.T) / float(temperature)
    logits = logits.masked_fill(torch.eye(2 * batch, dtype=torch.bool, device=logits.device), -1e9)
    positives = torch.cat([torch.arange(batch, 2 * batch, device=z1.device), torch.arange(0, batch, device=z1.device)])
    return F.cross_entropy(logits, positives)


def cluster_aware_contrastive_loss(z1: torch.Tensor, z2: torch.Tensor, labels: torch.Tensor, temperature: float) -> torch.Tensor:
    if z1.ndim != 2 or z2.ndim != 2 or z1.shape != z2.shape or labels.ndim != 1 or labels.shape[0] != z1.shape[0]:
        raise ValueError("z1/z2 must share [batch, dim] and labels must be [batch]")
    features = torch.cat([F.normalize(z1, dim=1), F.normalize(z2, dim=1)], dim=0)
    label_all = torch.cat([labels.long(), labels.long()], dim=0)
    logits = torch.matmul(features, features.T) / float(temperature)
    self_mask = torch.eye(logits.shape[0], dtype=torch.bool, device=logits.device)
    same_cluster = label_all.unsqueeze(0).eq(label_all.unsqueeze(1)) & ~self_mask
    logits_exp = torch.exp(logits.masked_fill(self_mask, -1e9))
    numerator = (logits_exp * same_cluster.to(dtype=logits_exp.dtype)).sum(dim=1).clamp_min(1e-8)
    denominator = logits_exp.sum(dim=1).clamp_min(1e-8)
    return (-torch.log(numerator / denominator)).mean()


def cicl_loss(
    outputs_clean: dict[str, torch.Tensor],
    outputs_aug1: dict[str, torch.Tensor],
    outputs_aug2: dict[str, torch.Tensor],
    clean: torch.Tensor,
    mask: torch.Tensor,
    centers: torch.Tensor,
    pseudo_labels: torch.Tensor,
    *,
    alpha: float,
    temperature: float,
    cluster_weight: float,
    reconstruction_weight: float,
    mask_weight: float,
    kl_weight: float,
) -> CICLLossParts:
    reconstruction = masked_reconstruction_loss(outputs_clean["reconstruction"], clean, mask)
    mask_loss = F.binary_cross_entropy_with_logits(outputs_clean["mask_logits"], mask)
    instance = instance_contrastive_loss(outputs_aug1["projection"], outputs_aug2["projection"], temperature)
    cluster = cluster_aware_contrastive_loss(outputs_aug1["projection"], outputs_aug2["projection"], pseudo_labels, temperature)
    q = student_t_distribution(outputs_clean["embedding"], centers.to(device=clean.device, dtype=clean.dtype), alpha)
    p = target_distribution(q).detach()
    cluster_kl = F.kl_div((q + 1e-8).log(), p, reduction="batchmean")
    total = (
        float(reconstruction_weight) * reconstruction
        + float(mask_weight) * mask_loss
        + instance
        + float(cluster_weight) * cluster
        + float(kl_weight) * cluster_kl
    )
    return CICLLossParts(total, reconstruction, mask_loss, instance, cluster, cluster_kl, mask.mean())
