from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
import torch.nn.functional as F


@dataclass
class JOAOLossParts:
    total: torch.Tensor
    contrastive: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor


def project_to_simplex(values: np.ndarray) -> np.ndarray:
    if values.ndim != 1:
        raise ValueError("values must be a vector")
    if values.size == 0:
        raise ValueError("values must be non-empty")
    u = np.sort(values)[::-1]
    cssv = np.cumsum(u) - 1.0
    ind = np.arange(1, values.size + 1)
    cond = u - cssv / ind > 0
    if not np.any(cond):
        return np.ones_like(values) / float(values.size)
    rho = ind[cond][-1]
    theta = cssv[cond][-1] / float(rho)
    projected = np.maximum(values - theta, 0.0)
    return projected / projected.sum()


def joao_update_probabilities(
    aug_prob: np.ndarray,
    aug_losses: np.ndarray,
    beta: float,
    gamma: float,
) -> np.ndarray:
    if aug_prob.shape != aug_losses.shape:
        raise ValueError("aug_prob and aug_losses must have the same shape")
    uniform = np.ones_like(aug_prob) / float(aug_prob.size)
    proposal = aug_prob + float(beta) * (aug_losses - float(gamma) * (aug_prob - uniform))
    return project_to_simplex(proposal)


def sample_augmentation_pair(aug_prob: np.ndarray, rng: np.random.Generator) -> tuple[int, int]:
    if aug_prob.ndim != 1:
        raise ValueError("aug_prob must be one-dimensional")
    first = int(rng.choice(aug_prob.size, p=aug_prob))
    second = int(rng.choice(aug_prob.size, p=aug_prob))
    return first, second


def apply_graph_augmentation(
    x: torch.Tensor,
    adj: torch.Tensor,
    aug_id: int,
    ratio: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    if x.ndim != 2 or adj.ndim != 2 or adj.shape[0] != adj.shape[1] or adj.shape[0] != x.shape[0]:
        raise ValueError("x must be [batch, genes] and adj must be [batch, batch]")
    if not 0.0 <= float(ratio) < 1.0:
        raise ValueError("ratio must be in [0, 1)")
    x_aug = x.clone()
    adj_aug = adj.clone()
    mask = torch.zeros_like(x)
    if aug_id == 0:
        edge_mask = (torch.rand_like(adj_aug) > float(ratio)).float()
        edge_mask = torch.triu(edge_mask, diagonal=1)
        edge_mask = edge_mask + edge_mask.t()
        edge_mask.fill_diagonal_(1.0)
        adj_aug = adj_aug * edge_mask
    elif aug_id == 1:
        mask = (torch.rand_like(x_aug) < float(ratio)).float()
        x_aug = x_aug * (1.0 - mask)
    elif aug_id == 2:
        mask = (torch.rand_like(x_aug) < float(ratio)).float()
        noise = torch.randn_like(x_aug) * x_aug.std(dim=0, keepdim=True).clamp_min(1e-3)
        x_aug = x_aug + mask * noise
    elif aug_id == 3:
        cell_mask = (torch.rand(x_aug.shape[0], device=x_aug.device) < float(ratio)).float()
        x_aug = x_aug * (1.0 - cell_mask.view(-1, 1))
        adj_aug = adj_aug * (1.0 - cell_mask.view(-1, 1)) * (1.0 - cell_mask.view(1, -1))
        adj_aug.fill_diagonal_(1.0)
        mask = cell_mask.view(-1, 1).expand_as(x_aug)
    elif aug_id == 4:
        keep = (torch.rand(x_aug.shape[0], device=x_aug.device) >= float(ratio)).float()
        if keep.sum() < 2:
            keep = torch.ones_like(keep)
        adj_aug = adj_aug * keep.view(-1, 1) * keep.view(1, -1)
        adj_aug.fill_diagonal_(1.0)
        x_aug = x_aug * keep.view(-1, 1)
        mask = (1.0 - keep).view(-1, 1).expand_as(x_aug)
    else:
        raise ValueError(f"unknown augmentation id {aug_id}")
    degree = adj_aug.sum(dim=1).clamp_min(1.0)
    adj_aug = degree.rsqrt().view(-1, 1) * (adj_aug > 0).float() * degree.rsqrt().view(1, -1)
    return x_aug, adj_aug, mask


def nt_xent_loss(z1: torch.Tensor, z2: torch.Tensor, temperature: float) -> torch.Tensor:
    if z1.shape != z2.shape or z1.ndim != 2:
        raise ValueError(f"z1 and z2 must share [batch, dim], got {tuple(z1.shape)} and {tuple(z2.shape)}")
    if z1.shape[0] < 2:
        raise ValueError("contrastive loss requires batch size at least 2")
    z1 = F.normalize(z1, dim=1)
    z2 = F.normalize(z2, dim=1)
    logits = z1 @ z2.t() / float(temperature)
    labels = torch.arange(z1.shape[0], device=z1.device)
    return 0.5 * (F.cross_entropy(logits, labels) + F.cross_entropy(logits.t(), labels))


def joao_scmae_loss(
    online_view1: dict[str, torch.Tensor],
    online_view2: dict[str, torch.Tensor],
    reconstruction_view: dict[str, torch.Tensor],
    clean: torch.Tensor,
    reconstruction_mask: torch.Tensor,
    *,
    temperature: float,
    reconstruction_weight: float,
    mask_weight: float,
) -> JOAOLossParts:
    if reconstruction_view["reconstruction"].shape != clean.shape:
        raise ValueError("reconstruction and clean expression must share [batch, genes] shape")
    if reconstruction_mask.shape != clean.shape:
        raise ValueError("reconstruction_mask must be [batch, genes]")
    denom = reconstruction_mask.sum().clamp_min(1.0)
    contrastive = nt_xent_loss(online_view1["projection"], online_view2["projection"], temperature)
    reconstruction = (F.smooth_l1_loss(reconstruction_view["reconstruction"], clean, reduction="none") * reconstruction_mask).sum() / denom
    mask = F.binary_cross_entropy_with_logits(reconstruction_view["mask_logits"], reconstruction_mask)
    total = contrastive + float(reconstruction_weight) * reconstruction + float(mask_weight) * mask
    return JOAOLossParts(total, contrastive, reconstruction, mask, reconstruction_mask.mean())
