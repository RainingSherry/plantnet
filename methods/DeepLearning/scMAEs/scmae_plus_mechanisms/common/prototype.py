from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans


class PrototypeHead(nn.Module):
    def __init__(self, n_clusters: int, dim: int) -> None:
        super().__init__()
        self.n_clusters = int(n_clusters)
        self.dim = int(dim)
        centers = torch.empty(self.n_clusters, self.dim)
        nn.init.xavier_uniform_(centers)
        self.centers = nn.Parameter(centers)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        z_norm = F.normalize(z, dim=1)
        c_norm = F.normalize(self.centers, dim=1)
        return z_norm @ c_norm.t()

    @torch.no_grad()
    def set_centers(self, centers) -> None:
        centers_t = torch.as_tensor(centers, dtype=self.centers.dtype, device=self.centers.device)
        if tuple(centers_t.shape) != tuple(self.centers.shape):
            raise ValueError(f"Expected centers {tuple(self.centers.shape)}, got {tuple(centers_t.shape)}")
        self.centers.copy_(centers_t)

    @torch.no_grad()
    def normalize_centers(self) -> None:
        self.centers.copy_(F.normalize(self.centers, dim=1))


def kmeans_centers(embedding: np.ndarray, n_clusters: int, seed: int) -> np.ndarray:
    km = KMeans(n_clusters=int(n_clusters), n_init=20, random_state=int(seed))
    km.fit(np.nan_to_num(embedding, nan=0.0, posinf=0.0, neginf=0.0))
    return km.cluster_centers_.astype(np.float32)


def student_t_assignments(z: torch.Tensor, centers: torch.Tensor, alpha: float = 1.0) -> torch.Tensor:
    dist = torch.cdist(z, centers, p=2.0).pow(2)
    q = (1.0 + dist / float(alpha)).pow(-(float(alpha) + 1.0) / 2.0)
    q = q / q.sum(dim=1, keepdim=True).clamp_min(1e-12)
    return q


def dec_target_distribution(q: torch.Tensor) -> torch.Tensor:
    weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-12)
    return (weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-12)).detach()


def dec_kl_loss(
    z: torch.Tensor,
    centers: torch.Tensor,
    alpha: float,
    confidence_threshold: float,
) -> tuple[torch.Tensor, dict]:
    q = student_t_assignments(z, centers, alpha=alpha)
    p = dec_target_distribution(q)
    confidence = q.max(dim=1).values
    per_sample = (p * (p.clamp_min(1e-12).log() - q.clamp_min(1e-12).log())).sum(dim=1)
    keep = confidence > float(confidence_threshold)
    if torch.any(keep):
        loss = per_sample[keep].mean()
        used_fraction = keep.float().mean()
    else:
        loss = per_sample.mean() * 0.0
        used_fraction = torch.zeros((), device=z.device)
    stats = {
        "prototype_confidence_mean": float(confidence.detach().mean().cpu()),
        "prototype_used_fraction": float(used_fraction.detach().cpu()),
    }
    return loss, stats


@torch.no_grad()
def sinkhorn_assignments(logits: torch.Tensor, temperature: float, iterations: int) -> torch.Tensor:
    q = torch.exp((logits / float(temperature)).float()).t()
    q = q / q.sum().clamp_min(1e-12)
    n_clusters, batch_size = q.shape
    for _ in range(int(iterations)):
        q = q / q.sum(dim=1, keepdim=True).clamp_min(1e-12)
        q = q / float(n_clusters)
        q = q / q.sum(dim=0, keepdim=True).clamp_min(1e-12)
        q = q / float(batch_size)
    q = q * float(batch_size)
    return q.t().detach()


def swapped_assignment_loss(
    logits_a: torch.Tensor,
    logits_b: torch.Tensor,
    assign_a: torch.Tensor,
    assign_b: torch.Tensor,
    pred_temperature: float,
) -> torch.Tensor:
    logp_a = F.log_softmax(logits_a / float(pred_temperature), dim=1)
    logp_b = F.log_softmax(logits_b / float(pred_temperature), dim=1)
    loss_ab = -(assign_a.detach() * logp_b).sum(dim=1).mean()
    loss_ba = -(assign_b.detach() * logp_a).sum(dim=1).mean()
    return 0.5 * (loss_ab + loss_ba)


@torch.no_grad()
def confidence_from_embedding(embedding: np.ndarray, centers: torch.Tensor, device: torch.device) -> float:
    z = torch.as_tensor(embedding, dtype=torch.float32, device=device)
    q = student_t_assignments(z, centers.detach(), alpha=1.0)
    return float(q.max(dim=1).values.mean().detach().cpu())

