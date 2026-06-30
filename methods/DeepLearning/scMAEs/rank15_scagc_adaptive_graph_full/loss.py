from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class ScAGCLossParts:
    total: torch.Tensor
    graph: torch.Tensor
    zinb: torch.Tensor
    contrastive: torch.Tensor
    clustering: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor
    adaptive_density: torch.Tensor


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
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    mask = (torch.rand_like(x) < float(mask_ratio)).float()
    replacement = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else torch.zeros_like(x)
    corrupted = torch.where(mask.bool(), replacement, x)
    effective_mask = (corrupted != x).float()
    return corrupted, effective_mask


def zinb_negative_log_likelihood(
    x: torch.Tensor,
    mu: torch.Tensor,
    theta: torch.Tensor,
    pi_logits: torch.Tensor,
    ridge_lambda: float = 0.0,
) -> torch.Tensor:
    if x.shape != mu.shape or x.shape != theta.shape or x.shape != pi_logits.shape:
        raise ValueError("x, mu, theta, and pi_logits must share [cells, genes] shape")
    eps = 1e-8
    x = x.clamp_min(0.0)
    mu = mu.clamp_min(eps)
    theta = theta.clamp_min(eps)
    softplus_pi = F.softplus(-pi_logits)
    log_theta_mu = torch.log(theta + mu + eps)
    nb_case = (
        torch.lgamma(theta + eps)
        + torch.lgamma(x + 1.0)
        - torch.lgamma(x + theta + eps)
        + (theta + x) * log_theta_mu
        - theta * torch.log(theta + eps)
        - x * torch.log(mu + eps)
    )
    zero_nb = torch.pow(theta / (theta + mu + eps), theta)
    zero_case = -torch.log(torch.sigmoid(pi_logits) + torch.sigmoid(-pi_logits) * zero_nb + eps)
    result = torch.where(x < eps, zero_case, softplus_pi + nb_case)
    if ridge_lambda > 0.0:
        result = result + float(ridge_lambda) * torch.sigmoid(pi_logits).square()
    return result.mean()


def graph_reconstruction_loss(predicted_adjacency: torch.Tensor, target_adjacency: torch.Tensor) -> torch.Tensor:
    if predicted_adjacency.shape != target_adjacency.shape or predicted_adjacency.ndim != 2:
        raise ValueError("predicted_adjacency and target_adjacency must share [cells, cells] shape")
    return F.mse_loss(predicted_adjacency, target_adjacency)


def graph_evolution_contrastive(z_previous: torch.Tensor, z_current: torch.Tensor, temperature: float) -> torch.Tensor:
    if z_previous.shape != z_current.shape or z_previous.ndim != 2:
        raise ValueError("z_previous and z_current must share [cells, hidden] shape")
    cells = z_previous.shape[0]
    if cells < 2:
        return z_previous.new_zeros(())
    features = torch.cat([F.normalize(z_previous, dim=1), F.normalize(z_current, dim=1)], dim=0)
    logits = torch.matmul(features, features.T) / float(temperature)
    logits = logits.masked_fill(torch.eye(2 * cells, dtype=torch.bool, device=logits.device), -1e9)
    labels = torch.cat([torch.arange(cells, 2 * cells, device=logits.device), torch.arange(0, cells, device=logits.device)])
    return F.cross_entropy(logits, labels)


def student_t_distribution(embedding: torch.Tensor, centers: torch.Tensor) -> torch.Tensor:
    if embedding.ndim != 2 or centers.ndim != 2 or embedding.shape[1] != centers.shape[1]:
        raise ValueError("embedding must be [cells, hidden] and centers [clusters, hidden]")
    dist_sq = torch.sum((embedding.unsqueeze(1) - centers.unsqueeze(0)) ** 2, dim=2)
    numerator = (1.0 + dist_sq).pow(-1.0)
    return numerator / numerator.sum(dim=1, keepdim=True).clamp_min(1e-8)


def target_distribution(q: torch.Tensor) -> torch.Tensor:
    if q.ndim != 2:
        raise ValueError(f"q must be [cells, clusters], got {tuple(q.shape)}")
    weight = q.square() / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
    return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)


def clustering_kl_loss(embedding: torch.Tensor, centers: torch.Tensor) -> torch.Tensor:
    q = student_t_distribution(embedding, centers)
    p = target_distribution(q).detach()
    return F.kl_div((q + 1e-8).log(), p, reduction="batchmean")


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and mask must share [cells, genes] shape")
    denom = mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum() / denom


def scagc_loss(
    outputs: dict[str, torch.Tensor],
    clean_scaled: torch.Tensor,
    zinb_target: torch.Tensor,
    mask: torch.Tensor,
    centers: torch.Tensor,
    *,
    graph_weight: float,
    zinb_weight: float,
    contrastive_weight: float,
    clustering_weight: float,
    reconstruction_weight: float,
    mask_weight: float,
    ridge_lambda: float,
    contrastive_temperature: float,
) -> ScAGCLossParts:
    graph = graph_reconstruction_loss(outputs["adjacency_reconstruction"], outputs["initial_adjacency"].detach())
    zinb = zinb_negative_log_likelihood(
        zinb_target,
        outputs["zinb_mu"],
        outputs["zinb_theta"],
        outputs["zinb_pi_logits"],
        ridge_lambda,
    )
    contrastive = graph_evolution_contrastive(outputs["previous_embedding"], outputs["embedding"], contrastive_temperature)
    clustering = clustering_kl_loss(outputs["embedding"], centers.to(device=clean_scaled.device, dtype=clean_scaled.dtype))
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean_scaled, mask)
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
    total = (
        float(graph_weight) * graph
        + float(zinb_weight) * zinb
        + float(contrastive_weight) * contrastive
        + float(clustering_weight) * clustering
        + float(reconstruction_weight) * reconstruction
        + float(mask_weight) * mask_loss
    )
    return ScAGCLossParts(
        total,
        graph,
        zinb,
        contrastive,
        clustering,
        reconstruction,
        mask_loss,
        mask.mean(),
        outputs["adaptive_adjacency"].detach().mean(),
    )
