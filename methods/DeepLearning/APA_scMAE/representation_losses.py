from __future__ import annotations

import torch
import torch.nn.functional as F


def invariance_loss(z_a: torch.Tensor, z_b: torch.Tensor) -> torch.Tensor:
    return F.mse_loss(z_a, z_b)


def variance_loss(z: torch.Tensor, *, margin: float = 1.0, eps: float = 1.0e-4) -> torch.Tensor:
    if z.shape[0] < 2:
        return z.sum() * 0.0
    std = torch.sqrt(z.var(dim=0, unbiased=False) + float(eps))
    return torch.relu(float(margin) - std).mean()


def covariance_loss(z: torch.Tensor) -> torch.Tensor:
    if z.shape[0] < 2 or z.shape[1] < 2:
        return z.sum() * 0.0
    centered = z - z.mean(dim=0, keepdim=True)
    cov = centered.T.matmul(centered) / float(z.shape[0] - 1)
    return off_diagonal(cov).pow(2).sum() / float(z.shape[1])


def off_diagonal(matrix: torch.Tensor) -> torch.Tensor:
    rows, cols = matrix.shape
    if rows != cols:
        raise ValueError(f"off_diagonal expects a square matrix, got {tuple(matrix.shape)}")
    return matrix.flatten()[:-1].view(rows - 1, rows + 1)[:, 1:].flatten()


def vicreg_losses(
    z_clean: torch.Tensor,
    z_masked: torch.Tensor,
    *,
    variance_margin: float = 1.0,
    eps: float = 1.0e-4,
) -> dict[str, torch.Tensor]:
    inv = invariance_loss(z_clean, z_masked)
    var = 0.5 * (
        variance_loss(z_clean, margin=variance_margin, eps=eps)
        + variance_loss(z_masked, margin=variance_margin, eps=eps)
    )
    cov = 0.5 * (covariance_loss(z_clean) + covariance_loss(z_masked))
    return {"loss_repr_invariance": inv, "loss_repr_variance": var, "loss_repr_covariance": cov}


def teacher_consistency_loss(z_student: torch.Tensor, z_teacher: torch.Tensor) -> torch.Tensor:
    return F.mse_loss(z_student, z_teacher.detach())


def soft_assignment(z: torch.Tensor, prototypes: torch.Tensor, *, temperature: float) -> torch.Tensor:
    if prototypes.ndim != 2:
        raise ValueError(f"prototypes must be 2D, got {tuple(prototypes.shape)}")
    logits = -torch.cdist(z, prototypes).pow(2) / max(float(temperature), 1.0e-8)
    return torch.softmax(logits, dim=1)


def prototype_kl_loss(q_clean: torch.Tensor, q_masked: torch.Tensor) -> torch.Tensor:
    return F.kl_div(torch.log(q_masked.clamp_min(1.0e-8)), q_clean.detach(), reduction="batchmean")


def balanced_assignment_loss(q: torch.Tensor) -> torch.Tensor:
    if q.shape[1] < 2:
        return q.sum() * 0.0
    mean_q = q.mean(dim=0)
    uniform = torch.full_like(mean_q, 1.0 / float(q.shape[1]))
    return F.mse_loss(mean_q, uniform)
