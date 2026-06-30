from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class ScVGAELossParts:
    total: torch.Tensor
    zinb: torch.Tensor
    reconstruction: torch.Tensor
    kl: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor


def apply_mask_corruption(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    if x.ndim != 2:
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    mask = (torch.rand_like(x) < float(mask_ratio)).float()
    if x.shape[0] > 1:
        replacement = x[torch.randperm(x.shape[0], device=x.device)]
    else:
        replacement = torch.zeros_like(x)
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
        raise ValueError("x, mu, theta, and pi_logits must share [batch, genes] shape")
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


def masked_reconstruction_loss(reconstruction: torch.Tensor, clean: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if reconstruction.shape != clean.shape or mask.shape != clean.shape:
        raise ValueError("reconstruction, clean, and mask must share [batch, genes] dimensions")
    denom = mask.sum().clamp_min(1.0)
    return (F.smooth_l1_loss(reconstruction, clean, reduction="none") * mask).sum() / denom


def kl_divergence_standard_normal(mean: torch.Tensor, logvar: torch.Tensor) -> torch.Tensor:
    if mean.shape != logvar.shape:
        raise ValueError("mean and logvar must have the same shape")
    return (-0.5 * (1.0 + logvar - mean.square() - logvar.exp()).sum(dim=1)).mean()


def scvgae_loss(
    outputs: dict[str, torch.Tensor],
    clean_scaled: torch.Tensor,
    zinb_target: torch.Tensor,
    mask: torch.Tensor,
    *,
    zinb_weight: float,
    reconstruction_weight: float,
    kl_weight: float,
    mask_weight: float,
    ridge_lambda: float,
) -> ScVGAELossParts:
    zinb = zinb_negative_log_likelihood(
        zinb_target,
        outputs["zinb_mu"],
        outputs["zinb_theta"],
        outputs["zinb_pi_logits"],
        ridge_lambda,
    )
    reconstruction = masked_reconstruction_loss(outputs["reconstruction"], clean_scaled, mask)
    kl = kl_divergence_standard_normal(outputs["z_mean"], outputs["z_logvar"])
    mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
    total = (
        float(zinb_weight) * zinb
        + float(reconstruction_weight) * reconstruction
        + float(kl_weight) * kl
        + float(mask_weight) * mask_loss
    )
    return ScVGAELossParts(total, zinb, reconstruction, kl, mask_loss, mask.mean())
