from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class ConsistencyLossParts:
    total: torch.Tensor
    consistency: torch.Tensor
    reconstruction: torch.Tensor
    mask: torch.Tensor
    mask_rate: torch.Tensor
    sigma_mean: torch.Tensor


def update_ema(target: torch.nn.Module, source: torch.nn.Module, rate: float) -> None:
    if not 0.0 <= float(rate) <= 1.0:
        raise ValueError("EMA rate must be in [0, 1]")
    with torch.no_grad():
        for target_param, source_param in zip(target.parameters(), source.parameters()):
            target_param.mul_(float(rate)).add_(source_param, alpha=1.0 - float(rate))
        for target_buffer, source_buffer in zip(target.buffers(), source.buffers()):
            target_buffer.copy_(source_buffer)


def karras_adjacent_sigmas(
    batch_size: int,
    num_scales: int,
    sigma_min: float,
    sigma_max: float,
    rho: float,
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor]:
    if num_scales < 2:
        raise ValueError("num_scales must be at least 2")
    indices = torch.randint(0, int(num_scales) - 1, (int(batch_size),), device=device)
    ramp = indices.float() / float(num_scales - 1)
    ramp_next = (indices.float() + 1.0) / float(num_scales - 1)
    min_inv = float(sigma_min) ** (1.0 / float(rho))
    max_inv = float(sigma_max) ** (1.0 / float(rho))
    sigmas = (max_inv + ramp * (min_inv - max_inv)) ** float(rho)
    next_sigmas = (max_inv + ramp_next * (min_inv - max_inv)) ** float(rho)
    return sigmas, next_sigmas


def masked_gaussian_corruption(
    clean: torch.Tensor,
    sigma: torch.Tensor,
    mask_ratio: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    if clean.ndim != 2:
        raise ValueError(f"clean must be [batch, genes], got {tuple(clean.shape)}")
    if sigma.ndim != 1 or sigma.shape[0] != clean.shape[0]:
        raise ValueError(f"sigma must be [batch], got {tuple(sigma.shape)} for clean {tuple(clean.shape)}")
    if not 0.0 < float(mask_ratio) <= 1.0:
        raise ValueError("mask_ratio must be in (0, 1]")
    mask = (torch.rand_like(clean) < float(mask_ratio)).float()
    noise = torch.randn_like(clean)
    corrupted = clean + mask * noise * sigma.view(-1, 1)
    return corrupted, mask, noise


def euler_step_to_next_sigma(
    x_t: torch.Tensor,
    sigma: torch.Tensor,
    next_sigma: torch.Tensor,
    clean_target: torch.Tensor,
    mask: torch.Tensor,
) -> torch.Tensor:
    if x_t.shape != clean_target.shape or mask.shape != x_t.shape:
        raise ValueError("x_t, clean_target, and mask must share [batch, genes] dimensions")
    direction = (x_t - clean_target) / sigma.view(-1, 1).clamp_min(1e-12)
    return x_t + mask * direction * (next_sigma - sigma).view(-1, 1)


def pseudo_huber(diff: torch.Tensor, delta: float) -> torch.Tensor:
    if delta <= 0:
        raise ValueError("delta must be positive")
    scale = float(delta)
    return scale * scale * (torch.sqrt(1.0 + (diff / scale).square()) - 1.0)


def masked_mean(values: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
    if values.shape != mask.shape:
        raise ValueError(f"values and mask must have the same shape, got {tuple(values.shape)} and {tuple(mask.shape)}")
    denom = mask.sum().clamp_min(1.0)
    return (values * mask).sum() / denom


def consistency_training_loss(
    model: torch.nn.Module,
    target_model: torch.nn.Module,
    clean: torch.Tensor,
    *,
    num_scales: int,
    sigma_min: float,
    sigma_max: float,
    rho: float,
    mask_ratio: float,
    pseudo_huber_delta: float,
    consistency_weight: float,
    reconstruction_weight: float,
    mask_weight: float,
) -> ConsistencyLossParts:
    sigma, next_sigma = karras_adjacent_sigmas(
        clean.shape[0],
        num_scales,
        sigma_min,
        sigma_max,
        rho,
        clean.device,
    )
    x_t, mask, _ = masked_gaussian_corruption(clean, sigma, mask_ratio)
    online = model(x_t, sigma)
    with torch.no_grad():
        x_next = euler_step_to_next_sigma(x_t, sigma, next_sigma, clean, mask)
        target = target_model(x_next, next_sigma)["denoised"]

    consistency = masked_mean(pseudo_huber(online["denoised"] - target, pseudo_huber_delta), mask)
    reconstruction = masked_mean(pseudo_huber(online["denoised"] - clean, pseudo_huber_delta), mask)
    mask_loss = F.binary_cross_entropy_with_logits(online["mask_logits"], mask)
    total = (
        float(consistency_weight) * consistency
        + float(reconstruction_weight) * reconstruction
        + float(mask_weight) * mask_loss
    )
    return ConsistencyLossParts(
        total=total,
        consistency=consistency,
        reconstruction=reconstruction,
        mask=mask_loss,
        mask_rate=mask.mean(),
        sigma_mean=sigma.mean(),
    )
