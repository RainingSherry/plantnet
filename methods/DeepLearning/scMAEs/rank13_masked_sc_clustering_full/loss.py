from __future__ import annotations

from dataclasses import dataclass

import torch
import torch.nn.functional as F


@dataclass
class MaskSCLossParts:
    total: torch.Tensor
    reconstruction: torch.Tensor
    token_variance: torch.Tensor
    mask_rate: torch.Tensor


def make_fixed_count_mask_indices(
    batch_size: int,
    num_patches: int,
    mask_ratio: float,
    device: torch.device,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    if batch_size <= 0 or num_patches <= 1:
        raise ValueError("batch_size must be positive and num_patches must be > 1")
    if not 0.0 < float(mask_ratio) < 1.0:
        raise ValueError("mask_ratio must be in (0, 1)")
    num_masked = int(round(num_patches * float(mask_ratio)))
    num_masked = max(1, min(num_patches - 1, num_masked))
    noise = torch.rand(batch_size, num_patches, device=device)
    order = torch.argsort(noise, dim=1)
    masked_indices = order[:, :num_masked]
    visible_indices = order[:, num_masked:]
    patch_mask = torch.zeros(batch_size, num_patches, device=device)
    patch_mask.scatter_(1, masked_indices, 1.0)
    return visible_indices, masked_indices, patch_mask


def gather_patch_targets(target_features: torch.Tensor, masked_indices: torch.Tensor) -> torch.Tensor:
    if target_features.ndim != 3 or masked_indices.ndim != 2:
        raise ValueError("target_features must be [batch, patches, dim] and masked_indices must be [batch, masked]")
    if target_features.shape[0] != masked_indices.shape[0]:
        raise ValueError("target_features and masked_indices batch sizes differ")
    gather = masked_indices.to(device=target_features.device).unsqueeze(-1).expand(-1, -1, target_features.shape[-1])
    return torch.gather(target_features, dim=1, index=gather)


def sequence_guided_reconstruction(prediction: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    if prediction.shape != target.shape:
        raise ValueError(f"prediction and target must share shape, got {tuple(prediction.shape)} vs {tuple(target.shape)}")
    denom = torch.tensor(prediction.numel(), dtype=prediction.dtype, device=prediction.device).clamp_min(1.0)
    return F.mse_loss(prediction, target, reduction="sum") / denom


def token_variance_regularizer(tokens: torch.Tensor) -> torch.Tensor:
    if tokens.ndim != 3:
        raise ValueError(f"tokens must be [batch, visible_patches, hidden], got {tuple(tokens.shape)}")
    if tokens.shape[0] < 2:
        return tokens.new_zeros(())
    std = torch.sqrt(tokens.flatten(0, 1).var(dim=0, unbiased=False) + 1e-4)
    return torch.relu(1.0 - std).mean()


def mask_sc_loss(
    outputs: dict[str, torch.Tensor],
    target_features: torch.Tensor,
    patch_mask: torch.Tensor,
    *,
    variance_weight: float,
) -> MaskSCLossParts:
    reconstruction = sequence_guided_reconstruction(outputs["masked_prediction"], target_features)
    token_variance = token_variance_regularizer(outputs["visible_tokens"])
    total = reconstruction + float(variance_weight) * token_variance
    return MaskSCLossParts(total, reconstruction, token_variance, patch_mask.mean())


def contrastive_sequence_loss(z1: torch.Tensor, z2: torch.Tensor, temperature: float) -> torch.Tensor:
    if z1.ndim != 2 or z2.ndim != 2 or z1.shape != z2.shape:
        raise ValueError("z1 and z2 must share [batch, dim] shape")
    logits = torch.matmul(F.normalize(z1, dim=1), F.normalize(z2, dim=1).T) / float(temperature)
    labels = torch.arange(z1.shape[0], device=z1.device)
    return 0.5 * (F.cross_entropy(logits, labels) + F.cross_entropy(logits.T, labels))
