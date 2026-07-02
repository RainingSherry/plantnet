from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


def instance_contrastive_loss(z1: torch.Tensor, z2: torch.Tensor) -> torch.Tensor:
    batch, steps = z1.shape[0], z1.shape[1]
    if batch <= 1:
        return z1.new_tensor(0.0)
    z = torch.cat([z1, z2], dim=0).transpose(0, 1)
    sim = torch.matmul(z, z.transpose(1, 2))
    logits = torch.tril(sim, diagonal=-1)[:, :, :-1] + torch.triu(sim, diagonal=1)[:, :, 1:]
    logits = -F.log_softmax(logits, dim=-1)
    idx = torch.arange(batch, device=z1.device)
    return (logits[:, idx, batch + idx - 1].mean() + logits[:, batch + idx, idx].mean()) * 0.5


def temporal_contrastive_loss(z1: torch.Tensor, z2: torch.Tensor) -> torch.Tensor:
    steps = z1.shape[1]
    if steps <= 1:
        return z1.new_tensor(0.0)
    z = torch.cat([z1, z2], dim=1)
    sim = torch.matmul(z, z.transpose(1, 2))
    logits = torch.tril(sim, diagonal=-1)[:, :, :-1] + torch.triu(sim, diagonal=1)[:, :, 1:]
    logits = -F.log_softmax(logits, dim=-1)
    idx = torch.arange(steps, device=z1.device)
    return (logits[:, idx, steps + idx - 1].mean() + logits[:, steps + idx, idx].mean()) * 0.5


def hierarchical_contrastive_loss(z1: torch.Tensor, z2: torch.Tensor, alpha: float = 0.5, temporal_unit: int = 0) -> torch.Tensor:
    z1 = F.normalize(z1, dim=-1)
    z2 = F.normalize(z2, dim=-1)
    loss = z1.new_tensor(0.0)
    depth = 0
    while z1.shape[1] > 1:
        if alpha > 0:
            loss = loss + alpha * instance_contrastive_loss(z1, z2)
        if depth >= temporal_unit and alpha < 1:
            loss = loss + (1.0 - alpha) * temporal_contrastive_loss(z1, z2)
        depth += 1
        z1 = F.max_pool1d(z1.transpose(1, 2), kernel_size=2).transpose(1, 2)
        z2 = F.max_pool1d(z2.transpose(1, 2), kernel_size=2).transpose(1, 2)
    if z1.shape[1] == 1 and alpha > 0:
        loss = loss + alpha * instance_contrastive_loss(z1, z2)
        depth += 1
    return loss / max(1, depth)


class TS2VecGeneModuleLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        module_weight: float = 0.08,
        contrast_weight: float = 0.04,
        variance_weight: float = 0.01,
        contrast_alpha: float = 0.5,
        temporal_unit: int = 0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.module_weight = float(module_weight)
        self.contrast_weight = float(contrast_weight)
        self.variance_weight = float(variance_weight)
        self.contrast_alpha = float(contrast_alpha)
        self.temporal_unit = int(temporal_unit)

    @staticmethod
    def variance_loss(z: torch.Tensor) -> torch.Tensor:
        std = torch.sqrt(z.var(dim=0) + 1e-4)
        return F.relu(0.5 - std).mean()

    def forward(
        self,
        out: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        module_target: torch.Tensor,
        view1: torch.Tensor,
        view2: torch.Tensor,
    ) -> tuple[torch.Tensor, dict[str, float]]:
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec = (weights * F.smooth_l1_loss(out["reconstruction"], target_expr, reduction="none")).mean()
        mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
        scmae = (1.0 - self.mask_weight) * rec + self.mask_weight * mask_loss
        module_loss = F.smooth_l1_loss(out["module_prediction"], module_target)
        contrast = hierarchical_contrastive_loss(view1, view2, self.contrast_alpha, self.temporal_unit)
        var_loss = self.variance_loss(out["latent"])
        total = scmae + self.module_weight * module_loss + self.contrast_weight * contrast + self.variance_weight * var_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "module_loss": float(module_loss.detach().cpu()),
            "contrast_loss": float(contrast.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "view_steps": float(view1.shape[1]),
            "module_repr_norm": float(out["module_repr"].norm(dim=-1).mean().detach().cpu()),
        }
