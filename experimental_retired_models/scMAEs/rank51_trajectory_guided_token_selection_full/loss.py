from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class TrajectoryGuidedLoss(nn.Module):
    """scMAE reconstruction/mask loss plus a light sampler regularizer."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        trajectory_weight: float = 0.04,
        sampler_target_weight: float = 0.01,
        sampler_entropy_weight: float = 0.002,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.trajectory_weight = float(trajectory_weight)
        self.sampler_target_weight = float(sampler_target_weight)
        self.sampler_entropy_weight = float(sampler_entropy_weight)
        self.variance_weight = float(variance_weight)

    def forward(
        self,
        outputs: dict[str, torch.Tensor],
        target_expr: torch.Tensor,
        mask: torch.Tensor,
        trajectory_target: torch.Tensor,
        sampler_logits: torch.Tensor | None,
    ) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.smooth_l1_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        traj_target = trajectory_target.detach().view(1, -1).expand_as(outputs["trajectory_pred"])
        trajectory_loss = F.smooth_l1_loss(outputs["trajectory_pred"], traj_target)
        if sampler_logits is None:
            entropy_loss = outputs["embedding"].new_tensor(0.0)
            sampler_target_loss = outputs["embedding"].new_tensor(0.0)
        else:
            probs = torch.softmax(sampler_logits, dim=0).clamp_min(1e-8)
            entropy = -(probs * probs.log()).sum() / torch.log(torch.tensor(float(probs.numel()), device=probs.device).clamp_min(2.0))
            entropy_loss = -entropy
            target = trajectory_target.detach()
            target = (target - target.min()) / (target.max() - target.min()).clamp_min(1e-6)
            sampler_target_loss = F.binary_cross_entropy_with_logits(sampler_logits, target)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = (
            scmae_loss
            + self.trajectory_weight * trajectory_loss
            + self.sampler_target_weight * sampler_target_loss
            + self.sampler_entropy_weight * entropy_loss
            + self.variance_weight * variance_loss
        )
        parts = {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "trajectory_loss": float(trajectory_loss.detach().cpu()),
            "sampler_target_loss": float(sampler_target_loss.detach().cpu()),
            "sampler_entropy_loss": float(entropy_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
        }
        return total, parts, rec_raw.detach().mean(dim=0)
