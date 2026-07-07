from __future__ import annotations

import torch
import torch.nn.functional as F
from torch import nn


class FuzzyRoughLoss(nn.Module):
    """scMAE loss plus fuzzy-rough core/boundary constraints."""

    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
        core_weight: float = 0.04,
        balance_weight: float = 0.01,
        separation_weight: float = 0.002,
        variance_weight: float = 0.0,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)
        self.core_weight = float(core_weight)
        self.balance_weight = float(balance_weight)
        self.separation_weight = float(separation_weight)
        self.variance_weight = float(variance_weight)

    @staticmethod
    def target_distribution(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def forward(self, model, outputs: dict[str, torch.Tensor], target_expr: torch.Tensor, mask: torch.Tensor, fuzzy_active: bool, core_threshold: float) -> tuple[torch.Tensor, dict[str, float], torch.Tensor]:
        rec_raw = F.mse_loss(outputs["reconstruction"], target_expr, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        rec_loss = (weights * rec_raw).mean()
        mask_loss = F.binary_cross_entropy_with_logits(outputs["mask_logits"], mask)
        q = outputs["membership"].clamp_min(1e-8)
        confidence = q.max(dim=1).values.detach()
        core_gate = (confidence >= core_threshold).float()
        if fuzzy_active and core_gate.sum() > 0:
            p = self.target_distribution(q).detach()
            kl_each = (p * (p.clamp_min(1e-8).log() - q.log())).sum(dim=1)
            core_loss = (kl_each * core_gate).sum() / core_gate.sum().clamp_min(1.0)
        else:
            core_loss = outputs["embedding"].new_tensor(0.0)
        mean_q = q.mean(dim=0)
        uniform = torch.full_like(mean_q, 1.0 / mean_q.numel())
        balance_loss = F.kl_div(mean_q.log(), uniform, reduction="batchmean") if fuzzy_active else outputs["embedding"].new_tensor(0.0)
        centers = F.normalize(model.cluster_centers, dim=1)
        sim = centers @ centers.T
        offdiag = sim[~torch.eye(sim.shape[0], device=sim.device, dtype=torch.bool)]
        separation_loss = F.relu(offdiag - 0.2).mean() if fuzzy_active and offdiag.numel() else outputs["embedding"].new_tensor(0.0)
        std = torch.sqrt(outputs["embedding"].var(dim=0, unbiased=False) + 1e-4)
        variance_loss = F.relu(0.5 - std).mean()
        scmae_loss = (1.0 - self.mask_loss_weight) * rec_loss + self.mask_loss_weight * mask_loss
        total = scmae_loss + self.core_weight * core_loss + self.balance_weight * balance_loss + self.separation_weight * separation_loss + self.variance_weight * variance_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae_loss.detach().cpu()),
            "recon_loss": float(rec_loss.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "core_loss": float(core_loss.detach().cpu()),
            "balance_loss": float(balance_loss.detach().cpu()),
            "separation_loss": float(separation_loss.detach().cpu()),
            "variance_loss": float(variance_loss.detach().cpu()),
            "core_fraction": float(core_gate.mean().detach().cpu()),
            "membership_entropy": float((-(q * q.log()).sum(dim=1).mean() / torch.log(torch.tensor(float(q.shape[1]), device=q.device))).detach().cpu()),
        }, rec_raw.detach().mean(dim=0)
