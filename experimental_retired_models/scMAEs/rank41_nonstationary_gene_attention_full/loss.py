from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class NonStationaryGeneAttentionLoss(nn.Module):
    def __init__(
        self,
        masked_data_weight: float = 0.75,
        mask_weight: float = 0.60,
        module_weight: float = 0.08,
        consistency_weight: float = 0.04,
        variance_weight: float = 0.01,
        factor_weight: float = 0.001,
    ):
        super().__init__()
        self.masked_data_weight = float(masked_data_weight)
        self.mask_weight = float(mask_weight)
        self.module_weight = float(module_weight)
        self.consistency_weight = float(consistency_weight)
        self.variance_weight = float(variance_weight)
        self.factor_weight = float(factor_weight)

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
        consistency = 0.5 * (F.smooth_l1_loss(view1, view2.detach()) + F.smooth_l1_loss(view2, view1.detach()))
        var_loss = self.variance_loss(out["latent"])
        factor_loss = (out["tau_mean"] - 1.0).pow(2) + 0.01 * out["delta_abs_mean"]
        total = scmae + self.module_weight * module_loss + self.consistency_weight * consistency + self.variance_weight * var_loss + self.factor_weight * factor_loss
        return total, {
            "loss": float(total.detach().cpu()),
            "scmae_loss": float(scmae.detach().cpu()),
            "reconstruction_loss": float(rec.detach().cpu()),
            "mask_loss": float(mask_loss.detach().cpu()),
            "module_loss": float(module_loss.detach().cpu()),
            "consistency_loss": float(consistency.detach().cpu()),
            "variance_loss": float(var_loss.detach().cpu()),
            "factor_loss": float(factor_loss.detach().cpu()),
            "tau_mean": float(out["tau_mean"].detach().cpu()),
            "delta_abs_mean": float(out["delta_abs_mean"].detach().cpu()),
            "module_repr_norm": float(out["module_repr"].norm(dim=-1).mean().detach().cpu()),
        }
