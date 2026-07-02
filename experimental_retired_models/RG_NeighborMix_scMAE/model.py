from __future__ import annotations

from typing import Dict, Tuple

import torch
import torch.nn.functional as F

from experimental_retired_models.NeighborMix_scMAE.model import AutoEncoder as _BaseAutoEncoder


class AutoEncoder(_BaseAutoEncoder):
    """scMAE AutoEncoder with gate-weighted per-sample loss."""

    def loss_mask_weighted(
        self,
        x: torch.Tensor,
        y: torch.Tensor,
        mask: torch.Tensor,
        sample_weight: torch.Tensor | None = None,
        mask_loss_scale: float = 1.0,
    ) -> Tuple[torch.Tensor, torch.Tensor, Dict[str, torch.Tensor]]:
        self._check_expression_shape(x, "x")
        self._check_expression_shape(y, "y")
        self._check_expression_shape(mask, "mask")
        if x.shape != y.shape or x.shape != mask.shape:
            raise ValueError("x, y, and mask must have identical shapes.")

        mask = mask.to(dtype=x.dtype, device=x.device)
        y = y.to(dtype=x.dtype, device=x.device)
        latent, mask_logits, reconstruction = self.forward_mask(x)
        raw_mse = F.mse_loss(reconstruction, y, reduction="none")
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        weighted_mse = weights * raw_mse
        if self.normalize_reconstruction_by_weight:
            rec_per = weighted_mse.sum(dim=1) / weights.sum(dim=1).clamp_min(1e-8)
        else:
            rec_per = weighted_mse.mean(dim=1)
        rec_per = (1.0 - self.mask_loss_weight) * rec_per
        mask_per = F.binary_cross_entropy_with_logits(mask_logits, mask, reduction="none").mean(dim=1)
        mask_per = self.mask_loss_weight * mask_per
        total_per = rec_per + float(mask_loss_scale) * mask_per

        if sample_weight is None:
            loss = total_per.mean()
        else:
            w = sample_weight.to(dtype=x.dtype, device=x.device).view(-1)
            loss = (total_per * w).sum() / w.sum().clamp_min(1e-8)
        parts = {
            "reconstruction_loss": rec_per.mean().detach(),
            "mask_loss": mask_per.mean().detach(),
            "total_loss": loss.detach(),
            "mask_positive_rate": mask.mean().detach(),
            "per_sample_loss": total_per.detach(),
        }
        return latent, loss, parts
