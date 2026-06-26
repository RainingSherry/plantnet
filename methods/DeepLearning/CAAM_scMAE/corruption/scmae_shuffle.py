from __future__ import annotations

import torch

from .base import CorruptionResult, hard_corruption, straight_through_corruption


class ScMAEShuffleCorruption:
    corruption_type = "scmae_shuffle"

    def corrupt(
        self,
        x: torch.Tensor,
        mask: torch.Tensor,
        replacement: torch.Tensor,
        eligibility: torch.Tensor,
        donor_indices: torch.Tensor | None = None,
        *,
        straight_through: bool = False,
        replacement_info: dict | None = None,
    ) -> CorruptionResult:
        x_tilde = straight_through_corruption(x, mask, replacement) if straight_through else hard_corruption(x, mask, replacement)
        delta = (x_tilde.detach() - x.detach()).abs()
        effective = (mask > 0).float() * (delta > 0).float()
        budget_deficit = torch.zeros(x.shape[0], dtype=torch.float32, device=x.device)
        info = {
            "corruption_type": self.corruption_type,
            "selected_mask_rate": float(mask.detach().mean().cpu()),
            "effective_changed_rate": float((effective > 0).float().mean().cpu()),
        }
        if replacement_info:
            info.update({k: v for k, v in replacement_info.items() if not torch.is_tensor(v)})
        return CorruptionResult(
            x_tilde=x_tilde,
            replacement=replacement,
            selected_mask=mask,
            effective_mask=effective,
            eligibility=eligibility.float(),
            donor_indices=donor_indices,
            budget_deficit=budget_deficit,
            info=info,
        )
