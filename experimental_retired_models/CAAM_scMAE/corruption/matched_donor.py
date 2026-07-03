from __future__ import annotations

import torch

from .base import CorruptionResult, hard_corruption, straight_through_corruption


class MatchedDonorCorruption:
    corruption_type = "matched_donor"

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
        effective = (mask > 0).float() * eligibility.float()
        if torch.any((mask > 0.5) & (~eligibility)):
            warning = "selected_mask differs from effective_mask because ineligible positions were selected"
        else:
            warning = ""
        budget_deficit = torch.zeros(x.shape[0], dtype=torch.float32, device=x.device)
        info = {
            "selected_mask_rate": float(mask.detach().mean().cpu()),
            "effective_changed_rate": float(((x_tilde.detach() - x).abs() > 0).float().mean().cpu()),
            "corruption_type": self.corruption_type,
            "warning": warning,
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
