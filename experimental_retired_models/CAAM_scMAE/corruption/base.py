from __future__ import annotations

from dataclasses import dataclass

import torch


@dataclass
class CorruptionResult:
    x_tilde: torch.Tensor
    replacement: torch.Tensor
    selected_mask: torch.Tensor
    effective_mask: torch.Tensor
    eligibility: torch.Tensor
    donor_indices: torch.Tensor | None
    budget_deficit: torch.Tensor
    info: dict


def hard_corruption(x: torch.Tensor, mask_hard: torch.Tensor, replacement: torch.Tensor) -> torch.Tensor:
    return x * (1.0 - mask_hard) + replacement.detach() * mask_hard


def straight_through_corruption(x: torch.Tensor, mask_st: torch.Tensor, replacement: torch.Tensor) -> torch.Tensor:
    return x * (1.0 - mask_st) + replacement.detach() * mask_st

