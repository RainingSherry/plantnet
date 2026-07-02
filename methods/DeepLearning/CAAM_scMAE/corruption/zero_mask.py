from __future__ import annotations

import torch


class ZeroMaskCorruption:
    corruption_type = "zero_mask"

    def corrupt(self, x: torch.Tensor, mask: torch.Tensor, metadata: dict | None = None):
        effective = mask.bool() & (x > 0)
        zeros = torch.zeros_like(x)
        x_tilde = torch.where(effective, zeros, x)
        info = {
            "selected_mask_rate": float(mask.float().mean().detach().cpu()),
            "effective_changed_rate": float((x_tilde != x).float().mean().detach().cpu()),
            "corruption_type": self.corruption_type,
        }
        return x_tilde, zeros, info

