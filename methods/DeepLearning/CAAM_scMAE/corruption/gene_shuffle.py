from __future__ import annotations

import torch


class GeneWiseShuffleCorruption:
    corruption_type = "gene_shuffle"

    def corrupt(self, x: torch.Tensor, mask: torch.Tensor, metadata: dict | None = None):
        perm = torch.argsort(torch.rand(x.shape, device=x.device), dim=0)
        shuffled = torch.gather(x, dim=0, index=perm)
        x_tilde = torch.where(mask.bool(), shuffled, x)
        value_used = torch.where(mask.bool(), shuffled, x)
        info = {
            "selected_mask_rate": float(mask.float().mean().detach().cpu()),
            "effective_changed_rate": float((x_tilde != x).float().mean().detach().cpu()),
            "corruption_type": self.corruption_type,
        }
        return x_tilde, value_used, info

