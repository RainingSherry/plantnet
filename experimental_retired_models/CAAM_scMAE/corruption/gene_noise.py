from __future__ import annotations

import torch


class GeneNoiseCorruption:
    corruption_type = "gene_noise"

    def __init__(self, lambda_noise: float = 0.1) -> None:
        self.lambda_noise = float(lambda_noise)

    def corrupt(self, x: torch.Tensor, mask: torch.Tensor, metadata: dict | None = None):
        sigma = self.lambda_noise * x.std(dim=0, unbiased=False).clamp_min(1.0e-8)
        noise = torch.randn_like(x) * sigma.view(1, -1)
        noisy = torch.clamp(x + noise, min=0.0)
        x_tilde = torch.where(mask.bool(), noisy, x)
        info = {
            "selected_mask_rate": float(mask.float().mean().detach().cpu()),
            "effective_changed_rate": float((x_tilde != x).float().mean().detach().cpu()),
            "corruption_type": self.corruption_type,
        }
        return x_tilde, noisy, info

