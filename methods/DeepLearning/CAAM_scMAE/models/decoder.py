from __future__ import annotations

import torch
import torch.nn as nn


class Decoder(nn.Module):
    def __init__(self, latent_dim: int, n_genes: int, conditioning: str = "pred_detached") -> None:
        super().__init__()
        if conditioning not in {"none", "pred_detached", "oracle_mask"}:
            raise ValueError(f"Unknown decoder mask conditioning: {conditioning}")
        self.conditioning = conditioning
        in_dim = latent_dim if conditioning == "none" else latent_dim + n_genes
        self.net = nn.Sequential(nn.Linear(in_dim, max(latent_dim, n_genes)), nn.Mish(), nn.Linear(max(latent_dim, n_genes), n_genes))

    def forward(self, z: torch.Tensor, mask_logits: torch.Tensor, oracle_mask: torch.Tensor | None = None) -> torch.Tensor:
        if self.conditioning == "none":
            return self.net(z)
        if self.conditioning == "oracle_mask":
            if oracle_mask is None:
                raise ValueError("oracle_mask conditioning requires oracle_mask")
            cond = oracle_mask.detach()
        else:
            cond = torch.sigmoid(mask_logits).detach()
        return self.net(torch.cat([z, cond], dim=1))

