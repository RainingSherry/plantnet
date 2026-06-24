from __future__ import annotations

import torch
import torch.nn as nn


class MLPEncoder(nn.Module):
    def __init__(self, n_genes: int, latent_dim: int, hidden_dim: int = 256, dropout: float = 0.0) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(n_genes, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.Mish(),
            nn.Linear(hidden_dim, latent_dim),
            nn.LayerNorm(latent_dim),
            nn.Mish(),
            nn.Linear(latent_dim, latent_dim),
            nn.LayerNorm(latent_dim),
        )

    def forward(self, x: torch.Tensor, **kwargs) -> dict[str, torch.Tensor | None]:
        return {"z": self.net(x), "module_tokens": None, "gene_attn": None, "cell_attn": None}

