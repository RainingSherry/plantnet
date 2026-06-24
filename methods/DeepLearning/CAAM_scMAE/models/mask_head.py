from __future__ import annotations

import torch
import torch.nn as nn


class MaskHead(nn.Module):
    def __init__(self, latent_dim: int, n_genes: int) -> None:
        super().__init__()
        self.net = nn.Linear(latent_dim, n_genes)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.net(z)

