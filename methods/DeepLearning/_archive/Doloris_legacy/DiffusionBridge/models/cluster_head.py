from __future__ import annotations

import torch
from torch import nn


class ClusterHead(nn.Module):
    def __init__(self, input_dim: int, n_clusters: int, hidden_dim: int = 256, alpha: float = 1.0):
        super().__init__()
        self.alpha = alpha
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.ReLU(),
        )
        self.prototypes = nn.Parameter(torch.randn(n_clusters, hidden_dim) * 0.02)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        h = self.encoder(z)
        dist = torch.cdist(h, self.prototypes, p=2).pow(2)
        q = 1.0 / (1.0 + dist / self.alpha)
        q = q.pow((self.alpha + 1.0) / 2.0)
        q = q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)
        return q

    @staticmethod
    def target_distribution(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def hard_assign(self, z: torch.Tensor) -> torch.Tensor:
        return self.forward(z).argmax(dim=1)
