from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class FuzzyRoughScMAE(nn.Module):
    """scMAE with fuzzy-rough core/boundary clustering adapter."""

    def __init__(self, input_dim: int, n_clusters: int, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.0):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.n_clusters = int(n_clusters)
        self.encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.Mish(inplace=True),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(nn.Linear(hidden_size + input_dim, decoder_hidden), nn.Mish(inplace=True), nn.Linear(decoder_hidden, input_dim))
        self.cluster_centers = nn.Parameter(torch.randn(n_clusters, hidden_size) * 0.02)
        self.membership_temperature = nn.Parameter(torch.tensor(1.0))

    def memberships(self, z: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(z, self.cluster_centers, p=2).pow(2)
        temp = F.softplus(self.membership_temperature).clamp(0.2, 5.0)
        q = (1.0 + dist / temp).pow(-1.0)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def set_centers(self, centers: torch.Tensor) -> None:
        with torch.no_grad():
            self.cluster_centers.copy_(centers.to(self.cluster_centers.device, dtype=self.cluster_centers.dtype))

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        q = self.memberships(z)
        return {"embedding": z, "reconstruction": recon, "mask_logits": mask_logits, "membership": q}
