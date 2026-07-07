from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class CICLClusterScMAE(nn.Module):
    """scMAE body with CICL-style projection and Student-t cluster head."""

    def __init__(
        self,
        num_genes: int,
        n_clusters: int,
        hidden_size: int = 128,
        projection_dim: int = 64,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(hidden_size + self.num_genes, self.num_genes)
        self.projector = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_size, projection_dim),
        )
        self.cluster_centers = nn.Parameter(torch.randn(self.n_clusters, hidden_size) * 0.02)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        projection = F.normalize(self.projector(latent), dim=1)
        q = self.student_q(latent)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "projection": projection,
            "cluster_q": q,
        }

    def student_q(self, latent: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(latent, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def random_mask_view(self, x: torch.Tensor, mask_prob: float, jitter_std: float = 0.0) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        view = x.masked_fill(mask.bool(), 0.0)
        if jitter_std > 0:
            view = view + torch.randn_like(view) * float(jitter_std) * (1.0 - mask)
        return view, mask

