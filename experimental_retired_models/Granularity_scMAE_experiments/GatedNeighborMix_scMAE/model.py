from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class GatedNeighborMixScMAE(nn.Module):
    """rank13's DEC-style scMAE (reproduced faithfully) as the backbone.

    Identical to rank13_masked_sc_cluster_target's MaskedScClusterScMAE so the
    proven DEC behavior (Quake ARI 0.91) is preserved exactly. The gated
    NeighborMix branch is added in the training loop / loss, not here.
    """

    def __init__(self, num_genes: int, n_clusters: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.hidden_size = int(hidden_size)
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
        self.cluster_centers = nn.Parameter(torch.randn(self.n_clusters, hidden_size) * 0.02)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        q = self.student_q(latent)
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction, "cluster_q": q}

    def student_q(self, latent: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(latent, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def random_mask(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        """Zero-masking (rank13). Returns (masked_x, mask)."""
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    @torch.no_grad()
    def initialize_centers(self, centers: torch.Tensor) -> None:
        if centers.shape != self.cluster_centers.shape:
            raise ValueError(f"center shape {tuple(centers.shape)} != {tuple(self.cluster_centers.shape)}")
        self.cluster_centers.copy_(centers)

    @staticmethod
    def target_distribution(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)
