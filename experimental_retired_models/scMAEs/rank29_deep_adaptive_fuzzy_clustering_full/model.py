from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class DeepAdaptiveFuzzyScMAE(nn.Module):
    """scMAE with adaptive fuzzy membership and core-cell clustering head."""

    def __init__(self, num_genes: int, n_clusters: int, hidden_size: int = 128, anchor_dim: int = 64, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.hidden_size = int(hidden_size)
        self.anchor_dim = int(anchor_dim)
        self.expr_encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.anchor_encoder = nn.Sequential(
            nn.Linear(anchor_dim, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.fusion = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.cluster_centers = nn.Parameter(torch.randn(n_clusters, hidden_size) * 0.02)
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, 256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, num_genes),
        )
        self.anchor_head = nn.Linear(hidden_size, anchor_dim)

    @torch.no_grad()
    def initialize_centers(self, centers: torch.Tensor) -> None:
        if centers.shape != self.cluster_centers.shape:
            raise ValueError(f"center shape {tuple(centers.shape)} != {tuple(self.cluster_centers.shape)}")
        self.cluster_centers.copy_(centers)

    def encode(self, x: torch.Tensor, anchor: torch.Tensor | None = None) -> torch.Tensor:
        expr = self.expr_encoder(x)
        if anchor is None:
            return expr
        anc = self.anchor_encoder(anchor)
        return self.fusion(torch.cat([expr, anc], dim=1)) + 0.25 * anc

    def soft_assign(self, z: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(z, self.cluster_centers).pow(2)
        q = (1.0 + dist).pow(-1.0)
        q = q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)
        return q

    def forward(self, x: torch.Tensor, anchor: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        latent = self.encode(x, anchor)
        q = self.soft_assign(latent)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        return {
            "latent": latent,
            "membership": q,
            "membership_confidence": q.max(dim=1).values,
            "membership_entropy": -(q * q.clamp_min(1e-8).log()).sum(dim=1),
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "anchor_pred": self.anchor_head(latent),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor, anchor: torch.Tensor | None = None) -> torch.Tensor:
        return self.encode(x, anchor)

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    @staticmethod
    def target_distribution(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return (weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)).detach()
