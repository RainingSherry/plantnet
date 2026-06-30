from __future__ import annotations

import torch
import torch.nn as nn


class FuzzyRoughBoundaryScMAE(nn.Module):
    """Expression autoencoder with a fuzzy rough pseudo-concept head."""

    def __init__(
        self,
        num_genes: int,
        n_clusters: int,
        hidden_size: int = 128,
        latent_size: int = 64,
        depth: int = 3,
        dropout: float = 0.1,
        student_alpha: float = 1.0,
    ) -> None:
        super().__init__()
        if num_genes <= 0 or hidden_size <= 0 or latent_size <= 0 or depth <= 0:
            raise ValueError("num_genes, hidden_size, latent_size, and depth must be positive")
        if n_clusters <= 1:
            raise ValueError("n_clusters must be greater than 1")
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.latent_size = int(latent_size)
        self.student_alpha = float(student_alpha)

        encoder_layers: list[nn.Module] = [nn.LayerNorm(num_genes), nn.Linear(num_genes, hidden_size), nn.GELU()]
        for _ in range(max(0, depth - 1)):
            encoder_layers.extend([nn.Dropout(dropout), nn.Linear(hidden_size, hidden_size), nn.GELU()])
        encoder_layers.extend([nn.Dropout(dropout), nn.Linear(hidden_size, latent_size), nn.LayerNorm(latent_size)])
        self.encoder = nn.Sequential(*encoder_layers)

        decoder_layers: list[nn.Module] = [nn.Linear(latent_size, hidden_size), nn.GELU()]
        for _ in range(max(0, depth - 1)):
            decoder_layers.extend([nn.Dropout(dropout), nn.Linear(hidden_size, hidden_size), nn.GELU()])
        decoder_layers.append(nn.Linear(hidden_size, num_genes))
        self.decoder = nn.Sequential(*decoder_layers)
        self.mask_head = nn.Linear(latent_size, num_genes)
        self.cluster_centers = nn.Parameter(torch.empty(n_clusters, latent_size))
        self.boundary_gate = nn.Sequential(
            nn.Linear(latent_size + n_clusters, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, n_clusters),
        )
        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.trunc_normal_(self.cluster_centers, std=0.02)

    def student_membership(self, embedding: torch.Tensor) -> torch.Tensor:
        if embedding.ndim != 2 or embedding.shape[1] != self.latent_size:
            raise ValueError(f"embedding must be [cells, {self.latent_size}], got {tuple(embedding.shape)}")
        dist_sq = torch.sum((embedding.unsqueeze(1) - self.cluster_centers.unsqueeze(0)) ** 2, dim=2)
        numerator = (1.0 + dist_sq / self.student_alpha).pow(-(self.student_alpha + 1.0) / 2.0)
        return numerator / numerator.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        embedding = self.encoder(x)
        membership = self.student_membership(embedding)
        boundary_logits = self.boundary_gate(torch.cat([embedding, membership], dim=1))
        boundary_membership = torch.softmax(boundary_logits, dim=1)
        return {
            "embedding": embedding,
            "reconstruction": self.decoder(embedding),
            "mask_logits": self.mask_head(embedding),
            "membership": membership,
            "boundary_membership": boundary_membership,
            "centers": self.cluster_centers,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        return self.encoder(x)

