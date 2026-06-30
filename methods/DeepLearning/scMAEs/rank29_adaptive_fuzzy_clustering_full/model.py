from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class FuzzyMembershipLayer(nn.Module):
    """Fuzzy c-means membership layer in latent bottleneck space."""

    def __init__(
        self,
        n_clusters: int,
        latent_size: int,
        initial_fuzzifier: float = 2.0,
        min_fuzzifier: float = 1.1,
        eps: float = 1e-6,
    ) -> None:
        super().__init__()
        if n_clusters <= 1 or latent_size <= 0:
            raise ValueError("n_clusters must be > 1 and latent_size must be positive")
        if initial_fuzzifier <= min_fuzzifier:
            raise ValueError("initial_fuzzifier must be larger than min_fuzzifier")
        self.n_clusters = int(n_clusters)
        self.latent_size = int(latent_size)
        self.min_fuzzifier = float(min_fuzzifier)
        self.eps = float(eps)
        self.centers = nn.Parameter(torch.empty(n_clusters, latent_size))
        raw = math.log(math.exp(float(initial_fuzzifier) - self.min_fuzzifier) - 1.0)
        self.raw_fuzzifier = nn.Parameter(torch.tensor(raw, dtype=torch.float32))
        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.trunc_normal_(self.centers, std=0.02)

    def fuzzifier(self) -> torch.Tensor:
        return F.softplus(self.raw_fuzzifier) + self.min_fuzzifier

    def forward(self, embedding: torch.Tensor) -> dict[str, torch.Tensor]:
        if embedding.ndim != 2 or embedding.shape[1] != self.latent_size:
            raise ValueError(f"embedding must be [cells, {self.latent_size}], got {tuple(embedding.shape)}")
        dist_sq = torch.sum((embedding.unsqueeze(1) - self.centers.unsqueeze(0)) ** 2, dim=2).clamp_min(self.eps)
        m = self.fuzzifier().to(device=embedding.device, dtype=embedding.dtype)
        logits = -torch.log(dist_sq) / (m - 1.0).clamp_min(1e-4)
        membership = torch.softmax(logits, dim=1)
        weighted_membership = membership.pow(m)
        return {
            "membership": membership,
            "weighted_membership": weighted_membership,
            "distance_sq": dist_sq,
            "fuzzifier": m,
            "centers": self.centers,
        }


class AdaptiveFuzzyClusteringScMAE(nn.Module):
    """DAFC-inspired autoencoder with embedded adaptive fuzzy clustering."""

    def __init__(
        self,
        num_genes: int,
        n_clusters: int,
        hidden_size: int = 128,
        latent_size: int = 64,
        depth: int = 3,
        dropout: float = 0.1,
        initial_fuzzifier: float = 2.0,
    ) -> None:
        super().__init__()
        if num_genes <= 0 or hidden_size <= 0 or latent_size <= 0 or depth <= 0:
            raise ValueError("num_genes, hidden_size, latent_size, and depth must be positive")
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.hidden_size = int(hidden_size)
        self.latent_size = int(latent_size)

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
        self.membership_layer = FuzzyMembershipLayer(
            n_clusters=n_clusters,
            latent_size=latent_size,
            initial_fuzzifier=initial_fuzzifier,
        )
        self.cluster_decoder = nn.Sequential(
            nn.Linear(latent_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        embedding = self.encoder(x)
        reconstruction = self.decoder(embedding)
        mask_logits = self.mask_head(embedding)
        fuzzy = self.membership_layer(embedding)
        fuzzy_latent = fuzzy["membership"] @ fuzzy["centers"].to(dtype=embedding.dtype)
        fuzzy_reconstruction = self.cluster_decoder(fuzzy_latent)
        return {
            "embedding": embedding,
            "reconstruction": reconstruction,
            "mask_logits": mask_logits,
            "membership": fuzzy["membership"],
            "weighted_membership": fuzzy["weighted_membership"],
            "distance_sq": fuzzy["distance_sq"],
            "fuzzifier": fuzzy["fuzzifier"],
            "centers": fuzzy["centers"],
            "fuzzy_latent": fuzzy_latent,
            "fuzzy_reconstruction": fuzzy_reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        return self.encoder(x)

