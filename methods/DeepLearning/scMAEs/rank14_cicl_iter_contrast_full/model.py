from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class GenePatchTransformer(nn.Module):
    """Transformer encoder for expression patches used by the CICL loop."""

    def __init__(self, num_genes: int, patch_size: int, hidden_size: int, depth: int, num_heads: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or patch_size <= 0 or hidden_size <= 0:
            raise ValueError("num_genes, patch_size, and hidden_size must be positive")
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.num_patches = int(math.ceil(self.num_genes / self.patch_size))
        self.pad_size = self.num_patches * self.patch_size - self.num_genes
        self.projection = nn.Linear(self.patch_size, hidden_size)
        self.position = nn.Parameter(torch.randn(1, self.num_patches, hidden_size) * 0.02)
        layer = nn.TransformerEncoderLayer(
            d_model=hidden_size,
            nhead=num_heads,
            dim_feedforward=hidden_size * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.transformer = nn.TransformerEncoder(layer, num_layers=depth)
        self.norm = nn.LayerNorm(hidden_size)
        self.dropout = nn.Dropout(dropout)

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        padded = F.pad(x, (0, self.pad_size)) if self.pad_size else x
        return padded.view(x.shape[0], self.num_patches, self.patch_size)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        tokens = self.projection(self.patchify(x)) + self.position
        tokens = self.transformer(self.dropout(tokens))
        return self.norm(tokens).mean(dim=1)


class CICLScMAE(nn.Module):
    """CICL encoder/projection/clustering model with a masked reconstruction branch."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        projection_size: int = 128,
        patch_size: int = 25,
        depth: int = 2,
        num_heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.encoder = GenePatchTransformer(num_genes, patch_size, hidden_size, depth, num_heads, dropout)
        self.projection_head = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.ReLU(inplace=True),
            nn.Linear(hidden_size, projection_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def project(self, embedding: torch.Tensor) -> torch.Tensor:
        if embedding.ndim != 2 or embedding.shape[1] != self.hidden_size:
            raise ValueError(f"embedding must be [batch, {self.hidden_size}], got {tuple(embedding.shape)}")
        return F.normalize(self.projection_head(embedding), dim=1)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        embedding = self.encode(x)
        mask_logits = self.mask_predictor(embedding)
        reconstruction = self.decoder(torch.cat([embedding, mask_logits], dim=1))
        return {
            "embedding": embedding,
            "projection": self.project(embedding),
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode(x)
