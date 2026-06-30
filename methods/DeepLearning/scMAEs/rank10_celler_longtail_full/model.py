from __future__ import annotations

import torch
import torch.nn as nn


class ResidualBlock(nn.Module):
    def __init__(self, hidden_size: int, dropout: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, hidden_size),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return x + self.net(x)


class CellerLongTailScMAE(nn.Module):
    """scMAE encoder-decoder with Celler-inspired long-tail prototype head."""

    def __init__(
        self,
        num_genes: int,
        n_prototypes: int,
        hidden_size: int = 128,
        depth: int = 2,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if num_genes <= 0 or n_prototypes <= 1:
            raise ValueError("num_genes must be positive and n_prototypes must be > 1")
        self.num_genes = int(num_genes)
        self.n_prototypes = int(n_prototypes)
        self.encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
        )
        self.blocks = nn.ModuleList([ResidualBlock(hidden_size, dropout) for _ in range(depth)])
        self.prototype_head = nn.Linear(hidden_size, n_prototypes)
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        h = self.encoder(x)
        for block in self.blocks:
            h = block(h)
        return h

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        embedding = self.encode(x)
        mask_logits = self.mask_predictor(embedding)
        reconstruction = self.decoder(torch.cat([embedding, mask_logits], dim=1))
        return {
            "embedding": embedding,
            "prototype_logits": self.prototype_head(embedding),
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode(x)
