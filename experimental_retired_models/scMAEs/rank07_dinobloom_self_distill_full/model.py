from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class DinoBloomScMAE(nn.Module):
    """scMAE body with a DINO-style projection/prototype head."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        projection_dim: int = 128,
        n_prototypes: int = 256,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
        )
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Linear(self.hidden_size + self.num_genes, self.num_genes)
        self.projector = nn.Sequential(
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.GELU(),
            nn.Linear(self.hidden_size, projection_dim),
        )
        self.prototype_head = nn.utils.parametrizations.weight_norm(
            nn.Linear(projection_dim, n_prototypes, bias=False)
        )

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        projection = F.normalize(self.projector(latent), dim=1)
        proto_logits = self.prototype_head(projection)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "projection": projection,
            "proto_logits": proto_logits,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def corrupt_swap(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        should_swap = torch.bernoulli(float(mask_prob) * torch.ones_like(x)).bool()
        replacement = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else x
        corrupted = torch.where(should_swap, replacement, x)
        mask = (corrupted != x).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            corrupted[empty, cols] = replacement[empty, cols]
            mask[empty, cols] = (corrupted[empty, cols] != x[empty, cols]).float()
        return corrupted, mask

