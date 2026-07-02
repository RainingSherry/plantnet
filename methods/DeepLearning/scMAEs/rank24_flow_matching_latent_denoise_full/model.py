from __future__ import annotations

import math

import torch
from torch import nn


def time_embedding(t: torch.Tensor, dim: int) -> torch.Tensor:
    half = dim // 2
    freqs = torch.exp(torch.arange(half, device=t.device, dtype=t.dtype) * (-math.log(10000.0) / max(1, half - 1)))
    args = t[:, None] * freqs[None, :]
    emb = torch.cat([torch.sin(args), torch.cos(args)], dim=1)
    if emb.shape[1] < dim:
        emb = torch.nn.functional.pad(emb, (0, dim - emb.shape[1]))
    return emb


class FlowMatchingLatentScMAE(nn.Module):
    """scMAE encoder with a lightweight latent vector-field head."""

    def __init__(self, num_genes: int, hidden_size: int = 128, time_dim: int = 32, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.time_dim = int(time_dim)
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
        )
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(self.hidden_size + self.num_genes, 256),
            nn.Mish(),
            nn.Dropout(dropout),
            nn.Linear(256, self.num_genes),
        )
        self.flow_head = nn.Sequential(
            nn.Linear(self.hidden_size + self.time_dim, self.hidden_size * 2),
            nn.SiLU(),
            nn.LayerNorm(self.hidden_size * 2),
            nn.Linear(self.hidden_size * 2, self.hidden_size),
            nn.SiLU(),
            nn.Linear(self.hidden_size, self.hidden_size),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encode(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction}

    def predict_flow(self, z_t: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        return self.flow_head(torch.cat([z_t, time_embedding(t, self.time_dim)], dim=1))

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode(x)

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask
