from __future__ import annotations

import math

import torch
from torch import nn


class SinusoidalTimeEmbedding(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.dim = int(dim)

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        half = self.dim // 2
        freq = torch.exp(torch.arange(half, device=t.device, dtype=t.dtype) * (-math.log(10000.0) / max(1, half - 1)))
        emb = t[:, None] * freq[None, :]
        emb = torch.cat([torch.sin(emb), torch.cos(emb)], dim=1)
        if emb.shape[1] < self.dim:
            emb = torch.nn.functional.pad(emb, (0, self.dim - emb.shape[1]))
        return emb


class PlanetGrnAttentionScMAE(nn.Module):
    """scMAE with Planet-inspired time-guided regulatory edge attention."""

    def __init__(
        self,
        input_dim: int,
        hidden_size: int = 128,
        decoder_hidden: int = 128,
        gene_dim: int = 48,
        time_dim: int = 48,
        edge_hidden: int = 128,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
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
        self.gene_embedding = nn.Embedding(input_dim, gene_dim)
        self.time_embedding = SinusoidalTimeEmbedding(time_dim)
        self.context_projection = nn.Sequential(nn.Linear(hidden_size, gene_dim), nn.LayerNorm(gene_dim), nn.Mish(inplace=True))
        self.edge_mlp = nn.Sequential(
            nn.Linear(gene_dim * 4 + time_dim, edge_hidden),
            nn.LayerNorm(edge_hidden),
            nn.Mish(inplace=True),
            nn.Linear(edge_hidden, edge_hidden),
            nn.Mish(inplace=True),
            nn.Linear(edge_hidden, 1),
        )

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {"embedding": z, "reconstruction": recon, "mask_logits": mask_logits}

    def edge_logits(self, z: torch.Tensor, src: torch.Tensor, dst: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        src_e = self.gene_embedding(src)
        dst_e = self.gene_embedding(dst)
        pair = torch.cat([src_e, dst_e, torch.abs(src_e - dst_e), src_e * dst_e], dim=-1)
        time = self.time_embedding(t)
        context = self.context_projection(z)
        conditioned = pair[None, :, :] + torch.cat([context[:, None, :].expand(-1, pair.shape[0], -1)] * 4, dim=-1) * 0.05
        edge_input = torch.cat([conditioned, time[:, None, :].expand(-1, pair.shape[0], -1)], dim=-1)
        return self.edge_mlp(edge_input).squeeze(-1)
