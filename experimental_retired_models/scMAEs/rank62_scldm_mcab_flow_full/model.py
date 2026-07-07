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


class GeneTokenCrossAttention(nn.Module):
    """Small MCAB-style pooling adapter for exchangeable gene-expression tokens."""

    def __init__(self, input_dim: int, token_dim: int, n_latents: int, n_heads: int, hidden_size: int):
        super().__init__()
        self.gene_embedding = nn.Embedding(input_dim, token_dim)
        self.count_projection = nn.Sequential(nn.Linear(1, token_dim), nn.Mish(inplace=True), nn.Linear(token_dim, token_dim))
        self.latent_queries = nn.Parameter(torch.randn(n_latents, token_dim) * 0.02)
        self.cross_attention = nn.MultiheadAttention(token_dim, n_heads, batch_first=True)
        self.norm_tokens = nn.LayerNorm(token_dim)
        self.norm_latents = nn.LayerNorm(token_dim)
        self.output = nn.Sequential(nn.Linear(n_latents * token_dim, hidden_size), nn.LayerNorm(hidden_size), nn.Mish(inplace=True), nn.Linear(hidden_size, hidden_size))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        b, g = x.shape
        genes = torch.arange(g, device=x.device)
        gene_tokens = self.gene_embedding(genes)[None, :, :].expand(b, -1, -1)
        count_tokens = self.count_projection(x.unsqueeze(-1))
        tokens = self.norm_tokens(gene_tokens + count_tokens)
        queries = self.latent_queries[None, :, :].expand(b, -1, -1)
        latents, _ = self.cross_attention(queries, tokens, tokens, need_weights=False)
        latents = self.norm_latents(latents + queries)
        return self.output(latents.flatten(1))


class ScLdmMcabFlowScMAE(nn.Module):
    """scMAE with MCAB-style gene-token pooling and linear-interpolant flow auxiliary."""

    def __init__(self, input_dim: int, hidden_size: int = 128, decoder_hidden: int = 128, token_dim: int = 32, n_latents: int = 8, n_heads: int = 4, time_dim: int = 64, dropout: float = 0.0):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.vector_encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.Mish(inplace=True),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
        )
        self.attention_encoder = GeneTokenCrossAttention(input_dim, token_dim, n_latents, n_heads, hidden_size)
        self.fuse = nn.Sequential(nn.Linear(hidden_size * 2, hidden_size), nn.LayerNorm(hidden_size), nn.Mish(inplace=True), nn.Linear(hidden_size, hidden_size))
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(nn.Linear(hidden_size + input_dim, decoder_hidden), nn.Mish(inplace=True), nn.Linear(decoder_hidden, input_dim))
        self.time_embedding = SinusoidalTimeEmbedding(time_dim)
        self.flow_head = nn.Sequential(nn.Linear(hidden_size + time_dim, decoder_hidden), nn.LayerNorm(decoder_hidden), nn.Mish(inplace=True), nn.Linear(decoder_hidden, hidden_size))

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        vec = self.vector_encoder(x)
        attn = self.attention_encoder(x)
        return self.fuse(torch.cat([vec, attn], dim=1))

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encode(x)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {"embedding": z, "reconstruction": recon, "mask_logits": mask_logits}

    def predict_flow(self, z_t: torch.Tensor, t: torch.Tensor) -> torch.Tensor:
        return self.flow_head(torch.cat([z_t, self.time_embedding(t)], dim=1))
