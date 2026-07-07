from __future__ import annotations

import math
from typing import Tuple

import torch
import torch.nn as nn


class ScDiVaDiscreteDiffusionMAE(nn.Module):
    """Time-conditioned masked discrete diffusion model for scRNA expression tokens.

    This is not a renamed MLP scMAE.  Each cell is represented as a length-G gene-token
    sequence.  A forward masked discrete diffusion process replaces expression tokens
    by a learned mask token at a time-dependent rate.  The transformer denoiser receives
    gene identity embeddings, corrupted expression-token embeddings, a diffusion-time
    embedding, and a CLS token.  It predicts the original expression token, the continuous
    expression value, and the corrupted positions.
    """

    def __init__(
        self,
        num_genes: int,
        n_bins: int = 32,
        hidden_size: int = 128,
        depth: int = 4,
        num_heads: int = 4,
        dropout: float = 0.1,
        max_steps: int = 100,
    ) -> None:
        super().__init__()
        if num_genes <= 0:
            raise ValueError("num_genes must be positive")
        if n_bins < 4:
            raise ValueError("n_bins must be at least 4")
        self.num_genes = int(num_genes)
        self.n_bins = int(n_bins)
        self.mask_token_id = int(n_bins)
        self.max_steps = int(max_steps)

        self.gene_embedding = nn.Embedding(num_genes, hidden_size)
        self.value_embedding = nn.Embedding(n_bins + 1, hidden_size)
        self.time_embedding = nn.Embedding(max_steps + 1, hidden_size)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.input_norm = nn.LayerNorm(hidden_size)

        layer = nn.TransformerEncoderLayer(
            d_model=hidden_size,
            nhead=num_heads,
            dim_feedforward=hidden_size * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, num_layers=depth)
        self.final_norm = nn.LayerNorm(hidden_size)
        self.token_head = nn.Linear(hidden_size, n_bins)
        self.value_head = nn.Linear(hidden_size, 1)
        self.mask_head = nn.Linear(hidden_size, 1)
        self.cls_projector = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
        )
        nn.init.normal_(self.cls_token, std=0.02)

    def _gene_ids(self, device: torch.device) -> torch.Tensor:
        return torch.arange(self.num_genes, device=device).view(1, self.num_genes)

    def encode(self, corrupted_tokens: torch.Tensor, t: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        if corrupted_tokens.ndim != 2 or corrupted_tokens.shape[1] != self.num_genes:
            raise ValueError(f"tokens must be [batch, {self.num_genes}], got {tuple(corrupted_tokens.shape)}")
        bsz = corrupted_tokens.shape[0]
        if t.ndim == 0:
            t = t.expand(bsz)
        t = t.to(device=corrupted_tokens.device, dtype=torch.long).clamp(0, self.max_steps)

        gene = self.gene_embedding(self._gene_ids(corrupted_tokens.device)).expand(bsz, -1, -1)
        value = self.value_embedding(corrupted_tokens.clamp(0, self.mask_token_id))
        time = self.time_embedding(t).view(bsz, 1, -1)
        tokens = self.input_norm(gene + value + time)
        cls = self.cls_token.expand(bsz, -1, -1) + time
        encoded = self.encoder(torch.cat([cls, tokens], dim=1))
        encoded = self.final_norm(encoded)
        return encoded[:, 0], encoded[:, 1:]

    def forward(self, corrupted_tokens: torch.Tensor, t: torch.Tensor):
        cls, gene_states = self.encode(corrupted_tokens, t)
        token_logits = self.token_head(gene_states)
        value_pred = self.value_head(gene_states).squeeze(-1)
        mask_logits = self.mask_head(gene_states).squeeze(-1)
        return self.cls_projector(cls), token_logits, value_pred, mask_logits

    @torch.no_grad()
    def feature(self, clean_tokens: torch.Tensor) -> torch.Tensor:
        t0 = torch.zeros(clean_tokens.shape[0], dtype=torch.long, device=clean_tokens.device)
        cls, _ = self.encode(clean_tokens, t0)
        return self.cls_projector(cls)


def cosine_mask_schedule(t: torch.Tensor, max_steps: int, max_mask_ratio: float) -> torch.Tensor:
    """Monotone masked-diffusion schedule from clean state to heavily masked state."""
    frac = t.float().clamp_min(0) / float(max(1, max_steps))
    return max_mask_ratio * torch.sin(0.5 * math.pi * frac).clamp(0.0, 1.0)
