from __future__ import annotations

import torch
from torch import nn


class MoleBERTContextScMAE(nn.Module):
    """scMAE with a Mole-BERT-style context tokenizer objective."""

    def __init__(self, input_dim: int, num_tokens: int = 128, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.05):
        super().__init__()
        self.input_dim = input_dim
        self.num_tokens = num_tokens
        self.cell_encoder = nn.Sequential(
            nn.Linear(input_dim, decoder_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(decoder_hidden, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.context_encoder = nn.Sequential(
            nn.Linear(input_dim, decoder_hidden),
            nn.GELU(),
            nn.Linear(decoder_hidden, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.context_gate = nn.Sequential(nn.Linear(hidden_size * 2, hidden_size), nn.Sigmoid())
        self.fusion = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size, decoder_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(decoder_hidden, input_dim),
        )
        self.mask_predictor = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, input_dim))
        self.token_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden_size, num_tokens))
        self.edge_head = nn.Sequential(nn.Linear(hidden_size * 3, hidden_size), nn.GELU(), nn.Linear(hidden_size, 1))

    def encode_base(self, x: torch.Tensor) -> torch.Tensor:
        return self.cell_encoder(x)

    def encode(self, x: torch.Tensor, neigh_x: torch.Tensor | None = None) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.cell_encoder(x)
        if neigh_x is None:
            ctx = torch.zeros_like(h)
        else:
            ctx = self.context_encoder(neigh_x.reshape(-1, neigh_x.shape[-1])).view(neigh_x.shape[0], neigh_x.shape[1], -1).mean(dim=1)
        gate = self.context_gate(torch.cat([h, ctx], dim=-1))
        z = self.fusion(torch.cat([h, gate * ctx], dim=-1)) + h
        return z, gate

    def forward(self, x: torch.Tensor, neigh_x: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        z, gate = self.encode(x, neigh_x)
        return {
            "embedding": z,
            "reconstruction": self.decoder(z),
            "mask_logits": self.mask_predictor(z),
            "token_logits": self.token_head(z),
            "context_gate": gate,
        }

    def edge_logits(self, z_a: torch.Tensor, z_b: torch.Tensor) -> torch.Tensor:
        return self.edge_head(torch.cat([z_a, z_b, torch.abs(z_a - z_b)], dim=-1)).squeeze(-1)

    @staticmethod
    def mask_view(x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = torch.rand_like(x) < mask_prob
        return x.masked_fill(mask, 0.0), mask.float()
