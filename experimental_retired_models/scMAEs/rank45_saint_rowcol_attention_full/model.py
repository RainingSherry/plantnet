from __future__ import annotations

import torch
from torch import nn


class SAINTRowColScMAE(nn.Module):
    """SAINT-style feature-token and row-attention scMAE."""

    def __init__(self, input_dim: int, token_count: int = 64, hidden_size: int = 128, decoder_hidden: int = 128, heads: int = 4, dropout: float = 0.05):
        super().__init__()
        self.input_dim = input_dim
        self.token_count = token_count
        self.hidden_size = hidden_size
        self.tokenizer = nn.Linear(input_dim, token_count * hidden_size)
        self.cls = nn.Parameter(torch.zeros(1, 1, hidden_size))
        enc_layer = nn.TransformerEncoderLayer(d_model=hidden_size, nhead=heads, dim_feedforward=hidden_size * 2, dropout=dropout, activation="gelu", batch_first=True, norm_first=True)
        self.column_encoder = nn.TransformerEncoder(enc_layer, num_layers=1)
        self.row_attention = nn.MultiheadAttention(hidden_size, heads, dropout=dropout, batch_first=True)
        self.row_norm = nn.LayerNorm(hidden_size)
        self.gate = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, input_dim), nn.Sigmoid())
        self.decoder = nn.Sequential(nn.Linear(hidden_size, decoder_hidden), nn.GELU(), nn.Dropout(dropout), nn.Linear(decoder_hidden, input_dim))
        self.mask_predictor = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, input_dim))
        self.proj_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, hidden_size))

    def encode(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        b = x.shape[0]
        tokens = self.tokenizer(x).view(b, self.token_count, self.hidden_size)
        tokens = torch.cat([self.cls.expand(b, -1, -1), tokens], dim=1)
        col_tokens = self.column_encoder(tokens)
        cls = col_tokens[:, 0]
        row_out, attn = self.row_attention(cls.unsqueeze(0), cls.unsqueeze(0), cls.unsqueeze(0), need_weights=True)
        z = self.row_norm(cls + row_out.squeeze(0))
        return z, attn.squeeze(0), self.gate(z)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z, attn, gate = self.encode(x)
        recon = self.decoder(z) * (0.5 + gate)
        return {"embedding": z, "reconstruction": recon, "mask_logits": self.mask_predictor(z), "projection": self.proj_head(z), "row_attention": attn, "feature_gate": gate}

    @staticmethod
    def mask_view(x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = torch.rand_like(x) < mask_prob
        return x.masked_fill(mask, 0.0), mask.float()
