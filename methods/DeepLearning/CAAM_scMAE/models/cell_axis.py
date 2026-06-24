from __future__ import annotations

import torch
import torch.nn as nn


class CellAxisContextAttention(nn.Module):
    def __init__(self, token_dim: int, heads: int, dropout: float) -> None:
        super().__init__()
        if token_dim % heads != 0:
            raise ValueError("token_dim must be divisible by heads")
        self.heads = int(heads)
        self.head_dim = int(token_dim // heads)
        self.scale = self.head_dim ** -0.5
        self.q = nn.Linear(token_dim, token_dim)
        self.k = nn.Linear(token_dim, token_dim)
        self.v = nn.Linear(token_dim, token_dim)
        self.out = nn.Linear(token_dim, token_dim)
        self.dropout = nn.Dropout(dropout)
        self.ln = nn.LayerNorm(token_dim)
        self.ffn_ln = nn.LayerNorm(token_dim)
        self.ffn = nn.Sequential(nn.Linear(token_dim, token_dim * 4), nn.Mish(), nn.Dropout(dropout), nn.Linear(token_dim * 4, token_dim))

    def forward(
        self,
        query_tokens: torch.Tensor,
        context_tokens: torch.Tensor,
        query_indices: torch.Tensor | None,
        context_indices: torch.Tensor | None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        b, m, d = query_tokens.shape
        c = int(context_tokens.shape[0])
        q = self.q(self.ln(query_tokens)).view(b, m, self.heads, self.head_dim)
        k = self.k(context_tokens).view(c, m, self.heads, self.head_dim)
        v = self.v(context_tokens).view(c, m, self.heads, self.head_dim)
        scores = torch.einsum("bmhd,cmhd->bmhc", q, k) * self.scale
        valid = torch.ones((b, c), dtype=torch.bool, device=query_tokens.device)
        if query_indices is not None and context_indices is not None:
            valid = query_indices.view(-1, 1).to(query_tokens.device) != context_indices.view(1, -1).to(query_tokens.device)
        scores = scores.masked_fill(~valid.view(b, 1, 1, c), -1.0e9)
        attn = torch.softmax(scores, dim=-1)
        attn = attn * valid.view(b, 1, 1, c).float()
        attn = attn / attn.sum(dim=-1, keepdim=True).clamp_min(1.0e-8)
        attn = self.dropout(attn)
        out = torch.einsum("bmhc,cmhd->bmhd", attn, v).reshape(b, m, d)
        x = query_tokens + self.out(out)
        x = x + self.ffn(self.ffn_ln(x))
        return x, attn

