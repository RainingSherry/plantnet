from __future__ import annotations

import torch
import torch.nn as nn


class GeneAxisBlock(nn.Module):
    def __init__(self, token_dim: int, heads: int, dropout: float) -> None:
        super().__init__()
        self.ln1 = nn.LayerNorm(token_dim)
        self.attn = nn.MultiheadAttention(token_dim, heads, dropout=dropout, batch_first=True)
        self.ln2 = nn.LayerNorm(token_dim)
        self.ffn = nn.Sequential(
            nn.Linear(token_dim, token_dim * 4),
            nn.Mish(),
            nn.Dropout(dropout),
            nn.Linear(token_dim * 4, token_dim),
        )

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.ln1(x)
        out, attn = self.attn(h, h, h, need_weights=True, average_attn_weights=False)
        x = x + out
        x = x + self.ffn(self.ln2(x))
        return x, attn


class GeneAxisEncoder(nn.Module):
    def __init__(self, token_dim: int, heads: int, layers: int, dropout: float) -> None:
        super().__init__()
        self.blocks = nn.ModuleList([GeneAxisBlock(token_dim, heads, dropout) for _ in range(layers)])

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        attn = None
        for block in self.blocks:
            x, attn = block(x)
        return x, attn

