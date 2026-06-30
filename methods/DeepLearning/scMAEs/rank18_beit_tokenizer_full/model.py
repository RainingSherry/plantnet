from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class ExpressionPatchEmbedding(nn.Module):
    def __init__(self, num_genes: int, patch_size: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or patch_size <= 0 or hidden_size <= 0:
            raise ValueError("num_genes, patch_size, and hidden_size must be positive")
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.num_patches = int(math.ceil(self.num_genes / self.patch_size))
        self.pad_size = self.num_patches * self.patch_size - self.num_genes
        self.projection = nn.Linear(self.patch_size, hidden_size)
        self.dropout = nn.Dropout(dropout)

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        padded = F.pad(x, (0, self.pad_size)) if self.pad_size else x
        return padded.view(x.shape[0], self.num_patches, self.patch_size)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.dropout(self.projection(self.patchify(x)))


class ExpressionBEiTScMAE(nn.Module):
    """BEiT-style masked discrete token prediction for expression patches."""

    def __init__(
        self,
        num_genes: int,
        vocab_size: int = 32,
        patch_size: int = 20,
        hidden_size: int = 128,
        depth: int = 2,
        num_heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if vocab_size <= 1:
            raise ValueError("vocab_size must be > 1")
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")
        self.num_genes = int(num_genes)
        self.vocab_size = int(vocab_size)
        self.patch_size = int(patch_size)
        self.patch_embed = ExpressionPatchEmbedding(num_genes, patch_size, hidden_size, dropout)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.mask_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.position = nn.Parameter(torch.randn(1, self.patch_embed.num_patches + 1, hidden_size) * 0.02)
        layer = nn.TransformerEncoderLayer(
            d_model=hidden_size,
            nhead=num_heads,
            dim_feedforward=hidden_size * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.blocks = nn.TransformerEncoder(layer, num_layers=depth)
        self.norm = nn.LayerNorm(hidden_size)
        self.lm_head = nn.Linear(hidden_size, vocab_size)
        self.patch_decoder = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, patch_size),
        )
        self.mask_decoder = nn.Linear(hidden_size, patch_size)
        self.embedding_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.Tanh())

    @property
    def num_patches(self) -> int:
        return self.patch_embed.num_patches

    def unpatchify(self, patch_values: torch.Tensor) -> torch.Tensor:
        if patch_values.ndim != 3 or patch_values.shape[1] != self.num_patches or patch_values.shape[2] != self.patch_size:
            raise ValueError(
                f"patch_values must be [batch, {self.num_patches}, {self.patch_size}], got {tuple(patch_values.shape)}"
            )
        return patch_values.reshape(patch_values.shape[0], -1)[:, : self.num_genes]

    def forward_features(self, x: torch.Tensor, patch_mask: torch.Tensor | None) -> torch.Tensor:
        tokens = self.patch_embed(x)
        if patch_mask is not None:
            if patch_mask.ndim != 2 or patch_mask.shape != tokens.shape[:2]:
                raise ValueError(f"patch_mask must be [batch, {self.num_patches}], got {tuple(patch_mask.shape)}")
            mask = patch_mask.to(dtype=tokens.dtype).unsqueeze(-1)
            tokens = tokens * (1.0 - mask) + self.mask_token.to(dtype=tokens.dtype) * mask
        cls = self.cls_token.expand(tokens.shape[0], -1, -1).to(dtype=tokens.dtype)
        tokens = torch.cat([cls, tokens], dim=1) + self.position.to(dtype=tokens.dtype)
        return self.norm(self.blocks(tokens))

    def forward(self, x: torch.Tensor, patch_mask: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        h = self.forward_features(x, patch_mask)
        patch_tokens = h[:, 1:]
        reconstruction = self.unpatchify(self.patch_decoder(patch_tokens))
        mask_logits = self.unpatchify(self.mask_decoder(patch_tokens))
        return {
            "embedding": self.embedding_head(h[:, 0]),
            "tokens": patch_tokens,
            "token_logits": self.lm_head(patch_tokens),
            "reconstruction": reconstruction,
            "mask_logits": mask_logits,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.embedding_head(self.forward_features(x, None)[:, 0])
