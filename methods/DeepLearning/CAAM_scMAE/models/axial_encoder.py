from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn

from .cell_axis import CellAxisContextAttention
from .gene_axis import GeneAxisEncoder
from .module_tokenizer import GeneModuleTokenizer


class AxialEncoder(nn.Module):
    def __init__(
        self,
        assignment: np.ndarray,
        token_dim: int,
        latent_dim: int,
        gene_attention_heads: int,
        gene_attention_layers: int,
        attention_dropout: float,
        cell_attention_heads: int | None = None,
    ) -> None:
        super().__init__()
        heads = int(cell_attention_heads or gene_attention_heads)
        self.tokenizer = GeneModuleTokenizer(assignment, token_dim)
        self.gene_axis = GeneAxisEncoder(token_dim, gene_attention_heads, gene_attention_layers, attention_dropout)
        self.cell_axis = CellAxisContextAttention(token_dim, heads, attention_dropout)
        self.proj = nn.Sequential(nn.Linear(token_dim, latent_dim), nn.LayerNorm(latent_dim))
        self.context_tokens: torch.Tensor | None = None
        self.context_indices: torch.Tensor | None = None

    def encode_gene(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        tokens = self.tokenizer(x)
        return self.gene_axis(tokens)

    def set_context_cache(self, context_tokens: torch.Tensor, context_indices: torch.Tensor) -> None:
        self.context_tokens = context_tokens.detach()
        self.context_indices = context_indices.detach().long()

    def refresh_context_cache(self, context_x: torch.Tensor, context_indices: torch.Tensor) -> None:
        was_training = self.training
        self.eval()
        with torch.no_grad():
            context_tokens, _ = self.encode_gene(context_x)
        self.set_context_cache(context_tokens, context_indices)
        self.train(was_training)

    def context_cache_checksum(self) -> float:
        if self.context_tokens is None:
            return 0.0
        return float(self.context_tokens.detach().sum().cpu())

    def forward(self, x: torch.Tensor, query_indices: torch.Tensor | None = None) -> dict[str, torch.Tensor | None]:
        gene_tokens, gene_attn = self.encode_gene(x)
        if self.context_tokens is not None:
            cell_tokens, cell_attn = self.cell_axis(
                gene_tokens,
                self.context_tokens.to(x.device),
                query_indices,
                self.context_indices.to(x.device) if self.context_indices is not None else None,
            )
        else:
            cell_tokens = gene_tokens
            b, m, _ = gene_tokens.shape
            cell_attn = torch.zeros((b, m, 1, 1), dtype=gene_tokens.dtype, device=gene_tokens.device)
        z = self.proj(cell_tokens.mean(dim=1))
        return {"z": z, "module_tokens": cell_tokens, "gene_attn": gene_attn, "cell_attn": cell_attn}

