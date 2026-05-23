from typing import Sequence, Tuple

import torch
import torch.nn as nn


class SupportPooling(nn.Module):
    def __init__(self, gene_dim: int, cell_dim: int, module_count: int):
        super().__init__()
        self.gene_dim = gene_dim
        self.cell_dim = cell_dim
        self.module_count = module_count
        self.proj = nn.Sequential(
            nn.Linear(gene_dim, cell_dim),
            nn.LayerNorm(cell_dim),
            nn.GELU(),
            nn.Linear(cell_dim, cell_dim),
        )
        self.module_proj = nn.Linear(gene_dim, module_count)

    def forward(
        self,
        gene_embeddings: torch.Tensor,
        support_indices: Sequence[torch.Tensor],
        support_weights: Sequence[torch.Tensor],
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        cell_vectors = []
        module_vectors = []
        for indices, weights in zip(support_indices, support_weights):
            selected = gene_embeddings.index_select(0, indices)
            w = weights.unsqueeze(-1)
            pooled = (selected * w).sum(dim=0)
            cell_vectors.append(self.proj(pooled))
            module_vectors.append(self.module_proj(pooled))
        return torch.stack(cell_vectors, dim=0), torch.stack(module_vectors, dim=0)
