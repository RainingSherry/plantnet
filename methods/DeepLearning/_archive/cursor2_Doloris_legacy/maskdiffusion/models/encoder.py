import torch
import torch.nn as nn
import torch.nn.functional as F


class SparseEncoder(nn.Module):
    """Encoder: gene expression vector X -> latent embedding z.

    Architecture:
        X (n_genes)
          -> Linear(n_genes, hidden_dim) + LayerNorm + GELU + Dropout
          -> Linear(hidden_dim, hidden_dim) + LayerNorm + GELU + Dropout
          -> Linear(hidden_dim, latent_dim)
        Output: z (latent_dim,)
    """

    def __init__(
        self,
        n_genes: int,
        latent_dim: int = 32,
        hidden_dim: int = 256,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, latent_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)
