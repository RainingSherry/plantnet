import torch
import torch.nn as nn


class GeneDecoder(nn.Module):
    """Decoder: latent embedding z -> reconstructed gene expression vector.

    Architecture:
        z (latent_dim)
          -> Linear(latent_dim, hidden_dim) + LayerNorm + GELU + Dropout
          -> Linear(hidden_dim, hidden_dim) + LayerNorm + GELU + Dropout
          -> Linear(hidden_dim, n_genes)
        Output: x_recon (n_genes,)
    """

    def __init__(
        self,
        latent_dim: int,
        n_genes: int,
        hidden_dim: int = 256,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_genes),
        )

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.net(z)
