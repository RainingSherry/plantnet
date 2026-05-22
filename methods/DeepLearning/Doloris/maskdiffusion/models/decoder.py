import torch.nn as nn


class GeneDecoder(nn.Module):
    def __init__(self, latent_dim: int, n_genes: int, hidden_dim: int = 256):
        super().__init__()
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, n_genes),
        )

    def forward(self, z):
        return self.decoder(z)
