import torch
import torch.nn as nn


class SupportMaskNet(nn.Module):
    def __init__(self, n_genes: int, hidden_dim: int = 256):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(n_genes, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, n_genes),
        )

    def forward(self, x):
        return self.net(x)

    def predict_proba(self, x):
        return torch.sigmoid(self.forward(x))
