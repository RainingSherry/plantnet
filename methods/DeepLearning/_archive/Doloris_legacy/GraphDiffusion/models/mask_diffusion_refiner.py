import torch
import torch.nn as nn


class MaskDiffusionRefiner(nn.Module):
    def __init__(self, cell_dim: int, steps: int = 3, dropout: float = 0.1):
        super().__init__()
        self.steps = steps
        self.noise_proj = nn.Linear(cell_dim, cell_dim)
        self.blocks = nn.ModuleList([
            nn.Sequential(
                nn.Linear(cell_dim, cell_dim),
                nn.LayerNorm(cell_dim),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(cell_dim, cell_dim),
            ) for _ in range(steps)
        ])
        self.alpha = nn.Parameter(torch.tensor(0.5))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        h = x
        for block in self.blocks:
            noise = self.noise_proj(torch.randn_like(h) * 0.01)
            h = h + torch.tanh(self.alpha) * noise + block(h)
        return h
