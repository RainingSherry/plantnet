from __future__ import annotations

import torch
from torch import nn


class CMTScMAE(nn.Module):
    """scMAE with collaborative masking and collaborative feature targets."""

    def __init__(self, input_dim: int, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.0):
        super().__init__()
        self.input_dim = input_dim
        self.hidden_size = hidden_size
        self.encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.Mish(inplace=True),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + input_dim, decoder_hidden),
            nn.Mish(inplace=True),
            nn.Linear(decoder_hidden, input_dim),
        )
        self.target_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.LayerNorm(hidden_size))

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {"embedding": z, "reconstruction": recon, "mask_logits": mask_logits, "target_pred": self.target_head(z)}


@torch.no_grad()
def update_ema(online: nn.Module, target: nn.Module, momentum: float) -> None:
    for p_online, p_target in zip(online.parameters(), target.parameters()):
        p_target.data.mul_(momentum).add_(p_online.data, alpha=1.0 - momentum)
