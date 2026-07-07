from __future__ import annotations

import torch
from torch import nn


class AsymmetricMAEScRNA(nn.Module):
    """Independent scRNA MAE with scMAE mask prediction and asymmetric decoder."""

    def __init__(self, input_dim: int, hidden_size: int = 128, decoder_hidden: int = 256, dropout: float = 0.0):
        super().__init__()
        self.input_dim = input_dim
        self.encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.Mish(inplace=True),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder_embed = nn.Sequential(nn.Linear(hidden_size, decoder_hidden), nn.LayerNorm(decoder_hidden), nn.Mish(inplace=True))
        self.mask_token = nn.Parameter(torch.zeros(1, input_dim))
        self.decoder = nn.Sequential(
            nn.Linear(decoder_hidden + input_dim, decoder_hidden),
            nn.Mish(inplace=True),
            nn.Linear(decoder_hidden, input_dim),
        )
        nn.init.normal_(self.mask_token, std=0.02)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        dec = self.decoder_embed(z)
        reconstruction = self.decoder(torch.cat([dec, mask_logits + self.mask_token], dim=1))
        return {"embedding": z, "reconstruction": reconstruction, "mask_logits": mask_logits}

    @staticmethod
    def swap_mask_view(x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < mask_prob).float()
        shuffled = x[torch.randperm(x.shape[0], device=x.device)]
        corrupted = torch.where(mask > 0, shuffled, x)
        changed = (corrupted != x).float() * mask
        return corrupted, changed
