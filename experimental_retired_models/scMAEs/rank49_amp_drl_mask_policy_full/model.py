from __future__ import annotations

import torch
from torch import nn


class DuelingMaskPolicy(nn.Module):
    """Tiny dueling-Q controller for adaptive masking actions."""

    def __init__(self, state_dim: int, n_actions: int, hidden_size: int = 64):
        super().__init__()
        self.trunk = nn.Sequential(
            nn.Linear(state_dim, hidden_size),
            nn.GELU(),
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
        )
        self.value = nn.Linear(hidden_size, 1)
        self.advantage = nn.Linear(hidden_size, n_actions)

    def forward(self, state: torch.Tensor) -> torch.Tensor:
        h = self.trunk(state)
        adv = self.advantage(h)
        return self.value(h) + adv - adv.mean(dim=-1, keepdim=True)


class AMPDRLMaskScMAE(nn.Module):
    """scMAE backbone with an adaptive-masking policy adapter."""

    def __init__(self, input_dim: int, n_actions: int, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.0):
        super().__init__()
        self.input_dim = input_dim
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
        self.policy = DuelingMaskPolicy(state_dim=6, n_actions=n_actions, hidden_size=64)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        reconstruction = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {"embedding": z, "reconstruction": reconstruction, "mask_logits": mask_logits}

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)
