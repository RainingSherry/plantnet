from __future__ import annotations

import torch
from torch import nn


class GraphWaveAdapter(nn.Module):
    """Wave-style latent adapter driven by a graph-propagated expression context."""

    def __init__(self, input_dim: int, latent_dim: int, hidden_dim: int = 256, wave_weight: float = 0.12):
        super().__init__()
        self.wave_weight = float(wave_weight)
        self.context_encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, latent_dim),
        )
        self.gate = nn.Sequential(
            nn.LayerNorm(latent_dim * 3),
            nn.Linear(latent_dim * 3, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, 1),
            nn.Sigmoid(),
        )
        self.refine = nn.Sequential(nn.LayerNorm(latent_dim), nn.Linear(latent_dim, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim))

    def forward(self, base_z: torch.Tensor, wave_context_x: torch.Tensor) -> dict:
        wave_z = self.context_encoder(wave_context_x)
        gate = self.gate(torch.cat([base_z, wave_z, wave_z - base_z], dim=1))
        delta = self.refine(wave_z - base_z)
        latent = base_z + self.wave_weight * gate * delta
        return {"latent": latent, "wave_z": wave_z, "wave_gate": gate}


class GraphWaveScMAE(nn.Module):
    """Independent scMAE with graph wave propagation auxiliary context."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        wave_weight: float = 0.12,
    ):
        super().__init__()
        self.mask_prob = float(mask_prob)
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, latent_dim),
        )
        self.wave = GraphWaveAdapter(input_dim, latent_dim, hidden_dim // 2, wave_weight)
        self.mask_predictor = nn.Sequential(nn.Linear(latent_dim, hidden_dim // 2), nn.GELU(), nn.Linear(hidden_dim // 2, input_dim))
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, input_dim),
        )

    def corrupt(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < self.mask_prob).float()
        return x * (1.0 - mask), mask

    def encode_with_wave(self, x: torch.Tensor, wave_context_x: torch.Tensor) -> dict:
        base = self.encoder(x)
        out = self.wave(base, wave_context_x)
        out["base_latent"] = base
        return out

    def forward(self, x: torch.Tensor, wave_context_x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        out = self.encode_with_wave(corrupted, wave_context_x)
        latent = out["latent"]
        return {
            "latent": latent,
            "base_latent": out["base_latent"],
            "wave_z": out["wave_z"],
            "wave_gate": out["wave_gate"],
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor, wave_context_x: torch.Tensor | None = None) -> torch.Tensor:
        if wave_context_x is None:
            return self.encoder(x)
        return self.encode_with_wave(x, wave_context_x)["latent"]
