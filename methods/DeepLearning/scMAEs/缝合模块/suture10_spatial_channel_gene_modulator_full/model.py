from __future__ import annotations

import torch
from torch import nn


class LatentSpatialChannelModulator(nn.Module):
    """Vector-form SCFM: channel gate plus cell/sample gate on latent features."""

    def __init__(self, latent_dim: int, hidden_dim: int = 128, modulator_weight: float = 0.05):
        super().__init__()
        self.modulator_weight = float(modulator_weight)
        bottleneck = max(8, latent_dim // 2)
        self.norm = nn.LayerNorm(latent_dim)
        self.channel_gate = nn.Sequential(
            nn.Linear(latent_dim, bottleneck),
            nn.GELU(),
            nn.Linear(bottleneck, latent_dim),
            nn.Sigmoid(),
        )
        self.cell_gate = nn.Sequential(
            nn.Linear(latent_dim * 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, 1),
            nn.Sigmoid(),
        )
        self.channel_proj = nn.Sequential(nn.Linear(latent_dim, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim))
        self.cell_proj = nn.Sequential(nn.Linear(latent_dim, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim))

    def forward(self, base_z: torch.Tensor) -> dict:
        z = self.norm(base_z)
        channel_gate = self.channel_gate(z)
        cell_descriptor = torch.cat([z, torch.abs(z)], dim=1)
        cell_gate = self.cell_gate(cell_descriptor)
        channel_delta = self.channel_proj(channel_gate * z)
        cell_delta = self.cell_proj(cell_gate * z)
        latent = base_z + self.modulator_weight * (channel_delta + cell_delta)
        return {
            "latent": latent,
            "channel_gate": channel_gate,
            "cell_gate": cell_gate,
            "channel_delta": channel_delta,
            "cell_delta": cell_delta,
        }


class SpatialChannelGeneModulatorScMAE(nn.Module):
    """Independent scMAE with latent spatial-channel feature modulation."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        modulator_weight: float = 0.05,
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
        self.modulator = LatentSpatialChannelModulator(latent_dim, hidden_dim // 4, modulator_weight)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, input_dim),
        )
        self.mask_predictor = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Linear(hidden_dim // 2, input_dim),
        )

    def corrupt(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < self.mask_prob).float()
        return x * (1.0 - mask), mask

    def encode_modulated(self, x: torch.Tensor) -> dict:
        base_z = self.encoder(x)
        out = self.modulator(base_z)
        out["base_latent"] = base_z
        return out

    def forward(self, x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        out = self.encode_modulated(corrupted)
        latent = out["latent"]
        return {
            "latent": latent,
            "base_latent": out["base_latent"],
            "channel_gate": out["channel_gate"],
            "cell_gate": out["cell_gate"],
            "channel_delta": out["channel_delta"],
            "cell_delta": out["cell_delta"],
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode_modulated(x)["latent"]
