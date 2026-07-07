from __future__ import annotations

import torch
from torch import nn


class LatentFeatureCorrection(nn.Module):
    """FCM-inspired masked/unmasked latent correction without 2D spatial convolutions."""

    def __init__(self, latent_dim: int, min_gate: float = 0.05):
        super().__init__()
        self.min_gate = float(min_gate)
        self.channel_gate = nn.Sequential(
            nn.LayerNorm(latent_dim * 4),
            nn.Linear(latent_dim * 4, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, latent_dim),
            nn.Sigmoid(),
        )
        self.global_weights = nn.Parameter(torch.ones(2))
        self.refine = nn.Sequential(nn.LayerNorm(latent_dim), nn.Linear(latent_dim, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim))

    def forward(self, masked_z: torch.Tensor, clean_z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        clean_teacher = clean_z.detach()
        stats = torch.cat([masked_z, clean_teacher, clean_teacher - masked_z, masked_z * clean_teacher], dim=1)
        raw_gate = self.channel_gate(stats)
        gate = self.min_gate + (1.0 - self.min_gate) * raw_gate
        weights = torch.relu(self.global_weights)
        weights = weights / (weights.sum() + 1e-6)
        corrected = masked_z + weights[0] * gate * (clean_teacher - masked_z)
        corrected = corrected + weights[1] * self.refine(corrected)
        return corrected, gate, weights


class FCMMaskUnmaskScMAE(nn.Module):
    """Independent scMAE with train-time masked/unmasked latent correction fusion."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        module_weight: float = 0.35,
        min_gate: float = 0.05,
    ):
        super().__init__()
        self.mask_prob = float(mask_prob)
        self.module_weight = float(module_weight)
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
        self.correction = LatentFeatureCorrection(latent_dim, min_gate=min_gate)
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

    def forward(self, x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        masked_z = self.encoder(corrupted)
        clean_z = self.encoder(x)
        corrected_z, gate, weights = self.correction(masked_z, clean_z)
        latent = masked_z + self.module_weight * (corrected_z - masked_z)
        return {
            "latent": latent,
            "masked_z": masked_z,
            "clean_z": clean_z,
            "correction_gate": gate,
            "fusion_weights": weights,
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)
