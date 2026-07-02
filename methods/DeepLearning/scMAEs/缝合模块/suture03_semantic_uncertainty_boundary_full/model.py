from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class SemanticUncertaintyGate(nn.Module):
    """Core / boundary / rare-risk decoupling gate for cell embeddings."""

    def __init__(self, latent_dim: int):
        super().__init__()
        self.gate = nn.Sequential(
            nn.LayerNorm(latent_dim),
            nn.Linear(latent_dim, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, 3),
        )
        self.core_refine = nn.Sequential(nn.Linear(latent_dim, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim))
        self.boundary_guard = nn.Sequential(nn.Linear(latent_dim, latent_dim), nn.Tanh())
        self.rare_refine = nn.Sequential(nn.Linear(latent_dim, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim))

    def forward(self, z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        gates = F.softmax(self.gate(z), dim=1)
        core = self.core_refine(z)
        boundary = self.boundary_guard(z)
        rare = self.rare_refine(z)
        delta = gates[:, 0:1] * core + gates[:, 2:3] * rare - 0.5 * gates[:, 1:2] * boundary
        return delta, gates


class SemanticUncertaintyScMAE(nn.Module):
    """Independent scMAE with SID-inspired semantic uncertainty gates."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        module_weight: float = 0.2,
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
        self.semantic_gate = SemanticUncertaintyGate(latent_dim)
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

    def encode_with_gate(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        base = self.encoder(x)
        delta, gates = self.semantic_gate(base)
        latent = base + self.module_weight * delta
        return latent, base, delta, gates

    def forward(self, x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        latent, base, delta, gates = self.encode_with_gate(corrupted)
        return {
            "latent": latent,
            "base_latent": base,
            "delta": delta,
            "gates": gates,
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
            "corrupted": corrupted,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        latent, _, _, _ = self.encode_with_gate(x)
        return latent
