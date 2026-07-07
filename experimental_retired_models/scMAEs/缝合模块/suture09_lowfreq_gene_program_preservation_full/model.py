from __future__ import annotations

import torch
from torch import nn


class LowFrequencyProgramAdapter(nn.Module):
    """Latent adapter guided by a coarse low-frequency gene-program target."""

    def __init__(self, program_dim: int, latent_dim: int, hidden_dim: int = 128, adapter_weight: float = 0.04):
        super().__init__()
        self.adapter_weight = float(adapter_weight)
        self.program_encoder = nn.Sequential(
            nn.Linear(program_dim, hidden_dim),
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
        self.delta = nn.Sequential(
            nn.LayerNorm(latent_dim),
            nn.Linear(latent_dim, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, latent_dim),
        )
        self.program_head = nn.Sequential(
            nn.LayerNorm(latent_dim),
            nn.Linear(latent_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, program_dim),
        )

    def forward(self, base_z: torch.Tensor, lowfreq_program: torch.Tensor) -> dict:
        program_z = self.program_encoder(lowfreq_program)
        gate = self.gate(torch.cat([base_z, program_z, program_z - base_z], dim=1))
        latent = base_z + self.adapter_weight * gate * self.delta(program_z - base_z)
        return {
            "latent": latent,
            "program_z": program_z,
            "program_gate": gate,
            "program_pred": self.program_head(latent),
        }


class LowFreqGeneProgramScMAE(nn.Module):
    """Independent scMAE with low-frequency gene-program preservation."""

    def __init__(
        self,
        input_dim: int,
        program_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        adapter_weight: float = 0.04,
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
        self.lowfreq = LowFrequencyProgramAdapter(program_dim, latent_dim, hidden_dim // 4, adapter_weight)
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

    def encode_with_program(self, x: torch.Tensor, lowfreq_program: torch.Tensor) -> dict:
        base_z = self.encoder(x)
        out = self.lowfreq(base_z, lowfreq_program)
        out["base_latent"] = base_z
        return out

    def forward(self, x: torch.Tensor, lowfreq_program: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        out = self.encode_with_program(corrupted, lowfreq_program)
        latent = out["latent"]
        return {
            "latent": latent,
            "base_latent": out["base_latent"],
            "program_z": out["program_z"],
            "program_gate": out["program_gate"],
            "program_pred": out["program_pred"],
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor, lowfreq_program: torch.Tensor | None = None) -> torch.Tensor:
        if lowfreq_program is None:
            return self.encoder(x)
        return self.encode_with_program(x, lowfreq_program)["latent"]
