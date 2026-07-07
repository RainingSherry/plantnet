from __future__ import annotations

import math

import torch
from torch import nn


class SinusoidalTimeEmbedding(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.dim = int(dim)

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        half = self.dim // 2
        freq = torch.exp(torch.arange(half, device=t.device, dtype=t.dtype) * (-math.log(10000.0) / max(1, half - 1)))
        emb = t[:, None] * freq[None, :]
        emb = torch.cat([torch.sin(emb), torch.cos(emb)], dim=1)
        if emb.shape[1] < self.dim:
            emb = torch.nn.functional.pad(emb, (0, self.dim - emb.shape[1]))
        return emb


class ScDiffusionPseudoCondScMAE(nn.Module):
    """scMAE backbone with pseudo-condition DDPM latent noise prediction."""

    def __init__(self, input_dim: int, n_conditions: int, hidden_size: int = 128, decoder_hidden: int = 128, time_dim: int = 64, cond_dim: int = 32, dropout: float = 0.0):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.n_conditions = int(n_conditions)
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
        self.decoder = nn.Sequential(nn.Linear(hidden_size + input_dim, decoder_hidden), nn.Mish(inplace=True), nn.Linear(decoder_hidden, input_dim))
        self.time_embedding = SinusoidalTimeEmbedding(time_dim)
        self.condition_embedding = nn.Embedding(n_conditions, cond_dim)
        self.noise_predictor = nn.Sequential(
            nn.Linear(hidden_size + time_dim + cond_dim, decoder_hidden),
            nn.LayerNorm(decoder_hidden),
            nn.Mish(inplace=True),
            nn.Linear(decoder_hidden, decoder_hidden),
            nn.Mish(inplace=True),
            nn.Linear(decoder_hidden, hidden_size),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return recon, mask_logits

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encode(x)
        recon, mask_logits = self.decode(z)
        return {"embedding": z, "reconstruction": recon, "mask_logits": mask_logits}

    def predict_noise(self, noisy_z: torch.Tensor, t: torch.Tensor, condition_id: torch.Tensor) -> torch.Tensor:
        time = self.time_embedding(t)
        cond = self.condition_embedding(condition_id)
        return self.noise_predictor(torch.cat([noisy_z, time, cond], dim=1))
