from __future__ import annotations

import math

import torch
from torch import nn


class NoiseEmbedding(nn.Module):
    def __init__(self, dim: int):
        super().__init__()
        self.dim = int(dim)
        self.mlp = nn.Sequential(nn.Linear(dim, dim), nn.Mish(), nn.Linear(dim, dim))

    def forward(self, sigma: torch.Tensor) -> torch.Tensor:
        if sigma.dim() == 2:
            sigma = sigma[:, 0]
        sigma = torch.log(sigma.float().clamp_min(1e-4))
        half = max(1, self.dim // 2)
        freqs = torch.exp(
            torch.arange(half, device=sigma.device, dtype=torch.float32)
            * (-math.log(10000.0) / max(1, half - 1))
        )
        args = sigma[:, None] * freqs[None, :]
        emb = torch.cat([torch.sin(args), torch.cos(args)], dim=1)
        if emb.shape[1] < self.dim:
            emb = torch.nn.functional.pad(emb, (0, self.dim - emb.shape[1]))
        return self.mlp(emb[:, : self.dim])


class ConsistencyEMAScMAE(nn.Module):
    """
    scMAE backbone with a noise-level embedding and latent projection head.

    The runner owns the EMA teacher; this module is intentionally a complete
    standalone student/teacher body rather than a shared FlexibleScMAE wrapper.
    """

    def __init__(self, num_genes: int, hidden_size: int = 128, projection_dim: int = 128, dropout: float = 0.1):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.mask_value = nn.Parameter(torch.zeros(self.num_genes))
        self.encoder = nn.Sequential(
            nn.Dropout(p=float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
        )
        self.noise_embedding = NoiseEmbedding(self.hidden_size)
        self.latent_norm = nn.LayerNorm(self.hidden_size)
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(self.hidden_size + self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.num_genes),
        )
        self.projector = nn.Sequential(
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, projection_dim),
        )

    def corrupt(
        self,
        x: torch.Tensor,
        mask_prob: float,
        sigma_min: float,
        sigma_max: float,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mask = torch.bernoulli(float(mask_prob) * torch.ones_like(x)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        log_min = math.log(float(sigma_min))
        log_max = math.log(float(sigma_max))
        sigma = torch.exp(torch.empty(x.shape[0], device=x.device).uniform_(log_min, log_max))
        sentinel = self.mask_value.to(dtype=x.dtype, device=x.device)[None, :]
        masked = torch.where(mask.bool(), sentinel.expand_as(x), x)
        noisy = masked + sigma[:, None] * torch.randn_like(x)
        return noisy, mask, sigma

    def forward(self, x: torch.Tensor, sigma: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        latent = self.latent_norm(latent + self.noise_embedding(sigma))
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        projection = self.projector(latent)
        return {
            "latent": latent,
            "projection": projection,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        sigma = torch.zeros(x.shape[0], device=x.device, dtype=x.dtype) + 1e-4
        return self.forward(x, sigma)["latent"]

