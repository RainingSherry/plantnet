from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class GeneTokenStatistics(nn.Module):
    """TSSA-inspired linear gene-token statistics auxiliary branch."""

    def __init__(self, n_genes: int, token_dim: int = 32, latent_dim: int = 32, dropout: float = 0.1):
        super().__init__()
        self.gene_embed = nn.Parameter(torch.randn(n_genes, token_dim) * 0.02)
        self.gene_bias = nn.Parameter(torch.zeros(n_genes, token_dim))
        self.token_proj = nn.Linear(token_dim, token_dim, bias=False)
        self.temperature = nn.Parameter(torch.ones(1))
        self.out = nn.Sequential(nn.LayerNorm(token_dim * 3), nn.Linear(token_dim * 3, latent_dim), nn.GELU(), nn.Dropout(dropout), nn.Linear(latent_dim, latent_dim))

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        tokens = x.unsqueeze(-1) * self.gene_embed.unsqueeze(0) + self.gene_bias.unsqueeze(0)
        w = self.token_proj(tokens)
        w_norm = F.normalize(w, dim=1)
        importance = F.softmax((w_norm.pow(2).sum(dim=-1) * self.temperature.clamp(0.1, 10.0)), dim=1)
        mean = torch.sum(importance.unsqueeze(-1) * w, dim=1)
        second = torch.sum(importance.unsqueeze(-1) * w.pow(2), dim=1)
        dispersion = torch.sqrt(torch.clamp(second - mean.pow(2), min=1e-6))
        max_token = torch.amax(w * importance.unsqueeze(-1), dim=1)
        aux = self.out(torch.cat([mean, dispersion, max_token], dim=1))
        return aux, importance


class TokenStatisticsAuxScMAE(nn.Module):
    """Independent scMAE with a gene-token statistics auxiliary latent branch."""

    def __init__(self, input_dim: int, hidden_dim: int = 512, latent_dim: int = 32, token_dim: int = 32, dropout: float = 0.1, mask_prob: float = 0.4, aux_weight: float = 0.15):
        super().__init__()
        self.mask_prob = float(mask_prob)
        self.aux_weight = float(aux_weight)
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
        self.token_stats = GeneTokenStatistics(input_dim, token_dim, latent_dim, dropout)
        self.mask_predictor = nn.Sequential(nn.Linear(latent_dim, hidden_dim // 2), nn.GELU(), nn.Linear(hidden_dim // 2, input_dim))
        self.decoder = nn.Sequential(nn.Linear(latent_dim, hidden_dim // 2), nn.LayerNorm(hidden_dim // 2), nn.GELU(), nn.Dropout(dropout), nn.Linear(hidden_dim // 2, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, input_dim))
        self.stat_head = nn.Linear(latent_dim, 3)

    def corrupt(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < self.mask_prob).float()
        return x * (1.0 - mask), mask

    def forward(self, x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        base = self.encoder(corrupted)
        aux, importance = self.token_stats(corrupted)
        latent = base + self.aux_weight * aux
        return {
            "latent": latent,
            "base_latent": base,
            "aux_latent": aux,
            "token_importance": importance,
            "stat_pred": self.stat_head(latent),
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        base = self.encoder(x)
        aux, _ = self.token_stats(x)
        return base + self.aux_weight * aux
