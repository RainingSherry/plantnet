from __future__ import annotations

import torch
from torch import nn


class SelfGuidedMaskScMAE(nn.Module):
    """scMAE with self-guided masking and gene-specific rank-token prediction."""

    def __init__(self, input_dim: int, num_bins: int = 8, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.05):
        super().__init__()
        self.input_dim = input_dim
        self.num_bins = num_bins
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, decoder_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(decoder_hidden, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.decoder = nn.Sequential(nn.Linear(hidden_size, decoder_hidden), nn.GELU(), nn.Dropout(dropout), nn.Linear(decoder_hidden, input_dim))
        self.mask_predictor = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, input_dim))
        self.rank_head = nn.Sequential(nn.Linear(hidden_size, decoder_hidden), nn.GELU(), nn.Dropout(dropout), nn.Linear(decoder_hidden, input_dim * num_bins))
        self.difficulty_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, input_dim), nn.Sigmoid())

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encode(x)
        return {
            "embedding": z,
            "reconstruction": self.decoder(z),
            "mask_logits": self.mask_predictor(z),
            "rank_logits": self.rank_head(z).view(x.shape[0], self.input_dim, self.num_bins),
            "difficulty": self.difficulty_head(z),
        }

    @staticmethod
    def guided_mask(x: torch.Tensor, gene_weights: torch.Tensor, mask_prob: float, uniform_mix: float = 0.25) -> tuple[torch.Tensor, torch.Tensor]:
        weights = gene_weights.to(x.device).float().clamp_min(1e-6)
        weights = weights / weights.mean().clamp_min(1e-6)
        probs = mask_prob * ((1.0 - uniform_mix) * weights.view(1, -1) + uniform_mix)
        probs = probs.clamp(0.0, 0.95)
        mask = torch.rand_like(x) < probs
        return x.masked_fill(mask, 0.0), mask.float()
