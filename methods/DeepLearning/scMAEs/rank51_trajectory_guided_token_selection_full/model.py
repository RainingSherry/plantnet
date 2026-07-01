from __future__ import annotations

import torch
from torch import nn


class TrajectoryTokenSampler(nn.Module):
    """Gene-token sampler driven by reconstruction-error trajectories."""

    def __init__(self, feature_dim: int = 4, hidden_dim: int = 32):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(feature_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, 1),
        )

    def forward(self, gene_features: torch.Tensor, temperature: float = 1.0) -> torch.Tensor:
        logits = self.net(gene_features).squeeze(-1)
        return logits / max(float(temperature), 1e-4)


class TrajectoryGuidedScMAE(nn.Module):
    """scMAE with an adaptive trajectory-guided gene-token masking adapter."""

    def __init__(
        self,
        input_dim: int,
        hidden_size: int = 128,
        decoder_hidden: int = 128,
        dropout: float = 0.0,
        sampler_hidden: int = 32,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + input_dim, decoder_hidden),
            nn.GELU(),
            nn.Linear(decoder_hidden, input_dim),
        )
        self.sampler = TrajectoryTokenSampler(feature_dim=4, hidden_dim=sampler_hidden)
        self.trajectory_head = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, input_dim),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def sampler_logits(self, gene_features: torch.Tensor, temperature: float = 1.0) -> torch.Tensor:
        return self.sampler(gene_features, temperature)

    def forward(self, x: torch.Tensor, gene_features: torch.Tensor | None = None, temperature: float = 1.0) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        out = {
            "embedding": z,
            "reconstruction": recon,
            "mask_logits": mask_logits,
            "trajectory_pred": self.trajectory_head(z),
        }
        if gene_features is not None:
            out["sampler_logits"] = self.sampler_logits(gene_features, temperature)
        return out
