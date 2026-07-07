from __future__ import annotations

import torch
from torch import nn


class GroundGANCausalMaskScMAE(nn.Module):
    """scMAE with a GRouNdGAN-inspired masked causal dependency adapter."""

    def __init__(
        self,
        input_dim: int,
        regulator_k: int,
        hidden_size: int = 128,
        decoder_hidden: int = 128,
        gene_dim: int = 48,
        causal_hidden: int = 128,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.input_dim = int(input_dim)
        self.regulator_k = int(regulator_k)
        self.hidden_size = int(hidden_size)
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
        self.target_embedding = nn.Embedding(input_dim, gene_dim)
        self.regulator_encoder = nn.Sequential(nn.Linear(regulator_k, gene_dim), nn.LayerNorm(gene_dim), nn.Mish(inplace=True))
        self.context_projection = nn.Sequential(nn.Linear(hidden_size, gene_dim), nn.LayerNorm(gene_dim), nn.Mish(inplace=True))
        self.causal_predictor = nn.Sequential(
            nn.Linear(gene_dim * 3, causal_hidden),
            nn.LayerNorm(causal_hidden),
            nn.Mish(inplace=True),
            nn.Linear(causal_hidden, causal_hidden),
            nn.Mish(inplace=True),
            nn.Linear(causal_hidden, 1),
        )

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {"embedding": z, "reconstruction": recon, "mask_logits": mask_logits}

    def causal_predict(self, z: torch.Tensor, target_idx: torch.Tensor, regulator_values: torch.Tensor, regulator_weights: torch.Tensor) -> torch.Tensor:
        weighted_values = regulator_values * regulator_weights[None, :, :]
        reg = self.regulator_encoder(weighted_values)
        tgt = self.target_embedding(target_idx)[None, :, :].expand(z.shape[0], -1, -1)
        ctx = self.context_projection(z)[:, None, :].expand(-1, target_idx.shape[0], -1)
        return self.causal_predictor(torch.cat([reg, tgt, ctx], dim=-1)).squeeze(-1)
