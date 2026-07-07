from __future__ import annotations

import torch
from torch import nn


class ConditionalMultiMaskScMAE(nn.Module):
    """scMAE with Conditional-MAE style input and latent-stage masking."""

    def __init__(self, input_dim: int, hidden_size: int = 128, decoder_hidden: int = 128, latent_dim: int = 256, dropout: float = 0.0):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.latent_dim = int(latent_dim)
        self.encoder_stage1 = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(input_dim, latent_dim),
            nn.LayerNorm(latent_dim),
            nn.GELU(),
        )
        self.latent_mask_logits = nn.Sequential(nn.Linear(latent_dim, latent_dim), nn.LayerNorm(latent_dim))
        self.latent_mask_token = nn.Parameter(torch.zeros(1, latent_dim))
        self.encoder_stage2 = nn.Sequential(
            nn.Linear(latent_dim, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(nn.Linear(hidden_size + input_dim, decoder_hidden), nn.GELU(), nn.Linear(decoder_hidden, input_dim))
        self.condition_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.GELU(), nn.Linear(hidden_size, latent_dim))
        nn.init.normal_(self.latent_mask_token, std=0.02)

    def encode(self, x: torch.Tensor, latent_mask_ratio: float = 0.0, force_no_latent_mask: bool = True) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h1 = self.encoder_stage1(x)
        logits = self.latent_mask_logits(h1)
        if force_no_latent_mask or latent_mask_ratio <= 0:
            latent_mask = torch.zeros_like(h1)
            h1_masked = h1
        else:
            score = torch.sigmoid(logits.detach())
            score = score / score.mean(dim=1, keepdim=True).clamp_min(1e-6)
            probs = (float(latent_mask_ratio) * score).clamp(0.02, 0.85)
            latent_mask = (torch.rand_like(h1) < probs).float()
            h1_masked = h1 * (1.0 - latent_mask) + self.latent_mask_token.expand_as(h1) * latent_mask
        z = self.encoder_stage2(h1_masked)
        return z, logits, latent_mask

    def forward(self, x: torch.Tensor, latent_mask_ratio: float = 0.0) -> dict[str, torch.Tensor]:
        z, latent_logits, latent_mask = self.encode(x, latent_mask_ratio, force_no_latent_mask=False)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {
            "embedding": z,
            "reconstruction": recon,
            "mask_logits": mask_logits,
            "latent_mask_logits": latent_logits,
            "latent_mask": latent_mask,
            "condition_pred": self.condition_head(z),
        }

    @torch.no_grad()
    def embed(self, x: torch.Tensor) -> torch.Tensor:
        z, _, _ = self.encode(x, latent_mask_ratio=0.0, force_no_latent_mask=True)
        return z
