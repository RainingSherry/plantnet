from __future__ import annotations

import math

import torch
from torch import nn


class SinusoidalTimeEmbedding(nn.Module):
    """Small continuous-time embedding for absorbing-mask diffusion steps."""

    def __init__(self, dim: int):
        super().__init__()
        self.dim = int(dim)
        self.proj = nn.Sequential(
            nn.Linear(self.dim, self.dim),
            nn.Mish(),
            nn.Linear(self.dim, self.dim),
        )

    def forward(self, t: torch.Tensor) -> torch.Tensor:
        if t.dim() == 2 and t.shape[1] == 1:
            t = t[:, 0]
        half = max(1, self.dim // 2)
        device = t.device
        freqs = torch.exp(
            torch.arange(half, device=device, dtype=torch.float32)
            * (-math.log(10000.0) / max(1, half - 1))
        )
        args = t.float().unsqueeze(1) * freqs.unsqueeze(0)
        emb = torch.cat([torch.sin(args), torch.cos(args)], dim=1)
        if emb.shape[1] < self.dim:
            emb = torch.nn.functional.pad(emb, (0, self.dim - emb.shape[1]))
        return self.proj(emb[:, : self.dim])


class ScDiVaAbsorbingScMAE(nn.Module):
    """
    ScDiVa-inspired independent scMAE variant.

    The model keeps scMAE's mask prediction and masked expression reconstruction,
    then adds a gene-specific quantile-token head trained only on absorbed
    positions. The absorbing state is represented by a learnable per-gene
    sentinel so the original expression value is removed rather than swapped
    from another cell.
    """

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        token_bins: int = 8,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.token_bins = int(token_bins)

        self.absorbing_value = nn.Parameter(torch.zeros(self.num_genes))
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
        self.time_embedding = SinusoidalTimeEmbedding(self.hidden_size)
        self.latent_norm = nn.LayerNorm(self.hidden_size)
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.expression_decoder = nn.Sequential(
            nn.Linear(self.hidden_size + self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.num_genes),
        )
        self.token_head = nn.Linear(self.hidden_size, self.num_genes * self.token_bins)

    def sample_absorbing_mask(
        self,
        x: torch.Tensor,
        t_min: float,
        t_max: float,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        batch = x.shape[0]
        t = torch.empty(batch, device=x.device).uniform_(float(t_min), float(t_max))
        mask = torch.bernoulli(t[:, None].expand_as(x)).float()

        empty_rows = mask.sum(dim=1) == 0
        if bool(empty_rows.any()):
            cols = torch.randint(0, x.shape[1], (int(empty_rows.sum()),), device=x.device)
            mask[empty_rows, cols] = 1.0

        sentinel = self.absorbing_value.to(dtype=x.dtype, device=x.device).unsqueeze(0)
        corrupted = torch.where(mask.bool(), sentinel.expand_as(x), x)
        return corrupted, mask, t

    def forward(self, x: torch.Tensor, t: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        latent = self.latent_norm(latent + self.time_embedding(t))
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.expression_decoder(torch.cat([latent, mask_logits], dim=1))
        token_logits = self.token_head(latent).view(
            x.shape[0], self.num_genes, self.token_bins
        )
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "token_logits": token_logits,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        t0 = torch.zeros(x.shape[0], device=x.device, dtype=x.dtype)
        return self.forward(x, t0)["latent"]

