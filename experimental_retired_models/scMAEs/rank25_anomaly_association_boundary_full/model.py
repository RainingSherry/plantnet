from __future__ import annotations

import math

import torch
from torch import nn
import torch.nn.functional as F


def _pad_to_patch(x: torch.Tensor, total_genes: int) -> torch.Tensor:
    if x.shape[1] == total_genes:
        return x
    return F.pad(x, (0, total_genes - x.shape[1]))


class AnomalyAssociationBlock(nn.Module):
    """Gene-module adaptation of Anomaly Transformer's series/prior attention."""

    def __init__(self, hidden_size: int, n_heads: int, dropout: float):
        super().__init__()
        if hidden_size % n_heads != 0:
            raise ValueError("hidden_size must be divisible by n_heads")
        self.hidden_size = int(hidden_size)
        self.n_heads = int(n_heads)
        self.head_dim = self.hidden_size // self.n_heads
        self.q_proj = nn.Linear(hidden_size, hidden_size)
        self.k_proj = nn.Linear(hidden_size, hidden_size)
        self.v_proj = nn.Linear(hidden_size, hidden_size)
        self.out_proj = nn.Linear(hidden_size, hidden_size)
        self.sigma_proj = nn.Linear(hidden_size, n_heads)
        self.norm1 = nn.LayerNorm(hidden_size)
        self.norm2 = nn.LayerNorm(hidden_size)
        self.ffn = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, hidden_size),
        )
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        bsz, n_tokens, _ = x.shape
        h = self.norm1(x)
        q = self.q_proj(h).view(bsz, n_tokens, self.n_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(h).view(bsz, n_tokens, self.n_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(h).view(bsz, n_tokens, self.n_heads, self.head_dim).transpose(1, 2)
        scores = torch.matmul(q, k.transpose(-2, -1)) / math.sqrt(self.head_dim)
        series = torch.softmax(scores, dim=-1).clamp_min(1e-8)

        pos = torch.arange(n_tokens, device=x.device, dtype=x.dtype)
        dist2 = (pos[None, :] - pos[:, None]).pow(2.0)
        sigma = F.softplus(self.sigma_proj(h)).transpose(1, 2).unsqueeze(-1) + 1e-3
        prior_logits = -dist2.view(1, 1, n_tokens, n_tokens) / (2.0 * sigma.pow(2.0))
        prior = torch.softmax(prior_logits, dim=-1).clamp_min(1e-8)

        context = torch.matmul(series, v).transpose(1, 2).contiguous().view(bsz, n_tokens, self.hidden_size)
        x = x + self.dropout(self.out_proj(context))
        x = x + self.dropout(self.ffn(self.norm2(x)))
        discrepancy = 0.5 * (
            (series * (series.log() - prior.log())).sum(dim=-1)
            + (prior * (prior.log() - series.log())).sum(dim=-1)
        ).mean(dim=1)
        return x, series, prior, discrepancy


class AnomalyAssociationScMAE(nn.Module):
    """scMAE with gene-module association discrepancy as a boundary-risk signal."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 25,
        hidden_size: int = 128,
        n_heads: int = 4,
        depth: int = 2,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.n_patches = int(math.ceil(self.num_genes / self.patch_size))
        self.total_genes = self.n_patches * self.patch_size
        self.hidden_size = int(hidden_size)
        self.patch_embed = nn.Sequential(
            nn.Linear(self.patch_size, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.pos_embed = nn.Parameter(torch.zeros(1, self.n_patches, hidden_size))
        self.blocks = nn.ModuleList([AnomalyAssociationBlock(hidden_size, n_heads, dropout) for _ in range(depth)])
        self.pool_norm = nn.LayerNorm(hidden_size)
        self.mask_predictor = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, self.num_genes),
        )
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + self.num_genes, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, self.num_genes),
        )
        self.boundary_head = nn.Sequential(
            nn.Linear(hidden_size + self.n_patches, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, 1),
        )
        nn.init.trunc_normal_(self.pos_embed, std=0.02)

    def encode_tokens(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        x_pad = _pad_to_patch(x, self.total_genes)
        patches = x_pad.view(x.shape[0], self.n_patches, self.patch_size)
        tokens = self.patch_embed(patches) + self.pos_embed
        discrepancies = []
        for block in self.blocks:
            tokens, _, _, disc = block(tokens)
            discrepancies.append(disc)
        assoc_profile = torch.stack(discrepancies, dim=0).mean(dim=0)
        latent = self.pool_norm(tokens.mean(dim=1))
        return latent, assoc_profile

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        latent, _ = self.encode_tokens(x)
        return latent

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent, assoc_profile = self.encode_tokens(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        boundary_logit = self.boundary_head(torch.cat([latent, assoc_profile], dim=1)).squeeze(-1)
        return {
            "latent": latent,
            "assoc_profile": assoc_profile,
            "association_risk": assoc_profile.mean(dim=1),
            "boundary_logit": boundary_logit,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode(x)

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask
