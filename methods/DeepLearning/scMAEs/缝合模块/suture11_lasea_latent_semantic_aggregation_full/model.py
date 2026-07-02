from __future__ import annotations

import torch
from torch import nn


class LatentSemanticAggregation(nn.Module):
    """1D latent rewrite of LaSEA multi-scale semantic extraction."""

    def __init__(self, latent_dim: int, semantic_weight: float = 0.05):
        super().__init__()
        self.semantic_weight = float(semantic_weight)
        self.norm = nn.LayerNorm(latent_dim)
        branch_channels = 4
        self.branches = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Conv1d(1, branch_channels, kernel_size=3, padding=d, dilation=d),
                    nn.GELU(),
                    nn.Conv1d(branch_channels, branch_channels, kernel_size=1),
                    nn.GELU(),
                )
                for d in (1, 2, 3, 4)
            ]
        )
        self.fusion = nn.Sequential(
            nn.Conv1d(branch_channels * 4, branch_channels, kernel_size=1),
            nn.GELU(),
            nn.Conv1d(branch_channels, 1, kernel_size=1),
        )
        self.semantic_attention = nn.Sequential(
            nn.Linear(latent_dim, max(8, latent_dim // 2)),
            nn.GELU(),
            nn.Linear(max(8, latent_dim // 2), latent_dim),
            nn.Sigmoid(),
        )
        self.semantic_gate = nn.Sequential(
            nn.Linear(latent_dim * 2, max(8, latent_dim // 2)),
            nn.GELU(),
            nn.Linear(max(8, latent_dim // 2), 1),
            nn.Sigmoid(),
        )
        self.summary = nn.Sequential(nn.LayerNorm(latent_dim), nn.Linear(latent_dim, latent_dim), nn.GELU(), nn.Linear(latent_dim, latent_dim))

    def semantic_summary(self, z: torch.Tensor) -> torch.Tensor:
        return self.summary(self.norm(z))

    def forward(self, base_z: torch.Tensor) -> dict:
        z = self.norm(base_z)
        seq = z[:, None, :]
        multi = torch.cat([branch(seq) for branch in self.branches], dim=1)
        fused = self.fusion(multi).squeeze(1)
        attention = self.semantic_attention(z)
        gate = self.semantic_gate(torch.cat([z, torch.abs(fused)], dim=1))
        delta = attention * fused
        latent = base_z + self.semantic_weight * gate * delta
        return {
            "latent": latent,
            "semantic_delta": delta,
            "semantic_attention": attention,
            "semantic_gate": gate,
            "semantic_summary": self.semantic_summary(latent),
        }


class LaSEALatentSemanticScMAE(nn.Module):
    """Independent scMAE with latent-aware semantic extraction and aggregation."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        semantic_weight: float = 0.05,
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
        self.semantic = LatentSemanticAggregation(latent_dim, semantic_weight)
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

    def encode_semantic(self, x: torch.Tensor) -> dict:
        base_z = self.encoder(x)
        out = self.semantic(base_z)
        out["base_latent"] = base_z
        return out

    def forward(self, x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        out = self.encode_semantic(corrupted)
        with torch.no_grad():
            clean_z = self.encoder(x)
            clean_summary = self.semantic.semantic_summary(clean_z)
        latent = out["latent"]
        return {
            "latent": latent,
            "base_latent": out["base_latent"],
            "semantic_delta": out["semantic_delta"],
            "semantic_attention": out["semantic_attention"],
            "semantic_gate": out["semantic_gate"],
            "semantic_summary": out["semantic_summary"],
            "clean_semantic_summary": clean_summary,
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode_semantic(x)["latent"]
