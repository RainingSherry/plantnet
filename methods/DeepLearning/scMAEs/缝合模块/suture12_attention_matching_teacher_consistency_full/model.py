from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class CleanMaskedAttentionMatcher(nn.Module):
    """Attention matching between masked latent and detached clean latent."""

    def __init__(self, latent_dim: int, hidden_dim: int = 128, match_weight: float = 0.035):
        super().__init__()
        self.match_weight = float(match_weight)
        self.masked_proj = nn.Sequential(nn.LayerNorm(latent_dim), nn.Linear(latent_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, latent_dim))
        self.clean_proj = nn.Sequential(nn.LayerNorm(latent_dim), nn.Linear(latent_dim, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, latent_dim))
        self.fusion = nn.Sequential(nn.LayerNorm(latent_dim * 2), nn.Linear(latent_dim * 2, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, latent_dim))
        self.gate = nn.Sequential(nn.LayerNorm(latent_dim * 3), nn.Linear(latent_dim * 3, hidden_dim), nn.GELU(), nn.Linear(hidden_dim, 1), nn.Sigmoid())

    def forward(self, masked_z: torch.Tensor, clean_z: torch.Tensor) -> dict:
        clean_z = clean_z.detach()
        masked_proj = self.masked_proj(masked_z)
        clean_proj = self.clean_proj(clean_z)
        similarity = F.cosine_similarity(masked_proj, clean_proj, dim=1, eps=1e-6).unsqueeze(1)
        attention = torch.sigmoid(similarity)
        stable_masked = attention * masked_proj
        stable_clean = attention * clean_proj
        fused = self.fusion(torch.cat([stable_masked, stable_clean], dim=1))
        gate = self.gate(torch.cat([masked_z, clean_z, clean_z - masked_z], dim=1))
        latent = masked_z + self.match_weight * gate * fused
        return {
            "latent": latent,
            "masked_proj": masked_proj,
            "clean_proj": clean_proj,
            "match_similarity": similarity,
            "match_attention": attention,
            "match_gate": gate,
            "matched_delta": fused,
        }


class AttentionMatchingTeacherScMAE(nn.Module):
    """Independent scMAE with weak clean/masked latent attention matching."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        match_weight: float = 0.035,
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
        self.matcher = CleanMaskedAttentionMatcher(latent_dim, hidden_dim // 4, match_weight)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, input_dim),
        )
        self.mask_predictor = nn.Sequential(nn.Linear(latent_dim, hidden_dim // 2), nn.GELU(), nn.Linear(hidden_dim // 2, input_dim))

    def corrupt(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < self.mask_prob).float()
        return x * (1.0 - mask), mask

    def encode_matched(self, x: torch.Tensor, clean_x: torch.Tensor | None = None) -> dict:
        masked_z = self.encoder(x)
        clean_z = self.encoder(clean_x) if clean_x is not None else masked_z.detach()
        out = self.matcher(masked_z, clean_z)
        out["base_latent"] = masked_z
        out["clean_latent"] = clean_z.detach()
        return out

    def forward(self, x: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        out = self.encode_matched(corrupted, x)
        latent = out["latent"]
        return {
            "latent": latent,
            "base_latent": out["base_latent"],
            "clean_latent": out["clean_latent"],
            "masked_proj": out["masked_proj"],
            "clean_proj": out["clean_proj"],
            "match_similarity": out["match_similarity"],
            "match_attention": out["match_attention"],
            "match_gate": out["match_gate"],
            "matched_delta": out["matched_delta"],
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode_matched(x, x)["latent"]
