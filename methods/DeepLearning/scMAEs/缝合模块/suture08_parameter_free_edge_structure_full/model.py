from __future__ import annotations

import torch
from torch import nn


class EdgeStructureReliabilityAdapter(nn.Module):
    """A light latent adapter controlled by a parameter-free graph reliability score."""

    def __init__(self, input_dim: int, latent_dim: int, hidden_dim: int = 256, adapter_weight: float = 0.06):
        super().__init__()
        self.adapter_weight = float(adapter_weight)
        self.context_encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, latent_dim),
        )
        self.delta = nn.Sequential(
            nn.LayerNorm(latent_dim),
            nn.Linear(latent_dim, latent_dim),
            nn.GELU(),
            nn.Linear(latent_dim, latent_dim),
        )

    def forward(self, base_z: torch.Tensor, structure_context_x: torch.Tensor, reliability: torch.Tensor) -> dict:
        context_z = self.context_encoder(structure_context_x)
        reliability = reliability.clamp(0.0, 1.0)
        if reliability.ndim == 1:
            reliability = reliability[:, None]
        delta = self.delta(context_z - base_z)
        latent = base_z + self.adapter_weight * reliability * delta
        return {
            "latent": latent,
            "context_z": context_z,
            "adapter_delta": delta,
            "reliability": reliability,
        }


class ParameterFreeEdgeStructureScMAE(nn.Module):
    """Independent scMAE with graph-derived parameter-free reliability control."""

    def __init__(
        self,
        input_dim: int,
        hidden_dim: int = 512,
        latent_dim: int = 32,
        dropout: float = 0.1,
        mask_prob: float = 0.4,
        adapter_weight: float = 0.06,
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
        self.reliability_adapter = EdgeStructureReliabilityAdapter(
            input_dim=input_dim,
            latent_dim=latent_dim,
            hidden_dim=hidden_dim // 2,
            adapter_weight=adapter_weight,
        )
        self.mask_predictor = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.GELU(),
            nn.Linear(hidden_dim // 2, input_dim),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hidden_dim // 2),
            nn.LayerNorm(hidden_dim // 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim // 2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, input_dim),
        )

    def corrupt(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < self.mask_prob).float()
        return x * (1.0 - mask), mask

    def encode_with_reliability(
        self,
        x: torch.Tensor,
        structure_context_x: torch.Tensor,
        reliability: torch.Tensor,
    ) -> dict:
        base_z = self.encoder(x)
        out = self.reliability_adapter(base_z, structure_context_x, reliability)
        out["base_latent"] = base_z
        return out

    def forward(self, x: torch.Tensor, structure_context_x: torch.Tensor, reliability: torch.Tensor) -> dict:
        corrupted, mask = self.corrupt(x)
        out = self.encode_with_reliability(corrupted, structure_context_x, reliability)
        latent = out["latent"]
        return {
            "latent": latent,
            "base_latent": out["base_latent"],
            "context_z": out["context_z"],
            "adapter_delta": out["adapter_delta"],
            "reliability": out["reliability"],
            "reconstruction": self.decoder(latent),
            "mask_logits": self.mask_predictor(latent),
            "mask": mask,
        }

    def feature(
        self,
        x: torch.Tensor,
        structure_context_x: torch.Tensor | None = None,
        reliability: torch.Tensor | None = None,
    ) -> torch.Tensor:
        if structure_context_x is None or reliability is None:
            return self.encoder(x)
        return self.encode_with_reliability(x, structure_context_x, reliability)["latent"]
