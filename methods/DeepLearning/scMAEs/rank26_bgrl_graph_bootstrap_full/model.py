from __future__ import annotations

import copy
import math

import torch
from torch import nn
import torch.nn.functional as F


class ShallowGraphEncoder(nn.Module):
    """Residual single-hop graph adapter for scRNA cell KNN context."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.cell_encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, 256),
            nn.LayerNorm(256),
            nn.PReLU(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.PReLU(),
        )
        self.neighbor_encoder = nn.Sequential(
            nn.Linear(num_genes, 256),
            nn.LayerNorm(256),
            nn.PReLU(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.PReLU(),
        )
        self.gate = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size),
            nn.Sigmoid(),
        )
        self.out_norm = nn.LayerNorm(hidden_size)

    def forward(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> torch.Tensor:
        base = self.cell_encoder(x)
        if neighbor_x is None:
            return self.out_norm(base)
        if neighbor_x.dim() == 3:
            neighbor_x = neighbor_x.mean(dim=1)
        neigh = self.neighbor_encoder(neighbor_x)
        gate = self.gate(torch.cat([base, neigh], dim=1))
        return self.out_norm(base + gate * (neigh - base))


class BGRLGraphScMAE(nn.Module):
    """scMAE with BGRL-style online/EMA target graph bootstrapping."""

    def __init__(self, num_genes: int, hidden_size: int = 128, predictor_hidden: int = 256, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.online_encoder = ShallowGraphEncoder(num_genes, hidden_size, dropout)
        self.target_encoder = copy.deepcopy(self.online_encoder)
        for param in self.target_encoder.parameters():
            param.requires_grad = False
        self.predictor = nn.Sequential(
            nn.Linear(hidden_size, predictor_hidden),
            nn.LayerNorm(predictor_hidden),
            nn.PReLU(),
            nn.Linear(predictor_hidden, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, 256),
            nn.PReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, num_genes),
        )
        self.edge_head = nn.Sequential(
            nn.Linear(hidden_size * 3, hidden_size),
            nn.PReLU(),
            nn.Linear(hidden_size, 1),
        )

    def encode_online(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> torch.Tensor:
        return self.online_encoder(x, neighbor_x)

    @torch.no_grad()
    def encode_target(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> torch.Tensor:
        return self.target_encoder(x, neighbor_x)

    def forward(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        latent = self.encode_online(x, neighbor_x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction}

    def bootstrap(self, x1: torch.Tensor, n1: torch.Tensor, x2: torch.Tensor, n2: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        online = self.predictor(self.encode_online(x1, n1))
        with torch.no_grad():
            target = self.encode_target(x2, n2).detach()
        return online, target

    def edge_logits(self, z: torch.Tensor, pos_z: torch.Tensor, neg_z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        pos = self.edge_head(torch.cat([z, pos_z, torch.abs(z - pos_z)], dim=1)).squeeze(-1)
        neg = self.edge_head(torch.cat([z, neg_z, torch.abs(z - neg_z)], dim=1)).squeeze(-1)
        return pos, neg

    @torch.no_grad()
    def update_target(self, momentum: float) -> None:
        if not 0.0 <= momentum <= 1.0:
            raise ValueError(f"EMA momentum must be in [0, 1], got {momentum}")
        for online, target in zip(self.online_encoder.parameters(), self.target_encoder.parameters()):
            target.data.mul_(momentum).add_(online.data, alpha=1.0 - momentum)

    @torch.no_grad()
    def feature(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> torch.Tensor:
        return self.encode_online(x, neighbor_x)

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    def feature_drop(self, x: torch.Tensor, drop_prob: float) -> torch.Tensor:
        if drop_prob <= 0:
            return x
        keep = (torch.rand(x.shape[-1], device=x.device, dtype=x.dtype) >= float(drop_prob)).float()
        if keep.sum() == 0:
            keep[torch.randint(0, x.shape[-1], (1,), device=x.device)] = 1.0
        view_shape = [1] * x.dim()
        view_shape[-1] = x.shape[-1]
        return x * keep.view(*view_shape)

    @staticmethod
    def cosine_momentum(base: float, step: int, total_steps: int) -> float:
        if total_steps <= 1:
            return 1.0
        return 1.0 - (1.0 - base) * (math.cos(math.pi * step / total_steps) + 1.0) / 2.0
