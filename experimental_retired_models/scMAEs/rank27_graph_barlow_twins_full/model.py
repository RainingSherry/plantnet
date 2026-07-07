from __future__ import annotations

import torch
from torch import nn


class ShallowGraphEncoder(nn.Module):
    """Symmetric shallow graph encoder for Graph Barlow Twins on cell KNN graphs."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
        self.cell_encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, 256),
            nn.BatchNorm1d(256, momentum=0.01),
            nn.PReLU(),
            nn.Linear(256, hidden_size),
        )
        self.neighbor_encoder = nn.Sequential(
            nn.Linear(num_genes, 256),
            nn.BatchNorm1d(256, momentum=0.01),
            nn.PReLU(),
            nn.Linear(256, hidden_size),
        )
        self.gate = nn.Sequential(nn.Linear(hidden_size * 2, hidden_size), nn.Sigmoid())

    def forward(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> torch.Tensor:
        base = self.cell_encoder(x)
        if neighbor_x is None:
            return base
        if neighbor_x.dim() == 3:
            neighbor_x = neighbor_x.mean(dim=1)
        neigh = self.neighbor_encoder(neighbor_x)
        gate = self.gate(torch.cat([base, neigh], dim=1))
        return base + gate * (neigh - base)


class GraphBarlowTwinsScMAE(nn.Module):
    """scMAE with symmetric Graph Barlow Twins redundancy-reduction views."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.encoder = ShallowGraphEncoder(num_genes, hidden_size, dropout)
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

    def encode(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> torch.Tensor:
        return self.encoder(x, neighbor_x)

    def forward(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        latent = self.encode(x, neighbor_x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction}

    def edge_logits(self, z: torch.Tensor, pos_z: torch.Tensor, neg_z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        pos = self.edge_head(torch.cat([z, pos_z, torch.abs(z - pos_z)], dim=1)).squeeze(-1)
        neg = self.edge_head(torch.cat([z, neg_z, torch.abs(z - neg_z)], dim=1)).squeeze(-1)
        return pos, neg

    @torch.no_grad()
    def feature(self, x: torch.Tensor, neighbor_x: torch.Tensor | None = None) -> torch.Tensor:
        return self.encode(x, neighbor_x)

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
