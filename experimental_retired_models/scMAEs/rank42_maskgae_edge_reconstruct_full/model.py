from __future__ import annotations

import torch
from torch import nn


class ShallowMaskedGraphEncoder(nn.Module):
    """Single-hop residual graph encoder for partially visible cell KNN edges."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
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
        self.edge_gate = nn.Sequential(nn.Linear(hidden_size * 2, hidden_size), nn.Sigmoid())
        self.out_norm = nn.LayerNorm(hidden_size)

    def forward(self, x: torch.Tensor, visible_neighbors: torch.Tensor | None) -> torch.Tensor:
        base = self.cell_encoder(x)
        if visible_neighbors is None:
            return self.out_norm(base)
        if visible_neighbors.dim() == 3:
            neighbor_x = visible_neighbors.mean(dim=1)
        else:
            neighbor_x = visible_neighbors
        neigh = self.neighbor_encoder(neighbor_x)
        gate = self.edge_gate(torch.cat([base, neigh], dim=1))
        return self.out_norm(base + gate * (neigh - base))


class MaskGAEEdgeScMAE(nn.Module):
    """scMAE with MaskGAE-style masked edge and degree reconstruction."""

    def __init__(self, num_genes: int, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.encoder = ShallowMaskedGraphEncoder(num_genes, hidden_size, dropout)
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, 256),
            nn.PReLU(),
            nn.Dropout(dropout),
            nn.Linear(256, num_genes),
        )
        self.edge_decoder = nn.Sequential(
            nn.Linear(hidden_size * 3, decoder_hidden),
            nn.PReLU(),
            nn.Dropout(dropout),
            nn.Linear(decoder_hidden, 1),
        )
        self.degree_decoder = nn.Sequential(
            nn.Linear(hidden_size, decoder_hidden),
            nn.PReLU(),
            nn.Dropout(dropout),
            nn.Linear(decoder_hidden, 1),
        )

    def encode(self, x: torch.Tensor, visible_neighbors: torch.Tensor | None = None) -> torch.Tensor:
        return self.encoder(x, visible_neighbors)

    def forward(self, x: torch.Tensor, visible_neighbors: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        latent = self.encode(x, visible_neighbors)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        degree = self.degree_decoder(latent).squeeze(-1)
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction, "degree_pred": degree}

    def edge_logits(self, z: torch.Tensor, pos_z: torch.Tensor, neg_z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        pos_feat = torch.cat([z, pos_z, torch.abs(z - pos_z)], dim=1)
        neg_feat = torch.cat([z, neg_z, torch.abs(z - neg_z)], dim=1)
        return self.edge_decoder(pos_feat).squeeze(-1), self.edge_decoder(neg_feat).squeeze(-1)

    @torch.no_grad()
    def feature(self, x: torch.Tensor, visible_neighbors: torch.Tensor | None = None) -> torch.Tensor:
        return self.encode(x, visible_neighbors)

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask
