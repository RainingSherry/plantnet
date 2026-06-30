from __future__ import annotations

import torch
import torch.nn as nn


def normalize_adjacency(adjacency: torch.Tensor) -> torch.Tensor:
    if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
        raise ValueError(f"adjacency must be [cells, cells], got {tuple(adjacency.shape)}")
    eye = torch.eye(adjacency.shape[0], device=adjacency.device, dtype=adjacency.dtype)
    a_hat = adjacency + eye
    degree = a_hat.sum(dim=1).clamp_min(1e-6)
    d_inv_sqrt = torch.rsqrt(degree)
    return d_inv_sqrt[:, None] * a_hat * d_inv_sqrt[None, :]


class DenseGraphConvLayer(nn.Module):
    def __init__(self, in_features: int, out_features: int, dropout: float) -> None:
        super().__init__()
        self.self_linear = nn.Linear(in_features, out_features, bias=False)
        self.neighbor_linear = nn.Linear(in_features, out_features, bias=False)
        self.norm = nn.LayerNorm(out_features)
        self.activation = nn.PReLU(1)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2:
            raise ValueError(f"x must be [cells, features], got {tuple(x.shape)}")
        if adjacency.shape != (x.shape[0], x.shape[0]):
            raise ValueError(f"adjacency must be [cells, cells], got {tuple(adjacency.shape)}")
        neighbor = normalize_adjacency(adjacency) @ x
        out = self.self_linear(x) + self.neighbor_linear(neighbor)
        return self.dropout(self.activation(self.norm(out)))


class GraphBarlowEncoder(nn.Module):
    def __init__(self, num_genes: int, hidden_size: int, latent_size: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or hidden_size <= 0 or latent_size <= 0:
            raise ValueError("num_genes, hidden_size, and latent_size must be positive")
        self.num_genes = int(num_genes)
        self.latent_size = int(latent_size)
        self.input_norm = nn.LayerNorm(num_genes)
        self.layer1 = DenseGraphConvLayer(num_genes, hidden_size, dropout)
        self.layer2 = DenseGraphConvLayer(hidden_size, hidden_size, dropout)
        self.layer3 = DenseGraphConvLayer(hidden_size, latent_size, dropout)
        self.skip1 = nn.Linear(num_genes, hidden_size, bias=False)
        self.skip2 = nn.Linear(num_genes, hidden_size, bias=False)

    def forward(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        x_norm = self.input_norm(x)
        h1 = self.layer1(x_norm, adjacency)
        h2 = self.layer2(h1 + self.skip1(x_norm), adjacency)
        return self.layer3(h1 + h2 + self.skip2(x_norm), adjacency)


class ProjectionHead(nn.Module):
    def __init__(self, latent_size: int, projection_size: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(latent_size, hidden_size),
            nn.BatchNorm1d(hidden_size),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, projection_size),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class GraphBarlowTwinsScMAE(nn.Module):
    """Graph Barlow Twins with a dense mini-batch cell-graph encoder."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        latent_size: int = 64,
        projection_size: int = 64,
        projector_hidden: int = 256,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.num_genes = int(num_genes)
        self.encoder = GraphBarlowEncoder(num_genes, hidden_size, latent_size, dropout)
        self.projector = ProjectionHead(latent_size, projection_size, projector_hidden, dropout)
        self.reconstruction_head = nn.Sequential(
            nn.Linear(latent_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )
        self.mask_head = nn.Linear(latent_size, num_genes)

    def forward(
        self,
        view1: torch.Tensor,
        adjacency1: torch.Tensor,
        view2: torch.Tensor,
        adjacency2: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        z1 = self.encoder(view1, adjacency1)
        z2 = self.encoder(view2, adjacency2)
        p1 = self.projector(z1)
        p2 = self.projector(z2)
        return {
            "embedding": 0.5 * (z1 + z2),
            "z1": z1,
            "z2": z2,
            "projection1": p1,
            "projection2": p2,
            "reconstruction1": self.reconstruction_head(z1),
            "reconstruction2": self.reconstruction_head(z2),
            "mask_logits1": self.mask_head(z1),
            "mask_logits2": self.mask_head(z2),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        return self.encoder(x, adjacency)
