from __future__ import annotations

import torch
import torch.nn as nn


def build_knn_adjacency(x: torch.Tensor, k: int, self_loop: bool = True) -> torch.Tensor:
    if x.ndim != 2:
        raise ValueError(f"x must be [batch, genes], got {tuple(x.shape)}")
    n = x.shape[0]
    if n == 0:
        raise ValueError("empty batch")
    k_eff = min(max(1, int(k)), max(1, n - 1))
    with torch.no_grad():
        x_norm = torch.nn.functional.normalize(x, dim=1)
        sim = x_norm @ x_norm.t()
        sim.fill_diagonal_(-float("inf"))
        idx = sim.topk(k_eff, dim=1).indices
        adj = torch.zeros((n, n), dtype=x.dtype, device=x.device)
        rows = torch.arange(n, device=x.device).view(-1, 1).expand_as(idx)
        adj[rows, idx] = 1.0
        adj = torch.maximum(adj, adj.t())
        if self_loop:
            adj.fill_diagonal_(1.0)
        degree = adj.sum(dim=1).clamp_min(1.0)
        inv_sqrt = degree.rsqrt()
        adj = inv_sqrt.view(-1, 1) * adj * inv_sqrt.view(1, -1)
    return adj


class GraphConvBlock(nn.Module):
    def __init__(self, hidden_size: int, dropout: float) -> None:
        super().__init__()
        self.norm = nn.LayerNorm(hidden_size)
        self.linear = nn.Linear(hidden_size, hidden_size)
        self.dropout = nn.Dropout(dropout)
        self.act = nn.GELU()

    def forward(self, x: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        if adj.ndim != 2 or adj.shape[0] != adj.shape[1] or adj.shape[0] != x.shape[0]:
            raise ValueError(f"adj must be [batch, batch], got {tuple(adj.shape)} for x {tuple(x.shape)}")
        h = adj @ self.norm(x)
        h = self.linear(h)
        h = self.dropout(self.act(h))
        return x + h


class JOAOScMAEGraphEncoder(nn.Module):
    """GraphCL/JOAO style encoder over mini-batch cell graphs.

    Each mini-batch is treated as a sparse cell graph built from expression
    similarities. Two augmented graph views are encoded by the same GCN stack and
    optimized with NT-Xent. A reconstruction branch keeps the scMAE masked
    expression objective active.
    """

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        depth: int = 3,
        projection_size: int = 128,
        dropout: float = 0.1,
        knn_k: int = 15,
    ) -> None:
        super().__init__()
        if num_genes <= 0:
            raise ValueError("num_genes must be positive")
        if depth <= 0:
            raise ValueError("depth must be positive")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.knn_k = int(knn_k)
        self.stem = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
        )
        self.layers = nn.ModuleList([GraphConvBlock(hidden_size, dropout) for _ in range(depth)])
        self.final_norm = nn.LayerNorm(hidden_size)
        self.projector = nn.Sequential(
            nn.Linear(hidden_size, projection_size),
            nn.BatchNorm1d(projection_size),
            nn.ReLU(inplace=True),
            nn.Linear(projection_size, projection_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, num_genes),
        )

    def encode(self, x: torch.Tensor, adj: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        h = self.stem(x)
        for layer in self.layers:
            h = layer(h, adj)
        return self.final_norm(h)

    def forward(self, x: torch.Tensor, adj: torch.Tensor) -> dict[str, torch.Tensor]:
        emb = self.encode(x, adj)
        mask_logits = self.mask_predictor(emb)
        reconstruction = self.decoder(torch.cat([emb, mask_logits], dim=1))
        return {
            "embedding": emb,
            "projection": self.projector(emb),
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        adj = build_knn_adjacency(x, self.knn_k)
        return self.encode(x, adj)
