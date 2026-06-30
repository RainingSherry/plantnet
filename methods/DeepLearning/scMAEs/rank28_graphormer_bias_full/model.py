from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F


def floyd_warshall_distance(adjacency: torch.Tensor, max_distance: int) -> torch.Tensor:
    if adjacency.ndim != 2 or adjacency.shape[0] != adjacency.shape[1]:
        raise ValueError(f"adjacency must be [cells, cells], got {tuple(adjacency.shape)}")
    cells = adjacency.shape[0]
    inf = int(max_distance) + 1
    dist = torch.full((cells, cells), inf, device=adjacency.device, dtype=torch.long)
    dist[adjacency.bool()] = 1
    dist.fill_diagonal_(0)
    for k in range(cells):
        dist = torch.minimum(dist, dist[:, k : k + 1] + dist[k : k + 1, :])
    return dist.clamp_max(inf)


class GraphormerSelfAttention(nn.Module):
    def __init__(self, hidden_size: int, num_heads: int, dropout: float) -> None:
        super().__init__()
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")
        self.hidden_size = int(hidden_size)
        self.num_heads = int(num_heads)
        self.head_dim = hidden_size // num_heads
        self.scale = self.head_dim ** -0.5
        self.q_proj = nn.Linear(hidden_size, hidden_size)
        self.k_proj = nn.Linear(hidden_size, hidden_size)
        self.v_proj = nn.Linear(hidden_size, hidden_size)
        self.out_proj = nn.Linear(hidden_size, hidden_size)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, attn_bias: torch.Tensor) -> torch.Tensor:
        if x.ndim != 3:
            raise ValueError(f"x must be [batch, tokens, hidden], got {tuple(x.shape)}")
        if attn_bias.shape != (x.shape[0], self.num_heads, x.shape[1], x.shape[1]):
            raise ValueError(f"attn_bias shape {tuple(attn_bias.shape)} incompatible with x {tuple(x.shape)}")
        batch, tokens, _ = x.shape
        q = self.q_proj(x).view(batch, tokens, self.num_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(batch, tokens, self.num_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(batch, tokens, self.num_heads, self.head_dim).transpose(1, 2)
        logits = torch.matmul(q, k.transpose(-2, -1)) * self.scale + attn_bias
        weights = torch.softmax(logits, dim=-1)
        weights = self.dropout(weights)
        out = torch.matmul(weights, v).transpose(1, 2).reshape(batch, tokens, self.hidden_size)
        return self.out_proj(out)


class GraphormerBlock(nn.Module):
    def __init__(self, hidden_size: int, num_heads: int, dropout: float, mlp_ratio: int) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_size)
        self.attn = GraphormerSelfAttention(hidden_size, num_heads, dropout)
        self.drop = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(hidden_size)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * mlp_ratio, hidden_size),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor, attn_bias: torch.Tensor) -> torch.Tensor:
        x = x + self.drop(self.attn(self.norm1(x), attn_bias))
        return x + self.mlp(self.norm2(x))


class GraphormerCellEncoder(nn.Module):
    """Graphormer-style graph-biased attention over a mini-batch cell graph."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        depth: int = 3,
        num_heads: int = 4,
        max_degree: int = 64,
        max_distance: int = 8,
        dropout: float = 0.1,
        mlp_ratio: int = 4,
    ) -> None:
        super().__init__()
        if num_genes <= 0 or hidden_size <= 0 or depth <= 0:
            raise ValueError("num_genes, hidden_size, and depth must be positive")
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.num_heads = int(num_heads)
        self.max_degree = int(max_degree)
        self.max_distance = int(max_distance)
        self.node_projection = nn.Sequential(
            nn.LayerNorm(num_genes),
            nn.Linear(num_genes, hidden_size),
        )
        self.in_degree_embedding = nn.Embedding(max_degree + 2, hidden_size, padding_idx=0)
        self.out_degree_embedding = nn.Embedding(max_degree + 2, hidden_size, padding_idx=0)
        self.spatial_pos_embedding = nn.Embedding(max_distance + 3, num_heads, padding_idx=0)
        self.edge_embedding = nn.Embedding(3, num_heads, padding_idx=0)
        self.graph_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.graph_token_distance = nn.Parameter(torch.zeros(num_heads))
        self.blocks = nn.ModuleList([GraphormerBlock(hidden_size, num_heads, dropout, mlp_ratio) for _ in range(depth)])
        self.norm = nn.LayerNorm(hidden_size)
        self.dropout = nn.Dropout(dropout)
        nn.init.trunc_normal_(self.graph_token, std=0.02)
        nn.init.normal_(self.graph_token_distance, std=0.02)

    def build_bias(self, adjacency: torch.Tensor) -> torch.Tensor:
        cells = adjacency.shape[0]
        degrees = adjacency.sum(dim=1).long().clamp(0, self.max_degree) + 1
        dist = floyd_warshall_distance(adjacency, self.max_distance)
        spatial_ids = (dist + 1).clamp(0, self.max_distance + 2)
        edge_ids = adjacency.long().clamp(0, 1) + 1
        spatial_bias = self.spatial_pos_embedding(spatial_ids).permute(2, 0, 1)
        edge_bias = self.edge_embedding(edge_ids).permute(2, 0, 1)
        bias = torch.zeros(self.num_heads, cells + 1, cells + 1, device=adjacency.device, dtype=spatial_bias.dtype)
        bias[:, 1:, 1:] = spatial_bias + edge_bias
        token_bias = self.graph_token_distance.view(self.num_heads, 1)
        bias[:, 1:, 0] = bias[:, 1:, 0] + token_bias
        bias[:, 0, :] = bias[:, 0, :] + token_bias
        return bias.unsqueeze(0), degrees

    def forward(self, x: torch.Tensor, adjacency: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [cells, {self.num_genes}], got {tuple(x.shape)}")
        if adjacency.shape != (x.shape[0], x.shape[0]):
            raise ValueError(f"adjacency must be [cells, cells], got {tuple(adjacency.shape)}")
        attn_bias, degrees = self.build_bias(adjacency)
        node = self.node_projection(x)
        node = node + self.in_degree_embedding(degrees) + self.out_degree_embedding(degrees)
        graph_token = self.graph_token.expand(1, -1, -1).to(dtype=node.dtype, device=node.device)
        h = torch.cat([graph_token, node.unsqueeze(0)], dim=1)
        h = self.dropout(h)
        for block in self.blocks:
            h = block(h, attn_bias.to(dtype=h.dtype))
        h = self.norm(h)
        return h[:, 0].squeeze(0), h[:, 1:].squeeze(0)


class GraphormerBiasScMAE(nn.Module):
    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        depth: int = 3,
        num_heads: int = 4,
        max_degree: int = 64,
        max_distance: int = 8,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.num_genes = int(num_genes)
        self.encoder = GraphormerCellEncoder(
            num_genes,
            hidden_size=hidden_size,
            depth=depth,
            num_heads=num_heads,
            max_degree=max_degree,
            max_distance=max_distance,
            dropout=dropout,
        )
        self.reconstruction_head = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, num_genes),
        )
        self.mask_head = nn.Linear(hidden_size, num_genes)
        self.graph_decoder_scale = hidden_size ** -0.5

    def forward(self, x: torch.Tensor, adjacency: torch.Tensor) -> dict[str, torch.Tensor]:
        graph_embedding, node_embedding = self.encoder(x, adjacency)
        reconstruction = self.reconstruction_head(node_embedding)
        mask_logits = self.mask_head(node_embedding)
        graph_logits = torch.sigmoid((node_embedding @ node_embedding.T) * self.graph_decoder_scale)
        return {
            "embedding": node_embedding,
            "graph_embedding": graph_embedding,
            "reconstruction": reconstruction,
            "mask_logits": mask_logits,
            "adjacency_reconstruction": graph_logits,
            "adjacency": adjacency,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor, adjacency: torch.Tensor) -> torch.Tensor:
        _, node_embedding = self.encoder(x, adjacency)
        return node_embedding
