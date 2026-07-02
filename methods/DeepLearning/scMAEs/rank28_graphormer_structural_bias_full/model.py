from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class GraphormerBiasAttention(nn.Module):
    """Multi-head self-attention with additive structural attention bias."""

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
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, bias: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        bsz, n_tokens, _ = x.shape
        q = self.q_proj(x).view(bsz, n_tokens, self.n_heads, self.head_dim).transpose(1, 2)
        k = self.k_proj(x).view(bsz, n_tokens, self.n_heads, self.head_dim).transpose(1, 2)
        v = self.v_proj(x).view(bsz, n_tokens, self.n_heads, self.head_dim).transpose(1, 2)
        logits = torch.matmul(q, k.transpose(-2, -1)) / (self.head_dim ** 0.5)
        logits = logits + bias
        attn = torch.softmax(logits, dim=-1)
        out = torch.matmul(self.dropout(attn), v).transpose(1, 2).contiguous().view(bsz, n_tokens, self.hidden_size)
        return self.out_proj(out), attn


class GraphormerLocalLayer(nn.Module):
    def __init__(self, hidden_size: int, n_heads: int, dropout: float):
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_size)
        self.attn = GraphormerBiasAttention(hidden_size, n_heads, dropout)
        self.norm2 = nn.LayerNorm(hidden_size)
        self.ffn = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, hidden_size),
        )
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, bias: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        attn_out, attn = self.attn(self.norm1(x), bias)
        x = x + self.dropout(attn_out)
        x = x + self.dropout(self.ffn(self.norm2(x)))
        return x, attn


class GraphormerStructuralScMAE(nn.Module):
    """scMAE with local Graphormer centrality and spatial attention bias."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        n_heads: int = 4,
        depth: int = 2,
        max_neighbors: int = 15,
        dropout: float = 0.05,
        svd_dim: int = 64,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.n_heads = int(n_heads)
        self.max_neighbors = int(max_neighbors)
        self.svd_dim = int(svd_dim)
        self.expr_embed = nn.Sequential(
            nn.Linear(num_genes, 256),
            nn.LayerNorm(256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, hidden_size),
        )
        self.anchor_embed = nn.Sequential(
            nn.Linear(svd_dim, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
        )
        self.centrality_embed = nn.Embedding(max_neighbors + 2, hidden_size)
        self.spatial_bias = nn.Embedding(max_neighbors + 2, n_heads)
        self.edge_bias = nn.Embedding(3, n_heads)
        self.layers = nn.ModuleList([GraphormerLocalLayer(hidden_size, n_heads, dropout) for _ in range(depth)])
        self.out_norm = nn.LayerNorm(hidden_size)
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size + num_genes, 256),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(256, num_genes),
        )
        self.edge_head = nn.Sequential(
            nn.Linear(hidden_size * 3, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, 1),
        )
        self.anchor_head = nn.Linear(hidden_size, svd_dim)

    def _structural_bias(self, ranks: torch.Tensor, keep_mask: torch.Tensor | None) -> torch.Tensor:
        bsz, n_tokens = ranks.shape
        dist = torch.maximum(ranks[:, :, None], ranks[:, None, :]).clamp(max=self.max_neighbors + 1)
        bias = self.spatial_bias(dist).permute(0, 3, 1, 2)
        edge_type = torch.ones_like(dist)
        edge_type[:, 0, :] = 2
        edge_type[:, :, 0] = 2
        if keep_mask is not None:
            dropped = keep_mask <= 0
            pair_dropped = dropped[:, :, None] | dropped[:, None, :]
            edge_type = edge_type.masked_fill(pair_dropped, 0)
        bias = bias + self.edge_bias(edge_type).permute(0, 3, 1, 2)
        return bias

    def encode(
        self,
        x: torch.Tensor,
        neighbor_x: torch.Tensor | None = None,
        anchor: torch.Tensor | None = None,
        neighbor_anchor: torch.Tensor | None = None,
        keep_mask: torch.Tensor | None = None,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if neighbor_x is None:
            neighbor_x = x[:, None, :].expand(-1, self.max_neighbors, -1)
        if neighbor_x.shape[1] > self.max_neighbors:
            neighbor_x = neighbor_x[:, : self.max_neighbors]
            if neighbor_anchor is not None:
                neighbor_anchor = neighbor_anchor[:, : self.max_neighbors]
            if keep_mask is not None:
                keep_mask = keep_mask[:, : self.max_neighbors]
        if neighbor_x.shape[1] < self.max_neighbors:
            pad = self.max_neighbors - neighbor_x.shape[1]
            neighbor_x = F.pad(neighbor_x, (0, 0, 0, pad))
            if neighbor_anchor is not None:
                neighbor_anchor = F.pad(neighbor_anchor, (0, 0, 0, pad))
            if keep_mask is not None:
                keep_mask = F.pad(keep_mask, (0, pad))
        tokens_x = torch.cat([x[:, None, :], neighbor_x], dim=1)
        tokens = self.expr_embed(tokens_x)
        if anchor is not None:
            if neighbor_anchor is None:
                neighbor_anchor = anchor[:, None, :].expand(-1, self.max_neighbors, -1)
            tokens_anchor = torch.cat([anchor[:, None, :], neighbor_anchor], dim=1)
            tokens = tokens + self.anchor_embed(tokens_anchor)
        bsz, n_tokens, _ = tokens.shape
        ranks = torch.arange(n_tokens, device=x.device).view(1, -1).expand(bsz, -1)
        centrality = torch.zeros_like(ranks)
        centrality[:, 0] = self.max_neighbors + 1
        centrality[:, 1:] = torch.clamp(self.max_neighbors + 1 - ranks[:, 1:], min=1)
        tokens = tokens + self.centrality_embed(centrality)
        if keep_mask is not None:
            keep_full = torch.cat([torch.ones(keep_mask.shape[0], 1, device=keep_mask.device, dtype=keep_mask.dtype), keep_mask], dim=1)
        else:
            keep_full = None
        bias = self._structural_bias(ranks, keep_full)
        attn_last = None
        for layer in self.layers:
            tokens, attn_last = layer(tokens, bias)
        latent = self.out_norm(tokens[:, 0])
        if anchor is not None:
            latent = latent + 0.2 * self.anchor_embed(anchor)
        if attn_last is None:
            attn_last = torch.zeros(x.shape[0], self.n_heads, self.max_neighbors + 1, self.max_neighbors + 1, device=x.device)
        return latent, attn_last

    def forward(
        self,
        x: torch.Tensor,
        neighbor_x: torch.Tensor | None = None,
        anchor: torch.Tensor | None = None,
        neighbor_anchor: torch.Tensor | None = None,
        keep_mask: torch.Tensor | None = None,
    ) -> dict[str, torch.Tensor]:
        latent, attn = self.encode(x, neighbor_x, anchor, neighbor_anchor, keep_mask)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        anchor_pred = self.anchor_head(latent)
        return {"latent": latent, "attention": attn, "mask_logits": mask_logits, "reconstruction": reconstruction, "anchor_pred": anchor_pred}

    def edge_logits(self, z: torch.Tensor, pos_z: torch.Tensor, neg_z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        pos = self.edge_head(torch.cat([z, pos_z, torch.abs(z - pos_z)], dim=1)).squeeze(-1)
        neg = self.edge_head(torch.cat([z, neg_z, torch.abs(z - neg_z)], dim=1)).squeeze(-1)
        return pos, neg

    @torch.no_grad()
    def feature(
        self,
        x: torch.Tensor,
        neighbor_x: torch.Tensor | None = None,
        anchor: torch.Tensor | None = None,
        neighbor_anchor: torch.Tensor | None = None,
    ) -> torch.Tensor:
        latent, _ = self.encode(x, neighbor_x, anchor, neighbor_anchor, None)
        return latent

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
        shape = [1] * x.dim()
        shape[-1] = x.shape[-1]
        return x * keep.view(*shape)
