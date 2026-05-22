from typing import Optional

import torch
import torch.nn as nn
import torch.nn.functional as F


class GeneGraphAttentionLayer(nn.Module):
    def __init__(self, in_dim: int, out_dim: int, dropout: float = 0.1):
        super().__init__()
        self.query = nn.Linear(in_dim, out_dim, bias=False)
        self.key = nn.Linear(in_dim, out_dim, bias=False)
        self.value = nn.Linear(in_dim, out_dim, bias=False)
        self.edge_gate = nn.Linear(1, out_dim, bias=False)
        self.dropout = dropout
        self.norm = nn.LayerNorm(out_dim)

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, edge_weight: Optional[torch.Tensor] = None) -> torch.Tensor:
        num_nodes = x.size(0)
        src, dst = edge_index
        q = self.query(x)
        k = self.key(x)
        v = self.value(x)

        attn_logits = (q[src] * k[dst]).sum(dim=-1) / (q.size(-1) ** 0.5)
        if edge_weight is None:
            edge_weight = torch.ones_like(attn_logits)
        attn_logits = attn_logits + torch.log(edge_weight.clamp_min(1e-8))
        attn = torch.zeros_like(attn_logits)
        for node in torch.unique(src):
            mask = src == node
            attn[mask] = F.softmax(attn_logits[mask], dim=0)
        attn = F.dropout(attn, p=self.dropout, training=self.training)

        messages = v[dst] * attn.unsqueeze(-1)
        out = torch.zeros_like(v)
        out.index_add_(0, src, messages)
        out = self.norm(out + v)
        return F.gelu(out)


class GeneGATEncoder(nn.Module):
    def __init__(
        self,
        input_dim: int,
        hidden_dim: int,
        output_dim: int,
        num_layers: int = 2,
        dropout: float = 0.1,
    ):
        super().__init__()
        layers = []
        dims = [input_dim] + [hidden_dim] * max(0, num_layers - 1) + [output_dim]
        for in_dim, out_dim in zip(dims[:-1], dims[1:]):
            layers.append(GeneGraphAttentionLayer(in_dim, out_dim, dropout=dropout))
        self.layers = nn.ModuleList(layers)
        self.input_proj = nn.Linear(input_dim, dims[1]) if input_dim != dims[1] else nn.Identity()
        self.dropout = dropout

    def forward(self, x: torch.Tensor, edge_index: torch.Tensor, edge_weight: Optional[torch.Tensor] = None) -> torch.Tensor:
        if x.dim() != 2:
            raise ValueError(f"Expected 2D gene feature tensor, got {tuple(x.shape)}")
        h = self.input_proj(x)
        h = F.dropout(h, p=self.dropout, training=self.training)
        for layer in self.layers:
            h = layer(h, edge_index=edge_index, edge_weight=edge_weight)
        return h
