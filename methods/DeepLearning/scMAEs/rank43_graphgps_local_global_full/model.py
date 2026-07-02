from __future__ import annotations

import torch
from torch import nn


class GPSAdapterLayer(nn.Module):
    """Shallow GraphGPS-style layer: local neighbor message + global attention."""

    def __init__(self, hidden_size: int, pe_dim: int, num_heads: int, dropout: float):
        super().__init__()
        self.pe_proj = nn.Sequential(
            nn.Linear(pe_dim, hidden_size),
            nn.GELU(),
            nn.LayerNorm(hidden_size),
        )
        self.local_msg = nn.Sequential(
            nn.Linear(hidden_size * 2, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, hidden_size),
        )
        self.local_gate = nn.Sequential(nn.Linear(hidden_size * 2, hidden_size), nn.Sigmoid())
        self.local_norm = nn.LayerNorm(hidden_size)
        self.attn = nn.MultiheadAttention(hidden_size, num_heads=num_heads, dropout=dropout, batch_first=True)
        self.attn_norm = nn.LayerNorm(hidden_size)
        self.ffn = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, hidden_size),
            nn.Dropout(dropout),
        )
        self.ffn_norm = nn.LayerNorm(hidden_size)

    def forward(self, h: torch.Tensor, neigh_h: torch.Tensor, pe: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h = h + self.pe_proj(pe)
        neigh_mean = neigh_h.mean(dim=1)
        local_delta = self.local_msg(torch.cat([h, neigh_mean], dim=-1))
        gate = self.local_gate(torch.cat([h, neigh_mean], dim=-1))
        h_local = self.local_norm(h + gate * local_delta)

        attn_out, attn_weights = self.attn(h.unsqueeze(0), h.unsqueeze(0), h.unsqueeze(0), need_weights=True)
        h_global = self.attn_norm(h + attn_out.squeeze(0))
        h_out = self.ffn_norm(h_local + h_global + self.ffn(h_local + h_global))
        return h_out, gate, attn_weights.squeeze(0)


class GraphGPSScMAE(nn.Module):
    """scMAE with a GraphGPS-inspired local-global graph adapter."""

    def __init__(
        self,
        input_dim: int,
        pe_dim: int = 8,
        hidden_size: int = 256,
        decoder_hidden: int = 512,
        num_heads: int = 4,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.input_dim = input_dim
        self.pe_dim = pe_dim
        self.input_encoder = nn.Sequential(
            nn.Linear(input_dim, decoder_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(decoder_hidden, hidden_size),
            nn.LayerNorm(hidden_size),
        )
        self.gps = GPSAdapterLayer(hidden_size, pe_dim, num_heads, dropout)
        self.mask_predictor = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, input_dim),
        )
        self.decoder = nn.Sequential(
            nn.Linear(hidden_size, decoder_hidden),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(decoder_hidden, input_dim),
        )
        self.pe_decoder = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, pe_dim),
        )
        self.edge_head = nn.Sequential(
            nn.Linear(hidden_size * 3, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, 1),
        )

    def encode_base(self, x: torch.Tensor) -> torch.Tensor:
        return self.input_encoder(x)

    def encode(self, x: torch.Tensor, neigh_x: torch.Tensor, pe: torch.Tensor, neigh_pe: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h = self.encode_base(x)
        neigh_h = self.encode_base(neigh_x.reshape(-1, neigh_x.shape[-1])).view(neigh_x.shape[0], neigh_x.shape[1], -1)
        neigh_h = neigh_h + self.gps.pe_proj(neigh_pe.reshape(-1, neigh_pe.shape[-1])).view(neigh_pe.shape[0], neigh_pe.shape[1], -1)
        return self.gps(h, neigh_h, pe)

    def forward(self, x: torch.Tensor, neigh_x: torch.Tensor, pe: torch.Tensor, neigh_pe: torch.Tensor) -> dict[str, torch.Tensor]:
        z, gate, attn = self.encode(x, neigh_x, pe, neigh_pe)
        return {
            "embedding": z,
            "reconstruction": self.decoder(z),
            "mask_logits": self.mask_predictor(z),
            "pe_reconstruction": self.pe_decoder(z),
            "local_gate": gate,
            "attention": attn,
        }

    def edge_logits(self, z_a: torch.Tensor, z_b: torch.Tensor) -> torch.Tensor:
        pair = torch.cat([z_a, z_b, torch.abs(z_a - z_b)], dim=-1)
        return self.edge_head(pair).squeeze(-1)

    @staticmethod
    def mask_view(x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = torch.rand_like(x) < mask_prob
        corrupted = x.masked_fill(mask, 0.0)
        return corrupted, mask.float()
