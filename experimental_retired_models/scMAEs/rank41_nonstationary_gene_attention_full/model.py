from __future__ import annotations

import math

import torch
from torch import nn
import torch.nn.functional as F


class Projector(nn.Module):
    """Learns de-stationary factors from raw module sequence and its statistics."""

    def __init__(self, module_size: int, n_modules: int, hidden_size: int, output_dim: int):
        super().__init__()
        self.series_conv = nn.Conv1d(n_modules, 1, kernel_size=3, padding=1, padding_mode="circular", bias=False)
        self.net = nn.Sequential(
            nn.Linear(2 * module_size, hidden_size),
            nn.ReLU(),
            nn.Linear(hidden_size, hidden_size),
            nn.ReLU(),
            nn.Linear(hidden_size, output_dim, bias=False),
        )

    def forward(self, raw_seq: torch.Tensor, stats: torch.Tensor) -> torch.Tensor:
        pooled = self.series_conv(raw_seq)
        feat = torch.cat([pooled, stats], dim=1).flatten(1)
        return self.net(feat)


class DeStationarySelfAttention(nn.Module):
    def __init__(self, d_model: int, n_heads: int = 4, dropout: float = 0.05):
        super().__init__()
        if d_model % n_heads != 0:
            raise ValueError("d_model must be divisible by n_heads")
        self.n_heads = int(n_heads)
        self.head_dim = int(d_model // n_heads)
        self.q = nn.Linear(d_model, d_model)
        self.k = nn.Linear(d_model, d_model)
        self.v = nn.Linear(d_model, d_model)
        self.out = nn.Linear(d_model, d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, tau: torch.Tensor, delta: torch.Tensor) -> torch.Tensor:
        batch, steps, d_model = x.shape
        q = self.q(x).view(batch, steps, self.n_heads, self.head_dim).transpose(1, 2)
        k = self.k(x).view(batch, steps, self.n_heads, self.head_dim).transpose(1, 2)
        v = self.v(x).view(batch, steps, self.n_heads, self.head_dim).transpose(1, 2)
        scores = torch.matmul(q, k.transpose(-2, -1))
        scores = scores * tau.view(batch, 1, 1, 1) + delta.view(batch, 1, 1, steps)
        attn = torch.softmax(scores / math.sqrt(self.head_dim), dim=-1)
        out = torch.matmul(self.dropout(attn), v).transpose(1, 2).contiguous().view(batch, steps, d_model)
        return self.out(out)


class DeStationaryEncoderLayer(nn.Module):
    def __init__(self, d_model: int, n_heads: int, d_ff: int, dropout: float):
        super().__init__()
        self.attn = DeStationarySelfAttention(d_model, n_heads, dropout)
        self.ffn = nn.Sequential(nn.Linear(d_model, d_ff), nn.GELU(), nn.Dropout(dropout), nn.Linear(d_ff, d_model))
        self.norm1 = nn.LayerNorm(d_model)
        self.norm2 = nn.LayerNorm(d_model)
        self.dropout = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, tau: torch.Tensor, delta: torch.Tensor) -> torch.Tensor:
        x = self.norm1(x + self.dropout(self.attn(x, tau, delta)))
        return self.norm2(x + self.dropout(self.ffn(x)))


class NonStationaryGeneAttentionScMAE(nn.Module):
    """scMAE with non-stationary Transformer attention over SVD-ordered gene modules."""

    def __init__(
        self,
        num_genes: int,
        gene_order: torch.Tensor,
        module_size: int = 25,
        hidden_size: int = 128,
        d_model: int = 64,
        d_ff: int = 128,
        e_layers: int = 2,
        n_heads: int = 4,
        projector_hidden: int = 64,
        dropout: float = 0.05,
        timestamp_mask_prob: float = 0.25,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.module_size = int(module_size)
        self.hidden_size = int(hidden_size)
        self.d_model = int(d_model)
        self.timestamp_mask_prob = float(timestamp_mask_prob)
        order = torch.as_tensor(gene_order, dtype=torch.long)
        if order.numel() != self.num_genes:
            raise ValueError("gene_order length must equal num_genes")
        self.register_buffer("gene_order", order, persistent=True)
        self.n_modules = int((self.num_genes + self.module_size - 1) // self.module_size)
        self.padded_genes = int(self.n_modules * self.module_size)

        self.cell_encoder = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
        )
        self.module_input = nn.Sequential(nn.Linear(self.module_size, self.d_model), nn.LayerNorm(self.d_model), nn.GELU())
        self.layers = nn.ModuleList([DeStationaryEncoderLayer(self.d_model, n_heads, d_ff, dropout) for _ in range(max(1, int(e_layers)))])
        self.tau_learner = Projector(self.module_size, self.n_modules, projector_hidden, 1)
        self.delta_learner = Projector(self.module_size, self.n_modules, projector_hidden, self.n_modules)
        self.module_to_cell = nn.Sequential(nn.Linear(self.d_model, self.hidden_size), nn.LayerNorm(self.hidden_size), nn.GELU())
        self.fusion_norm = nn.LayerNorm(self.hidden_size)
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Sequential(nn.Linear(self.hidden_size + self.num_genes, 256), nn.Mish(), nn.Dropout(dropout), nn.Linear(256, self.num_genes))
        self.module_head = nn.Sequential(nn.Linear(self.d_model, self.d_model), nn.GELU(), nn.Linear(self.d_model, 1))

    def to_module_sequence(self, x: torch.Tensor) -> torch.Tensor:
        ordered = x.index_select(1, self.gene_order)
        if self.padded_genes > self.num_genes:
            ordered = F.pad(ordered, (0, self.padded_genes - self.num_genes))
        return ordered.view(x.shape[0], self.n_modules, self.module_size)

    def module_targets(self, x: torch.Tensor) -> torch.Tensor:
        return self.to_module_sequence(x).mean(dim=2)

    def encode_module_sequence(self, seq: torch.Tensor, apply_timestamp_mask: bool) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        raw_seq = seq
        mean = raw_seq.mean(dim=1, keepdim=True).detach()
        centered = raw_seq - mean
        std = torch.sqrt(centered.var(dim=1, keepdim=True, unbiased=False) + 1e-5).detach()
        z_seq = centered / std
        tau = self.tau_learner(raw_seq.detach(), std).exp().clamp(0.1, 10.0)
        delta = self.delta_learner(raw_seq.detach(), mean).clamp(-5.0, 5.0)
        z = self.module_input(z_seq)
        if apply_timestamp_mask and self.training and self.timestamp_mask_prob > 0:
            keep = torch.rand(z.shape[:2], device=z.device) > self.timestamp_mask_prob
            all_masked = keep.sum(dim=1) == 0
            if bool(all_masked.any()):
                cols = torch.randint(0, z.shape[1], (int(all_masked.sum()),), device=z.device)
                keep[all_masked, cols] = True
            z = z.masked_fill(~keep.unsqueeze(-1), 0.0)
        for layer in self.layers:
            z = layer(z, tau, delta)
        return z, tau, delta

    def encode_modules(self, x: torch.Tensor, apply_timestamp_mask: bool = False) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        return self.encode_module_sequence(self.to_module_sequence(x), apply_timestamp_mask)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        module_repr, _, _ = self.encode_modules(x, apply_timestamp_mask=False)
        module_pool = module_repr.max(dim=1).values
        return self.fusion_norm(self.cell_encoder(x) + self.module_to_cell(module_pool))

    def consistency_views(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        seq = self.to_module_sequence(x)
        return self.encode_module_sequence(seq, True)[0], self.encode_module_sequence(seq, True)[0]

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        module_repr, tau, delta = self.encode_modules(x, apply_timestamp_mask=False)
        module_pool = module_repr.max(dim=1).values
        latent = self.fusion_norm(self.cell_encoder(x) + self.module_to_cell(module_pool))
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "module_prediction": self.module_head(module_repr).squeeze(-1),
            "module_repr": module_repr,
            "tau_mean": tau.mean(),
            "delta_abs_mean": delta.abs().mean(),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode(x)

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask
