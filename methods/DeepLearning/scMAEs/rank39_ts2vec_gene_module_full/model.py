from __future__ import annotations

import random

import torch
from torch import nn
import torch.nn.functional as F


class DilatedResidualBlock(nn.Module):
    def __init__(self, channels: int, dilation: int, dropout: float):
        super().__init__()
        pad = dilation
        self.net = nn.Sequential(
            nn.Conv1d(channels, channels, kernel_size=3, padding=pad, dilation=dilation),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Conv1d(channels, channels, kernel_size=3, padding=pad, dilation=dilation),
            nn.GELU(),
        )
        self.norm = nn.GroupNorm(1, channels)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.norm(x + self.net(x))


class TS2VecGeneModuleScMAE(nn.Module):
    """scMAE with a TS2Vec-style contextual encoder over data-driven gene modules."""

    def __init__(
        self,
        num_genes: int,
        gene_order: torch.Tensor,
        module_size: int = 25,
        hidden_size: int = 128,
        module_dim: int = 96,
        depth: int = 4,
        dropout: float = 0.05,
        timestamp_mask_prob: float = 0.25,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.module_size = int(module_size)
        self.hidden_size = int(hidden_size)
        self.module_dim = int(module_dim)
        self.timestamp_mask_prob = float(timestamp_mask_prob)
        order = torch.as_tensor(gene_order, dtype=torch.long)
        if order.numel() != self.num_genes:
            raise ValueError("gene_order length must equal num_genes")
        self.register_buffer("gene_order", order, persistent=True)
        self.n_modules = int((self.num_genes + self.module_size - 1) // self.module_size)
        self.padded_genes = int(self.n_modules * self.module_size)

        self.cell_encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
        )
        self.module_input = nn.Sequential(
            nn.Linear(self.module_size, self.module_dim),
            nn.LayerNorm(self.module_dim),
            nn.GELU(),
        )
        dilations = [2**i for i in range(max(1, int(depth)))]
        self.module_encoder = nn.Sequential(*[DilatedResidualBlock(self.module_dim, d, dropout) for d in dilations])
        self.module_to_cell = nn.Sequential(
            nn.Linear(self.module_dim, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.GELU(),
        )
        self.fusion_norm = nn.LayerNorm(self.hidden_size)
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(self.hidden_size + self.num_genes, 256),
            nn.Mish(),
            nn.Dropout(dropout),
            nn.Linear(256, self.num_genes),
        )
        self.module_head = nn.Sequential(
            nn.Linear(self.module_dim, self.module_dim),
            nn.GELU(),
            nn.Linear(self.module_dim, 1),
        )

    def to_module_sequence(self, x: torch.Tensor) -> torch.Tensor:
        ordered = x.index_select(1, self.gene_order)
        if self.padded_genes > self.num_genes:
            ordered = F.pad(ordered, (0, self.padded_genes - self.num_genes))
        return ordered.view(x.shape[0], self.n_modules, self.module_size)

    def module_targets(self, x: torch.Tensor) -> torch.Tensor:
        seq = self.to_module_sequence(x)
        return seq.mean(dim=2)

    def encode_module_sequence(self, seq: torch.Tensor, apply_timestamp_mask: bool) -> torch.Tensor:
        z = self.module_input(seq)
        if apply_timestamp_mask and self.training and self.timestamp_mask_prob > 0:
            keep = torch.rand(z.shape[:2], device=z.device) > self.timestamp_mask_prob
            all_masked = keep.sum(dim=1) == 0
            if bool(all_masked.any()):
                cols = torch.randint(0, z.shape[1], (int(all_masked.sum()),), device=z.device)
                keep[all_masked, cols] = True
            z = z.masked_fill(~keep.unsqueeze(-1), 0.0)
        z = z.transpose(1, 2)
        z = self.module_encoder(z)
        return z.transpose(1, 2)

    def encode_modules(self, x: torch.Tensor, apply_timestamp_mask: bool = False) -> torch.Tensor:
        return self.encode_module_sequence(self.to_module_sequence(x), apply_timestamp_mask)

    def sample_context_views(self, x: torch.Tensor, temporal_unit: int = 0) -> tuple[torch.Tensor, torch.Tensor]:
        seq = self.to_module_sequence(x)
        ts_l = seq.shape[1]
        if ts_l <= 2:
            out = self.encode_module_sequence(seq, apply_timestamp_mask=True)
            return out, out
        min_l = min(ts_l, max(2, 2 ** (int(temporal_unit) + 1)))
        crop_l = random.randint(min_l, ts_l)
        crop_left = random.randint(0, ts_l - crop_l)
        crop_right = crop_left + crop_l
        crop_eleft = random.randint(0, crop_left)
        crop_eright = random.randint(crop_right, ts_l)
        out1 = self.encode_module_sequence(seq[:, crop_eleft:crop_right], apply_timestamp_mask=True)
        out2 = self.encode_module_sequence(seq[:, crop_left:crop_eright], apply_timestamp_mask=True)
        return out1[:, -crop_l:], out2[:, :crop_l]

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        cell_latent = self.cell_encoder(x)
        module_repr = self.encode_modules(x, apply_timestamp_mask=False)
        module_pool = module_repr.max(dim=1).values
        return self.fusion_norm(cell_latent + self.module_to_cell(module_pool))

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        module_repr = self.encode_modules(x, apply_timestamp_mask=False)
        module_pool = module_repr.max(dim=1).values
        cell_latent = self.cell_encoder(x)
        latent = self.fusion_norm(cell_latent + self.module_to_cell(module_pool))
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        module_prediction = self.module_head(module_repr).squeeze(-1)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "module_prediction": module_prediction,
            "module_repr": module_repr,
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
