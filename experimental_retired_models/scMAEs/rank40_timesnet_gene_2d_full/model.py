from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class Inception2DBlock(nn.Module):
    """Parameter-efficient 2D inception block adapted from TimesNet."""

    def __init__(self, in_channels: int, out_channels: int, num_kernels: int = 3):
        super().__init__()
        self.kernels = nn.ModuleList(
            [nn.Conv2d(in_channels, out_channels, kernel_size=2 * i + 1, padding=i) for i in range(max(1, int(num_kernels)))]
        )
        for module in self.modules():
            if isinstance(module, nn.Conv2d):
                nn.init.kaiming_normal_(module.weight, mode="fan_out", nonlinearity="relu")
                if module.bias is not None:
                    nn.init.zeros_(module.bias)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return torch.stack([kernel(x) for kernel in self.kernels], dim=-1).mean(dim=-1)


def fft_for_period(x: torch.Tensor, top_k: int) -> tuple[list[int], torch.Tensor]:
    # x: [B, T, C]. Frequency 0 is removed to avoid selecting the mean component.
    xf = torch.fft.rfft(x, dim=1)
    frequency_score = xf.abs().mean(dim=0).mean(dim=-1)
    if frequency_score.numel() <= 1:
        weights = x.new_ones((x.shape[0], 1))
        return [max(1, x.shape[1])], weights
    frequency_score = frequency_score.clone()
    frequency_score[0] = 0
    k = min(max(1, int(top_k)), frequency_score.numel() - 1)
    _, top_idx = torch.topk(frequency_score, k)
    periods = [max(1, int(x.shape[1] // max(1, int(idx)))) for idx in top_idx.detach().cpu().tolist()]
    weights = xf.abs().mean(dim=-1)[:, top_idx]
    return periods, weights


class TimesBlock(nn.Module):
    """FFT period discovery + 2D variation modeling on gene-module sequences."""

    def __init__(self, d_model: int, d_ff: int, top_k: int = 3, num_kernels: int = 3, dropout: float = 0.05):
        super().__init__()
        self.top_k = int(top_k)
        self.conv = nn.Sequential(
            Inception2DBlock(d_model, d_ff, num_kernels),
            nn.GELU(),
            nn.Dropout(dropout),
            Inception2DBlock(d_ff, d_model, num_kernels),
        )
        self.norm = nn.LayerNorm(d_model)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        batch, steps, channels = x.shape
        periods, period_weight = fft_for_period(x, self.top_k)
        outputs = []
        for period in periods:
            length = steps if steps % period == 0 else ((steps // period) + 1) * period
            padded = x if length == steps else torch.cat([x, x.new_zeros(batch, length - steps, channels)], dim=1)
            view = padded.reshape(batch, length // period, period, channels).permute(0, 3, 1, 2).contiguous()
            view = self.conv(view)
            view = view.permute(0, 2, 3, 1).reshape(batch, length, channels)
            outputs.append(view[:, :steps])
        stacked = torch.stack(outputs, dim=-1)
        weight = F.softmax(period_weight, dim=1).view(batch, 1, 1, -1)
        return self.norm(torch.sum(stacked * weight, dim=-1) + x)


class TimesNetGene2DScMAE(nn.Module):
    """scMAE with TimesNet-style 2D variation modeling over SVD-ordered gene modules."""

    def __init__(
        self,
        num_genes: int,
        gene_order: torch.Tensor,
        module_size: int = 25,
        hidden_size: int = 128,
        d_model: int = 64,
        d_ff: int = 96,
        e_layers: int = 2,
        top_k: int = 3,
        num_kernels: int = 3,
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
        self.module_input = nn.Sequential(
            nn.Linear(self.module_size, self.d_model),
            nn.LayerNorm(self.d_model),
            nn.GELU(),
        )
        self.times_blocks = nn.ModuleList([TimesBlock(self.d_model, d_ff, top_k, num_kernels, dropout) for _ in range(max(1, int(e_layers)))])
        self.module_to_cell = nn.Sequential(nn.Linear(self.d_model, self.hidden_size), nn.LayerNorm(self.hidden_size), nn.GELU())
        self.fusion_norm = nn.LayerNorm(self.hidden_size)
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Sequential(
            nn.Linear(self.hidden_size + self.num_genes, 256),
            nn.Mish(),
            nn.Dropout(dropout),
            nn.Linear(256, self.num_genes),
        )
        self.module_head = nn.Sequential(nn.Linear(self.d_model, self.d_model), nn.GELU(), nn.Linear(self.d_model, 1))

    def to_module_sequence(self, x: torch.Tensor) -> torch.Tensor:
        ordered = x.index_select(1, self.gene_order)
        if self.padded_genes > self.num_genes:
            ordered = F.pad(ordered, (0, self.padded_genes - self.num_genes))
        return ordered.view(x.shape[0], self.n_modules, self.module_size)

    def module_targets(self, x: torch.Tensor) -> torch.Tensor:
        return self.to_module_sequence(x).mean(dim=2)

    def encode_module_sequence(self, seq: torch.Tensor, apply_timestamp_mask: bool) -> torch.Tensor:
        z = self.module_input(seq)
        if apply_timestamp_mask and self.training and self.timestamp_mask_prob > 0:
            keep = torch.rand(z.shape[:2], device=z.device) > self.timestamp_mask_prob
            all_masked = keep.sum(dim=1) == 0
            if bool(all_masked.any()):
                cols = torch.randint(0, z.shape[1], (int(all_masked.sum()),), device=z.device)
                keep[all_masked, cols] = True
            z = z.masked_fill(~keep.unsqueeze(-1), 0.0)
        for block in self.times_blocks:
            z = block(z)
        return z

    def encode_modules(self, x: torch.Tensor, apply_timestamp_mask: bool = False) -> torch.Tensor:
        return self.encode_module_sequence(self.to_module_sequence(x), apply_timestamp_mask)

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        cell_latent = self.cell_encoder(x)
        module_repr = self.encode_modules(x, apply_timestamp_mask=False)
        module_pool = module_repr.max(dim=1).values
        return self.fusion_norm(cell_latent + self.module_to_cell(module_pool))

    def consistency_views(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        seq = self.to_module_sequence(x)
        return self.encode_module_sequence(seq, True), self.encode_module_sequence(seq, True)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        module_repr = self.encode_modules(x, apply_timestamp_mask=False)
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
