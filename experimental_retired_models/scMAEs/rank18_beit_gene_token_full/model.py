from __future__ import annotations

import torch
from torch import nn


class BEiTGeneTokenScMAE(nn.Module):
    """scMAE body with BEiT-style discrete gene-token prediction."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        token_bins: int = 8,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.token_bins = int(token_bins)
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
        )
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Linear(self.hidden_size + self.num_genes, self.num_genes)
        self.token_head = nn.Sequential(
            nn.LayerNorm(self.hidden_size),
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.GELU(),
            nn.Linear(self.hidden_size, self.num_genes * self.token_bins),
        )
        self.replaced_head = nn.Linear(self.hidden_size, self.num_genes)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        token_logits = self.token_head(latent).view(x.shape[0], self.num_genes, self.token_bins)
        replaced_logits = self.replaced_head(latent)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "token_logits": token_logits,
            "replaced_logits": replaced_logits,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def beit_corrupt(
        self,
        x: torch.Tensor,
        mask_prob: float,
        replace_prob: float,
        mask_value: float = 0.0,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        replacement = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else x
        replace_gate = ((torch.rand_like(x) < float(replace_prob)).float() * mask).float()
        corrupted = torch.where(replace_gate.bool(), replacement, x)
        corrupted = torch.where((mask - replace_gate).bool(), torch.full_like(x, float(mask_value)), corrupted)
        return corrupted, mask, replace_gate

