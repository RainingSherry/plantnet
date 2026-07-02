from __future__ import annotations

import torch
from torch import nn


class MultiMAEGeneScMAE(nn.Module):
    """Shared scMAE encoder with lightweight task-specific decoders."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        token_bins: int = 8,
        module_size: int = 25,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.token_bins = int(token_bins)
        self.module_size = int(module_size)
        self.n_modules = (self.num_genes + self.module_size - 1) // self.module_size
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
        self.expr_decoder = nn.Sequential(
            nn.Linear(self.hidden_size + self.num_genes, 256),
            nn.Mish(),
            nn.Linear(256, self.num_genes),
        )
        self.token_decoder = nn.Sequential(
            nn.LayerNorm(self.hidden_size),
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.GELU(),
            nn.Linear(self.hidden_size, self.num_genes * self.token_bins),
        )
        self.module_decoder = nn.Sequential(
            nn.LayerNorm(self.hidden_size),
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.n_modules),
        )

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.expr_decoder(torch.cat([latent, mask_logits], dim=1))
        token_logits = self.token_decoder(latent).view(x.shape[0], self.num_genes, self.token_bins)
        module_reconstruction = self.module_decoder(latent)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "token_logits": token_logits,
            "module_reconstruction": module_reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def multitask_mask(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

