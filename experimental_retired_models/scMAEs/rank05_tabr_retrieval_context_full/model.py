from __future__ import annotations

import torch
from torch import nn


class TabRContextScMAE(nn.Module):
    """Retrieval-context scMAE with train/inference-consistent context input."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.0):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.encoder = nn.Sequential(
            nn.Dropout(p=float(dropout)),
            nn.Linear(self.num_genes * 3, 384),
            nn.LayerNorm(384),
            nn.Mish(),
            nn.Linear(384, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.hidden_size),
        )
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Linear(self.hidden_size + self.num_genes, self.num_genes)
        self.context_projector = nn.Sequential(
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, self.hidden_size),
        )

    def encode_input(self, x: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        return torch.cat([x, context, x - context], dim=1)

    def forward(self, x: torch.Tensor, context: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(self.encode_input(x, context))
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        context_latent = self.context_projector(context)
        return {
            "latent": latent,
            "context_latent": context_latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor, context: torch.Tensor) -> torch.Tensor:
        return self.forward(x, context)["latent"]

    def corrupt_swap(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        should_swap = torch.bernoulli(float(mask_prob) * torch.ones_like(x)).bool()
        repl = x[torch.randperm(x.shape[0], device=x.device)] if x.shape[0] > 1 else x
        corrupted = torch.where(should_swap, repl, x)
        mask = (corrupted != x).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            corrupted[empty, cols] = repl[empty, cols]
            mask[empty, cols] = (corrupted[empty, cols] != x[empty, cols]).float()
        return corrupted, mask

