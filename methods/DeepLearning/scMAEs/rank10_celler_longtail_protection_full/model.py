from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class CellerLongTailScMAE(nn.Module):
    """scMAE body with Celler-inspired gene tokens and prototype logits."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        token_bins: int = 16,
        n_prototypes: int = 64,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.token_bins = int(token_bins)
        self.n_prototypes = int(n_prototypes)
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
        self.prototype_head = nn.Sequential(
            nn.LayerNorm(self.hidden_size),
            nn.Linear(self.hidden_size, self.n_prototypes),
        )
        self.token_head = nn.Sequential(
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.Mish(),
            nn.Linear(self.hidden_size, self.num_genes * self.token_bins),
        )

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        proto_logits = self.prototype_head(latent)
        token_logits = self.token_head(latent).view(x.shape[0], self.num_genes, self.token_bins)
        return {
            "latent": latent,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "proto_logits": proto_logits,
            "token_logits": token_logits,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def mask_nonzero_tokens(self, x: torch.Tensor, log_expr: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        nonzero = log_expr > 1e-8
        mask = (torch.rand_like(x) < float(mask_prob)) & nonzero
        fallback = mask.sum(dim=1) == 0
        if bool(fallback.any()):
            cols = torch.randint(0, x.shape[1], (int(fallback.sum()),), device=x.device)
            mask[fallback, cols] = True
        corrupted = x.masked_fill(mask, 0.0)
        return corrupted, mask.float()

    def corrupt_boundary_view(self, x: torch.Tensor, log_expr: torch.Tensor, mask_prob: float, jitter_std: float) -> tuple[torch.Tensor, torch.Tensor]:
        corrupted, mask = self.mask_nonzero_tokens(x, log_expr, mask_prob)
        if jitter_std > 0:
            noise = torch.randn_like(corrupted) * float(jitter_std)
            corrupted = corrupted + noise * (1.0 - mask)
        return corrupted, mask

    @staticmethod
    def normalized_entropy(logits: torch.Tensor) -> torch.Tensor:
        probs = F.softmax(logits, dim=1)
        ent = -(probs * torch.log(probs.clamp_min(1e-8))).sum(dim=1)
        return ent / max(1.0, float(torch.log(torch.tensor(logits.shape[1], dtype=logits.dtype, device=logits.device))))

