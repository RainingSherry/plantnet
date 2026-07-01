from __future__ import annotations

import torch
from torch import nn
import torch.nn.functional as F


class ScAGCAdaptiveGraphScMAE(nn.Module):
    """scMAE body with a shallow adaptive graph adapter and edge decoder."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.graph_adapter = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, hidden_size),
            nn.Mish(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(hidden_size + self.num_genes, self.num_genes)
        self.edge_decoder = nn.Sequential(
            nn.Linear(hidden_size * 4, hidden_size),
            nn.Mish(),
            nn.Linear(hidden_size, 1),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        base = self.encoder(x)
        return base + 0.25 * self.graph_adapter(base)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encode(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction}

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encode(x)

    def edge_logits(self, z1: torch.Tensor, z2: torch.Tensor) -> torch.Tensor:
        pair = torch.cat([z1, z2, torch.abs(z1 - z2), z1 * z2], dim=1)
        return self.edge_decoder(pair).squeeze(1)

    def random_mask(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask


def confidence_gated_mix(
    x: torch.Tensor,
    neighbor_x: torch.Tensor,
    edge_conf: torch.Tensor,
    rare_risk: torch.Tensor,
    edge_threshold: float,
    rare_threshold: float,
    alpha: float,
    keep_prob: float,
) -> tuple[torch.Tensor, torch.Tensor]:
    gate = ((edge_conf >= float(edge_threshold)) & (rare_risk <= float(rare_threshold))).float()
    if keep_prob < 1.0:
        gate = gate * (torch.rand_like(gate) < float(keep_prob)).float()
    mixed = gate[:, None] * (float(alpha) * x + (1.0 - float(alpha)) * neighbor_x) + (1.0 - gate[:, None]) * x
    return mixed, gate

