from __future__ import annotations

import torch
from torch import nn


class SoftGraphEdgeAdapter(nn.Module):
    """Shallow edge-confidence adapter for soft graph regularization."""

    def __init__(self, hidden_size: int):
        super().__init__()
        self.edge_mlp = nn.Sequential(
            nn.Linear(hidden_size * 4, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
            nn.Linear(hidden_size, 1),
        )

    def forward(self, z: torch.Tensor, neighbor_z: torch.Tensor) -> torch.Tensor:
        z_expand = z[:, None, :].expand_as(neighbor_z)
        pair = torch.cat([z_expand, neighbor_z, torch.abs(z_expand - neighbor_z), z_expand * neighbor_z], dim=-1)
        return self.edge_mlp(pair).squeeze(-1)


class SoftGraphScMAE(nn.Module):
    """scMAE backbone with a shallow soft-graph edge confidence adapter."""

    def __init__(self, input_dim: int, hidden_size: int = 128, decoder_hidden: int = 128, dropout: float = 0.0):
        super().__init__()
        self.input_dim = int(input_dim)
        self.hidden_size = int(hidden_size)
        self.encoder = nn.Sequential(
            nn.Dropout(p=dropout),
            nn.Linear(input_dim, 256),
            nn.LayerNorm(256),
            nn.Mish(inplace=True),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, input_dim)
        self.decoder = nn.Sequential(nn.Linear(hidden_size + input_dim, decoder_hidden), nn.Mish(inplace=True), nn.Linear(decoder_hidden, input_dim))
        self.edge_adapter = SoftGraphEdgeAdapter(hidden_size)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z = self.encoder(x)
        mask_logits = self.mask_predictor(z)
        recon = self.decoder(torch.cat([z, mask_logits], dim=1))
        return {"embedding": z, "reconstruction": recon, "mask_logits": mask_logits}

    def edge_logits(self, z: torch.Tensor, neighbor_z: torch.Tensor) -> torch.Tensor:
        return self.edge_adapter(z, neighbor_z)
