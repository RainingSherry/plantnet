from __future__ import annotations

import torch
from torch import nn


class Data2VecScMAE(nn.Module):
    """scMAE body exposing layerwise features for data2vec targets."""

    def __init__(self, num_genes: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.drop = nn.Dropout(float(dropout))
        self.fc1 = nn.Linear(self.num_genes, 256)
        self.norm1 = nn.LayerNorm(256)
        self.fc2 = nn.Linear(256, self.hidden_size)
        self.norm2 = nn.LayerNorm(self.hidden_size)
        self.fc3 = nn.Linear(self.hidden_size, self.hidden_size)
        self.norm3 = nn.LayerNorm(self.hidden_size)
        self.act = nn.Mish()
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.decoder = nn.Linear(self.hidden_size + self.num_genes, self.num_genes)
        self.predictor = nn.Sequential(
            nn.LayerNorm(self.hidden_size),
            nn.Linear(self.hidden_size, self.hidden_size),
            nn.GELU(),
            nn.Linear(self.hidden_size, self.hidden_size),
        )

    def encode_layers(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h1 = self.act(self.norm1(self.fc1(self.drop(x))))
        h2 = self.act(self.norm2(self.fc2(h1)))
        h3 = self.norm3(self.fc3(h2))
        return h1[:, : self.hidden_size], h2, h3

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        _, h2, latent = self.encode_layers(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        prediction = self.predictor(latent)
        return {
            "latent": latent,
            "teacher_layer": 0.5 * (h2 + latent),
            "prediction": prediction,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.forward(x)["latent"]

    def mask_view(self, x: torch.Tensor, mask_prob: float, jitter_std: float = 0.0) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        view = x.masked_fill(mask.bool(), 0.0)
        if jitter_std > 0:
            view = view + torch.randn_like(view) * float(jitter_std) * (1.0 - mask)
        return view, mask

