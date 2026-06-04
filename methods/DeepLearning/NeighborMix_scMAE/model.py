from __future__ import annotations

import torch
import torch.nn as nn
from torch.nn.functional import binary_cross_entropy_with_logits as bce_logits
from torch.nn.functional import mse_loss as mse


class AutoEncoder(torch.nn.Module):
    """scMAE masked autoencoder kept compatible with OtherMode/scMAE-main."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        dropout: float = 0.0,
        masked_data_weight: float = 0.75,
        mask_loss_weight: float = 0.7,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.masked_data_weight = float(masked_data_weight)
        self.mask_loss_weight = float(mask_loss_weight)

        self.encoder = nn.Sequential(
            nn.Dropout(p=float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(inplace=True),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(inplace=True),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(in_features=hidden_size + self.num_genes, out_features=self.num_genes)

    def forward_mask(self, x: torch.Tensor):
        latent = self.encoder(x)
        predicted_mask = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, predicted_mask], dim=1))
        return latent, predicted_mask, reconstruction

    def loss_mask(self, x: torch.Tensor, y: torch.Tensor, mask: torch.Tensor):
        latent, predicted_mask, reconstruction = self.forward_mask(x)
        weights = mask * self.masked_data_weight + (1.0 - mask) * (1.0 - self.masked_data_weight)
        reconstruction_loss = (1.0 - self.mask_loss_weight) * torch.mul(
            weights,
            mse(reconstruction, y, reduction="none"),
        )
        mask_loss = self.mask_loss_weight * bce_logits(predicted_mask, mask, reduction="mean")
        loss = reconstruction_loss.mean() + mask_loss
        return latent, loss

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

