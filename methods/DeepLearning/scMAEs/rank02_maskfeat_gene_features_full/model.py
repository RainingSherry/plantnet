from __future__ import annotations

import math

import torch
import torch.nn as nn


class GenePatchMaskFeatMAE(nn.Module):
    """MaskFeat-style masked feature prediction for scRNA gene patches.

    MaskFeat predicts hand-crafted features of masked visual patches rather than
    reconstructing pixels only. This scRNA adaptation partitions genes into
    fixed-size patches, replaces masked patch embeddings with a learned mask
    token, and predicts deterministic gene-patch features for masked patches.
    """

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 16,
        hidden_size: int = 128,
        depth: int = 4,
        num_heads: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        if num_genes <= 0:
            raise ValueError("num_genes must be positive")
        if patch_size <= 1:
            raise ValueError("patch_size must be greater than 1")
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")

        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.num_patches = int(math.ceil(num_genes / patch_size))
        self.padded_genes = self.num_patches * self.patch_size
        self.feature_dim = self.patch_size * 3

        self.patch_embed = nn.Linear(self.patch_size, hidden_size)
        self.mask_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.cls_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.pos_embed = nn.Parameter(torch.zeros(1, self.num_patches + 1, hidden_size))

        layer = nn.TransformerEncoderLayer(
            d_model=hidden_size,
            nhead=num_heads,
            dim_feedforward=hidden_size * 4,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(layer, num_layers=depth)
        self.norm = nn.LayerNorm(hidden_size)
        self.feature_head = nn.Linear(hidden_size, self.feature_dim)
        self.reconstruction_head = nn.Linear(hidden_size, self.patch_size)
        self.mask_head = nn.Linear(hidden_size, 1)
        self.projector = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
        )
        nn.init.normal_(self.mask_token, std=0.02)
        nn.init.normal_(self.cls_token, std=0.02)
        nn.init.normal_(self.pos_embed, std=0.02)

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        if self.padded_genes > self.num_genes:
            pad = x.new_zeros(x.shape[0], self.padded_genes - self.num_genes)
            x = torch.cat([x, pad], dim=1)
        return x.view(x.shape[0], self.num_patches, self.patch_size)

    def random_mask(self, batch_size: int, mask_ratio: float, device: torch.device) -> torch.Tensor:
        if not 0.0 < float(mask_ratio) < 1.0:
            raise ValueError(f"mask_ratio must be in (0,1), got {mask_ratio}")
        len_keep = max(1, int(self.num_patches * (1.0 - float(mask_ratio))))
        noise = torch.rand(batch_size, self.num_patches, device=device)
        ids_shuffle = torch.argsort(noise, dim=1)
        ids_restore = torch.argsort(ids_shuffle, dim=1)
        mask = torch.ones(batch_size, self.num_patches, device=device)
        mask[:, :len_keep] = 0.0
        return torch.gather(mask, dim=1, index=ids_restore)

    def encode_patches(self, patches: torch.Tensor, patch_mask: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        if patches.ndim != 3 or patches.shape[1:] != (self.num_patches, self.patch_size):
            raise ValueError(
                f"patches must be [batch, {self.num_patches}, {self.patch_size}], got {tuple(patches.shape)}"
            )
        if patch_mask.shape != patches.shape[:2]:
            raise ValueError(f"patch_mask must be {tuple(patches.shape[:2])}, got {tuple(patch_mask.shape)}")
        x = self.patch_embed(patches)
        mask = patch_mask.to(dtype=x.dtype, device=x.device).unsqueeze(-1)
        x = x * (1.0 - mask) + self.mask_token.expand(x.shape[0], x.shape[1], -1) * mask
        cls = self.cls_token.expand(x.shape[0], -1, -1)
        x = torch.cat([cls, x], dim=1) + self.pos_embed
        encoded = self.norm(self.encoder(x))
        return encoded[:, 0], encoded[:, 1:]

    def forward(self, x: torch.Tensor, mask_ratio: float = 0.75, patch_mask: torch.Tensor | None = None):
        patches = self.patchify(x)
        if patch_mask is None:
            patch_mask = self.random_mask(x.shape[0], mask_ratio, x.device)
        cls, states = self.encode_patches(patches, patch_mask)
        return {
            "embedding": self.projector(cls),
            "patches": patches,
            "patch_mask": patch_mask,
            "feature_pred": self.feature_head(states),
            "reconstruction_pred": self.reconstruction_head(states),
            "mask_logits": self.mask_head(states).squeeze(-1),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        patches = self.patchify(x)
        patch_mask = x.new_zeros(x.shape[0], self.num_patches)
        cls, _ = self.encode_patches(patches, patch_mask)
        return self.projector(cls)

