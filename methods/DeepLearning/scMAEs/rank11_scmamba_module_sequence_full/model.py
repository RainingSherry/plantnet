from __future__ import annotations

import torch
from torch import nn


class ModuleMixerBlock(nn.Module):
    """Lightweight SSM-style gated depthwise sequence mixer."""

    def __init__(self, hidden_size: int, kernel_size: int = 7, dropout: float = 0.05):
        super().__init__()
        self.norm = nn.LayerNorm(hidden_size)
        self.depthwise = nn.Conv1d(hidden_size, hidden_size, kernel_size, padding=kernel_size // 2, groups=hidden_size)
        self.gate = nn.Linear(hidden_size, hidden_size)
        self.mlp = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, hidden_size),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        residual = x
        z = self.norm(x)
        mixed = self.depthwise(z.transpose(1, 2)).transpose(1, 2)
        x = residual + mixed * torch.sigmoid(self.gate(z))
        return x + self.mlp(x)


class ScMambaModuleScMAE(nn.Module):
    """Patch-tokenized scMAE with module-sequence mixer."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 25,
        hidden_size: int = 128,
        n_layers: int = 3,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.hidden_size = int(hidden_size)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.padded_genes = self.n_patches * self.patch_size
        self.patch_embed = nn.Sequential(
            nn.Linear(self.patch_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.GELU(),
        )
        self.pos_embed = nn.Parameter(torch.zeros(1, self.n_patches, self.hidden_size))
        self.blocks = nn.ModuleList([ModuleMixerBlock(self.hidden_size, dropout=dropout) for _ in range(n_layers)])
        self.pool_norm = nn.LayerNorm(self.hidden_size)
        self.decoder = nn.Sequential(
            nn.Linear(self.hidden_size, 256),
            nn.Mish(),
            nn.Linear(256, self.num_genes),
        )
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.module_decoder = nn.Linear(self.hidden_size, self.n_patches)

    def _pad(self, x: torch.Tensor) -> torch.Tensor:
        if self.padded_genes == self.num_genes:
            return x
        return torch.nn.functional.pad(x, (0, self.padded_genes - self.num_genes))

    def encode_tokens(self, x: torch.Tensor) -> torch.Tensor:
        patches = self._pad(x).view(x.shape[0], self.n_patches, self.patch_size)
        tokens = self.patch_embed(patches) + self.pos_embed
        for block in self.blocks:
            tokens = block(tokens)
        return tokens

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        tokens = self.encode_tokens(x)
        latent = self.pool_norm(tokens.mean(dim=1))
        return {
            "tokens": tokens,
            "latent": latent,
            "mask_logits": self.mask_predictor(latent),
            "reconstruction": self.decoder(latent),
            "module_reconstruction": self.module_decoder(latent),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.forward(x)["latent"]

    def block_mask(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        b = x.shape[0]
        patch_mask = (torch.rand(b, self.n_patches, device=x.device) < float(mask_prob))
        empty = patch_mask.sum(dim=1) == 0
        if bool(empty.any()):
            patch_mask[empty, torch.randint(0, self.n_patches, (int(empty.sum()),), device=x.device)] = True
        gene_mask = patch_mask.repeat_interleave(self.patch_size, dim=1)[:, : self.num_genes].float()
        return x.masked_fill(gene_mask.bool(), 0.0), gene_mask

