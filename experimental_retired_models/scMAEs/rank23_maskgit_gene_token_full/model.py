from __future__ import annotations

import math

import torch
from torch import nn


def sinusoidal_positional_encoding(length: int, dim: int) -> torch.Tensor:
    pe = torch.zeros(length, dim)
    position = torch.arange(0, length, dtype=torch.float32).unsqueeze(1)
    div_term = torch.exp(torch.arange(0, dim, 2, dtype=torch.float32) * (-math.log(10000.0) / max(1, dim)))
    pe[:, 0::2] = torch.sin(position * div_term)
    if dim > 1:
        pe[:, 1::2] = torch.cos(position * div_term[: pe[:, 1::2].shape[1]])
    return pe


class MaskGITGeneTokenScMAE(nn.Module):
    """Bidirectional module-token transformer for MaskGIT-style scMAE."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        token_bins: int = 8,
        n_layers: int = 3,
        num_heads: int = 4,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.hidden_size = int(hidden_size)
        self.token_bins = int(token_bins)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.padded_genes = self.n_patches * self.patch_size
        self.patch_embed = nn.Sequential(
            nn.Linear(self.patch_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        pe = sinusoidal_positional_encoding(self.n_patches, self.hidden_size)
        self.register_buffer("pos_embed", pe.unsqueeze(0), persistent=False)
        layer = nn.TransformerEncoderLayer(
            d_model=self.hidden_size,
            nhead=num_heads,
            dim_feedforward=self.hidden_size * 3,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.blocks = nn.TransformerEncoder(layer, num_layers=n_layers)
        self.norm = nn.LayerNorm(self.hidden_size)
        self.pool_norm = nn.LayerNorm(self.hidden_size)
        self.expr_head = nn.Linear(self.hidden_size, self.patch_size)
        self.mask_head = nn.Linear(self.hidden_size, self.patch_size)
        self.replaced_head = nn.Linear(self.hidden_size, self.patch_size)
        self.token_head = nn.Linear(self.hidden_size, self.patch_size * self.token_bins)
        self.apply(self._init_weights)

    def _init_weights(self, module: nn.Module) -> None:
        if isinstance(module, nn.Linear):
            nn.init.xavier_uniform_(module.weight)
            if module.bias is not None:
                nn.init.zeros_(module.bias)
        elif isinstance(module, nn.LayerNorm):
            nn.init.ones_(module.weight)
            nn.init.zeros_(module.bias)

    def _pad(self, x: torch.Tensor) -> torch.Tensor:
        if self.padded_genes == self.num_genes:
            return x
        return torch.nn.functional.pad(x, (0, self.padded_genes - self.num_genes))

    def _unpad(self, x: torch.Tensor) -> torch.Tensor:
        return x.reshape(x.shape[0], self.padded_genes)[:, : self.num_genes]

    def encode_tokens(self, x: torch.Tensor) -> torch.Tensor:
        patches = self._pad(x).view(x.shape[0], self.n_patches, self.patch_size)
        tokens = self.patch_embed(patches) + self.pos_embed.to(x.device)
        return self.norm(self.blocks(tokens))

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        tokens = self.encode_tokens(x)
        latent = self.pool_norm(tokens.mean(dim=1))
        reconstruction = self._unpad(self.expr_head(tokens))
        mask_logits = self._unpad(self.mask_head(tokens))
        replaced_logits = self._unpad(self.replaced_head(tokens))
        token_logits = self.token_head(tokens).view(x.shape[0], self.n_patches, self.patch_size, self.token_bins)
        token_logits = token_logits.reshape(x.shape[0], self.padded_genes, self.token_bins)[:, : self.num_genes, :]
        confidence = torch.softmax(token_logits.detach(), dim=-1).amax(dim=-1)
        return {
            "latent": latent,
            "module_tokens": tokens,
            "reconstruction": reconstruction,
            "mask_logits": mask_logits,
            "replaced_logits": replaced_logits,
            "token_logits": token_logits,
            "token_confidence": confidence,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.forward(x)["latent"]

    def maskgit_corrupt(
        self,
        x: torch.Tensor,
        mask_ratio: float,
        replace_prob: float,
        mask_value: float = 0.0,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        bsz, genes = x.shape
        n_mask = max(1, min(genes - 1, int(round(genes * float(mask_ratio)))))
        noise = torch.rand(bsz, genes, device=x.device)
        ids = torch.argsort(noise, dim=1)[:, :n_mask]
        mask = torch.zeros(bsz, genes, device=x.device)
        mask.scatter_(1, ids, 1.0)
        replacement = x[torch.randperm(bsz, device=x.device)] if bsz > 1 else x
        replace_gate = ((torch.rand_like(x) < float(replace_prob)).float() * mask).float()
        corrupted = torch.where(replace_gate.bool(), replacement, x)
        corrupted = torch.where((mask - replace_gate).bool(), torch.full_like(x, float(mask_value)), corrupted)
        return corrupted, mask, replace_gate
