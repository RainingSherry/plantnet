from __future__ import annotations

import copy
import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class GenePatchEmbedding(nn.Module):
    def __init__(self, num_genes: int, patch_size: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or patch_size <= 0 or hidden_size <= 0:
            raise ValueError("num_genes, patch_size, and hidden_size must be positive")
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.num_patches = int(math.ceil(self.num_genes / self.patch_size))
        self.pad_size = self.num_patches * self.patch_size - self.num_genes
        self.proj = nn.Linear(self.patch_size, hidden_size)
        self.dropout = nn.Dropout(dropout)

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        padded = F.pad(x, (0, self.pad_size)) if self.pad_size else x
        return padded.view(x.shape[0], self.num_patches, self.patch_size)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.dropout(self.proj(self.patchify(x)))


class TransformerBlock(nn.Module):
    def __init__(self, hidden_size: int, num_heads: int, dropout: float, mlp_ratio: int) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_size)
        self.attn = nn.MultiheadAttention(hidden_size, num_heads, dropout=dropout, batch_first=True)
        self.drop = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(hidden_size)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * mlp_ratio, hidden_size),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        q = self.norm1(x)
        attn, _ = self.attn(q, q, q, need_weights=False)
        x = x + self.drop(attn)
        return x + self.mlp(self.norm2(x))


class GenePatchEncoder(nn.Module):
    def __init__(
        self,
        num_genes: int,
        patch_size: int,
        hidden_size: int,
        depth: int,
        num_heads: int,
        dropout: float,
        mlp_ratio: int = 4,
    ) -> None:
        super().__init__()
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")
        self.patch_embed = GenePatchEmbedding(num_genes, patch_size, hidden_size, dropout)
        self.position = nn.Parameter(torch.randn(1, self.patch_embed.num_patches, hidden_size) * 0.02)
        self.blocks = nn.ModuleList([TransformerBlock(hidden_size, num_heads, dropout, mlp_ratio) for _ in range(depth)])
        self.norm = nn.LayerNorm(hidden_size)

    @property
    def num_genes(self) -> int:
        return self.patch_embed.num_genes

    @property
    def patch_size(self) -> int:
        return self.patch_embed.patch_size

    @property
    def num_patches(self) -> int:
        return self.patch_embed.num_patches

    def forward(self, x: torch.Tensor, keep_mask: torch.Tensor | None = None) -> torch.Tensor:
        tokens = self.patch_embed(x) + self.position.to(dtype=x.dtype)
        if keep_mask is not None:
            if keep_mask.ndim != 2 or keep_mask.shape != tokens.shape[:2]:
                raise ValueError(f"keep_mask must be [batch, {self.num_patches}], got {tuple(keep_mask.shape)}")
            tokens = tokens * keep_mask.to(dtype=tokens.dtype).unsqueeze(-1)
        for block in self.blocks:
            tokens = block(tokens)
        return self.norm(tokens)


class IjepaPredictor(nn.Module):
    def __init__(
        self,
        num_patches: int,
        hidden_size: int,
        predictor_size: int,
        depth: int,
        num_heads: int,
        dropout: float,
        mlp_ratio: int = 4,
    ) -> None:
        super().__init__()
        if predictor_size % num_heads != 0:
            raise ValueError("predictor_size must be divisible by num_heads")
        self.num_patches = int(num_patches)
        self.context_proj = nn.Linear(hidden_size, predictor_size)
        self.target_token = nn.Parameter(torch.zeros(1, 1, predictor_size))
        self.position = nn.Parameter(torch.randn(1, num_patches, predictor_size) * 0.02)
        self.blocks = nn.ModuleList([TransformerBlock(predictor_size, num_heads, dropout, mlp_ratio) for _ in range(depth)])
        self.norm = nn.LayerNorm(predictor_size)
        self.out = nn.Linear(predictor_size, hidden_size)
        nn.init.trunc_normal_(self.target_token, std=0.02)

    def forward(self, context_tokens: torch.Tensor, context_mask: torch.Tensor, target_mask: torch.Tensor) -> torch.Tensor:
        if context_mask.shape != target_mask.shape or context_mask.shape != context_tokens.shape[:2]:
            raise ValueError("context_tokens, context_mask, and target_mask must agree on [batch, patches]")
        context = self.context_proj(context_tokens)
        pos = self.position.to(dtype=context.dtype)
        context = context + pos
        pred_tokens = self.target_token.to(dtype=context.dtype).expand(context.shape[0], self.num_patches, -1) + pos
        tokens = context * context_mask.to(dtype=context.dtype).unsqueeze(-1)
        tokens = tokens + pred_tokens * target_mask.to(dtype=context.dtype).unsqueeze(-1)
        for block in self.blocks:
            tokens = block(tokens)
        return self.out(self.norm(tokens))


class IJEPAGeneContextScMAE(nn.Module):
    """I-JEPA-style joint-embedding prediction over contiguous gene-patch blocks."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        depth: int = 3,
        num_heads: int = 4,
        predictor_size: int = 64,
        predictor_depth: int = 2,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.context_encoder = GenePatchEncoder(num_genes, patch_size, hidden_size, depth, num_heads, dropout)
        self.target_encoder = copy.deepcopy(self.context_encoder)
        self.target_encoder.requires_grad_(False)
        self.predictor = IjepaPredictor(
            self.context_encoder.num_patches,
            hidden_size,
            predictor_size,
            predictor_depth,
            num_heads,
            dropout,
        )
        self.reconstruction_head = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, patch_size),
        )
        self.mask_classifier = nn.Linear(hidden_size, num_genes)

    @property
    def num_genes(self) -> int:
        return self.context_encoder.num_genes

    @property
    def patch_size(self) -> int:
        return self.context_encoder.patch_size

    @property
    def num_patches(self) -> int:
        return self.context_encoder.num_patches

    def unpatchify(self, patch_values: torch.Tensor) -> torch.Tensor:
        if patch_values.ndim != 3 or patch_values.shape[1:] != (self.num_patches, self.patch_size):
            raise ValueError(
                f"patch_values must be [batch, {self.num_patches}, {self.patch_size}], got {tuple(patch_values.shape)}"
            )
        return patch_values.reshape(patch_values.shape[0], -1)[:, : self.num_genes]

    @torch.no_grad()
    def target_features(self, x: torch.Tensor) -> torch.Tensor:
        self.target_encoder.eval()
        return F.layer_norm(self.target_encoder(x, keep_mask=None), (self.context_encoder.blocks[0].norm1.normalized_shape[0],))

    @torch.no_grad()
    def update_target(self, momentum: float) -> None:
        if not 0.0 <= float(momentum) <= 1.0:
            raise ValueError("momentum must be in [0, 1]")
        online = self.context_encoder.state_dict()
        target = self.target_encoder.state_dict()
        for key, target_value in target.items():
            online_value = online[key].detach()
            if torch.is_floating_point(target_value):
                target_value.mul_(float(momentum)).add_(online_value.to(dtype=target_value.dtype), alpha=1.0 - float(momentum))
            else:
                target_value.copy_(online_value)

    def forward(self, x: torch.Tensor, context_mask: torch.Tensor, target_mask: torch.Tensor) -> dict[str, torch.Tensor]:
        context_tokens = self.context_encoder(x, keep_mask=context_mask)
        predicted_targets = self.predictor(context_tokens, context_mask, target_mask)
        reconstruction = self.unpatchify(self.reconstruction_head(predicted_targets))
        pooled = context_tokens.mean(dim=1)
        return {
            "embedding": pooled,
            "predicted_targets": predicted_targets,
            "reconstruction": reconstruction,
            "mask_logits": self.mask_classifier(pooled),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        tokens = self.context_encoder(x, keep_mask=None)
        return tokens.mean(dim=1)
