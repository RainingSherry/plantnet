from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class ExpressionMatrixPatcher(nn.Module):
    """2D input matrix embedding for grouped gene expression sequences."""

    def __init__(self, num_genes: int, matrix_patch_size: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or matrix_patch_size <= 0 or hidden_size <= 0:
            raise ValueError("num_genes, matrix_patch_size, and hidden_size must be positive")
        self.num_genes = int(num_genes)
        self.matrix_patch_size = int(matrix_patch_size)
        raw_side = int(math.ceil(math.sqrt(self.num_genes)))
        self.matrix_side = int(math.ceil(raw_side / self.matrix_patch_size) * self.matrix_patch_size)
        self.padded_genes = self.matrix_side * self.matrix_side
        self.num_patch_rows = self.matrix_side // self.matrix_patch_size
        self.num_patches = self.num_patch_rows * self.num_patch_rows
        self.raw_patch_dim = self.matrix_patch_size * self.matrix_patch_size
        self.projection = nn.Conv2d(
            1,
            hidden_size,
            kernel_size=self.matrix_patch_size,
            stride=self.matrix_patch_size,
            bias=True,
        )
        self.position = nn.Parameter(torch.randn(1, self.num_patches, hidden_size) * 0.02)
        self.dropout = nn.Dropout(dropout)

    def to_matrix(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        padded = F.pad(x, (0, self.padded_genes - self.num_genes))
        return padded.view(x.shape[0], 1, self.matrix_side, self.matrix_side)

    def raw_patches(self, x: torch.Tensor) -> torch.Tensor:
        matrix = self.to_matrix(x)
        patches = F.unfold(matrix, kernel_size=self.matrix_patch_size, stride=self.matrix_patch_size)
        return patches.transpose(1, 2).contiguous()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        tokens = self.projection(self.to_matrix(x)).flatten(2).transpose(1, 2)
        return self.dropout(tokens + self.position.to(dtype=tokens.dtype))


class SequenceLevelEncoder(nn.Module):
    """Contrastive-sc style short expression sequence encoder used as a frozen target network."""

    def __init__(self, raw_patch_dim: int, target_dim: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        self.raw_patch_dim = int(raw_patch_dim)
        self.target_dim = int(target_dim)
        self.net = nn.Sequential(
            nn.LayerNorm(raw_patch_dim),
            nn.Linear(raw_patch_dim, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, target_dim),
        )

    def forward(self, patches: torch.Tensor) -> torch.Tensor:
        if patches.ndim != 3 or patches.shape[2] != self.raw_patch_dim:
            raise ValueError(f"patches must be [batch, patches, {self.raw_patch_dim}], got {tuple(patches.shape)}")
        return F.normalize(self.net(patches), dim=-1)


def _build_transformer_encoder(hidden_size: int, depth: int, num_heads: int, dropout: float) -> nn.TransformerEncoder:
    if hidden_size % num_heads != 0:
        raise ValueError("hidden_size must be divisible by num_heads")
    layer = nn.TransformerEncoderLayer(
        d_model=hidden_size,
        nhead=num_heads,
        dim_feedforward=hidden_size * 4,
        dropout=dropout,
        activation="gelu",
        batch_first=True,
        norm_first=True,
    )
    return nn.TransformerEncoder(layer, num_layers=depth)


class MaskSCClusteringModel(nn.Module):
    """mask-sc encoder/decoder with expression matrix patches and sequence-guided targets."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        decoder_size: int = 128,
        target_dim: int = 64,
        encoder_depth: int = 2,
        decoder_depth: int = 2,
        num_heads: int = 4,
        matrix_patch_size: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.decoder_size = int(decoder_size)
        self.target_dim = int(target_dim)
        self.patcher = ExpressionMatrixPatcher(num_genes, matrix_patch_size, hidden_size, dropout)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.encoder = _build_transformer_encoder(hidden_size, encoder_depth, num_heads, dropout)
        self.enc_norm = nn.LayerNorm(hidden_size)
        self.encoder_to_decoder = nn.Linear(hidden_size, decoder_size)
        self.decoder_position = nn.Parameter(torch.randn(1, self.patcher.num_patches, decoder_size) * 0.02)
        self.decoder_cls_position = nn.Parameter(torch.zeros(1, 1, decoder_size))
        self.mask_token = nn.Parameter(torch.zeros(1, 1, decoder_size))
        self.decoder = _build_transformer_encoder(decoder_size, decoder_depth, num_heads, dropout)
        self.pred_head = nn.Sequential(
            nn.LayerNorm(decoder_size),
            nn.Linear(decoder_size, decoder_size),
            nn.GELU(),
            nn.Linear(decoder_size, target_dim),
        )
        self.embedding_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.Tanh())

    @property
    def num_patches(self) -> int:
        return self.patcher.num_patches

    def _gather_tokens(self, tokens: torch.Tensor, indices: torch.Tensor) -> torch.Tensor:
        if indices.ndim != 2 or indices.shape[0] != tokens.shape[0]:
            raise ValueError("indices must be [batch, selected_patches]")
        gather = indices.to(device=tokens.device).unsqueeze(-1).expand(-1, -1, tokens.shape[-1])
        return torch.gather(tokens, dim=1, index=gather)

    def encode_visible(self, x: torch.Tensor, visible_indices: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        tokens = self.patcher(x)
        visible = self._gather_tokens(tokens, visible_indices)
        cls = self.cls_token.expand(x.shape[0], -1, -1).to(dtype=tokens.dtype)
        encoded = self.encoder(torch.cat([cls, visible], dim=1))
        encoded = self.enc_norm(encoded)
        return encoded[:, 0], encoded[:, 1:]

    def forward(
        self,
        x: torch.Tensor,
        visible_indices: torch.Tensor,
        masked_indices: torch.Tensor,
    ) -> dict[str, torch.Tensor]:
        cls_encoded, visible_encoded = self.encode_visible(x, visible_indices)
        dec_visible = self.encoder_to_decoder(visible_encoded)
        dec_cls = self.encoder_to_decoder(cls_encoded).unsqueeze(1) + self.decoder_cls_position
        visible_pos = self._gather_tokens(self.decoder_position.expand(x.shape[0], -1, -1), visible_indices)
        masked_pos = self._gather_tokens(self.decoder_position.expand(x.shape[0], -1, -1), masked_indices)
        dec_visible = dec_visible + visible_pos
        dec_masked = self.mask_token.expand(x.shape[0], masked_indices.shape[1], -1) + masked_pos
        decoded = self.decoder(torch.cat([dec_cls, dec_visible, dec_masked], dim=1))
        masked_decoded = decoded[:, 1 + visible_indices.shape[1] :]
        prediction = self.pred_head(masked_decoded)
        return {
            "embedding": self.embedding_head(visible_encoded.mean(dim=1)),
            "masked_prediction": F.normalize(prediction, dim=-1),
            "visible_tokens": visible_encoded,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        full_indices = torch.arange(self.num_patches, device=x.device).view(1, -1).expand(x.shape[0], -1)
        cls_encoded, tokens = self.encode_visible(x, full_indices)
        return self.embedding_head(tokens.mean(dim=1) + 0.0 * cls_encoded)
