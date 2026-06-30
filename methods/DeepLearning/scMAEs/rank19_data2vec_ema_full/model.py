from __future__ import annotations

import copy
import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class ExpressionPatchEmbedding(nn.Module):
    def __init__(self, num_genes: int, patch_size: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or patch_size <= 0 or hidden_size <= 0:
            raise ValueError("num_genes, patch_size, and hidden_size must be positive")
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.num_patches = int(math.ceil(self.num_genes / self.patch_size))
        self.pad_size = self.num_patches * self.patch_size - self.num_genes
        self.projection = nn.Linear(self.patch_size, hidden_size)
        self.dropout = nn.Dropout(dropout)

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        padded = F.pad(x, (0, self.pad_size)) if self.pad_size else x
        return padded.view(x.shape[0], self.num_patches, self.patch_size)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.dropout(self.projection(self.patchify(x)))


class TransformerBlock(nn.Module):
    def __init__(self, hidden_size: int, num_heads: int, dropout: float, mlp_ratio: int) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_size)
        self.attn = nn.MultiheadAttention(hidden_size, num_heads, dropout=dropout, batch_first=True)
        self.dropout = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(hidden_size)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * mlp_ratio, hidden_size),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        attn_out, _ = self.attn(self.norm1(x), self.norm1(x), self.norm1(x), need_weights=False)
        x = x + self.dropout(attn_out)
        return x + self.mlp(self.norm2(x))


class Data2VecEncoder(nn.Module):
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
        self.patch_embed = ExpressionPatchEmbedding(num_genes, patch_size, hidden_size, dropout)
        self.cls_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.mask_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.position = nn.Parameter(torch.randn(1, self.patch_embed.num_patches + 1, hidden_size) * 0.02)
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

    def forward(
        self,
        x: torch.Tensor,
        patch_mask: torch.Tensor | None = None,
        return_layers: bool = False,
    ) -> tuple[torch.Tensor, list[torch.Tensor]]:
        tokens = self.patch_embed(x)
        if patch_mask is not None:
            if patch_mask.ndim != 2 or patch_mask.shape != tokens.shape[:2]:
                raise ValueError(f"patch_mask must be [batch, {self.num_patches}], got {tuple(patch_mask.shape)}")
            mask = patch_mask.to(dtype=tokens.dtype).unsqueeze(-1)
            tokens = tokens * (1.0 - mask) + self.mask_token.to(dtype=tokens.dtype) * mask
        cls = self.cls_token.expand(tokens.shape[0], -1, -1).to(dtype=tokens.dtype)
        h = torch.cat([cls, tokens], dim=1) + self.position.to(dtype=tokens.dtype)
        layers: list[torch.Tensor] = []
        for block in self.blocks:
            h = block(h)
            if return_layers:
                layers.append(self.norm(h))
        h = self.norm(h)
        if not return_layers:
            layers = [h]
        return h, layers


class Data2VecEMAExpressionScMAE(nn.Module):
    """data2vec-style EMA teacher/student for expression patches."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        depth: int = 3,
        num_heads: int = 4,
        dropout: float = 0.1,
        average_top_k_layers: int = 2,
    ) -> None:
        super().__init__()
        if average_top_k_layers <= 0 or average_top_k_layers > depth:
            raise ValueError("average_top_k_layers must be in [1, depth]")
        self.student = Data2VecEncoder(num_genes, patch_size, hidden_size, depth, num_heads, dropout)
        self.teacher = copy.deepcopy(self.student)
        self.teacher.requires_grad_(False)
        self.average_top_k_layers = int(average_top_k_layers)
        self.final_proj = nn.Linear(hidden_size, hidden_size)
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.patch_decoder = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, patch_size),
        )

    @property
    def num_genes(self) -> int:
        return self.student.num_genes

    @property
    def patch_size(self) -> int:
        return self.student.patch_size

    @property
    def num_patches(self) -> int:
        return self.student.num_patches

    def unpatchify(self, patch_values: torch.Tensor) -> torch.Tensor:
        if patch_values.ndim != 3 or patch_values.shape[1] != self.num_patches or patch_values.shape[2] != self.patch_size:
            raise ValueError(
                f"patch_values must be [batch, {self.num_patches}, {self.patch_size}], got {tuple(patch_values.shape)}"
            )
        return patch_values.reshape(patch_values.shape[0], -1)[:, : self.num_genes]

    @torch.no_grad()
    def teacher_targets(self, x: torch.Tensor) -> torch.Tensor:
        self.teacher.eval()
        _, layers = self.teacher(x, patch_mask=None, return_layers=True)
        top_layers = layers[-self.average_top_k_layers :]
        normalized = [F.layer_norm(layer[:, 1:], layer[:, 1:].shape[-1:]) for layer in top_layers]
        return torch.stack(normalized, dim=0).mean(dim=0)

    @torch.no_grad()
    def update_teacher(self, decay: float) -> None:
        if not 0.0 <= float(decay) <= 1.0:
            raise ValueError("decay must be in [0, 1]")
        student_state = self.student.state_dict()
        teacher_state = self.teacher.state_dict()
        for key, teacher_value in teacher_state.items():
            student_value = student_state[key].detach()
            if torch.is_floating_point(teacher_value):
                teacher_value.mul_(float(decay)).add_(student_value.to(dtype=teacher_value.dtype), alpha=1.0 - float(decay))
            else:
                teacher_value.copy_(student_value)

    def forward(self, x: torch.Tensor, patch_mask: torch.Tensor) -> dict[str, torch.Tensor]:
        student_h, _ = self.student(x, patch_mask=patch_mask, return_layers=False)
        patch_tokens = student_h[:, 1:]
        prediction = self.final_proj(patch_tokens)
        mask_logits = self.mask_predictor(student_h[:, 0])
        reconstruction = self.unpatchify(self.patch_decoder(patch_tokens))
        return {
            "embedding": student_h[:, 0],
            "patch_prediction": prediction,
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        h, _ = self.student(x, patch_mask=None, return_layers=False)
        return h[:, 0]
