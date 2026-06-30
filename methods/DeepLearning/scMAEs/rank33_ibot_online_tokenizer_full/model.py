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
        q = self.norm1(x)
        attn_out, _ = self.attn(q, q, q, need_weights=False)
        x = x + self.dropout(attn_out)
        return x + self.mlp(self.norm2(x))


class ExpressionIBOTEncoder(nn.Module):
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

    def forward(self, x: torch.Tensor, patch_mask: torch.Tensor | None = None) -> torch.Tensor:
        tokens = self.patch_embed(x)
        if patch_mask is not None:
            if patch_mask.ndim != 2 or patch_mask.shape != tokens.shape[:2]:
                raise ValueError(f"patch_mask must be [batch, {self.num_patches}], got {tuple(patch_mask.shape)}")
            mask = patch_mask.to(dtype=tokens.dtype).unsqueeze(-1)
            tokens = tokens * (1.0 - mask) + self.mask_token.to(dtype=tokens.dtype) * mask
        cls = self.cls_token.expand(tokens.shape[0], -1, -1).to(dtype=tokens.dtype)
        h = torch.cat([cls, tokens], dim=1) + self.position.to(dtype=tokens.dtype)
        for block in self.blocks:
            h = block(h)
        return self.norm(h)


class IBOTProjectionHead(nn.Module):
    """Shared projection head used for both class and expression patch tokens."""

    def __init__(self, hidden_size: int, bottleneck_size: int, out_dim: int, dropout: float) -> None:
        super().__init__()
        if bottleneck_size <= 0 or out_dim <= 1:
            raise ValueError("bottleneck_size must be positive and out_dim must be > 1")
        self.mlp = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, bottleneck_size),
        )
        self.last_layer = nn.utils.parametrizations.weight_norm(nn.Linear(bottleneck_size, out_dim, bias=False))

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        z = self.mlp(x)
        z = F.normalize(z, dim=-1)
        return self.last_layer(z)


class IBOTOnlineTokenizerScMAE(nn.Module):
    """iBOT-style online tokenizer and self-distillation for expression patches."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        depth: int = 3,
        num_heads: int = 4,
        dropout: float = 0.1,
        out_dim: int = 256,
        bottleneck_size: int = 64,
    ) -> None:
        super().__init__()
        self.student = ExpressionIBOTEncoder(num_genes, patch_size, hidden_size, depth, num_heads, dropout)
        self.teacher = copy.deepcopy(self.student)
        self.teacher.requires_grad_(False)
        self.projection_head = IBOTProjectionHead(hidden_size, bottleneck_size, out_dim, dropout)
        self.teacher_projection_head = copy.deepcopy(self.projection_head)
        self.teacher_projection_head.requires_grad_(False)
        self.patch_decoder = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, patch_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, num_genes)
        self.embedding_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.Tanh())
        self.register_buffer("cls_center", torch.zeros(1, out_dim))
        self.register_buffer("patch_center", torch.zeros(1, 1, out_dim))

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

    def _student_view(self, x: torch.Tensor, patch_mask: torch.Tensor) -> dict[str, torch.Tensor]:
        h = self.student(x, patch_mask)
        patch_tokens = h[:, 1:]
        patch_values = self.patch_decoder(patch_tokens)
        return {
            "embedding": self.embedding_head(h[:, 0]),
            "cls_logits": self.projection_head(h[:, 0]),
            "patch_logits": self.projection_head(patch_tokens),
            "patch_tokens": patch_tokens,
            "reconstruction": self.unpatchify(patch_values),
            "mask_logits": self.mask_predictor(h[:, 0]),
        }

    @torch.no_grad()
    def _teacher_view(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        self.teacher.eval()
        self.teacher_projection_head.eval()
        h = self.teacher(x, patch_mask=None)
        patch_tokens = h[:, 1:]
        return {
            "cls_logits": self.teacher_projection_head(h[:, 0]),
            "patch_logits": self.teacher_projection_head(patch_tokens),
        }

    def forward(
        self,
        view1: torch.Tensor,
        mask1: torch.Tensor,
        view2: torch.Tensor,
        mask2: torch.Tensor,
    ) -> dict[str, dict[str, torch.Tensor]]:
        return {
            "student1": self._student_view(view1, mask1),
            "student2": self._student_view(view2, mask2),
            "teacher1": self._teacher_view(view1),
            "teacher2": self._teacher_view(view2),
        }

    @torch.no_grad()
    def update_teacher(self, decay: float) -> None:
        if not 0.0 <= float(decay) <= 1.0:
            raise ValueError("decay must be in [0, 1]")
        for student_module, teacher_module in (
            (self.student, self.teacher),
            (self.projection_head, self.teacher_projection_head),
        ):
            student_state = student_module.state_dict()
            teacher_state = teacher_module.state_dict()
            for key, teacher_value in teacher_state.items():
                student_value = student_state[key].detach()
                if torch.is_floating_point(teacher_value):
                    teacher_value.mul_(float(decay)).add_(student_value.to(dtype=teacher_value.dtype), alpha=1.0 - float(decay))
                else:
                    teacher_value.copy_(student_value)

    @torch.no_grad()
    def update_centers(self, teacher1: dict[str, torch.Tensor], teacher2: dict[str, torch.Tensor], momentum: float) -> None:
        if not 0.0 <= float(momentum) <= 1.0:
            raise ValueError("momentum must be in [0, 1]")
        cls_batch_center = torch.cat([teacher1["cls_logits"], teacher2["cls_logits"]], dim=0).mean(dim=0, keepdim=True)
        patch_batch_center = torch.cat([teacher1["patch_logits"], teacher2["patch_logits"]], dim=0).mean(dim=(0, 1), keepdim=True)
        self.cls_center.mul_(float(momentum)).add_(cls_batch_center.to(dtype=self.cls_center.dtype), alpha=1.0 - float(momentum))
        self.patch_center.mul_(float(momentum)).add_(patch_batch_center.to(dtype=self.patch_center.dtype), alpha=1.0 - float(momentum))

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        h = self.student(x, patch_mask=None)
        return self.embedding_head(h[:, 0])
