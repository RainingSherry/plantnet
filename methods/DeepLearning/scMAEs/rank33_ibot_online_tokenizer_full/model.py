from __future__ import annotations

import math

import torch
from torch import nn
import torch.nn.functional as F


def sinusoidal_positional_encoding(length: int, dim: int) -> torch.Tensor:
    pe = torch.zeros(length, dim)
    position = torch.arange(0, length, dtype=torch.float32).unsqueeze(1)
    div_term = torch.exp(torch.arange(0, dim, 2, dtype=torch.float32) * (-math.log(10000.0) / max(1, dim)))
    pe[:, 0::2] = torch.sin(position * div_term)
    if dim > 1:
        pe[:, 1::2] = torch.cos(position * div_term[: pe[:, 1::2].shape[1]])
    return pe


class GeneModuleBackbone(nn.Module):
    """Small token encoder over contiguous gene modules."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        n_layers: int = 2,
        num_heads: int = 4,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.hidden_size = int(hidden_size)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.padded_genes = self.n_patches * self.patch_size
        self.patch_embed = nn.Sequential(
            nn.Linear(self.patch_size, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        self.cls_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        pe = sinusoidal_positional_encoding(self.n_patches + 1, hidden_size)
        self.register_buffer("pos_embed", pe.unsqueeze(0), persistent=False)
        layer = nn.TransformerEncoderLayer(
            d_model=hidden_size,
            nhead=num_heads,
            dim_feedforward=hidden_size * 3,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.blocks = nn.TransformerEncoder(layer, num_layers=n_layers)
        self.norm = nn.LayerNorm(hidden_size)

    def _pad(self, x: torch.Tensor) -> torch.Tensor:
        if self.padded_genes == self.num_genes:
            return x
        return F.pad(x, (0, self.padded_genes - self.num_genes))

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        patches = self._pad(x).view(x.shape[0], self.n_patches, self.patch_size)
        patch_tokens = self.patch_embed(patches)
        cls = self.cls_token.expand(x.shape[0], -1, -1)
        tokens = torch.cat([cls, patch_tokens], dim=1) + self.pos_embed.to(x.device)
        tokens = self.norm(self.blocks(tokens))
        return tokens[:, 0], tokens[:, 1:]


class ProjectionHead(nn.Module):
    def __init__(self, hidden_size: int, out_dim: int, dropout: float = 0.05):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, hidden_size),
            nn.GELU(),
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, out_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class iBOTBranch(nn.Module):
    def __init__(
        self,
        num_genes: int,
        patch_size: int,
        hidden_size: int,
        cls_out_dim: int,
        patch_out_dim: int,
        n_layers: int,
        num_heads: int,
        dropout: float,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.padded_genes = self.n_patches * self.patch_size
        self.backbone = GeneModuleBackbone(num_genes, patch_size, hidden_size, n_layers, num_heads, dropout)
        self.cls_head = ProjectionHead(hidden_size, cls_out_dim, dropout)
        self.patch_head = ProjectionHead(hidden_size, patch_out_dim, dropout)
        self.expr_head = nn.Linear(hidden_size, patch_size)
        self.mask_head = nn.Linear(hidden_size, patch_size)
        self.latent_norm = nn.LayerNorm(hidden_size)

    def _unpad(self, x: torch.Tensor) -> torch.Tensor:
        return x.reshape(x.shape[0], self.padded_genes)[:, : self.num_genes]

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        cls, patch_tokens = self.backbone(x)
        reconstruction = self._unpad(self.expr_head(patch_tokens))
        mask_logits = self._unpad(self.mask_head(patch_tokens))
        return {
            "latent": self.latent_norm(cls),
            "module_tokens": patch_tokens,
            "class_logits": self.cls_head(cls),
            "patch_logits": self.patch_head(patch_tokens),
            "reconstruction": reconstruction,
            "mask_logits": mask_logits,
        }


class iBOTOnlineTokenizerScMAE(nn.Module):
    """EMA-teacher online tokenizer combined with scMAE reconstruction heads."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        cls_out_dim: int = 128,
        patch_out_dim: int = 128,
        n_layers: int = 2,
        num_heads: int = 4,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.padded_genes = self.n_patches * self.patch_size
        self.student = iBOTBranch(num_genes, patch_size, hidden_size, cls_out_dim, patch_out_dim, n_layers, num_heads, dropout)
        self.teacher = iBOTBranch(num_genes, patch_size, hidden_size, cls_out_dim, patch_out_dim, n_layers, num_heads, dropout)
        self.reset_teacher()
        self.apply(self._init_weights)
        self.reset_teacher()

    def _init_weights(self, module: nn.Module) -> None:
        if isinstance(module, nn.Linear):
            nn.init.xavier_uniform_(module.weight)
            if module.bias is not None:
                nn.init.zeros_(module.bias)
        elif isinstance(module, nn.LayerNorm):
            nn.init.ones_(module.weight)
            nn.init.zeros_(module.bias)
        elif isinstance(module, GeneModuleBackbone):
            nn.init.normal_(module.cls_token, std=0.02)

    @torch.no_grad()
    def reset_teacher(self) -> None:
        self.teacher.load_state_dict(self.student.state_dict(), strict=True)
        for param in self.teacher.parameters():
            param.requires_grad_(False)

    @torch.no_grad()
    def update_teacher(self, momentum: float) -> None:
        for student_param, teacher_param in zip(self.student.parameters(), self.teacher.parameters()):
            teacher_param.data.mul_(float(momentum)).add_(student_param.detach().data, alpha=1.0 - float(momentum))

    def forward_student(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        return self.student(x)

    @torch.no_grad()
    def forward_teacher(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        self.teacher.eval()
        return self.teacher(x)

    @torch.no_grad()
    def feature(self, x: torch.Tensor, use_teacher: bool = False) -> torch.Tensor:
        branch = self.teacher if use_teacher else self.student
        return branch(x)["latent"]

    def mask_view(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask

    def weak_view(self, x: torch.Tensor, noise_std: float = 0.03, dropout_prob: float = 0.05) -> torch.Tensor:
        y = x
        if noise_std > 0:
            y = y + torch.randn_like(y) * float(noise_std)
        if dropout_prob > 0:
            keep = (torch.rand_like(y) >= float(dropout_prob)).float()
            y = y * keep
        return y

    def gene_mask_to_module_mask(self, mask: torch.Tensor) -> torch.Tensor:
        if self.padded_genes != self.num_genes:
            mask = F.pad(mask, (0, self.padded_genes - self.num_genes))
        return mask.view(mask.shape[0], self.n_patches, self.patch_size).mean(dim=-1)
