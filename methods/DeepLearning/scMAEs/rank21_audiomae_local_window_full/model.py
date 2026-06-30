from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


def compute_gene_grid(num_genes: int, patch_size: int) -> tuple[int, int, int, int]:
    if num_genes <= 0 or patch_size <= 0:
        raise ValueError("num_genes and patch_size must be positive")
    base_h = int(math.ceil(math.sqrt(num_genes)))
    base_w = int(math.ceil(num_genes / base_h))
    grid_h = int(math.ceil(base_h / patch_size) * patch_size)
    grid_w = int(math.ceil(base_w / patch_size) * patch_size)
    return grid_h, grid_w, grid_h // patch_size, grid_w // patch_size


class TransformerBlock(nn.Module):
    def __init__(self, dim: int, num_heads: int, dropout: float, mlp_ratio: int = 4) -> None:
        super().__init__()
        if dim % num_heads != 0:
            raise ValueError("dim must be divisible by num_heads")
        self.norm1 = nn.LayerNorm(dim)
        self.attn = nn.MultiheadAttention(dim, num_heads, dropout=dropout, batch_first=True)
        self.drop = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, dim * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim * mlp_ratio, dim),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        q = self.norm1(x)
        attn, _ = self.attn(q, q, q, need_weights=False)
        x = x + self.drop(attn)
        return x + self.mlp(self.norm2(x))


def window_partition(x: torch.Tensor, window_size: tuple[int, int]) -> torch.Tensor:
    bsz, height, width, channels = x.shape
    wh, ww = window_size
    if height % wh != 0 or width % ww != 0:
        raise ValueError("height and width must be divisible by window_size")
    x = x.view(bsz, height // wh, wh, width // ww, ww, channels)
    return x.permute(0, 1, 3, 2, 4, 5).reshape(-1, wh * ww, channels)


def window_reverse(windows: torch.Tensor, window_size: tuple[int, int], height: int, width: int, batch_size: int) -> torch.Tensor:
    wh, ww = window_size
    x = windows.view(batch_size, height // wh, width // ww, wh, ww, -1)
    return x.permute(0, 1, 3, 2, 4, 5).reshape(batch_size, height, width, -1)


def shifted_window_mask(
    height: int,
    width: int,
    window_size: tuple[int, int],
    shift_size: tuple[int, int],
    device: torch.device,
) -> torch.Tensor | None:
    sh, sw = shift_size
    if sh == 0 and sw == 0:
        return None
    wh, ww = window_size
    img_mask = torch.zeros((1, height, width, 1), device=device)
    h_slices = (slice(0, -wh), slice(-wh, -sh), slice(-sh, None)) if sh > 0 else (slice(0, height),)
    w_slices = (slice(0, -ww), slice(-ww, -sw), slice(-sw, None)) if sw > 0 else (slice(0, width),)
    count = 0
    for h_slice in h_slices:
        for w_slice in w_slices:
            img_mask[:, h_slice, w_slice, :] = count
            count += 1
    mask_windows = window_partition(img_mask, window_size).squeeze(-1)
    attn_mask = mask_windows.unsqueeze(1) - mask_windows.unsqueeze(2)
    return attn_mask.masked_fill(attn_mask != 0, -100.0).masked_fill(attn_mask == 0, 0.0)


class WindowAttention(nn.Module):
    def __init__(self, dim: int, num_heads: int, dropout: float) -> None:
        super().__init__()
        if dim % num_heads != 0:
            raise ValueError("dim must be divisible by num_heads")
        self.num_heads = int(num_heads)
        self.head_dim = dim // self.num_heads
        self.scale = self.head_dim ** -0.5
        self.qkv = nn.Linear(dim, dim * 3)
        self.proj = nn.Linear(dim, dim)
        self.attn_drop = nn.Dropout(dropout)
        self.proj_drop = nn.Dropout(dropout)

    def forward(self, x: torch.Tensor, attn_mask: torch.Tensor | None) -> torch.Tensor:
        batch_windows, tokens, dim = x.shape
        qkv = self.qkv(x).view(batch_windows, tokens, 3, self.num_heads, self.head_dim)
        qkv = qkv.permute(2, 0, 3, 1, 4)
        q, k, v = qkv[0], qkv[1], qkv[2]
        attn = (q * self.scale) @ k.transpose(-2, -1)
        if attn_mask is not None:
            num_windows = attn_mask.shape[0]
            attn = attn.view(batch_windows // num_windows, num_windows, self.num_heads, tokens, tokens)
            attn = attn + attn_mask.unsqueeze(0).unsqueeze(2)
            attn = attn.view(-1, self.num_heads, tokens, tokens)
        attn = self.attn_drop(attn.softmax(dim=-1))
        out = (attn @ v).transpose(1, 2).reshape(batch_windows, tokens, dim)
        return self.proj_drop(self.proj(out))


class LocalWindowBlock(nn.Module):
    """Swin-style local window decoder block for AudioMAE-like reconstruction."""

    def __init__(
        self,
        dim: int,
        num_heads: int,
        window_size: tuple[int, int],
        shift_size: tuple[int, int],
        dropout: float,
        mlp_ratio: int = 4,
    ) -> None:
        super().__init__()
        self.window_size = window_size
        self.shift_size = shift_size
        self.norm1 = nn.LayerNorm(dim)
        self.attn = WindowAttention(dim, num_heads, dropout)
        self.drop = nn.Dropout(dropout)
        self.norm2 = nn.LayerNorm(dim)
        self.mlp = nn.Sequential(
            nn.Linear(dim, dim * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(dim * mlp_ratio, dim),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor, grid_hw: tuple[int, int]) -> torch.Tensor:
        bsz, tokens, channels = x.shape
        height, width = grid_hw
        if tokens != height * width:
            raise ValueError(f"tokens={tokens} does not match grid {grid_hw}")
        shortcut = x
        x_grid = self.norm1(x).view(bsz, height, width, channels)
        wh = min(self.window_size[0], height)
        ww = min(self.window_size[1], width)
        pad_h = (wh - height % wh) % wh
        pad_w = (ww - width % ww) % ww
        x_grid = F.pad(x_grid, (0, 0, 0, pad_w, 0, pad_h))
        padded_h, padded_w = height + pad_h, width + pad_w
        sh = min(self.shift_size[0], max(0, wh - 1))
        sw = min(self.shift_size[1], max(0, ww - 1))
        if sh > 0 or sw > 0:
            shifted = torch.roll(x_grid, shifts=(-sh, -sw), dims=(1, 2))
            attn_mask = shifted_window_mask(padded_h, padded_w, (wh, ww), (sh, sw), x.device)
        else:
            shifted = x_grid
            attn_mask = None
        windows = window_partition(shifted, (wh, ww))
        attended = self.attn(windows, attn_mask)
        shifted = window_reverse(attended, (wh, ww), padded_h, padded_w, bsz)
        if sh > 0 or sw > 0:
            x_grid = torch.roll(shifted, shifts=(sh, sw), dims=(1, 2))
        else:
            x_grid = shifted
        x_grid = x_grid[:, :height, :width, :].reshape(bsz, tokens, channels)
        x = shortcut + self.drop(x_grid)
        return x + self.mlp(self.norm2(x))


class AudioMAELocalWindowScMAE(nn.Module):
    """AudioMAE-inspired 2D masked autoencoder over a gene grid."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 4,
        hidden_size: int = 128,
        depth: int = 3,
        num_heads: int = 4,
        decoder_size: int = 128,
        decoder_depth: int = 4,
        decoder_heads: int = 4,
        window_size: int = 4,
        dropout: float = 0.1,
    ) -> None:
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.grid_h, self.grid_w, self.patch_h, self.patch_w = compute_gene_grid(num_genes, patch_size)
        self.num_patches = self.patch_h * self.patch_w
        self.patch_embed = nn.Conv2d(1, hidden_size, kernel_size=patch_size, stride=patch_size)
        self.cls_token = nn.Parameter(torch.randn(1, 1, hidden_size) * 0.02)
        self.pos_embed = nn.Parameter(torch.randn(1, self.num_patches + 1, hidden_size) * 0.02)
        self.blocks = nn.ModuleList([TransformerBlock(hidden_size, num_heads, dropout) for _ in range(depth)])
        self.norm = nn.LayerNorm(hidden_size)
        self.decoder_embed = nn.Linear(hidden_size, decoder_size)
        self.mask_token = nn.Parameter(torch.randn(1, 1, decoder_size) * 0.02)
        self.decoder_pos_embed = nn.Parameter(torch.randn(1, self.num_patches + 1, decoder_size) * 0.02)
        self.decoder_blocks = nn.ModuleList()
        for idx in range(decoder_depth):
            shift = (0, 0) if idx % 2 == 0 else (max(1, window_size // 2), 0)
            self.decoder_blocks.append(
                LocalWindowBlock(
                    decoder_size,
                    decoder_heads,
                    (int(window_size), int(window_size)),
                    shift,
                    dropout,
                )
            )
        self.decoder_norm = nn.LayerNorm(decoder_size)
        self.decoder_pred = nn.Linear(decoder_size, patch_size * patch_size)

    def vector_to_grid(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        total = self.grid_h * self.grid_w
        padded = F.pad(x, (0, total - self.num_genes)) if total > self.num_genes else x
        return padded.view(x.shape[0], 1, self.grid_h, self.grid_w)

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        grid = self.vector_to_grid(x)
        patches = grid.unfold(2, self.patch_size, self.patch_size).unfold(3, self.patch_size, self.patch_size)
        patches = patches.permute(0, 2, 3, 1, 4, 5).reshape(x.shape[0], self.num_patches, -1)
        return patches

    def valid_patch_weights(self, device: torch.device) -> torch.Tensor:
        valid = torch.zeros(1, self.grid_h * self.grid_w, device=device)
        valid[:, : self.num_genes] = 1.0
        valid = valid.view(1, 1, self.grid_h, self.grid_w)
        patches = valid.unfold(2, self.patch_size, self.patch_size).unfold(3, self.patch_size, self.patch_size)
        return patches.permute(0, 2, 3, 1, 4, 5).reshape(1, self.num_patches, -1)

    def random_masking(self, tokens: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        batch_size, length, dim = tokens.shape
        len_keep = max(1, min(length - 1, int(length * (1.0 - float(mask_ratio)))))
        noise = torch.rand(batch_size, length, device=tokens.device)
        ids_shuffle = torch.argsort(noise, dim=1)
        ids_restore = torch.argsort(ids_shuffle, dim=1)
        ids_keep = ids_shuffle[:, :len_keep]
        x_masked = torch.gather(tokens, dim=1, index=ids_keep.unsqueeze(-1).expand(-1, -1, dim))
        mask = torch.ones(batch_size, length, device=tokens.device)
        mask[:, :len_keep] = 0
        mask = torch.gather(mask, dim=1, index=ids_restore)
        return x_masked, mask, ids_restore

    def random_masking_2d(
        self,
        tokens: torch.Tensor,
        mask_t_prob: float,
        mask_f_prob: float,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        batch_size, length, dim = tokens.shape
        if length != self.num_patches:
            raise ValueError("token length does not match patch grid")
        keep_t = max(1, min(self.patch_h, int(self.patch_h * (1.0 - float(mask_t_prob)))))
        keep_f = max(1, min(self.patch_w, int(self.patch_w * (1.0 - float(mask_f_prob)))))
        noise_t = torch.rand(batch_size, self.patch_h, device=tokens.device)
        ids_restore_t = torch.argsort(torch.argsort(noise_t, dim=1), dim=1)
        noise_f = torch.rand(batch_size, self.patch_w, device=tokens.device)
        ids_restore_f = torch.argsort(torch.argsort(noise_f, dim=1), dim=1)
        mask_t = torch.ones(batch_size, self.patch_h, device=tokens.device)
        mask_t[:, :keep_t] = 0
        mask_t = torch.gather(mask_t, dim=1, index=ids_restore_t).unsqueeze(2).expand(-1, -1, self.patch_w)
        mask_f = torch.ones(batch_size, self.patch_w, device=tokens.device)
        mask_f[:, :keep_f] = 0
        mask_f = torch.gather(mask_f, dim=1, index=ids_restore_f).unsqueeze(1).expand(-1, self.patch_h, -1)
        mask_grid = 1.0 - (1.0 - mask_t) * (1.0 - mask_f)
        positions = torch.arange(length, device=tokens.device, dtype=tokens.dtype).view(1, self.patch_h, self.patch_w)
        ordered = positions + mask_grid * float(length + 1)
        ids_shuffle = torch.argsort(ordered.flatten(start_dim=1), dim=1)
        ids_restore = torch.argsort(ids_shuffle, dim=1)
        ids_keep = ids_shuffle[:, : keep_t * keep_f]
        x_masked = torch.gather(tokens, dim=1, index=ids_keep.unsqueeze(-1).expand(-1, -1, dim))
        return x_masked, mask_grid.flatten(start_dim=1), ids_restore

    def forward_encoder(
        self,
        x: torch.Tensor,
        mask_ratio: float,
        mask_2d: bool,
        mask_t_prob: float,
        mask_f_prob: float,
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        grid = self.vector_to_grid(x)
        tokens = self.patch_embed(grid).flatten(2).transpose(1, 2)
        tokens = tokens + self.pos_embed[:, 1:, :].to(dtype=tokens.dtype, device=tokens.device)
        if mask_2d:
            tokens, mask, ids_restore = self.random_masking_2d(tokens, mask_t_prob, mask_f_prob)
        else:
            tokens, mask, ids_restore = self.random_masking(tokens, mask_ratio)
        cls = self.cls_token.expand(tokens.shape[0], -1, -1).to(dtype=tokens.dtype, device=tokens.device)
        cls = cls + self.pos_embed[:, :1, :].to(dtype=tokens.dtype, device=tokens.device)
        h = torch.cat([cls, tokens], dim=1)
        for block in self.blocks:
            h = block(h)
        return self.norm(h), mask, ids_restore

    def forward_decoder(self, encoded: torch.Tensor, ids_restore: torch.Tensor) -> torch.Tensor:
        x = self.decoder_embed(encoded)
        mask_tokens = self.mask_token.expand(x.shape[0], ids_restore.shape[1] + 1 - x.shape[1], -1)
        x_without_cls = torch.cat([x[:, 1:, :], mask_tokens], dim=1)
        x_without_cls = torch.gather(
            x_without_cls,
            dim=1,
            index=ids_restore.unsqueeze(-1).expand(-1, -1, x.shape[-1]),
        )
        x = torch.cat([x[:, :1, :], x_without_cls], dim=1)
        x = x + self.decoder_pos_embed.to(dtype=x.dtype, device=x.device)
        x_patch = x[:, 1:, :]
        for block in self.decoder_blocks:
            x_patch = block(x_patch, (self.patch_h, self.patch_w))
        x_patch = self.decoder_norm(x_patch)
        return self.decoder_pred(x_patch)

    def forward(
        self,
        x: torch.Tensor,
        mask_ratio: float,
        mask_2d: bool,
        mask_t_prob: float,
        mask_f_prob: float,
    ) -> dict[str, torch.Tensor]:
        encoded, mask, ids_restore = self.forward_encoder(x, mask_ratio, mask_2d, mask_t_prob, mask_f_prob)
        pred = self.forward_decoder(encoded, ids_restore)
        return {"embedding": encoded[:, 0], "prediction": pred, "mask": mask}

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        grid = self.vector_to_grid(x)
        tokens = self.patch_embed(grid).flatten(2).transpose(1, 2)
        tokens = tokens + self.pos_embed[:, 1:, :].to(dtype=tokens.dtype, device=tokens.device)
        cls = self.cls_token.expand(tokens.shape[0], -1, -1).to(dtype=tokens.dtype, device=tokens.device)
        cls = cls + self.pos_embed[:, :1, :].to(dtype=tokens.dtype, device=tokens.device)
        h = torch.cat([cls, tokens], dim=1)
        for block in self.blocks:
            h = block(h)
        return self.norm(h)[:, 0]
