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


class ShiftedLocalDecoderBlock(nn.Module):
    """AudioMAE-style local window decoder block for ordered gene-module patches."""

    def __init__(self, hidden_size: int, num_heads: int = 4, window_size: int = 5, shift: int = 0, dropout: float = 0.05):
        super().__init__()
        self.window_size = max(1, int(window_size))
        self.shift = int(shift)
        self.norm1 = nn.LayerNorm(hidden_size)
        self.attn = nn.MultiheadAttention(hidden_size, num_heads, dropout=dropout, batch_first=True)
        self.norm2 = nn.LayerNorm(hidden_size)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * 2),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 2, hidden_size),
        )

    def _window_attention(self, x: torch.Tensor) -> torch.Tensor:
        bsz, length, dim = x.shape
        pad = (self.window_size - (length % self.window_size)) % self.window_size
        if pad:
            x = torch.nn.functional.pad(x, (0, 0, 0, pad))
        padded_len = x.shape[1]
        windows = x.view(bsz, padded_len // self.window_size, self.window_size, dim).reshape(-1, self.window_size, dim)
        attended, _ = self.attn(windows, windows, windows, need_weights=False)
        attended = attended.reshape(bsz, padded_len // self.window_size, self.window_size, dim).reshape(bsz, padded_len, dim)
        return attended[:, :length, :]

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        z = self.norm1(x)
        if self.shift:
            z_shift = torch.roll(z, shifts=-self.shift, dims=1)
            z_attn = torch.roll(self._window_attention(z_shift), shifts=self.shift, dims=1)
        else:
            z_attn = self._window_attention(z)
        x = x + z_attn
        return x + self.mlp(self.norm2(x))


class AudioMAEModulePatchScMAE(nn.Module):
    """Independent scMAE variant adapted from AudioMAE to coexpression-ordered gene patches."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        encoder_layers: int = 2,
        decoder_layers: int = 3,
        num_heads: int = 4,
        window_size: int = 5,
        dropout: float = 0.05,
        module_target_dim: int = 3,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.hidden_size = int(hidden_size)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.padded_genes = self.n_patches * self.patch_size
        self.module_target_dim = int(module_target_dim)

        self.patch_embed = nn.Sequential(
            nn.Linear(self.patch_size, self.hidden_size),
            nn.LayerNorm(self.hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        pe = sinusoidal_positional_encoding(self.n_patches, self.hidden_size)
        self.register_buffer("pos_embed", pe.unsqueeze(0), persistent=False)

        encoder_layer = nn.TransformerEncoderLayer(
            d_model=self.hidden_size,
            nhead=num_heads,
            dim_feedforward=self.hidden_size * 3,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(encoder_layer, num_layers=encoder_layers)
        self.encoder_norm = nn.LayerNorm(self.hidden_size)

        self.decoder_embed = nn.Linear(self.hidden_size, self.hidden_size)
        self.mask_token = nn.Parameter(torch.zeros(1, 1, self.hidden_size))
        decoder_blocks = []
        for layer_id in range(decoder_layers):
            shift = 0 if layer_id % 2 == 0 else max(1, window_size // 2)
            decoder_blocks.append(ShiftedLocalDecoderBlock(self.hidden_size, num_heads, window_size, shift, dropout))
        self.decoder_blocks = nn.ModuleList(decoder_blocks)
        self.decoder_norm = nn.LayerNorm(self.hidden_size)
        self.patch_decoder = nn.Linear(self.hidden_size, self.patch_size)

        self.pool_norm = nn.LayerNorm(self.hidden_size)
        self.mask_predictor = nn.Linear(self.hidden_size, self.num_genes)
        self.module_head = nn.Linear(self.hidden_size, self.n_patches * self.module_target_dim)

        self.apply(self._init_weights)
        nn.init.normal_(self.mask_token, std=0.02)

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

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        return self._pad(x).view(x.shape[0], self.n_patches, self.patch_size)

    def unpatchify(self, patches: torch.Tensor) -> torch.Tensor:
        return patches.reshape(patches.shape[0], self.padded_genes)[:, : self.num_genes]

    def make_audio_mask(self, batch_size: int, device: torch.device, mask_ratio: float, structured_prob: float = 0.35) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        len_keep = max(1, min(self.n_patches - 1, int(round(self.n_patches * (1.0 - float(mask_ratio))))))
        use_structured = torch.rand(batch_size, device=device) < float(structured_prob)
        noise = torch.rand(batch_size, self.n_patches, device=device)
        if bool(use_structured.any()):
            starts = torch.randint(0, self.n_patches, (batch_size,), device=device)
            block_len = max(1, self.n_patches - len_keep)
            positions = torch.arange(self.n_patches, device=device).unsqueeze(0)
            distance = (positions - starts.unsqueeze(1)).remainder(self.n_patches)
            structured_noise = torch.where(distance < block_len, noise + 1.0, noise)
            noise = torch.where(use_structured.unsqueeze(1), structured_noise, noise)
        ids_shuffle = torch.argsort(noise, dim=1)
        ids_restore = torch.argsort(ids_shuffle, dim=1)
        ids_keep = ids_shuffle[:, :len_keep]
        patch_mask = torch.ones(batch_size, self.n_patches, device=device)
        patch_mask[:, :len_keep] = 0.0
        patch_mask = torch.gather(patch_mask, dim=1, index=ids_restore)
        return patch_mask, ids_keep, ids_restore

    def encode_visible(self, patches: torch.Tensor, ids_keep: torch.Tensor) -> torch.Tensor:
        tokens = self.patch_embed(patches) + self.pos_embed.to(patches.device)
        visible = torch.gather(tokens, dim=1, index=ids_keep.unsqueeze(-1).expand(-1, -1, tokens.shape[-1]))
        return self.encoder_norm(self.encoder(visible))

    def decode_full(self, visible_encoded: torch.Tensor, ids_restore: torch.Tensor) -> torch.Tensor:
        dec_visible = self.decoder_embed(visible_encoded)
        mask_tokens = self.mask_token.expand(dec_visible.shape[0], ids_restore.shape[1] - dec_visible.shape[1], -1)
        x = torch.cat([dec_visible, mask_tokens], dim=1)
        x = torch.gather(x, dim=1, index=ids_restore.unsqueeze(-1).expand(-1, -1, x.shape[-1]))
        x = x + self.pos_embed.to(x.device)
        for block in self.decoder_blocks:
            x = block(x)
        return self.decoder_norm(x)

    def forward(self, x: torch.Tensor, mask_ratio: float = 0.65, structured_prob: float = 0.35) -> dict[str, torch.Tensor]:
        patches = self.patchify(x)
        patch_mask, ids_keep, ids_restore = self.make_audio_mask(x.shape[0], x.device, mask_ratio, structured_prob)
        visible = self.encode_visible(patches, ids_keep)
        latent = self.pool_norm(visible.mean(dim=1))
        decoded = self.decode_full(visible, ids_restore)
        patch_recon = self.patch_decoder(decoded)
        gene_mask = patch_mask.repeat_interleave(self.patch_size, dim=1)[:, : self.num_genes]
        return {
            "latent": latent,
            "decoded_tokens": decoded,
            "patch_mask": patch_mask,
            "gene_mask": gene_mask,
            "mask_logits": self.mask_predictor(latent),
            "reconstruction": self.unpatchify(patch_recon),
            "module_reconstruction": self.module_head(latent).view(x.shape[0], self.n_patches, self.module_target_dim),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        patches = self.patchify(x)
        tokens = self.patch_embed(patches) + self.pos_embed.to(x.device)
        encoded = self.encoder_norm(self.encoder(tokens))
        return self.pool_norm(encoded.mean(dim=1))
