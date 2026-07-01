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


class VideoMAEGeneModuleScMAE(nn.Module):
    """Asymmetric VideoMAE-style autoencoder over data-ordered gene modules."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        tube_frames: int = 4,
        hidden_size: int = 128,
        decoder_size: int = 96,
        n_layers: int = 2,
        decoder_layers: int = 2,
        num_heads: int = 4,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.tube_frames = int(max(1, tube_frames))
        self.hidden_size = int(hidden_size)
        self.decoder_size = int(decoder_size)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.padded_genes = self.n_patches * self.patch_size
        self.tube_frames = min(self.tube_frames, self.n_patches)
        self.patches_per_frame = (self.n_patches + self.tube_frames - 1) // self.tube_frames
        self.tube_slots = self.tube_frames * self.patches_per_frame

        self.patch_embed = nn.Sequential(
            nn.Linear(self.patch_size, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
        )
        pe = sinusoidal_positional_encoding(self.tube_slots, hidden_size)
        self.register_buffer("encoder_pos", pe.unsqueeze(0), persistent=False)
        enc_layer = nn.TransformerEncoderLayer(
            d_model=hidden_size,
            nhead=num_heads,
            dim_feedforward=hidden_size * 3,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(enc_layer, num_layers=n_layers)
        self.encoder_norm = nn.LayerNorm(hidden_size)
        self.to_decoder = nn.Linear(hidden_size, decoder_size, bias=False)
        self.mask_token = nn.Parameter(torch.zeros(1, 1, decoder_size))
        dec_pe = sinusoidal_positional_encoding(self.tube_slots, decoder_size)
        self.register_buffer("decoder_pos", dec_pe.unsqueeze(0), persistent=False)
        dec_heads = max(1, min(num_heads, decoder_size // 16))
        dec_layer = nn.TransformerEncoderLayer(
            d_model=decoder_size,
            nhead=dec_heads,
            dim_feedforward=decoder_size * 3,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.decoder = nn.TransformerEncoder(dec_layer, num_layers=decoder_layers)
        self.decoder_norm = nn.LayerNorm(decoder_size)
        self.expr_head = nn.Linear(decoder_size, patch_size)
        self.mask_head = nn.Linear(decoder_size, patch_size)
        self.module_head = nn.Linear(decoder_size, 1)
        self.latent_norm = nn.LayerNorm(hidden_size)
        self.apply(self._init_weights)

    def _init_weights(self, module: nn.Module) -> None:
        if isinstance(module, nn.Linear):
            nn.init.xavier_uniform_(module.weight)
            if module.bias is not None:
                nn.init.zeros_(module.bias)
        elif isinstance(module, nn.LayerNorm):
            nn.init.ones_(module.weight)
            nn.init.zeros_(module.bias)
        nn.init.normal_(self.mask_token, std=0.02)

    def _pad_genes(self, x: torch.Tensor) -> torch.Tensor:
        if self.padded_genes == self.num_genes:
            return x
        return F.pad(x, (0, self.padded_genes - self.num_genes))

    def _pad_patches(self, patches: torch.Tensor) -> torch.Tensor:
        if self.tube_slots == self.n_patches:
            return patches
        pad = patches.new_zeros(patches.shape[0], self.tube_slots - self.n_patches, patches.shape[2])
        return torch.cat([patches, pad], dim=1)

    def _unpad_genes(self, patches: torch.Tensor) -> torch.Tensor:
        return patches[:, : self.n_patches, :].reshape(patches.shape[0], self.padded_genes)[:, : self.num_genes]

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        return self._pad_genes(x).view(x.shape[0], self.n_patches, self.patch_size)

    def tube_mask(self, batch_size: int, mask_ratio: float, device: torch.device) -> torch.Tensor:
        masks = []
        for _ in range(batch_size):
            sample_mask = torch.ones(self.tube_slots, dtype=torch.bool, device=device)
            for frame in range(self.tube_frames):
                start = frame * self.patches_per_frame
                end = min(start + self.patches_per_frame, self.n_patches)
                n_valid = max(0, end - start)
                if n_valid == 0:
                    continue
                n_mask = max(1, min(n_valid - 1, int(round(n_valid * float(mask_ratio))))) if n_valid > 1 else 1
                frame_mask = torch.zeros(n_valid, dtype=torch.bool, device=device)
                idx = torch.argsort(torch.rand(n_valid, device=device))[:n_mask]
                frame_mask[idx] = True
                sample_mask[start:end] = frame_mask
            masks.append(sample_mask)
        mask = torch.stack(masks, dim=0)
        return mask

    def gene_mask_from_patch_mask(self, patch_mask: torch.Tensor) -> torch.Tensor:
        gene_mask = patch_mask[:, : self.n_patches].float().unsqueeze(-1).expand(-1, -1, self.patch_size)
        return gene_mask.reshape(patch_mask.shape[0], self.padded_genes)[:, : self.num_genes]

    def forward(self, x: torch.Tensor, patch_mask: torch.Tensor) -> dict[str, torch.Tensor]:
        patches = self.patchify(x)
        patch_targets_mean = patches.mean(dim=-1)
        patches_full = self._pad_patches(patches)
        tokens = self.patch_embed(patches_full) + self.encoder_pos[:, : self.tube_slots].to(x.device)
        visible = (~patch_mask).bool()
        x_vis = tokens[visible].reshape(x.shape[0], -1, self.hidden_size)
        x_vis = self.encoder_norm(self.encoder(x_vis))
        latent = self.latent_norm(x_vis.mean(dim=1))
        dec_vis = self.to_decoder(x_vis)
        full = self.mask_token.expand(x.shape[0], self.tube_slots, -1).clone()
        full[visible] = dec_vis.reshape(-1, self.decoder_size)
        full = full + self.decoder_pos[:, : self.tube_slots].to(x.device)
        decoded = self.decoder_norm(self.decoder(full))
        rec_patches = self.expr_head(decoded)
        reconstruction = self._unpad_genes(rec_patches)
        mask_logits = self._unpad_genes(self.mask_head(decoded))
        module_pred = self.module_head(decoded[:, : self.n_patches]).squeeze(-1)
        return {
            "latent": latent,
            "tokens": decoded[:, : self.n_patches],
            "reconstruction": reconstruction,
            "mask_logits": mask_logits,
            "module_pred": module_pred,
            "module_target": patch_targets_mean,
            "patch_mask": patch_mask[:, : self.n_patches].float(),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        patch_mask = torch.zeros(x.shape[0], self.tube_slots, dtype=torch.bool, device=x.device)
        return self.forward(x, patch_mask)["latent"]
