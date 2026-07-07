from __future__ import annotations

import copy
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


class ModulePatchEncoder(nn.Module):
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
        self.encoder = nn.TransformerEncoder(layer, num_layers=n_layers)
        self.norm = nn.LayerNorm(self.hidden_size)
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

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        return self._pad(x).view(x.shape[0], self.n_patches, self.patch_size)

    def forward(self, x: torch.Tensor, ids: torch.Tensor | None = None) -> torch.Tensor:
        patches = self.patchify(x)
        tokens = self.patch_embed(patches) + self.pos_embed.to(x.device)
        if ids is not None:
            tokens = torch.gather(tokens, dim=1, index=ids.unsqueeze(-1).expand(-1, -1, tokens.shape[-1]))
        return self.norm(self.encoder(tokens))


class TargetConditionedPredictor(nn.Module):
    def __init__(
        self,
        n_patches: int,
        hidden_size: int = 128,
        predictor_size: int = 96,
        depth: int = 3,
        num_heads: int = 4,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.n_patches = int(n_patches)
        self.hidden_size = int(hidden_size)
        self.context_proj = nn.Linear(hidden_size, predictor_size)
        self.mask_token = nn.Parameter(torch.zeros(1, 1, predictor_size))
        pe = sinusoidal_positional_encoding(self.n_patches, predictor_size)
        self.register_buffer("pos_embed", pe.unsqueeze(0), persistent=False)
        layer = nn.TransformerEncoderLayer(
            d_model=predictor_size,
            nhead=num_heads,
            dim_feedforward=predictor_size * 3,
            dropout=dropout,
            activation="gelu",
            batch_first=True,
            norm_first=True,
        )
        self.blocks = nn.TransformerEncoder(layer, num_layers=depth)
        self.norm = nn.LayerNorm(predictor_size)
        self.proj = nn.Linear(predictor_size, hidden_size)
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

    def forward(self, context_tokens: torch.Tensor, context_ids: torch.Tensor, target_ids: torch.Tensor) -> torch.Tensor:
        ctx = self.context_proj(context_tokens)
        pe = self.pos_embed.to(context_tokens.device)
        ctx = ctx + torch.gather(pe.expand(ctx.shape[0], -1, -1), 1, context_ids.unsqueeze(-1).expand(-1, -1, ctx.shape[-1]))
        tgt_pos = torch.gather(pe.expand(ctx.shape[0], -1, -1), 1, target_ids.unsqueeze(-1).expand(-1, -1, ctx.shape[-1]))
        tgt = self.mask_token.expand(ctx.shape[0], target_ids.shape[1], -1) + tgt_pos
        out = self.blocks(torch.cat([ctx, tgt], dim=1))
        out = self.norm(out[:, ctx.shape[1] :, :])
        return self.proj(out)


class IJEPAScMAE(nn.Module):
    """Context-target JEPA with scMAE reconstruction and mask prediction."""

    def __init__(
        self,
        num_genes: int,
        patch_size: int = 20,
        hidden_size: int = 128,
        encoder_layers: int = 2,
        predictor_size: int = 96,
        predictor_depth: int = 3,
        num_heads: int = 4,
        dropout: float = 0.05,
    ):
        super().__init__()
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.n_patches = (self.num_genes + self.patch_size - 1) // self.patch_size
        self.online_encoder = ModulePatchEncoder(num_genes, patch_size, hidden_size, encoder_layers, num_heads, dropout)
        self.target_encoder = copy.deepcopy(self.online_encoder)
        for param in self.target_encoder.parameters():
            param.requires_grad_(False)
        self.predictor = TargetConditionedPredictor(self.n_patches, hidden_size, predictor_size, predictor_depth, num_heads, dropout)
        self.pool_norm = nn.LayerNorm(hidden_size)
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.expr_decoder = nn.Sequential(
            nn.Linear(hidden_size + self.num_genes, 256),
            nn.Mish(),
            nn.Dropout(dropout),
            nn.Linear(256, self.num_genes),
        )

    def sample_masks(
        self,
        batch_size: int,
        device: torch.device,
        target_fraction: float = 0.18,
        n_targets: int = 2,
        context_fraction: float = 0.72,
    ) -> dict[str, torch.Tensor]:
        target_len = max(1, min(self.n_patches - 1, int(round(self.n_patches * float(target_fraction)))))
        n_targets = max(1, int(n_targets))
        target_mask = torch.zeros(batch_size, self.n_patches, device=device, dtype=torch.bool)
        target_ids_list = []
        positions = torch.arange(self.n_patches, device=device)
        for _ in range(n_targets):
            starts = torch.randint(0, self.n_patches, (batch_size,), device=device)
            block = (positions.unsqueeze(0) - starts.unsqueeze(1)).remainder(self.n_patches) < target_len
            ids = torch.topk(block.float(), k=target_len, dim=1).indices
            target_ids_list.append(ids)
            target_mask |= block
        all_target_ids = torch.cat(target_ids_list, dim=1)
        context_len = max(1, min(self.n_patches - target_len, int(round(self.n_patches * float(context_fraction)))))
        noise = torch.rand(batch_size, self.n_patches, device=device)
        noise = noise.masked_fill(target_mask, 2.0)
        context_ids = torch.argsort(noise, dim=1)[:, :context_len]
        gene_mask = target_mask.repeat_interleave(self.patch_size, dim=1)[:, : self.num_genes].float()
        return {
            "context_ids": context_ids,
            "target_ids": all_target_ids,
            "patch_target_mask": target_mask.float(),
            "gene_mask": gene_mask,
        }

    def forward(self, x: torch.Tensor, masks: dict[str, torch.Tensor]) -> dict[str, torch.Tensor]:
        context_tokens = self.online_encoder(x, masks["context_ids"])
        latent = self.pool_norm(context_tokens.mean(dim=1))
        pred_tokens = self.predictor(context_tokens, masks["context_ids"], masks["target_ids"])
        with torch.no_grad():
            target_tokens_all = self.target_encoder(x, None)
            target_tokens = torch.gather(target_tokens_all, 1, masks["target_ids"].unsqueeze(-1).expand(-1, -1, target_tokens_all.shape[-1]))
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.expr_decoder(torch.cat([latent, mask_logits], dim=1))
        return {
            "latent": latent,
            "prediction": pred_tokens,
            "target": target_tokens.detach(),
            "mask_logits": mask_logits,
            "reconstruction": reconstruction,
            "gene_mask": masks["gene_mask"],
            "patch_target_mask": masks["patch_target_mask"],
        }

    @torch.no_grad()
    def update_target_encoder(self, momentum: float) -> None:
        for online, target in zip(self.online_encoder.parameters(), self.target_encoder.parameters()):
            target.data.mul_(momentum).add_(online.data, alpha=1.0 - momentum)

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        tokens = self.online_encoder(x, None)
        return self.pool_norm(tokens.mean(dim=1))
