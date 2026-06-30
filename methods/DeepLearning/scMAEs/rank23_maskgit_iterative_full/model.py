from __future__ import annotations

import torch
import torch.nn as nn


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
        attn, _ = self.attn(q, q, q, need_weights=False)
        x = x + self.dropout(attn)
        return x + self.mlp(self.norm2(x))


class MaskGITExpressionTransformer(nn.Module):
    """Bidirectional MaskGIT token transformer for quantized expression patches."""

    def __init__(
        self,
        num_patches: int,
        patch_size: int,
        vocab_size: int = 64,
        hidden_size: int = 128,
        depth: int = 4,
        num_heads: int = 4,
        dropout: float = 0.1,
        mlp_ratio: int = 4,
    ) -> None:
        super().__init__()
        if num_patches <= 0 or patch_size <= 0 or vocab_size <= 1:
            raise ValueError("num_patches, patch_size, and vocab_size must be valid positive values")
        if hidden_size % num_heads != 0:
            raise ValueError("hidden_size must be divisible by num_heads")
        self.num_patches = int(num_patches)
        self.patch_size = int(patch_size)
        self.vocab_size = int(vocab_size)
        self.mask_token_id = int(vocab_size)
        self.token_embed = nn.Embedding(vocab_size + 1, hidden_size)
        self.position = nn.Parameter(torch.randn(1, num_patches, hidden_size) * 0.02)
        self.blocks = nn.ModuleList([TransformerBlock(hidden_size, num_heads, dropout, mlp_ratio) for _ in range(depth)])
        self.norm = nn.LayerNorm(hidden_size)
        self.mlm_dense = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.LayerNorm(hidden_size),
        )
        self.output_bias = nn.Parameter(torch.zeros(vocab_size))
        self.patch_decoder = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, patch_size),
        )

    def encode(self, token_ids: torch.Tensor) -> torch.Tensor:
        if token_ids.ndim != 2 or token_ids.shape[1] != self.num_patches:
            raise ValueError(f"token_ids must be [batch, {self.num_patches}], got {tuple(token_ids.shape)}")
        if token_ids.min().item() < 0 or token_ids.max().item() > self.mask_token_id:
            raise ValueError("token_ids contain values outside [0, mask_token_id]")
        h = self.token_embed(token_ids.long()) + self.position.to(dtype=self.token_embed.weight.dtype)
        for block in self.blocks:
            h = block(h)
        return self.norm(h)

    def logits_from_hidden(self, hidden: torch.Tensor) -> torch.Tensor:
        mlm = self.mlm_dense(hidden)
        logits = torch.matmul(mlm, self.token_embed.weight[: self.vocab_size].T)
        return logits + self.output_bias

    def forward(self, token_ids: torch.Tensor) -> dict[str, torch.Tensor]:
        hidden = self.encode(token_ids)
        logits = self.logits_from_hidden(hidden)
        patches = self.patch_decoder(hidden)
        return {
            "embedding": hidden.mean(dim=1),
            "hidden": hidden,
            "token_logits": logits,
            "patch_reconstruction": patches,
        }

    @torch.no_grad()
    def iterative_decode(
        self,
        token_ids: torch.Tensor,
        initial_mask: torch.Tensor,
        num_steps: int,
        temperature: float = 1.0,
    ) -> tuple[torch.Tensor, torch.Tensor]:
        if initial_mask.shape != token_ids.shape:
            raise ValueError("initial_mask must match token_ids")
        current = token_ids.clone()
        current[initial_mask.bool()] = self.mask_token_id
        original_unknown = initial_mask.bool()
        confidences = torch.zeros_like(current, dtype=torch.float32)
        for step in range(max(1, int(num_steps))):
            out = self(current)
            probs = torch.softmax(out["token_logits"], dim=-1)
            chosen_prob, chosen = probs.max(dim=-1)
            unknown = current.eq(self.mask_token_id)
            current = torch.where(unknown, chosen, current)
            confidences = torch.where(unknown, chosen_prob, confidences)
            if step == int(num_steps) - 1:
                break
            remaining = original_unknown & current.ne(self.mask_token_id)
            ratio = float(step + 1) / float(max(1, int(num_steps)))
            mask_ratio = max(1e-6, torch.cos(torch.tensor(ratio * 1.57079632679, device=current.device)).item())
            for row in range(current.shape[0]):
                candidates = torch.nonzero(remaining[row], as_tuple=False).flatten()
                if candidates.numel() <= 1:
                    continue
                k = max(1, min(candidates.numel() - 1, int(candidates.numel() * mask_ratio)))
                score = torch.log(confidences[row, candidates].clamp_min(1e-8))
                if temperature > 0:
                    score = score + (float(temperature) * (1.0 - ratio)) * (-torch.empty_like(score).exponential_().log())
                remask = candidates[torch.argsort(score)[:k]]
                current[row, remask] = self.mask_token_id
        return current.clamp_max(self.vocab_size - 1), confidences
