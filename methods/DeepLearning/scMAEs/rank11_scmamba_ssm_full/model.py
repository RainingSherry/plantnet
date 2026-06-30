from __future__ import annotations

import math

import torch
import torch.nn as nn
import torch.nn.functional as F


class PatchTokenizer(nn.Module):
    """scMamba-style expression patch tokenizer with explicit mask token."""

    def __init__(self, num_genes: int, patch_size: int, hidden_size: int, dropout: float) -> None:
        super().__init__()
        if num_genes <= 0 or patch_size <= 0 or hidden_size <= 0:
            raise ValueError("num_genes, patch_size, and hidden_size must be positive")
        self.num_genes = int(num_genes)
        self.patch_size = int(patch_size)
        self.num_patches = int(math.ceil(self.num_genes / self.patch_size))
        self.pad_size = self.num_patches * self.patch_size - self.num_genes
        self.projection = nn.Linear(self.patch_size, hidden_size)
        self.mask_token = nn.Parameter(torch.zeros(1, 1, hidden_size))
        self.position = nn.Parameter(torch.randn(1, self.num_patches, hidden_size) * 0.02)
        self.dropout = nn.Dropout(dropout)

    def patchify(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 2 or x.shape[1] != self.num_genes:
            raise ValueError(f"x must be [batch, {self.num_genes}], got {tuple(x.shape)}")
        padded = F.pad(x, (0, self.pad_size)) if self.pad_size else x
        return padded.view(x.shape[0], self.num_patches, self.patch_size)

    def forward(self, x: torch.Tensor, patch_mask: torch.Tensor | None = None) -> torch.Tensor:
        patches = self.patchify(x)
        tokens = self.projection(patches)
        if patch_mask is not None:
            if patch_mask.ndim != 2 or patch_mask.shape != tokens.shape[:2]:
                raise ValueError(f"patch_mask must be [batch, {self.num_patches}], got {tuple(patch_mask.shape)}")
            mask = patch_mask.to(dtype=tokens.dtype).unsqueeze(-1)
            tokens = tokens * (1.0 - mask) + self.mask_token.to(dtype=tokens.dtype) * mask
        return self.dropout(tokens + self.position.to(dtype=tokens.dtype))


class SelectiveSSMLayer(nn.Module):
    """Dependency-free selective SSM scan over gene patches.

    This mirrors the essential Mamba/SSM ingredients instead of replacing them
    with a feed-forward gate: learned negative state dynamics A, input B,
    output C, skip D, positive delta, causal depthwise convolution, and a
    recurrent selective scan along the patch sequence.
    """

    def __init__(
        self,
        hidden_size: int,
        state_size: int,
        conv_kernel: int,
        dropout: float,
        expand: int = 2,
    ) -> None:
        super().__init__()
        if hidden_size <= 0 or state_size <= 0 or conv_kernel <= 0:
            raise ValueError("hidden_size, state_size, and conv_kernel must be positive")
        inner = int(hidden_size * expand)
        self.hidden_size = int(hidden_size)
        self.state_size = int(state_size)
        self.inner_size = inner
        self.in_proj = nn.Linear(hidden_size, inner * 2)
        self.depthwise_conv = nn.Conv1d(
            inner,
            inner,
            kernel_size=conv_kernel,
            padding=conv_kernel - 1,
            groups=inner,
        )
        self.delta_proj = nn.Linear(inner, inner)
        self.bc_proj = nn.Linear(inner, inner * state_size * 2)
        self.A_log = nn.Parameter(torch.empty(inner, state_size))
        self.D = nn.Parameter(torch.ones(inner))
        self.delta_bias = nn.Parameter(torch.zeros(inner))
        self.out_proj = nn.Linear(inner, hidden_size)
        self.dropout = nn.Dropout(dropout)
        nn.init.uniform_(self.A_log, -3.0, -1.0)

    def selective_scan(self, u: torch.Tensor) -> torch.Tensor:
        if u.ndim != 3 or u.shape[2] != self.inner_size:
            raise ValueError(f"u must be [batch, patches, {self.inner_size}], got {tuple(u.shape)}")
        batch, seq_len, inner = u.shape
        dtype = u.dtype
        A = -torch.exp(self.A_log.float()).to(device=u.device, dtype=dtype)
        D = self.D.to(device=u.device, dtype=dtype)
        delta = F.softplus(self.delta_proj(u) + self.delta_bias.to(dtype=dtype))
        bc = self.bc_proj(u).view(batch, seq_len, inner, 2, self.state_size)
        B_t = bc[:, :, :, 0, :]
        C_t = bc[:, :, :, 1, :]
        state = u.new_zeros(batch, inner, self.state_size)
        outputs = []
        for step in range(seq_len):
            dt = delta[:, step].unsqueeze(-1)
            dA = torch.exp(dt * A.unsqueeze(0))
            dB_u = dt * B_t[:, step] * u[:, step].unsqueeze(-1)
            state = dA * state + dB_u
            y = (state * C_t[:, step]).sum(dim=-1) + D.unsqueeze(0) * u[:, step]
            outputs.append(y)
        return torch.stack(outputs, dim=1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if x.ndim != 3 or x.shape[2] != self.hidden_size:
            raise ValueError(f"x must be [batch, patches, {self.hidden_size}], got {tuple(x.shape)}")
        projected, gate = self.in_proj(x).chunk(2, dim=-1)
        conv_in = projected.transpose(1, 2)
        conv_out = self.depthwise_conv(conv_in)[..., : x.shape[1]].transpose(1, 2)
        u = F.silu(conv_out)
        scanned = self.selective_scan(u)
        mixed = scanned * F.silu(gate)
        return self.dropout(self.out_proj(mixed))


class SSMBlock(nn.Module):
    def __init__(self, hidden_size: int, state_size: int, conv_kernel: int, dropout: float, mlp_ratio: int) -> None:
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden_size)
        self.ssm = SelectiveSSMLayer(hidden_size, state_size, conv_kernel, dropout)
        self.norm2 = nn.LayerNorm(hidden_size)
        self.mlp = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * mlp_ratio),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * mlp_ratio, hidden_size),
            nn.Dropout(dropout),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = x + self.ssm(self.norm1(x))
        return x + self.mlp(self.norm2(x))


class ScMambaSSMScMAE(nn.Module):
    """Patch-tokenized scMAE whose encoder is a selective state-space model."""

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        depth: int = 2,
        patch_size: int = 32,
        state_size: int = 16,
        conv_kernel: int = 4,
        dropout: float = 0.1,
        mlp_ratio: int = 2,
        pool: str = "mean",
    ) -> None:
        super().__init__()
        if pool not in {"mean", "last"}:
            raise ValueError("pool must be 'mean' or 'last'")
        self.num_genes = int(num_genes)
        self.hidden_size = int(hidden_size)
        self.patch_size = int(patch_size)
        self.pool = pool
        self.tokenizer = PatchTokenizer(num_genes, patch_size, hidden_size, dropout)
        self.blocks = nn.ModuleList(
            [SSMBlock(hidden_size, state_size, conv_kernel, dropout, mlp_ratio) for _ in range(depth)]
        )
        self.norm = nn.LayerNorm(hidden_size)
        self.patch_decoder = nn.Sequential(
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size, patch_size),
        )
        self.mask_decoder = nn.Linear(hidden_size, patch_size)
        self.embedding_head = nn.Sequential(nn.Linear(hidden_size, hidden_size), nn.Tanh())

    @property
    def num_patches(self) -> int:
        return self.tokenizer.num_patches

    def encode_tokens(self, x: torch.Tensor, patch_mask: torch.Tensor | None = None) -> torch.Tensor:
        h = self.tokenizer(x, patch_mask)
        for block in self.blocks:
            h = block(h)
        return self.norm(h)

    def pool_tokens(self, h: torch.Tensor) -> torch.Tensor:
        pooled = h.mean(dim=1) if self.pool == "mean" else h[:, -1]
        return self.embedding_head(pooled)

    def unpatchify(self, patch_values: torch.Tensor) -> torch.Tensor:
        if patch_values.ndim != 3 or patch_values.shape[1] != self.num_patches or patch_values.shape[2] != self.patch_size:
            raise ValueError(
                f"patch_values must be [batch, {self.num_patches}, {self.patch_size}], got {tuple(patch_values.shape)}"
            )
        return patch_values.reshape(patch_values.shape[0], -1)[:, : self.num_genes]

    def forward(self, x: torch.Tensor, patch_mask: torch.Tensor | None = None) -> dict[str, torch.Tensor]:
        h = self.encode_tokens(x, patch_mask)
        reconstruction = self.unpatchify(self.patch_decoder(h))
        mask_logits = self.unpatchify(self.mask_decoder(h))
        return {
            "embedding": self.pool_tokens(h),
            "tokens": h,
            "reconstruction": reconstruction,
            "mask_logits": mask_logits,
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.pool_tokens(self.encode_tokens(x, None))
