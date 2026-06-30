from __future__ import annotations

import math

import torch
import torch.nn as nn


def timestep_embedding(values: torch.Tensor, dim: int, max_period: int = 10000) -> torch.Tensor:
    if values.ndim != 1:
        raise ValueError(f"values must be [batch], got {tuple(values.shape)}")
    half = dim // 2
    freqs = torch.exp(
        -math.log(max_period) * torch.arange(0, half, dtype=torch.float32, device=values.device) / max(1, half)
    )
    args = values.float().view(-1, 1) * freqs.view(1, -1)
    emb = torch.cat([torch.cos(args), torch.sin(args)], dim=1)
    if dim % 2:
        emb = torch.cat([emb, torch.zeros_like(emb[:, :1])], dim=1)
    return emb


class FiLMResidualBlock(nn.Module):
    def __init__(self, hidden_size: int, time_size: int, dropout: float) -> None:
        super().__init__()
        self.norm = nn.LayerNorm(hidden_size)
        self.film = nn.Sequential(
            nn.SiLU(),
            nn.Linear(time_size, hidden_size * 2),
        )
        self.net = nn.Sequential(
            nn.Linear(hidden_size, hidden_size * 4),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_size * 4, hidden_size),
        )

    def forward(self, x: torch.Tensor, time_emb: torch.Tensor) -> torch.Tensor:
        scale, shift = self.film(time_emb).chunk(2, dim=1)
        h = self.norm(x)
        h = h * (1.0 + scale) + shift
        return x + self.net(h)


class ConsistencyExpressionDenoiser(nn.Module):
    """Karras-style consistency denoiser for scMAE expression vectors.

    The network predicts the model-output term used in consistency models. The
    clean-expression estimate is assembled with the boundary-condition scaling
    from OpenAI's consistency model code:

        denoised = c_skip(sigma) * x_t + c_out(sigma) * model(x_t, sigma)

    `feature()` returns the hidden representation produced at the clean boundary
    (`sigma_min`) and is used for the KMeans benchmark.
    """

    def __init__(
        self,
        num_genes: int,
        hidden_size: int = 128,
        depth: int = 4,
        time_size: int = 128,
        dropout: float = 0.1,
        sigma_data: float = 0.5,
        sigma_min: float = 0.002,
    ) -> None:
        super().__init__()
        if num_genes <= 0:
            raise ValueError("num_genes must be positive")
        if hidden_size <= 0 or depth <= 0:
            raise ValueError("hidden_size and depth must be positive")
        self.num_genes = int(num_genes)
        self.sigma_data = float(sigma_data)
        self.sigma_min = float(sigma_min)
        self.time_size = int(time_size)

        self.time_mlp = nn.Sequential(
            nn.Linear(time_size, time_size),
            nn.SiLU(),
            nn.Linear(time_size, time_size),
        )
        self.input = nn.Sequential(
            nn.Dropout(dropout),
            nn.Linear(num_genes, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
        )
        self.blocks = nn.ModuleList(
            [FiLMResidualBlock(hidden_size, time_size, dropout) for _ in range(depth)]
        )
        self.final_norm = nn.LayerNorm(hidden_size)
        self.model_output = nn.Linear(hidden_size, num_genes)
        self.mask_head = nn.Linear(hidden_size, num_genes)
        self.projector = nn.Sequential(
            nn.LayerNorm(hidden_size),
            nn.Linear(hidden_size, hidden_size),
            nn.GELU(),
            nn.Linear(hidden_size, hidden_size),
        )

    def scalings_for_boundary_condition(self, sigma: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        sigma = sigma.float()
        sigma_data = torch.as_tensor(self.sigma_data, dtype=sigma.dtype, device=sigma.device)
        sigma_min = torch.as_tensor(self.sigma_min, dtype=sigma.dtype, device=sigma.device)
        c_skip = sigma_data.square() / ((sigma - sigma_min).square() + sigma_data.square())
        c_out = (sigma - sigma_min) * sigma_data / (sigma.square() + sigma_data.square()).sqrt()
        c_in = 1.0 / (sigma.square() + sigma_data.square()).sqrt()
        return c_skip, c_out, c_in

    def encode_noisy(self, x_t: torch.Tensor, sigma: torch.Tensor) -> torch.Tensor:
        if x_t.ndim != 2 or x_t.shape[1] != self.num_genes:
            raise ValueError(f"x_t must be [batch, {self.num_genes}], got {tuple(x_t.shape)}")
        if sigma.ndim != 1 or sigma.shape[0] != x_t.shape[0]:
            raise ValueError(f"sigma must be [batch], got {tuple(sigma.shape)} for x_t {tuple(x_t.shape)}")
        _, _, c_in = self.scalings_for_boundary_condition(sigma)
        h = self.input(x_t * c_in.view(-1, 1))
        rescaled_t = 250.0 * torch.log(sigma.float().clamp_min(1e-44))
        time_emb = self.time_mlp(timestep_embedding(rescaled_t, self.time_size))
        for block in self.blocks:
            h = block(h, time_emb)
        return self.final_norm(h)

    def forward(self, x_t: torch.Tensor, sigma: torch.Tensor) -> dict[str, torch.Tensor]:
        h = self.encode_noisy(x_t, sigma)
        raw = self.model_output(h)
        c_skip, c_out, _ = self.scalings_for_boundary_condition(sigma)
        denoised = c_skip.view(-1, 1) * x_t + c_out.view(-1, 1) * raw
        return {
            "embedding": self.projector(h),
            "model_output": raw,
            "denoised": denoised,
            "mask_logits": self.mask_head(h),
        }

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        sigma = torch.full((x.shape[0],), self.sigma_min, dtype=x.dtype, device=x.device)
        h = self.encode_noisy(x, sigma)
        return self.projector(h)
