"""
Neural network utilities for Dual-Domain Diffusion Bridge.
Adapted from DOLORIS codebase.
"""
from __future__ import annotations

import math

import torch as th
import torch.nn as nn


def timestep_embedding(timesteps: th.Tensor, dim: int, max_period: float = 10000.0) -> th.Tensor:
    """
    Create sinusoidal timestep embeddings.

    :param timesteps: a 1-D Tensor of N indices, one per batch element.
    :param dim: the dimension of the output.
    :param max_period: controls the minimum frequency of the embeddings.
    :return: an [N x dim] Tensor of positional embeddings.
    """
    half = dim // 2
    freqs = th.exp(
        -math.log(max_period) * th.arange(start=0, end=half, dtype=th.float32) / half
    ).to(device=timesteps.device)
    args = timesteps[:, None].float() * freqs[None]
    embedding = th.cat([th.cos(args), th.sin(args)], dim=-1)
    if dim % 2:
        embedding = th.cat([embedding, th.zeros_like(embedding[:, :1])], dim=-1)
    return embedding


def update_ema(target_params, source_params, rate: float = 0.99):
    """
    Update target parameters to be closer to those of source parameters using
    an exponential moving average.

    :param rate: the EMA rate (closer to 1 means slower).
    """
    for targ, src in zip(target_params, source_params):
        targ.detach().mul_(rate).add_(src, alpha=1 - rate)


def mean_flat(tensor: th.Tensor) -> th.Tensor:
    """Take the mean over all non-batch dimensions."""
    return tensor.mean(dim=list(range(1, len(tensor.shape))))


class SiLU(nn.Module):
    def forward(self, x):
        return x * th.sigmoid(x)
