from __future__ import annotations

import torch
from torch import nn

from .source_diffusion import DiffusionBridge


class BridgeSampler(nn.Module):
    def __init__(self, bridge: DiffusionBridge):
        super().__init__()
        self.bridge = bridge

    # Reverse path: raw sparse observation -> shared Gaussian latent.
    def ddim_reverse_sample_loop(
        self,
        model: nn.Module | None = None,
        image: torch.Tensor | None = None,
        raw_mask: torch.Tensor | None = None,
        **kwargs,
    ) -> torch.Tensor:
        if image is None:
            raise ValueError('image is required')
        return self.bridge.ddim_reverse_sample_loop(image, raw_mask=raw_mask, **kwargs)

    # Forward path: shared Gaussian latent -> cluster-friendly denoised embedding.
    def ddim_sample_loop(
        self,
        model: nn.Module | None = None,
        shape: tuple[int, int] | None = None,
        noise: torch.Tensor | None = None,
        raw_x: torch.Tensor | None = None,
        raw_mask: torch.Tensor | None = None,
        condition: torch.Tensor | None = None,
        **kwargs,
    ) -> torch.Tensor:
        if noise is None:
            if shape is None:
                raise ValueError('noise or shape is required')
            device = next(self.bridge.parameters()).device
            noise = torch.randn(*shape, device=device)
        return self.bridge.ddim_sample_loop(noise, raw_x=raw_x, raw_mask=raw_mask, condition=condition, **kwargs)
