from __future__ import annotations

import torch
from torch import nn


def build_support_mask(x: torch.Tensor, hvg_mask: torch.Tensor | None = None, topk: int | None = None) -> torch.Tensor:
    # Support mask marks the observed / trustworthy coordinates that should stay on-manifold.
    mask = (x > 0).float()
    if hvg_mask is not None:
        if hvg_mask.dim() == 1:
            hvg_mask = hvg_mask.unsqueeze(0)
        mask = mask * hvg_mask.float()
    if topk is not None and topk > 0 and topk < mask.shape[-1]:
        _, indices = torch.topk(x.abs(), topk, dim=-1)
        sparse = torch.zeros_like(mask)
        sparse.scatter_(dim=-1, index=indices, value=1.0)
        mask = mask * sparse
    return mask


def apply_support_projection(x: torch.Tensor, mask: torch.Tensor, blend: float = 0.2) -> torch.Tensor:
    if mask is None:
        return x
    # Unobserved coordinates are softly pulled toward the row-wise reference instead of being hard-clamped.
    reference = x.mean(dim=-1, keepdim=True)
    return x * mask + blend * reference * (1.0 - mask)


class GeneSupportMask(nn.Module):
    def __init__(self, mask: torch.Tensor | None = None, blend: float = 0.2):
        super().__init__()
        self.blend = blend
        if mask is None:
            self.register_buffer('gene_mask', torch.tensor([]), persistent=False)
        else:
            self.register_buffer('gene_mask', mask.float(), persistent=False)

    def forward(self, x: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
        gene_mask = mask
        if gene_mask is None and self.gene_mask.numel() > 0:
            gene_mask = self.gene_mask
        if gene_mask is None:
            return x
        if gene_mask.dim() == 1:
            gene_mask = gene_mask.unsqueeze(0)
        if gene_mask.shape[0] != x.shape[0]:
            gene_mask = gene_mask.expand(x.shape[0], -1)
        return apply_support_projection(x, gene_mask, blend=self.blend)
