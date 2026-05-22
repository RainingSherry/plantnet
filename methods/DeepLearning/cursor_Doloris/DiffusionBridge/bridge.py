"""
Bridge module: connects source and target diffusion models.

The bridge orchestrates the full pipeline:
  1. Source diffusion: raw sparse expression -> shared Gaussian latent (DDIM reverse)
  2. Target diffusion: shared Gaussian latent -> cluster-friendly embedding (DDIM sample)

This is isomorphic to DOLORIS's DDIB framework, adapted for clustering:
  - DOLORIS: control cells -> shared latent -> perturbed cells
  - Bridge: raw sparse counts -> shared latent -> denoised cluster embeddings

The bridge also integrates a support anchor mechanism that uses the raw cell
expression structure to guide the denoised embedding.
"""
from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F

from source_diffusion import LatentDomainDiffusion, SourceDiffusion, TargetDiffusion


class BridgeSampler(nn.Module):
    """
    High-level bridge sampler that wraps the DiffusionBridge.

    Provides clean interfaces for:
    - ddim_reverse_sample_loop: encode raw cells to latent
    - ddim_sample_loop: decode latent to target embeddings
    """

    def __init__(self, bridge: "DiffusionBridge"):
        super().__init__()
        self.bridge = bridge

    def encode(
        self,
        raw_x: torch.Tensor,
        raw_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """
        Encode raw sparse expression to shared Gaussian latent via DDIM inversion.
        """
        return self.bridge.ddim_reverse_sample_loop(raw_x, raw_mask=raw_mask)

    def decode(
        self,
        latent: torch.Tensor,
        raw_x: torch.Tensor | None = None,
        raw_mask: torch.Tensor | None = None,
        condition: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """
        Decode shared latent to cluster-friendly embedding via DDIM sampling.
        """
        return self.bridge.ddim_sample_loop(
            latent, raw_x=raw_x, raw_mask=raw_mask, condition=condition
        )

    def forward(
        self,
        raw_x: torch.Tensor,
        raw_mask: torch.Tensor | None = None,
    ):
        """
        Full bridge forward: encode then decode.

        Returns:
            shared_latent: encoded representation in Gaussian space
            target_embedding: decoded cluster-friendly embedding
            support_anchor: support anchor from raw expression
        """
        return self.bridge(raw_x, raw_mask=raw_mask)


class DiffusionBridge(nn.Module):
    """
    Dual-Domain Diffusion Bridge.

    Architecture:
        SourceDiffusion: learns raw sparse expression distribution
                         (domain_dim = n_genes)
        TargetDiffusion: learns cluster-friendly denoised representation
                         (domain_dim = latent_dim)
        SupportAnchor:   encodes raw expression structure to guide target decoding

    The bridge connects these through the shared Gaussian latent space:

        raw_x --[source.encode]--> shared_latent --[target.ddim_sample]--> target_embedding
                   |                                                        ^
                   |                                                        |
                   +--[support_encoder]--> support_anchor (condition) -------+
    """

    def __init__(
        self,
        source: SourceDiffusion,
        target: TargetDiffusion,
        support_mask: nn.Module | None = None,
        support_hidden_dim: int = 256,
    ):
        super().__init__()
        self.source = source
        self.target = target
        self.support_mask = support_mask

        # Encode raw expression structure into a support anchor
        # that conditions the target diffusion sampling
        self.support_encoder = nn.Sequential(
            nn.Linear(source.domain_dim, support_hidden_dim),
            nn.LayerNorm(support_hidden_dim),
            nn.SiLU(),
            nn.Linear(support_hidden_dim, support_hidden_dim),
            nn.LayerNorm(support_hidden_dim),
            nn.SiLU(),
            nn.Linear(support_hidden_dim, target.shared_dim),
        )

    def encode_support(
        self, raw_x: torch.Tensor, raw_mask: torch.Tensor | None = None
    ) -> torch.Tensor:
        """
        Encode raw expression to a support anchor that captures expression structure.

        The support mask optionally blends observed expression with row-wise reference
        to handle zeros/dropouts gracefully.
        """
        masked_x = raw_x
        if self.support_mask is not None:
            masked_x = self.support_mask(masked_x, mask=raw_mask)
        return self.support_encoder(masked_x)

    def ddim_reverse_sample_loop(
        self, x: torch.Tensor, raw_mask: torch.Tensor | None = None, *args, **kwargs
    ) -> torch.Tensor:
        """DDIM reverse: encode raw domain to shared Gaussian latent."""
        if self.support_mask is not None:
            x = self.support_mask(x, mask=raw_mask)
        return self.source.ddim_reverse_sample_loop(x)

    def ddim_sample_loop(
        self,
        z: torch.Tensor,
        raw_x: torch.Tensor | None = None,
        raw_mask: torch.Tensor | None = None,
        condition: torch.Tensor | None = None,
        *args,
        **kwargs,
    ) -> torch.Tensor:
        """DDIM sample: decode shared latent to target domain."""
        if condition is None and raw_x is not None:
            condition = self.encode_support(raw_x, raw_mask=raw_mask)
        return self.target.ddim_sample_loop(z, cond=condition)

    def forward(
        self, raw_x: torch.Tensor, raw_mask: torch.Tensor | None = None
    ) -> "BridgeOutput":
        """
        Full bridge forward pass.

        Args:
            raw_x: raw sparse expression matrix (batch_size, n_genes)
            raw_mask: support mask for handling sparsity

        Returns:
            BridgeOutput with:
                - shared_latent: Gaussian latent from source encoding
                - target_embedding: cluster-friendly embedding from target decoding
                - support_anchor: expression structure conditioning signal
        """
        masked_x = raw_x
        if self.support_mask is not None:
            masked_x = self.support_mask(masked_x, mask=raw_mask)

        # Source: raw -> shared Gaussian latent (DDIM reverse)
        shared_latent = self.source.ddim_reverse_sample_loop(masked_x)

        # Support anchor: captures raw expression structure
        support_anchor = self.encode_support(raw_x, raw_mask=raw_mask)

        # Target: shared latent -> cluster-friendly embedding (DDIM sample)
        target_embedding = self.target.ddim_sample_loop(shared_latent, cond=support_anchor)

        return BridgeOutput(
            shared_latent=shared_latent,
            target_embedding=target_embedding,
            support_anchor=support_anchor,
        )


class BridgeOutput:
    """Container for bridge forward pass outputs."""

    def __init__(
        self,
        shared_latent: torch.Tensor,
        target_embedding: torch.Tensor,
        support_anchor: torch.Tensor,
    ):
        self.shared_latent = shared_latent
        self.target_embedding = target_embedding
        self.support_anchor = support_anchor


class GaussianBridgePrior(nn.Module):
    """
    Optional Gaussian prior module for regularizing the bridge latent space.

    This can be used during training to enforce that the encoded latent
    remains close to an isotropic Gaussian, complementing the implicit
    regularization from the DDIM reverse process.
    """

    def __init__(self, dim: int):
        super().__init__()
        self.proj = nn.Linear(dim, dim)

    def forward(self, z: torch.Tensor) -> torch.Tensor:
        return self.proj(z)

    def prior_loss(self, z: torch.Tensor) -> torch.Tensor:
        """KL(q(z) || N(0,I)) approximation via MSE to prevent mode collapse."""
        return z.pow(2).mean()
