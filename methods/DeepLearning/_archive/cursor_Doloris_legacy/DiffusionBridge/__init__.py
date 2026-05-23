"""
Dual-Domain Diffusion Bridge for Single-Cell Clustering

A bridge-based representation learning framework for scRNA-seq clustering.

Core Architecture:
    Raw Sparse Count Domain
           |
           v  (Source Diffusion: DDIM reverse)
    Shared Gaussian Latent
           |
           v  (Target Diffusion: DDIM sample + support anchor)
    Cluster-Separable Denoised Domain
           |
           v  (DEC Cluster Head)
    Soft Cluster Assignments

Mathematical Narrative:
    Traditional dimensionality reduction is a lossy compression that projects
    all cells into a shared manifold without distinguishing cluster boundaries.
    Our bridge-based representation learning constructs a path through the
    high-dimensional space: raw sparse observations are first encoded into a
    Gaussian latent space (preserving distributional information), then decoded
    into a cluster-separable denoised space (separating cell types).

Reference:
    Inspired by DOLORIS (ICLR 2026): Dual Conditional Diffusion Implicit Bridges
    with Sparsity Masking Strategy for Unpaired Single-Cell Perturbation Estimation.
"""

from source_diffusion import (
    LatentDomainDiffusion,
    SourceDiffusion,
    TargetDiffusion,
)
from bridge import (
    BridgeOutput,
    BridgeSampler,
    DiffusionBridge,
    GaussianBridgePrior,
)
from support_mask import (
    GeneSupportMask,
    SparsityPredictor,
    apply_support_projection,
    build_support_mask,
)
from cluster_head import ClusterHead

__all__ = [
    # Diffusion models
    "LatentDomainDiffusion",
    "SourceDiffusion",
    "TargetDiffusion",
    # Bridge
    "BridgeOutput",
    "BridgeSampler",
    "DiffusionBridge",
    "GaussianBridgePrior",
    # Support
    "GeneSupportMask",
    "SparsityPredictor",
    "apply_support_projection",
    "build_support_mask",
    # Clustering
    "ClusterHead",
]
