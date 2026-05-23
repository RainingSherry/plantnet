"""
Target diffusion model for Dual-Domain Diffusion Bridge.

The target diffusion learns a cluster-separable denoised representation.
This is re-exported here for convenience.
"""
from source_diffusion import TargetDiffusion

__all__ = ["TargetDiffusion"]
