# -*- coding: utf-8 -*-
"""
ScSpade: Support-Masked Diffusion Autoencoder for scRNA-seq Clustering

A generative deep learning framework for single-cell RNA-seq clustering
inspired by DOLORIS's sparsity masking strategy.

Core idea:
    X → M = 1(X > 0) → z → z_diff → clustering

Key innovations:
    1. SupportMaskNet: Learns gene activation patterns (expressed vs. zero)
    2. LatentDiffusionAE: Denoises latent representations via diffusion
    3. Joint training: Mask guides diffusion, diffusion improves mask prediction
"""

__version__ = '0.1.0'

from .models import *
from .data import load_and_preprocess, create_dataloader
from .eval import evaluate_clustering, run_all_evaluations

__all__ = [
    'load_and_preprocess',
    'create_dataloader',
    'evaluate_clustering',
    'run_all_evaluations',
]
