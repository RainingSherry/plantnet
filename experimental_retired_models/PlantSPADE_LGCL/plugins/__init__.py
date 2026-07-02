"""
PlantSPADE-LGCL Plug-in Modules
================================
Plugins adapted from /data/luolie/缝合模块:

1. FCR (Frequency Contrastive Regularization, ECCV 2024)
   - fcr_loss.py: FFT-based contrastive loss for cell embeddings
2. PolaContrastiveLoss (Polarity-aware Contrastive Loss, adapted from PolaFormer ICLR 2025)
   - pola_linear_attention.py: Polarity decomposition for gene co-expression modeling
3. BiSSM1D (Bidirectional State Space Model, based on Mamba AAAI 2025)
   - bimamba2_1d.py: Bidirectional SSM for long-range dependencies
4. CTR-GC (Channel-wise Topology Refinement Graph Convolution, CVPR 2021)
   - ctr_gc.py: Dynamic gene co-expression topology learning
"""
from .fcr_loss import FCRLoss
from .pola_linear_attention import PolaContrastiveLoss
from .bimamba2_1d import BiSSM1D
from .ctr_gc import CTRGC

__all__ = [
    "FCRLoss",
    "PolaContrastiveLoss",
    "BiSSM1D",
    "CTRGC",
]
