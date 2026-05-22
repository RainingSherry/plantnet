"""
Support mask module: handles sparsity in raw single-cell expression.

The support mask strategy is inspired by DOLORIS's sparsity masking approach.
It marks observed / trustworthy expression values and softly blends
unobserved positions (zeros / dropouts) with a row-wise reference.

This prevents the model from being dominated by zeros while preserving
biological signal in the observed values.
"""
from __future__ import annotations

import torch
from torch import nn


def build_support_mask(
    x: torch.Tensor,
    hvg_mask: torch.Tensor | None = None,
    topk: int | None = None,
    sparsity_threshold: float = 0.0,
) -> torch.Tensor:
    """
    Build a support mask marking observed / trustworthy expression values.

    Strategy:
    1. Start with observed entries (x > sparsity_threshold)
    2. Optionally intersect with HVG mask
    3. Optionally keep only top-k per row by absolute expression

    Args:
        x: expression matrix (batch_size, n_genes)
        hvg_mask: binary mask for highly variable genes
        topk: if set, keep only top-k entries per row
        sparsity_threshold: threshold for considering a value as "observed"

    Returns:
        Binary mask tensor (batch_size, n_genes)
    """
    mask = (x > sparsity_threshold).float()

    if hvg_mask is not None:
        if hvg_mask.dim() == 1:
            hvg_mask = hvg_mask.unsqueeze(0)
        mask = mask * hvg_mask.float()

    if topk is not None and topk > 0 and topk < mask.shape[-1]:
        _, indices = torch.topk(x.abs(), topk, dim=-1)
        sparse_mask = torch.zeros_like(mask)
        sparse_mask.scatter_(dim=-1, index=indices, value=1.0)
        mask = mask * sparse_mask

    return mask


def apply_support_projection(
    x: torch.Tensor, mask: torch.Tensor, blend: float = 0.2
) -> torch.Tensor:
    """
    Apply support projection: softly blend unobserved positions with row-wise reference.

    Unobserved entries (mask=0) are not hard-clamped to zero. Instead, they are
    softly pulled toward the row-wise mean. This preserves gradient flow for
    dropout positions while respecting the observed structure.

    Args:
        x: expression matrix (batch_size, n_genes)
        mask: support mask (1=observed, 0=unobserved)
        blend: blending strength for unobserved positions

    Returns:
        Projected expression matrix
    """
    reference = x.mean(dim=-1, keepdim=True)
    return x * mask + blend * reference * (1.0 - mask)


class GeneSupportMask(nn.Module):
    """
    Learnable support mask module for gene expression.

    Applied to raw input before entering the bridge encoder.
    Keeps observed values as-is and softly blends unobserved positions.
    """

    def __init__(
        self,
        gene_mask: torch.Tensor | None = None,
        blend: float = 0.2,
        learnable_blend: bool = False,
    ):
        super().__init__()
        self.blend = blend
        self.learnable_blend = learnable_blend
        if gene_mask is None:
            self.register_buffer("gene_mask", torch.tensor([]), persistent=False)
        else:
            self.register_buffer("gene_mask", gene_mask.float(), persistent=False)
        if learnable_blend:
            self.blend_param = nn.Parameter(torch.tensor(blend))

    @property
    def blend_value(self) -> float:
        if self.learnable_blend:
            return torch.sigmoid(self.blend_param).item()
        return self.blend

    def forward(self, x: torch.Tensor, mask: torch.Tensor | None = None) -> torch.Tensor:
        """
        Apply support projection to input.

        Args:
            x: expression matrix (batch_size, n_genes)
            mask: optional per-sample support mask

        Returns:
            Projected matrix
        """
        gene_mask = mask
        if gene_mask is None and self.gene_mask.numel() > 0:
            gene_mask = self.gene_mask
        if gene_mask is None:
            return x
        if gene_mask.dim() == 1:
            gene_mask = gene_mask.unsqueeze(0)
        if gene_mask.shape[0] != x.shape[0]:
            gene_mask = gene_mask.expand(x.shape[0], -1)
        return apply_support_projection(x, gene_mask, blend=self.blend_value)


class SparsityPredictor(nn.Module):
    """
    Predicts gene activation (zero vs. non-zero) for each gene in each cell.

    This is inspired by DOLORIS's Mask Model. It learns to predict which genes
    are silenced (zero) under the target representation, helping the diffusion
    model focus on expressed genes.

    The predictor takes the cluster-friendly embedding as input and outputs
    a probability for each gene being active (non-zero).
    """

    def __init__(
        self,
        embedding_dim: int,
        n_genes: int,
        hidden_dim: int = 128,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(embedding_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.SiLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, n_genes),
        )

    def forward(self, embedding: torch.Tensor) -> torch.Tensor:
        """
        Predict gene activation probabilities.

        Args:
            embedding: cluster-friendly embedding (batch_size, embedding_dim)

        Returns:
            Activation probabilities (batch_size, n_genes) in [0, 1]
        """
        return torch.sigmoid(self.net(embedding))

    def mask_loss(
        self, pred_prob: torch.Tensor, target_binary: torch.Tensor
    ) -> torch.Tensor:
        """
        BCE loss between predicted activation probabilities and ground truth.

        Args:
            pred_prob: predicted probabilities (batch_size, n_genes)
            target_binary: ground truth binary mask (1=active, 0=silent)

        Returns:
            BCE loss (scalar)
        """
        return torch.nn.functional.binary_cross_entropy(
            pred_prob, target_binary, reduction="mean"
        )
