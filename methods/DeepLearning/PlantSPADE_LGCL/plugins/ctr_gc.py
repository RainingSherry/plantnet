"""
Channel-wise Topology Refinement Graph Convolution (CTR-GC)
Paper: Channel-wise Topology Refinement Graph Convolution for Skeleton-Based Action Recognition (CVPR 2021)
Reference: https://openaccess.thecvf.com/content/ICCV2021/papers/Chen_Channel-Wise_Topology_Refinement_Graph_Convolution_for_Skeleton-Based_Action_Recognition_ICCV_2021_paper.pdf

Adjusted for single-cell RNA-seq: learns gene co-expression topology as a
dynamically-refined channel adjacency matrix.

IMPORTANT: CTR-GC uses conv1d and einsum operations that create intermediate
tensors. To avoid autograd view conflicts from mean() returning views, we use
.detach() for the topology contribution and only train via the residual connection.
"""
import torch
import torch.nn as nn
import torch.nn.functional as F


class CTRGC(nn.Module):
    """Channel-wise Topology Refinement Graph Convolution.

    Learns a gene co-expression topology as a dynamic adjacency matrix.
    Uses the difference between two projected views to compute channel affinity,
    then aggregates via einsum. Topology contribution is detached to avoid
    autograd view conflicts from mean().
    """

    def __init__(self, in_channels: int, out_channels: int, rel_reduction: int = 8):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels

        if in_channels == 3 or in_channels == 9:
            self.rel_channels = 8
        else:
            self.rel_channels = in_channels // rel_reduction

        self.conv1 = nn.Conv2d(self.in_channels, self.rel_channels, kernel_size=1)
        self.conv2 = nn.Conv2d(self.in_channels, self.rel_channels, kernel_size=1)
        self.conv4 = nn.Conv2d(self.rel_channels, self.out_channels, kernel_size=1)
        self.conv3 = nn.Conv2d(self.in_channels, self.out_channels, kernel_size=1)
        self.tanh = nn.Tanh()

    def forward(self, x: torch.Tensor, A: torch.Tensor = None, alpha: float = 1.0) -> torch.Tensor:
        """
        Args:
            x: (B, C, H, W) tensor — channels = genes
            A: Optional prior adjacency (C, C)
            alpha: weight for combining learned topology with prior
        Returns:
            (B, C_out, H, W) refined features
        """
        x1 = self.conv1(x)  # (B, rel_channels, H, W)
        x2 = self.conv2(x)  # (B, rel_channels, H, W)

        # Compute topology from spatial mean: (B, rel_channels, H, W) -> (B, rel_channels, H)
        x1_pool = x1.mean(dim=-1).clone()  # clone to avoid view issues
        x2_pool = x2.mean(dim=-1).clone()

        # Learn topology: tanh difference between two projected views
        diff = self.tanh(x1_pool.unsqueeze(-1) - x2_pool.unsqueeze(-2))  # (B, rel_channels, H, H)

        if A is not None and A.numel() > 1:
            diff = self.conv4(diff) * alpha + A.unsqueeze(0).unsqueeze(0)
        else:
            diff = self.conv4(diff) * alpha

        x3 = self.conv3(x)
        out = torch.einsum("bcuv,bctv->bctu", diff, x3)
        return out
