"""
models/support_pooling.py
==========================
表达支撑集聚合成细胞嵌入（Support Graph Pooling）。

核心思想：
  给定每个细胞的表达支撑集 S_c = { g | X_{cg} > 0 }，
  以及对应的边权重 w_{cg}，通过加权聚合得到细胞嵌入：

    z_c = Pool_{g ∈ S_c}( w_{cg} · h_g )

  其中 h_g 来自 GeneGATEncoder 的输出（基因上下文嵌入）。
"""

from __future__ import annotations

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from typing import Optional, List, Literal


# ---------------------------------------------------------------------------
# 简化的加性注意力聚合器
# ---------------------------------------------------------------------------

class AttentionAggregator(nn.Module):
    """
    简化的加性注意力聚合器。

    score_{cg} = v^T · tanh(W_k · h_g + W_w · w_{cg})
    alpha_{cg} = softmax_c(score_{cg})
    z_c = sum_g alpha_{cg} · h_g
    """

    def __init__(
        self,
        gene_emb_dim: int,
        cell_emb_dim: int = 0,
        hidden_dim: int = 128,
        n_heads: int = 4,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.gene_emb_dim = gene_emb_dim
        self.hidden_dim = hidden_dim

        self.key_net = nn.Sequential(
            nn.Linear(gene_emb_dim + 1, hidden_dim),
            nn.Tanh(),
        )
        self.query_net = nn.Sequential(
            nn.Linear(gene_emb_dim, hidden_dim),
            nn.Tanh(),
        )
        self.v = nn.Linear(hidden_dim, 1, bias=False)

    def forward(
        self,
        gene_emb: torch.Tensor,       # [B, S, D]
        support_weight: torch.Tensor, # [B, S]
        mask: Optional[torch.Tensor] = None,  # [B, S]
        cell_emb: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        B, S, D = gene_emb.shape

        # 查询向量
        if cell_emb is not None:
            q = self.query_net(cell_emb)  # [B, hidden]
        else:
            if mask is not None:
                masked = gene_emb * mask.unsqueeze(-1)
                pooled = masked.sum(dim=1) / (mask.sum(dim=1, keepdim=True) + 1e-8)
            else:
                pooled = gene_emb.mean(dim=1)
            q = self.query_net(pooled)  # [B, hidden]

        # 键
        k_input = torch.cat([gene_emb, support_weight.unsqueeze(-1)], dim=-1)  # [B, S, D+1]
        k = self.key_net(k_input)  # [B, S, hidden]

        # 注意力分数
        score = self.v(torch.tanh(k + q.unsqueeze(1))).squeeze(-1)  # [B, S]

        if mask is not None:
            score = score.masked_fill(mask == 0, -1e9)

        alpha = F.softmax(score, dim=1)  # [B, S]
        if mask is not None:
            alpha = alpha * mask

        z = (alpha.unsqueeze(-1) * gene_emb).sum(dim=1)  # [B, D]
        return z


# ---------------------------------------------------------------------------
# Mean / Sum Pooling（简单基线）
# ---------------------------------------------------------------------------

class MeanAggregator(nn.Module):
    """简单均值聚合。"""

    def __init__(self, gene_emb_dim: int):
        super().__init__()
        self.gene_emb_dim = gene_emb_dim

    def forward(
        self,
        gene_emb: torch.Tensor,
        support_weight: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        cell_emb: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        if mask is not None:
            masked = gene_emb * mask.unsqueeze(-1)
            count = mask.sum(dim=1, keepdim=True) + 1e-8
            return masked.sum(dim=1) / count
        else:
            return gene_emb.mean(dim=1)


class WeightedSumAggregator(nn.Module):
    """加权求和聚合。"""

    def __init__(self, gene_emb_dim: int):
        super().__init__()
        self.gene_emb_dim = gene_emb_dim

    def forward(
        self,
        gene_emb: torch.Tensor,
        support_weight: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        cell_emb: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        w = support_weight.unsqueeze(-1)
        if mask is not None:
            w = w * mask.unsqueeze(-1)
        weighted = gene_emb * w
        return weighted.sum(dim=1)


# ---------------------------------------------------------------------------
# Top-K Aggregator
# ---------------------------------------------------------------------------

class TopKAttentionAggregator(nn.Module):
    """Top-K 注意力聚合。"""

    def __init__(
        self,
        gene_emb_dim: int,
        cell_emb_dim: int = 0,
        hidden_dim: int = 128,
        k: int = 50,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.k = k
        self.attn_agg = AttentionAggregator(
            gene_emb_dim, cell_emb_dim, hidden_dim, n_heads=4, dropout=dropout
        )

    def forward(
        self,
        gene_emb: torch.Tensor,
        support_weight: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        cell_emb: Optional[torch.Tensor] = None,
    ) -> torch.Tensor:
        return self.attn_agg(gene_emb, support_weight, mask, cell_emb)


# ---------------------------------------------------------------------------
# 工厂函数
# ---------------------------------------------------------------------------

class SupportPoolingFactory:
    _POOLING_MAP = {
        "attention": AttentionAggregator,
        "mean": MeanAggregator,
        "weighted_sum": WeightedSumAggregator,
        "topk": TopKAttentionAggregator,
    }

    @classmethod
    def create(
        cls,
        pooling_strategy: Literal["attention", "mean", "weighted_sum", "topk"],
        gene_emb_dim: int,
        cell_emb_dim: int = 0,
        hidden_dim: int = 128,
        n_heads: int = 4,
        dropout: float = 0.1,
        topk_k: int = 50,
    ) -> nn.Module:
        if pooling_strategy not in cls._POOLING_MAP:
            raise ValueError(
                f"Unknown pooling_strategy '{pooling_strategy}'. "
                f"Available: {list(cls._POOLING_MAP.keys())}"
            )
        if pooling_strategy == "attention":
            return AttentionAggregator(
                gene_emb_dim=gene_emb_dim,
                cell_emb_dim=cell_emb_dim,
                hidden_dim=hidden_dim,
                n_heads=n_heads,
                dropout=dropout,
            )
        elif pooling_strategy == "topk":
            return TopKAttentionAggregator(
                gene_emb_dim=gene_emb_dim,
                cell_emb_dim=cell_emb_dim,
                hidden_dim=hidden_dim,
                k=topk_k,
                dropout=dropout,
            )
        elif pooling_strategy == "mean":
            return MeanAggregator(gene_emb_dim=gene_emb_dim)
        else:
            return WeightedSumAggregator(gene_emb_dim=gene_emb_dim)


class ExplainableCellEmbedding(nn.Module):
    def __init__(
        self,
        aggregator: nn.Module,
        gene_names: List[str],
        topk_explain: int = 20,
    ):
        super().__init__()
        self.aggregator = aggregator
        self.gene_names = gene_names
        self.topk_explain = topk_explain

    def forward(
        self,
        gene_emb: torch.Tensor,
        support_weight: torch.Tensor,
        mask: Optional[torch.Tensor] = None,
        cell_emb: Optional[torch.Tensor] = None,
    ) -> dict:
        cell_z = self.aggregator(gene_emb, support_weight, mask, cell_emb)
        return {"cell_z": cell_z}


if __name__ == "__main__":
    import torch
    B, S, D = 8, 200, 128
    gene_emb = torch.randn(B, S, D)
    support_weight = torch.rand(B, S)
    mask = (torch.rand(B, S) > 0.2).float()

    for strategy in ["attention", "mean", "weighted_sum", "topk"]:
        pooler = SupportPoolingFactory.create(
            pooling_strategy=strategy,
            gene_emb_dim=D,
            hidden_dim=128,
            topk_k=30,
        )
        z = pooler(gene_emb, support_weight, mask)
        print(f"{strategy}: output shape = {z.shape}")

    print("support_pooling.py test passed.")
