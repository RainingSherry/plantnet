"""
graphs/build_cell_gene_graph.py
================================
构建细胞-基因二部支撑图（cell-gene bipartite support graph）。

核心思想：
  每个细胞不是通过全量 20,000 个基因进入模型，
  而是通过它"真正表达的基因子集"进入图神经网络做消息传递。

  S_c = { g | X_{cg} > 0 }          ← 细胞 c 的表达支撑集
  w_{cg} ∈ {log1p_count, rank, tfidf, norm_count}  ← 边权重

  z_c = Pool_{g ∈ S_c}( w_{cg} · h_g )   ← 聚合得到细胞嵌入

这样，零值基因不会引入噪声，而表达证据通过基因图结构扩散，
形成稳定的细胞状态表示——天然可解释（每个细胞嵌入可追溯到高贡献基因）。

图结构（二部图）：
  - 左边：n_cells 个细胞节点
  - 右边：n_hvg 个 HVG 基因节点
  - 边：细胞 c 表达基因 g（X_{cg} > 0），权重为 w_{cg}

支持三种权重策略：
  1. log1p_count  — log(1 + X_{cg})
  2. rank_weight  — 基于细胞内排名的归一化权重
  3. tfidf_style — TF-IDF 风格（gene rarity × cell expression）
  4. norm_count   — 归一化计数
"""

from __future__ import annotations

from typing import Literal, Optional, Tuple
import warnings

import numpy as np
import scipy.sparse as sp
import torch


# ---------------------------------------------------------------------------
# 权重计算策略
# ---------------------------------------------------------------------------

def compute_log1p_weight(X_expr: np.ndarray) -> np.ndarray:
    """
    边权重 = log(1 + X_{cg})，压缩大值、保留相对尺度。
    返回 (n_cells, n_hvg) 的密集权重矩阵。
    """
    return np.log1p(X_expr).astype(np.float32)


def compute_rank_weight(X_expr: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """
    边权重 = 基于细胞内排名的归一化权重。
    在每个细胞内，对基因表达值排名，然后归一化到 [0, 1]。
    排名的倒数——表达量越高，权重越高。
    """
    n_cells, n_genes = X_expr.shape
    rank_weight = np.zeros_like(X_expr, dtype=np.float32)

    for c in range(n_cells):
        row = X_expr[c]
        nonzero_idx = np.where(row > 0)[0]
        if len(nonzero_idx) == 0:
            continue
        # 排名：最高的 expr -> 最大的 rank（最前面）
        nonzero_vals = row[nonzero_idx]
        sorted_order = np.argsort(nonzero_vals)
        ranks = np.zeros(len(nonzero_idx), dtype=np.float32)
        ranks[sorted_order] = np.arange(len(nonzero_idx), 0, -1, dtype=np.float32)  # 最高表达 = 最大 rank
        # 归一化到 [0, 1]
        r_min, r_max = ranks.min(), ranks.max()
        if r_max > r_min:
            norm = (ranks - r_min) / (r_max - r_min + eps)
        else:
            norm = np.ones_like(ranks) * 0.5
        rank_weight[c, nonzero_idx] = norm

    return rank_weight


def compute_tfidf_style_weight(
    X_expr: np.ndarray,
    gene_freq: Optional[np.ndarray] = None,
    smooth_idf: float = 1.0,
) -> np.ndarray:
    """
    TF-IDF 风格权重：
      TF = log(1 + X_{cg})
      IDF = log( (N + 1) / (df(g) + 1) ) + 1
      权重 = TF × IDF
    其中 df(g) = 表达基因 g 的细胞数。
    """
    n_cells, n_genes = X_expr.shape
    tf = np.log1p(X_expr).astype(np.float32)

    if gene_freq is None:
        # df(g) = 非零细胞数
        gene_freq = (X_expr > 0).sum(axis=0)  # shape (n_genes,)

    # IDF 平滑
    idf = np.log((n_cells + smooth_idf) / (gene_freq + smooth_idf)) + smooth_idf
    idf = np.asarray(idf, dtype=np.float32).reshape(1, -1)

    weight = tf * idf
    # 再次归一化
    row_max = weight.max(axis=1, keepdims=True)
    row_max = np.where(row_max > 0, row_max, 1.0)
    weight = weight / (row_max + 1e-8)

    return weight


def compute_norm_count(X_expr: np.ndarray) -> np.ndarray:
    """边权重 = L2 归一化后的计数向量。"""
    norms = np.linalg.norm(X_expr, axis=1, keepdims=True)
    norms = np.where(norms > 0, norms, 1.0)
    return (X_expr / norms).astype(np.float32)


WEIGHT_STRATEGIES = {
    "log1p": compute_log1p_weight,
    "rank": compute_rank_weight,
    "tfidf": compute_tfidf_style_weight,
    "norm": compute_norm_count,
}


# ---------------------------------------------------------------------------
# 构建细胞-基因二部图（稀疏表示）
# ---------------------------------------------------------------------------

def build_cell_gene_bipartite_graph(
    X: np.ndarray,
    gene_names: np.ndarray,
    support_strategy: Literal["log1p", "rank", "tfidf", "norm"] = "log1p",
    dropout_rate: float = 0.0,
    random_seed: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    构建细胞-基因二部支撑图（稀疏表示，供后续 SupportPooling 使用）。

    参数
    ----
    X : (n_cells, n_hvg)  仅包含 HVG 基因的表达矩阵
        （通常来自 build_gene_graph 选出的 HVG）
    gene_names : (n_hvg,)  对应的基因名
    support_strategy : str  权重策略
        - "log1p"  : log(1 + count)，默认推荐
        - "rank"   : 排名归一化权重
        - "tfidf"  : TF-IDF 风格（突出稀有但重要的基因）
        - "norm"   : L2 归一化计数
    dropout_rate : float  随机丢弃表达边的概率（数据增强，防止过拟合）
    random_seed : int

    返回
    ----
    cell_gene_edge_src : [E]  源节点（细胞节点，0 ~ n_cells-1）
    cell_gene_edge_dst : [E]  目标节点（基因节点，0 ~ n_hvg-1）
    cell_gene_edge_w   : [E]  边权重（归一化后的权重）
    cell_gene_edge_mask : [E]  边掩码（0 = 被 dropout，1 = 保留）
    nonzero_gene_mask_per_cell : (n_cells, n_hvg) bool  原始非零掩码

    重要：
      返回边索引是稀疏 COO 格式（已去除权重为0的边）。
      细胞节点索引范围 [0, n_cells)
      基因节点索引范围 [n_cells, n_cells + n_hvg)
    """
    rng = np.random.default_rng(random_seed)
    n_cells, n_hvg = X.shape

    # ---- Step 1: 计算边权重 ----
    if support_strategy not in WEIGHT_STRATEGIES:
        warnings.warn(f"Unknown support_strategy '{support_strategy}', using 'log1p'")
        support_strategy = "log1p"

    weight_matrix = WEIGHT_STRATEGIES[support_strategy](X)
    nonzero_mask = (X > 0)

    # ---- Step 2: Dropout 增强 ----
    if dropout_rate > 0:
        keep_mask = rng.uniform(0, 1, size=X.shape) > dropout_rate
        effective_mask = nonzero_mask & keep_mask
    else:
        effective_mask = nonzero_mask.copy()

    # ---- Step 3: 构建稀疏 COO 格式边 ----
    row_idx, col_idx = np.where(effective_mask)
    edge_weights = weight_matrix[effective_mask]

    # 归一化边权重（按细胞内归一化，使每个细胞的出边权重和为 1）
    row_sum = np.zeros(n_cells, dtype=np.float32)
    np.add.at(row_sum, row_idx, edge_weights)
    row_sum = np.where(row_sum > 0, row_sum, 1.0)
    edge_weights_norm = edge_weights / row_sum[row_idx]

    # ---- Step 4: 建立节点索引映射 ----
    # 细胞节点: 0 ~ n_cells-1
    # 基因节点: n_cells ~ n_cells+n_hvg-1
    gene_node_offset = n_cells  # 基因节点在全局图中的起始索引

    cell_node_ids = row_idx.astype(np.int64)                                    # [E]
    gene_node_ids = (col_idx + gene_node_offset).astype(np.int64)             # [E]

    # 二部边：从细胞指向基因（cell → gene）
    edge_src = np.concatenate([cell_node_ids, gene_node_ids])                 # [2E]
    edge_dst = np.concatenate([gene_node_ids, cell_node_ids])                   # [2E]
    edge_w = np.concatenate([edge_weights_norm, edge_weights_norm])              # [2E]（无向化）

    # 原始边掩码（用于计算 Mask Loss）
    original_nonzero_flat = nonzero_mask[row_idx]  # [E]
    all_nonzero_flat = np.ones(len(edge_weights_norm), dtype=np.float32)       # [E]

    print(f"  [build_cell_gene_graph] {n_cells} cells × {n_hvg} HVG genes")
    print(f"  [build_cell_gene_graph] {len(row_idx)} active edges, dropout_rate={dropout_rate}")

    return (
        edge_src.astype(np.int64),
        edge_dst.astype(np.int64),
        edge_w.astype(np.float32),
        all_nonzero_flat.astype(np.float32),  # edge_mask
        nonzero_mask.astype(np.float32),       # (n_cells, n_hvg) nonzero per cell
    )


# ---------------------------------------------------------------------------
# SparseTensor 格式（供 torch_geometric 直接使用）
# ---------------------------------------------------------------------------

def build_pyg_bipartite_edges(
    X: np.ndarray,
    support_strategy: str = "log1p",
    dropout_rate: float = 0.0,
    random_seed: int = 42,
) -> dict:
    """
    构建 PyG 格式的二部图边，返回可直接用于 SparseAdj 或 Data 的字典。

    返回结构
    --------
    {
        "edge_index": [2, E],           # torch.long，边索引
        "edge_attr": [E, 1],           # torch.float32，边权重
        "n_cells": int,
        "n_genes": int,
        "n_edges": int,
    }
    """
    n_cells, n_genes = X.shape

    # 构建边
    if support_strategy not in WEIGHT_STRATEGIES:
        support_strategy = "log1p"

    weight_matrix = WEIGHT_STRATEGIES[support_strategy](X)
    nonzero_mask = (X > 0)

    rng = np.random.default_rng(random_seed)
    if dropout_rate > 0:
        keep_mask = rng.uniform(0, 1, size=X.shape) > dropout_rate
        effective_mask = nonzero_mask & keep_mask
    else:
        effective_mask = nonzero_mask

    row_idx, col_idx = np.where(effective_mask)
    edge_weights = weight_matrix[effective_mask]

    # 归一化
    row_sum = np.zeros(n_cells, dtype=np.float32)
    np.add.at(row_sum, row_idx, edge_weights)
    row_sum = np.where(row_sum > 0, row_sum, 1.0)
    edge_weights_norm = edge_weights / row_sum[row_idx]

    gene_offset = n_cells
    src = np.concatenate([row_idx, col_idx + gene_offset])
    dst = np.concatenate([col_idx + gene_offset, row_idx])
    w = np.concatenate([edge_weights_norm, edge_weights_norm])

    return {
        "edge_index": np.stack([src, dst], axis=0).astype(np.int64),
        "edge_attr": w.astype(np.float32),
        "n_cells": n_cells,
        "n_genes": n_genes,
        "n_edges": len(src),
    }


# ---------------------------------------------------------------------------
# 支持度感知的掩码预测（类似 DOLORIS 的 Sparsity Masking）
# ---------------------------------------------------------------------------

def compute_support_mask_labels(
    X_expr: np.ndarray,
    threshold: float = 0.0,
) -> np.ndarray:
    """
    计算每个细胞-基因对的二元掩码标签（用于 Mask Loss）。

    M_{cg} = 1  如果 X_{cg} > threshold（即基因被激活）
             0  如果 X_{cg} <= threshold（即基因沉默）

    返回 (n_cells, n_hvg) 的二元掩码。
    """
    return (X_expr > threshold).astype(np.float32)


if __name__ == "__main__":
    # 简单测试
    n_cells, n_genes = 100, 50
    rng = np.random.default_rng(42)
    X_fake = rng.exponential(2, size=(n_cells, n_genes)).astype(np.float32)
    X_fake[X_fake < 0.5] = 0  # 制造稀疏性

    gene_names = np.array([f"Gene_{i}" for i in range(n_genes)])

    for strategy in ["log1p", "rank", "tfidf", "norm"]:
        src, dst, w, mask, nonzero = build_cell_gene_bipartite_graph(
            X_fake, gene_names, support_strategy=strategy, dropout_rate=0.0
        )
        print(f"{strategy}: {len(src)} edges, weight_range=[{w.min():.4f}, {w.max():.4f}]")

    print("build_cell_gene_graph.py test passed.")
