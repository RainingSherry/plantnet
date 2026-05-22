"""
models/gene_gat_encoder.py
===========================
基于图注意力网络（GAT）的基因编码器。

核心思想（源自 DOLORIS 的 GRN_conditional_network）：
  h_g^{(l+1)} = GAT^{(l)}( h_g^{(l)}, {h_{neighbor}}^{(l)} )
  h_g^{(0)}  = 可学习的基因嵌入（从随机初始化或预训练嵌入）

改进点：
  - 不再依赖 scGPT 预训练嵌入，改用可学习的基因嵌入 + 数据驱动的共表达图
  - 支持边权重（相关系数）作为注意力偏置
  - 支持三层 GAT 以捕获多跳基因关系
  - 支持细胞类型条件嵌入（注入批次/组织信息）
"""

from __future__ import annotations

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import LayerNorm, Linear
from typing import Optional, Tuple, List


# ---------------------------------------------------------------------------
# GAT Layer（参考 DOLORIS 的 MultiLayerGAT + 边权重支持）
# ---------------------------------------------------------------------------

class GATConv(nn.Module):
    """
    简化版 GATConv：支持边权重的多头图注意力卷积。

    参数
    ----
    in_channels : int  输入特征维度
    out_channels : int  输出特征维度（单头输出）
    heads : int         注意力头数
    concat : bool       多头输出是否拼接（True）还是平均（False）
    dropout : float     Dropout 概率
    add_self_loops : bool 是否在图中添加自环
    edge_dim : int      边特征维度（如果有边权重）
    """

    def __init__(
        self,
        in_channels: int,
        out_channels: int,
        heads: int = 4,
        concat: bool = True,
        dropout: float = 0.1,
        add_self_loops: bool = True,
        edge_dim: Optional[int] = None,
        negative_slope: float = 0.2,
    ):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.heads = heads
        self.concat = concat
        self.dropout = dropout
        self.negative_slope = negative_slope
        self.add_self_loops = add_self_loops
        self.edge_dim = edge_dim

        head_dim = out_channels // heads if concat else out_channels
        assert concat and (out_channels % heads == 0) or (not concat), \
            "out_channels must be divisible by heads"

        self.head_dim = head_dim

        # 注意力参数
        self.lin_src = Linear(in_channels, heads * head_dim, bias=False)
        self.lin_dst = Linear(in_channels, heads * head_dim, bias=False)
        self.lin_edge = Linear(edge_dim, heads * head_dim, bias=False) if edge_dim else None

        # 注意力偏置（leaky_relu 的负斜率）
        self.att = nn.Parameter(torch.Tensor(1, heads, head_dim))

        # 输出变换
        if concat:
            self.out_lin = Linear(heads * head_dim, out_channels)
        else:
            self.out_lin = Linear(heads * head_dim, out_channels)

        self.reset_parameters()

    def reset_parameters(self):
        nn.init.xavier_uniform_(self.lin_src.weight)
        nn.init.xavier_uniform_(self.lin_dst.weight)
        if self.lin_edge is not None:
            nn.init.xavier_uniform_(self.lin_edge.weight)
        nn.init.xavier_uniform_(self.att)
        nn.init.zeros_(self.out_lin.bias)

    def forward(
        self,
        x_src: torch.Tensor,       # [N, in_channels]
        x_dst: Optional[torch.Tensor],  # [N, in_channels] or None (= x_src)
        edge_index: torch.Tensor,  # [2, E]
        edge_attr: Optional[torch.Tensor] = None,  # [E, edge_dim]
    ) -> torch.Tensor:
        """
        参数
        ----
        x_src : [N, D_in]  源节点特征
        x_dst : [N, D_in]  目标节点特征（用于非对称消息传递，默认同 x_src）
        edge_index : [2, E]  边索引
        edge_attr : [E, D_edge]  边特征（可学习边权重）

        返回
        ----
        out : [N, out_channels * heads] (concat=True) 或 [N, out_channels] (concat=False)
        """
        N = x_src.size(0)
        x_dst_ = x_src if x_dst is None else x_dst

        # ---- 线性变换 ----
        h_src = self.lin_src(x_src).view(N, self.heads, self.head_dim)     # [N, H, D]
        h_dst = self.lin_dst(x_dst_).view(N, self.heads, self.head_dim)   # [N, H, D]

        # ---- 边权重融入注意力 ----
        if self.lin_edge is not None and edge_attr is not None:
            # edge_attr: [E] → [E, 1] 以兼容 Linear 层
            edge_in = edge_attr.unsqueeze(-1) if edge_attr.dim() == 1 else edge_attr
            edge_emb = self.lin_edge(edge_in).view(-1, self.heads, self.head_dim)  # [E, H, D]
            h_src_edge = h_src[edge_index[0]] + edge_emb                                   # [E, H, D]
        else:
            h_src_edge = h_src[edge_index[0]]                                   # [E, H, D]

        h_dst_edge = h_dst[edge_index[1]]                                       # [E, H, D]

        # ---- 计算注意力分数 ----
        # alpha = LeakyReLU( (h_src_edge || h_dst_edge) · att )
        alpha = torch.einsum("ehd,hd->eh", h_src_edge + h_dst_edge, self.att.squeeze(0))
        alpha = F.leaky_relu(alpha, negative_slope=self.negative_slope)
        alpha = F.softmax(alpha, dim=0)  # 按目标节点归一化

        # ---- 消息传递 ----
        alpha = alpha.unsqueeze(-1)                                     # [E, H, 1]
        if self.training and self.dropout > 0:
            alpha = F.dropout(alpha, p=self.dropout, training=True)

        # 消息 = 源节点嵌入 × 注意力权重
        messages = h_src_edge * alpha                                     # [E, H, D]
        E, H, D_h = messages.size()

        # 用 index_add_ 按目标节点聚合消息（更简洁）
        messages_flat = messages.view(E, H * D_h)   # [E, H*D]
        out = torch.zeros(N, H * D_h, device=x_src.device, dtype=x_src.dtype)
        out.index_add_(0, edge_index[1], messages_flat)              # [N, H*D]

        if not self.concat:
            out = out.mean(dim=1, keepdim=False)                       # [N, D]

        out = self.out_lin(out)                                       # [N, out_channels]
        return out


class MultiLayerGAT(nn.Module):
    """
    多层 GAT（继承自 DOLORIS 的 MultiLayerGAT），支持：
      - 边权重（相关系数）
      - 残差连接
      - LayerNorm
      - 可配置头数和层数

    典型配置：
      dim = [gene_init_dim, hidden_dim, hidden_dim]
      heads = [4, 4]
    """

    def __init__(
        self,
        dims: List[int],
        heads: List[int] = [4, 4],
        dropout: float = 0.1,
        use_residual: bool = True,
        use_norm: bool = True,
        edge_dim: int = 1,
    ):
        """
        参数
        ----
        dims : [D_0, D_1, D_2, ...]  每层的输入/输出维度列表
        heads : 每层的注意力头数（长度 = len(dims) - 1）
        dropout : Dropout 概率
        use_residual : 是否使用残差连接
        use_norm : 是否使用 LayerNorm
        edge_dim : 边特征维度（默认 1，表示边权重）
        """
        super().__init__()
        if len(heads) != len(dims) - 1:
            raise ValueError(f"len(heads)={len(heads)} must equal len(dims)-1={len(dims)-1}")

        self.layers = nn.ModuleList()
        self.norms = nn.ModuleList()
        self.heads = heads
        self.use_residual = use_residual
        self.use_norm = use_norm

        for i in range(len(dims) - 1):
            in_dim = dims[i]
            out_dim = dims[i + 1]
            n_heads = heads[i]
            concat = True
            layer = GATConv(
                in_channels=in_dim,
                out_channels=out_dim,
                heads=n_heads,
                concat=concat,
                dropout=dropout,
                edge_dim=edge_dim,
                add_self_loops=True,
            )
            self.layers.append(layer)
            if use_norm:
                self.norms.append(LayerNorm(out_dim))

    def forward(
        self,
        x: torch.Tensor,            # [N, D_0]  节点特征
        edge_index: torch.Tensor,   # [2, E]
        edge_attr: Optional[torch.Tensor] = None,  # [E, D_edge]
    ) -> torch.Tensor:
        """
        参数
        ----
        x : [N, D_in]  节点特征矩阵
        edge_index : [2, E]  边索引
        edge_attr : [E, edge_dim]  边特征

        返回
        ----
        x : [N, D_out]  最后一层输出
        """
        for i, layer in enumerate(self.layers):
            h_in = x
            x = layer(x, x, edge_index, edge_attr)  # [N, D_i+1]

            if self.use_residual and x.shape == h_in.shape:
                x = x + h_in

            if self.use_norm:
                x = self.norms[i](x)

            x = F.gelu(x)
            if self.training:
                x = F.dropout(x, p=0.1, training=True)

        return x


# ---------------------------------------------------------------------------
# Gene Embedding + GAT 编码器
# ---------------------------------------------------------------------------

class GeneGATEncoder(nn.Module):
    """
    基因图 GAT 编码器。

    核心流程：
      1. 可学习的基因嵌入初始化（可加载预训练向量）
      2. 多层 GAT 在基因图上传播信息
      3. 输出每个基因的上下文感知嵌入 h_g

    数学表达：
      h_g^{(0)} = GeneEmbedding(g)          ← 可学习嵌入
      h_g^{(L)} = GATLayer^{(L)}(h_g^{(L-1)}, A)  ← L 层图注意力
      输出 = h_g^{(L)}
    """

    def __init__(
        self,
        n_genes: int,
        gene_dim: int = 64,
        hidden_dim: int = 256,
        n_layers: int = 2,
        heads: List[int] = [4, 4],
        dropout: float = 0.1,
        pretrained_gene_emb: Optional[torch.Tensor] = None,
        freeze_embedding: bool = False,
    ):
        """
        参数
        ----
        n_genes : 基因数量（HVG 数量）
        gene_dim : 基因嵌入维度（输入 GAT 维度）
        hidden_dim : GAT 隐藏层维度
        n_layers : GAT 层数
        heads : 每层注意力头数（长度 = n_layers）
        dropout : Dropout 概率
        pretrained_gene_emb : (n_genes, gene_dim) 预训练嵌入（可选）
        freeze_embedding : 是否冻结基因嵌入
        """
        super().__init__()
        self.n_genes = n_genes
        self.gene_dim = gene_dim
        self.hidden_dim = hidden_dim

        # ---- 基因嵌入层 ----
        self.gene_embedding = nn.Embedding(n_genes, gene_dim)
        if pretrained_gene_emb is not None:
            assert pretrained_gene_emb.shape == (n_genes, gene_dim), \
                f"Pretrained embedding shape {pretrained_gene_emb.shape} != ({n_genes}, {gene_dim})"
            self.gene_embedding.weight.data.copy_(pretrained_gene_emb)
            print(f"  [GeneGATEncoder] Loaded pretrained gene embeddings: {n_genes} × {gene_dim}")
        if freeze_embedding:
            self.gene_embedding.weight.requires_grad = False

        # ---- LayerNorm + 位置编码 ----
        self.ln_gene = LayerNorm(gene_dim)

        # ---- 多层 GAT ----
        if len(heads) != n_layers:
            heads = [heads[0]] * n_layers if isinstance(heads, list) and len(heads) == 1 else [4] * n_layers

        # GAT dims: input=gene_dim, hidden layers=hidden_dim, output=gene_dim
        # Output = [n_genes, gene_dim] so it matches embed_projector input dim
        gat_dims = [gene_dim] + [hidden_dim] * max(0, n_layers - 1) + [gene_dim]
        self.gnn = MultiLayerGAT(
            dims=gat_dims,
            heads=heads,
            dropout=dropout,
            use_residual=True,
            use_norm=True,
            edge_dim=1,
        )

        self.output_dim = gene_dim

    def forward(
        self,
        gene_ids: torch.Tensor,       # [N_genes] or [B, N_genes]  基因索引（通常用全部基因）
        edge_index: torch.Tensor,     # [2, E_gene]
        edge_weight: Optional[torch.Tensor] = None,  # [E_gene] 边权重
        cell_type_emb: Optional[torch.Tensor] = None,  # [B, cell_type_emb_dim] 细胞类型嵌入
    ) -> torch.Tensor:
        """
        参数
        ----
        gene_ids : 基因索引，默认使用全部基因（0 ~ n_genes-1）
        edge_index : 基因图的边索引 [2, E]
        edge_weight : 边权重 [E]（相关系数，归一化到 0.5~1.0）
        cell_type_emb : 细胞类型嵌入 [B, D]（可选，用于注入细胞类型条件）

        返回
        ----
        h_g : [n_genes, hidden_dim]  基因上下文嵌入
        """
        B = 1
        if gene_ids.dim() == 2:
            B = gene_ids.size(0)
            gene_ids = gene_ids.view(-1)

        # 嵌入基因
        h_g = self.gene_embedding(gene_ids)  # [B*n_genes, gene_dim] 或 [n_genes, gene_dim]

        if B > 1:
            h_g = h_g.view(B, self.n_genes, self.gene_dim)

        # LayerNorm
        h_g = self.ln_gene(h_g)

        # 多层 GAT（批处理模式：每个样本用同一个基因图）
        if B > 1:
            all_h = []
            for b in range(B):
                h_b = self.gnn(h_g[b], edge_index, edge_weight)
                all_h.append(h_b)
            h_g = torch.stack(all_h, dim=0)  # [B, n_genes, hidden_dim]
        else:
            h_g = self.gnn(h_g, edge_index, edge_weight)  # [n_genes, output_dim]
            h_g = h_g.unsqueeze(0)                              # [1, n_genes, output_dim]

        return h_g  # [B, n_genes, hidden_dim] 或 [1, n_genes, hidden_dim]


# ---------------------------------------------------------------------------
# 基因嵌入的稀疏读取（用于 SupportPooling）
# ---------------------------------------------------------------------------

class GeneLookupWithMask(nn.Module):
    """
    给定每个细胞的支撑集（expressed gene indices），从基因嵌入中索引对应的嵌入向量。
    输出 (n_cells, max_support_size, gene_emb_dim) 的稀疏嵌入矩阵。

    与 DOLORIS 的 knockout 索引类似，但这里我们索引的是"正向表达"的基因。
    """

    def __init__(
        self,
        n_genes: int,
        gene_emb_dim: int = 64,
        max_support_size: Optional[int] = None,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.gene_emb_dim = gene_emb_dim
        self.max_support_size = max_support_size

        # 可学习的基因嵌入（与 GeneGATEncoder 共享）
        self.gene_embedding = nn.Embedding(n_genes, gene_emb_dim, padding_idx=-1)

    def forward(
        self,
        gene_indices: torch.Tensor,   # [B, max_support]  每个细胞的表达基因索引
        mask: Optional[torch.Tensor] = None,  # [B, max_support]  有效位置掩码
    ) -> torch.Tensor:
        """
        参数
        ----
        gene_indices : [B, max_support]  每个细胞的表达基因索引（不足处用 -1 填充）
        mask : [B, max_support]  有效位掩码（1=有效，0=填充）

        返回
        ----
        gene_emb_batch : [B, max_support, gene_emb_dim]  每个细胞的支撑基因嵌入
        """
        B, S = gene_indices.shape
        gene_emb = self.gene_embedding(gene_indices)  # [B, S, D]

        if mask is not None:
            # 将无效位置置零
            gene_emb = gene_emb * mask.unsqueeze(-1)

        return gene_emb


if __name__ == "__main__":
    # 简单测试
    import torch

    n_genes = 500
    gene_dim = 64
    hidden_dim = 128
    n_cells = 32

    # 模拟基因图
    n_edges = 3000
    edge_index = torch.randint(0, n_genes, (2, n_edges))
    edge_weight = torch.rand(n_edges).clamp_(min=0.5, max=1.0)

    # 模拟支撑集（每个细胞表达 ~100 个基因）
    max_support = 200
    support_indices = torch.randint(0, n_genes, (n_cells, max_support))
    support_mask = (support_indices >= 0).float()

    # 测试 GeneGATEncoder
    encoder = GeneGATEncoder(
        n_genes=n_genes,
        gene_dim=gene_dim,
        hidden_dim=hidden_dim,
        n_layers=2,
        heads=[4, 4],
    )

    gene_emb = encoder(
        gene_ids=torch.arange(n_genes),
        edge_index=edge_index,
        edge_weight=edge_weight,
    )
    print(f"Gene embeddings: {gene_emb.shape}")  # [1, n_genes, hidden_dim]

    # 测试 GeneLookupWithMask
    lookup = GeneLookupWithMask(n_genes=n_genes, gene_emb_dim=gene_dim)
    support_emb = lookup(support_indices, support_mask)
    print(f"Support embeddings: {support_emb.shape}")  # [n_cells, max_support, gene_dim]

    print("gene_gat_encoder.py test passed.")
