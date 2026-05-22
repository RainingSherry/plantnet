# -*- coding: utf-8 -*-
"""
SupportMaskNet: 预测每个细胞-基因是否表达的掩码网络

灵感来自 DOLORIS 的稀疏度掩码策略，但应用于聚类任务：
- 输入：细胞的归一化表达向量
- 输出：每个基因被激活的概率 p(gene_active)
- 目标：学习"表达支撑结构"（observed support），而非判断零值的真假

关键设计：
1. 共享参数 encoder 学习细胞嵌入
2. 独立 mask head 预测基因激活概率
3. 使用 BCE 损失，目标为 M = 1(X > 0)
"""

import torch
import torch.nn as nn
import torch.nn.functional as F


class SupportMaskNet(nn.Module):
    """
    Support Mask Network for predicting gene activation probabilities.

    Architecture:
        X (n_cells, n_genes)
            │
            ▼
        Dropout
            │
            ▼
        Linear(n_genes → hidden1) + LayerNorm + Mish
            │
            ▼
        Linear(hidden1 → hidden2) + LayerNorm + Mish
            │
            ▼
        Linear(hidden2 → hidden3)
            │
            ├──────────────────────┐
            ▼                      ▼
        gene_activation_head    cell_embedding_head
        (n_cells, n_genes)     (n_cells, latent_dim)

    Losses:
        - BCE: 预测每个基因是否被激活
        - Weighted BCE: 高表达基因权重更高
    """

    def __init__(
        self,
        num_genes: int,
        hidden_dims: list = None,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256, 128]

        self.num_genes = num_genes

        # Shared encoder
        encoder_layers = []
        in_dim = num_genes
        for hidden_dim in hidden_dims:
            encoder_layers.append(nn.Dropout(p=dropout))
            encoder_layers.append(nn.Linear(in_dim, hidden_dim))
            encoder_layers.append(nn.LayerNorm(hidden_dim))
            encoder_layers.append(nn.Mish(inplace=True))
            in_dim = hidden_dim

        self.encoder = nn.Sequential(*encoder_layers)
        self.encoder_output_dim = hidden_dims[-1]

        # Gene activation prediction head
        self.gene_activation_head = nn.Linear(self.encoder_output_dim, num_genes)

        # Cell embedding head (for downstream clustering)
        self.cell_embedding_head = nn.Linear(self.encoder_output_dim, hidden_dims[-1])

    def forward(self, x: torch.Tensor) -> dict:
        """
        Forward pass.

        Args:
            x: Input gene expression (n_cells, n_genes)

        Returns:
            dict with keys:
                - gene_activation_prob: probability of gene activation (n_cells, n_genes)
                - cell_embedding: cell embedding (n_cells, latent_dim)
                - encoder_output: raw encoder output (n_cells, hidden_dim)
        """
        h = self.encoder(x)

        gene_activation_prob = torch.sigmoid(self.gene_activation_head(h))
        cell_embedding = self.cell_embedding_head(h)

        return {
            'gene_activation_prob': gene_activation_prob,
            'cell_embedding': cell_embedding,
            'encoder_output': h,
        }

    def get_gene_activation_loss(
        self,
        x: torch.Tensor,
        mask: torch.Tensor = None,
        pos_weight: float = 1.0,
    ) -> tuple:
        """
        Compute BCE loss for gene activation prediction.

        Args:
            x: Input gene expression (n_cells, n_genes)
            mask: Binary mask (n_cells, n_genes), 1=expressed, 0=zero
                  If None, computed as (x > 0).float()
            pos_weight: Weight for positive class (expressed genes)

        Returns:
            (loss, gene_activation_prob)
        """
        if mask is None:
            mask = (x > 0).float()

        output = self.forward(x)
        gene_activation_prob = output['gene_activation_prob']

        # BCE loss with optional positive weight
        # Higher weight for expressed genes (they carry more biological signal)
        loss = F.binary_cross_entropy(
            gene_activation_prob,
            mask,
            reduction='mean',
            pos_weight=torch.tensor(pos_weight).to(x.device) if pos_weight != 1.0 else None,
        )

        return loss, gene_activation_prob

    def get_support_aware_embedding(
        self,
        x: torch.Tensor,
        threshold: float = 0.5,
    ) -> torch.Tensor:
        """
        Get cell embedding weighted by predicted gene activation.

        Args:
            x: Input gene expression (n_cells, n_genes)
            threshold: Threshold for gene activation

        Returns:
            Support-weighted cell embedding (n_cells, latent_dim)
        """
        with torch.no_grad():
            output = self.forward(x)
            gene_activation_prob = output['gene_activation_prob']

            # Create support mask
            support_mask = (gene_activation_prob > threshold).float()

            # Apply mask to input
            masked_x = x * support_mask

        # Get embedding from masked input
        h = self.encoder(masked_x)
        embedding = self.cell_embedding_head(h)

        return embedding


class SupportMaskNetWithAttention(nn.Module):
    """
    Enhanced SupportMaskNet with cross-attention mechanism for gene-gene relationships.
    """

    def __init__(
        self,
        num_genes: int,
        hidden_dims: list = None,
        num_heads: int = 4,
        dropout: float = 0.1,
    ):
        super().__init__()
        if hidden_dims is None:
            hidden_dims = [512, 256, 128]

        self.num_genes = num_genes

        # Gene embedding layer
        self.gene_embedding = nn.Linear(1, hidden_dims[0] // 4)

        # Encoder
        encoder_layers = []
        in_dim = num_genes * (hidden_dims[0] // 4) + hidden_dims[0]
        for hidden_dim in hidden_dims[1:]:
            encoder_layers.append(nn.Dropout(p=dropout))
            encoder_layers.append(nn.Linear(in_dim, hidden_dim))
            encoder_layers.append(nn.LayerNorm(hidden_dim))
            encoder_layers.append(nn.Mish(inplace=True))
            in_dim = hidden_dim

        self.encoder = nn.Sequential(*encoder_layers)
        self.encoder_output_dim = hidden_dims[-1]

        # Multi-head attention for gene interactions
        self.gene_attention = nn.MultiheadAttention(
            embed_dim=hidden_dims[0],
            num_heads=num_heads,
            dropout=dropout,
            batch_first=True,
        )

        # Heads
        self.gene_activation_head = nn.Linear(self.encoder_output_dim, num_genes)
        self.cell_embedding_head = nn.Linear(self.encoder_output_dim, hidden_dims[-1])

    def forward(self, x: torch.Tensor) -> dict:
        """
        Forward pass with attention mechanism.

        Args:
            x: Input gene expression (n_cells, n_genes)

        Returns:
            dict with keys:
                - gene_activation_prob: probability of gene activation (n_cells, n_genes)
                - cell_embedding: cell embedding (n_cells, latent_dim)
        """
        batch_size, num_genes = x.shape

        # Gene-level embedding
        x_expanded = x.unsqueeze(-1)  # (n_cells, n_genes, 1)
        gene_emb = self.gene_embedding(x_expanded)  # (n_cells, n_genes, gene_emb_dim)

        # Self-attention across genes
        gene_emb_transposed = gene_emb.transpose(0, 1)  # (n_genes, n_cells, gene_emb_dim)
        attn_output, _ = self.gene_attention(
            gene_emb_transposed, gene_emb_transposed, gene_emb_transposed
        )
        attn_output = attn_output.transpose(0, 1)  # (n_cells, n_genes, gene_emb_dim)

        # Global cell representation
        cell_repr = attn_output.mean(dim=1)  # (n_cells, gene_emb_dim)

        # Encode cell representation and concatenate with attention output
        cell_repr_expanded = cell_repr.unsqueeze(1).expand(-1, num_genes, -1)
        combined = torch.cat([attn_output, cell_repr_expanded], dim=-1)

        # Flatten and encode
        combined_flat = combined.reshape(batch_size, -1)
        h = self.encoder(combined_flat)

        # Predictions
        gene_activation_prob = torch.sigmoid(self.gene_activation_head(h))
        cell_embedding = self.cell_embedding_head(h)

        return {
            'gene_activation_prob': gene_activation_prob,
            'cell_embedding': cell_embedding,
        }
