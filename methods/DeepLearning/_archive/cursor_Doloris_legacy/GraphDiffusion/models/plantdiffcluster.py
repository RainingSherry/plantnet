"""
models/plantdiffcluster.py
==========================
PlantDiffCluster 主模型类。

整合所有子模块：
  1. GeneGATEncoder   — 基因图注意力编码
  2. SupportPooling   — 稀疏感知基因→细胞聚合
  3. MaskDiffusionRefiner  — DDPM/DDIM 去噪精炼
  4. ClusterHead      — GMM / DEC / Contrastive 聚类头
  5. GeneDecoder       — 可选的重构解码器

前向流程：
  X → GeneGATEncoder → gene_emb
  gene_emb + support → SupportPooling → cell_emb_raw
  cell_emb_raw → MaskDiffusionRefiner → cell_emb_final
  cell_emb_final → ClusterHead → cluster_probs + loss
"""

from __future__ import annotations

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Module
from typing import Optional, Dict, List, Any
import numpy as np
import json
from pathlib import Path

from .gene_gat_encoder import GeneGATEncoder
from .support_pooling import SupportPoolingFactory
from .mask_diffusion_refiner import MaskDiffusionRefiner
from .cluster_head import ClusterHeadFactory


# ---------------------------------------------------------------------------
# 主模型
# ---------------------------------------------------------------------------

class PlantDiffCluster(Module):
    """
    PlantDiffCluster：用于植物单细胞 RNA-seq 聚类的扩散图神经网络。

    参数
    ----
    n_genes : int  高变基因数量
    gene_names : List[str]  基因名列表
    graph_dict : dict  基因图数据，包含：
        - edge_index: [2, n_edges] 边索引
        - edge_weight: [n_edges] 边权重（可选）
    config : dict  模型配置，包含所有超参数
    """

    def __init__(
        self,
        n_genes: int,
        gene_names: List[str],
        graph_dict: dict,
        config: dict,
    ):
        super().__init__()
        self.n_genes = n_genes
        self.gene_names = gene_names
        self.graph_dict = graph_dict
        self.config = config

        gene_dim = config.get("gene_dim", 64)
        hidden_dim = config.get("hidden_dim", 256)
        embed_dim = config.get("embed_dim", 128)
        time_embed_dim = config.get("time_embed_dim", 128)
        n_layers = config.get("n_layers", 2)
        heads = config.get("heads", [4] * n_layers)
        pooling_strategy = config.get("pooling_strategy", "attention")
        pooling_topk = config.get("pooling_topk", 50)
        n_clusters = config.get("n_clusters", 15)
        cluster_strategy = config.get("cluster_strategy", "gmm")
        use_diffusion = config.get("use_diffusion", True)
        use_mask_predictor = config.get("use_mask_predictor", True)
        num_timesteps = config.get("num_timesteps", 500)
        ddim_steps = config.get("ddim_steps", 20)
        beta_schedule = config.get("beta_schedule", "cosine")
        refiner_depth = config.get("refiner_depth", 3)
        refiner_hidden_dim = config.get("refiner_hidden_dim", hidden_dim)
        use_decoder = config.get("use_decoder", True)
        decoder_hidden_dim = config.get("decoder_hidden_dim", 512)
        dropout = config.get("dropout", 0.1)

        self.embed_dim = embed_dim
        self.n_clusters = n_clusters
        self.use_diffusion = use_diffusion

        # ---- 1. 基因图 GAT 编码器 ----
        self.gene_encoder = GeneGATEncoder(
            n_genes=n_genes,
            gene_dim=gene_dim,
            hidden_dim=hidden_dim,
            n_layers=n_layers,
            heads=heads,
            dropout=dropout,
        )

        # ---- 2. 支撑集聚合器 ----
        self.pooling = SupportPoolingFactory.create(
            pooling_strategy=pooling_strategy,
            gene_emb_dim=gene_dim,
            cell_emb_dim=0,
            hidden_dim=hidden_dim,
            n_heads=4,
            dropout=dropout,
            topk_k=pooling_topk,
        )

        # ---- 3. 嵌入投影层 ----
        self.gene_to_embed = nn.Linear(gene_dim, embed_dim)

        # ---- 4. 扩散精炼器 ----
        if use_diffusion:
            refiner_depth = config.get("refiner_depth", 3)
            self.diffusion_refiner = MaskDiffusionRefiner(
                embed_dim=embed_dim,
                n_genes=n_genes,
                hidden_dim=refiner_hidden_dim,
                time_embed_dim=time_embed_dim,
                cond_embed_dim=0,
                refiner_depth=refiner_depth,
                dropout=dropout,
                num_timesteps=num_timesteps,
                beta_schedule=beta_schedule,
                use_mask_predictor=use_mask_predictor,
                use_layer_norm=True,
            )
        else:
            self.diffusion_refiner = None

        # ---- 5. 聚类头 ----
        self.cluster_head = ClusterHeadFactory.create(
            strategy=cluster_strategy,
            embed_dim=embed_dim,
            n_clusters=config.get("cell_type_num", n_clusters),
            hidden_dim=hidden_dim,
        )

        # ---- 6. 解码器（可选） ----
        if use_decoder:
            self.decoder = nn.Sequential(
                nn.Linear(embed_dim, decoder_hidden_dim),
                nn.ReLU(),
                nn.Dropout(dropout),
                nn.Linear(decoder_hidden_dim, n_genes),
            )
        else:
            self.decoder = None

        # ---- 注册基因图到 buffer ----
        self._register_graph_buffers()

    def _register_graph_buffers(self):
        """将基因图数据注册为 buffer（自动跟随 .to() 迁移到正确设备）。"""
        edge_index = self.graph_dict["edge_index"]
        edge_weight = self.graph_dict.get("edge_weight")

        self.register_buffer(
            "gene_edge_index",
            torch.tensor(edge_index, dtype=torch.long, device=next(self.parameters()).device)
        )
        if edge_weight is not None and len(edge_weight) > 0:
            self.register_buffer(
                "gene_edge_weight",
                torch.tensor(edge_weight, dtype=torch.float32, device=next(self.parameters()).device)
            )
        else:
            self.register_buffer("gene_edge_weight", None)

    def initialize_clusters(self, embeddings: torch.Tensor):
        """用训练数据的嵌入初始化聚类中心。"""
        self.cluster_head.initialize(embeddings.detach().cpu().numpy(), method="kmeans")

    def forward(
        self,
        X: torch.Tensor,               # [B, n_genes]  原始基因表达
        cell_type: Optional[torch.Tensor] = None,  # [B]  细胞类型标签
        support_weight: Optional[torch.Tensor] = None,  # [B, max_support]  支撑集权重
        support_mask: Optional[torch.Tensor] = None,   # [B, max_support]  有效位掩码
        support_idx: Optional[torch.Tensor] = None,      # [B, max_support]  支撑基因索引
        t: Optional[torch.Tensor] = None,               # 扩散时间步（训练时）
        return_refined: bool = False,                   # 是否返回精炼嵌入
    ) -> Dict[str, torch.Tensor]:
        """
        前向传播。

        流程：
          1. 基因 GAT 编码
          2. 支撑集聚合 → cell_emb_raw
          3. 扩散精炼 → cell_emb_refined
          4. 聚类头 → cluster_probs + cluster_loss
          5. 解码器 → reconstruction（可选）
        """
        device = X.device
        B = X.size(0)

        # ---- Step 1: 基因图 GAT 编码 ----
        gene_emb = self.gene_encoder(
            gene_ids=torch.arange(self.n_genes, device=device),
            edge_index=self.gene_edge_index,
            edge_weight=self.gene_edge_weight,
        )  # [1, n_genes, hidden_dim]

        gene_emb = gene_emb.squeeze(0)  # [n_genes, hidden_dim]

        # ---- Step 2: 支撑集聚合 ----
        # Use support_idx to gather relevant gene embeddings for each cell
        if support_idx is not None:
            # support_idx: [B, S], support_weight: [B, S]
            # gene_emb: [n_genes, hidden_dim] -> [B, S, hidden_dim]
            gene_emb_batch = gene_emb.unsqueeze(0).expand(B, -1, -1)  # [B, n_genes, hidden_dim]
            # Ensure all tensors are on the same device
            if support_idx.device != gene_emb.device:
                support_idx = support_idx.to(gene_emb.device)
            if support_weight.device != gene_emb.device:
                support_weight = support_weight.to(gene_emb.device)
            if support_mask is not None and support_mask.device != gene_emb.device:
                support_mask = support_mask.to(gene_emb.device)
            # Handle -1 indices by replacing with 0 (will be masked anyway)
            safe_support_idx = support_idx.clone()
            safe_support_idx[safe_support_idx < 0] = 0
            gathered_gene_emb = torch.gather(
                gene_emb_batch, 1, safe_support_idx.unsqueeze(-1).expand(-1, -1, gene_emb.shape[-1])
            )  # [B, S, hidden_dim]
            # Mask out invalid positions
            valid_mask = (support_idx >= 0).unsqueeze(-1).float()  # [B, S, 1]
            gathered_gene_emb = gathered_gene_emb * valid_mask
            # Pool over the support set dimension
            cell_emb_raw = self.pooling(
                gene_emb=gathered_gene_emb,
                support_weight=support_weight,
                mask=support_mask,
            )  # [B, gene_dim]
        else:
            if support_weight.device != gene_emb.device:
                support_weight = support_weight.to(gene_emb.device)
            if support_mask is not None and support_mask.device != gene_emb.device:
                support_mask = support_mask.to(gene_emb.device)
            cell_emb_raw = self.pooling(
                gene_emb=gene_emb.unsqueeze(0).expand(B, -1, -1),
                support_weight=support_weight,
                mask=support_mask,
            )  # [B, gene_dim]


        # ---- Step 3: 投影到嵌入空间 ----
        cell_emb = self.gene_to_embed(cell_emb_raw)  # [B, embed_dim]

        # ---- Step 4: 扩散精炼 ----
        cell_emb_refined = cell_emb
        losses = {}

        if self.diffusion_refiner is not None:
            if t is not None:
                # Training mode: compute losses
                diff_output = self.diffusion_refiner.training_losses(
                    z0=cell_emb,
                    t=t,
                    cond_emb=None,
                    true_mask=None,
                )
                losses["mse"] = diff_output["mse"]
                losses["loss"] = diff_output["loss"]
                # Use the predicted clean embedding for clustering
                # The refiner predicts x0; we use it for downstream tasks
                cell_emb_refined = cell_emb  # stay with pre-refinement for loss computation
            else:
                # Inference mode: sample refined embedding
                refined = self.diffusion_refiner.sample(
                    cell_emb.shape,
                    cond_emb=None,
                    n_steps=20,
                    use_ddim=True,
                )
                cell_emb_refined = refined

        # ---- Step 5: 聚类 ----
        cluster_loss = torch.tensor(0.0, device=device)
        cluster_probs = None
        cluster_assign = None

        if self.cluster_head is not None:
            # cluster_head.forward expects z: [B, hidden_dim]
            cluster_assign, cell_emb_proj = self.cluster_head(cell_emb_refined)
            cluster_loss = torch.tensor(0.0, device=device)
            cluster_probs = None

        losses["cluster"] = cluster_loss

        # ---- Step 6: 重构 ----
        recon_loss = torch.tensor(0.0, device=device)
        reconstruction = None

        if self.decoder is not None:
            reconstruction = self.decoder(cell_emb_refined)  # [B, n_genes]
            if X.shape == reconstruction.shape:
                recon_loss = F.mse_loss(reconstruction, X)
            losses["recon"] = recon_loss

        losses["loss"] = (
            losses.get("recon", torch.tensor(0.0, device=device))
            + 0.1 * losses.get("mask", torch.tensor(0.0, device=device))
            + 0.1 * losses.get("cluster", torch.tensor(0.0, device=device))
        )

        return {
            "cell_z": cell_emb_refined,
            "cell_z_raw": cell_emb,
            "reconstruction": reconstruction,
            "cluster_probs": cluster_probs,
            "cluster_assign": cluster_assign,
            "losses": losses,
            "gene_emb": gene_emb,
        }

    @torch.no_grad()
    def encode(self, X: torch.Tensor, **kwargs) -> torch.Tensor:
        """推理时的编码接口。"""
        output = self.forward(X, t=None, **kwargs)
        return output["cell_z"]

    @torch.no_grad()
    def predict(self, X: torch.Tensor, **kwargs) -> torch.Tensor:
        """推理时的预测接口。"""
        output = self.forward(X, t=None, **kwargs)
        return output["cluster_assign"]


# ---------------------------------------------------------------------------
# Checkpoint utilities
# ---------------------------------------------------------------------------

def save_checkpoint(
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    epoch: int,
    path: str,
    best_metrics: Optional[dict] = None,
    embeddings: Optional[np.ndarray] = None,
    pred_labels: Optional[np.ndarray] = None,
    scheduler=None,
):
    """保存模型检查点。"""
    ckpt = {
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "epoch": epoch,
        "best_metrics": best_metrics or {},
    }
    if scheduler is not None:
        ckpt["scheduler_state_dict"] = scheduler.state_dict()
    if embeddings is not None:
        ckpt["embeddings"] = embeddings
    if pred_labels is not None:
        ckpt["pred_labels"] = pred_labels

    torch.save(ckpt, path)

    # 保存 config
    config_path = str(Path(path).parent / "config.json")
    with open(config_path, "w") as f:
        json.dump(model.config if hasattr(model, "config") else {}, f, indent=2, default=str)


def load_checkpoint(model: nn.Module, path: str, device: str = "cuda") -> dict:
    """加载模型检查点。"""
    ckpt = torch.load(path, map_location=device, weights_only=False)
    model.load_state_dict(ckpt["model_state_dict"])
    return ckpt
