"""
models/cluster_head.py
=======================
聚类头模块。

支持三种聚类策略：
  1. GMM（高斯混合模型）—— 基于 PhytoCluster 的 VAE+GMM 联合优化
  2. Contrastive Clustering —— 对比学习风格（SCAN / DCC / SimCSE）
  3. DEC-style Cluster Assignment —— DEC（Deep Embedded Clustering）风格软分配

结合 DOLORIS 的 DDIB 隐空间聚类思想：
  - 在扩散精炼后的隐空间中直接做聚类
  - 使用 GMM 建模每个聚类的分布
  - 支持软分配（soft assignment）和硬分配（hard assignment）
  - 提供聚类可解释性：每个 cluster 可追溯到 top contributing genes

损失函数组合：
  L_total = L_reconstruction + λ_mask · L_mask + λ_cluster · L_cluster

其中 L_cluster 可以是：
  - GMM 变分目标（ELBO）
  - 对比损失（NT-Xent）
  - DEC 交叉熵软分配
"""

from __future__ import annotations

import math
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import Module
from typing import Optional, Tuple, Literal, List
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
import numpy as np


# ---------------------------------------------------------------------------
# 辅助函数
# ---------------------------------------------------------------------------

def compute_snr(z: torch.Tensor) -> float:
    """计算隐空间信噪比（Signal-to-Noise Ratio）。"""
    between_var = z.var(dim=0).mean().item()
    within_var = ((z - z.mean(dim=0)) ** 2).mean().item()
    return between_var / (within_var + 1e-8)


# ---------------------------------------------------------------------------
# 1. GMM Cluster Head（参考 PhytoCluster）
# ---------------------------------------------------------------------------

class GMMClusterHead(Module):
    """
    GMM 聚类头（参考 PhytoCluster 的 VAE+GMM 联合优化）。

    核心思想：
      - 在隐空间 z 上建模 K 个高斯分布
      - 使用 EM 算法初始化 GMM 参数（以 K-Means 结果为起点）
      - 变分推断：q(k|z) = p(k|z) · π_k / Σ_j p(z|j)·π_j
      - 目标：最大化 ELBO，同时最小化聚类熵

    数学表达：
      f_gmm(h; θ) = Σ_{k=1}^K π_k · φ(h; μ_k, Σ_k)

    其中 φ 是高斯概率密度函数，π_k 是混合权重。
    """

    def __init__(
        self,
        embed_dim: int,
        n_clusters: int,
        hidden_dim: Optional[int] = None,
        use_diagonal_cov: bool = True,
        eps: float = 1e-8,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.n_clusters = n_clusters
        self.eps = eps

        # 可学习的聚类参数（用 MLE 初始化）
        self.register_buffer("cluster_centers", torch.zeros(n_clusters, embed_dim))
        self.register_buffer("cluster_cov", torch.eye(embed_dim).unsqueeze(0).repeat(n_clusters, 1, 1))
        self.register_buffer("cluster_weights", torch.ones(n_clusters) / n_clusters)

        # 可学习的非线性变换（将 z 投影到更适合聚类的空间）
        if hidden_dim:
            self.projection = nn.Sequential(
                nn.Linear(embed_dim, hidden_dim),
                nn.BatchNorm1d(hidden_dim),
                nn.ReLU(),
                nn.Dropout(0.1),
                nn.Linear(hidden_dim, embed_dim),
            )
        else:
            self.projection = nn.Identity()

    def initialize(self, z: torch.Tensor, method: str = "kmeans"):
        """
        用 z 的初始化来设置 GMM 参数。

        参数
        ----
        z : [N, D]  隐嵌入（通常是全部训练数据的嵌入）
        method : str  初始化方法
          - "kmeans": K-Means 聚类中心作为 GMM 均值
          - "random": 随机初始化
        """
        with torch.no_grad():
            z_proj = self.projection(z)

            if method == "kmeans":
                kmeans = KMeans(n_clusters=self.n_clusters, n_init=20, random_state=42)
                labels = kmeans.fit_predict(z_proj.cpu().numpy())
                centers = torch.from_numpy(kmeans.cluster_centers_).float()

                # 初始化聚类中心
                self.cluster_centers.copy_(centers.to(z.device))

                # 初始化协方差（用每个簇的样本协方差）
                for k in range(self.n_clusters):
                    mask = torch.from_numpy(labels == k).float()
                    if mask.sum() > 1:
                        cluster_z = z_proj[mask.bool()]
                        cov = torch.cov(cluster_z.T) + self.eps * torch.eye(self.embed_dim, device=z.device)
                        self.cluster_cov[k] = cov

                # 初始化权重
                for k in range(self.n_clusters):
                    self.cluster_weights[k] = mask.sum() / len(labels)

            elif method == "random":
                idx = torch.randperm(len(z_proj))[: self.n_clusters]
                self.cluster_centers.copy_(z_proj[idx])

    def forward(
        self,
        z: torch.Tensor,
        return_probs: bool = False,
    ) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        """
        前向传播：计算每个样本属于各聚类的概率。

        参数
        ----
        z : [B, D]  隐嵌入
        return_probs : bool  是否返回概率（True = 软分配，False = 硬分配）

        返回
        ----
        assignments : [B,]  聚类标签（硬分配）或
        assignments, probs : ([B,], [B, K])  （软分配）
        """
        z_proj = self.projection(z)

        # 计算每个簇的概率密度
        log_probs = self._gaussian_log_prob(z_proj)  # [B, K]
        probs = F.softmax(log_probs, dim=-1)          # [B, K]

        # 软分配（贝叶斯）
        assignments = (probs * self.cluster_weights[None, :]).sum(dim=-1)  # [B]
        assignments = assignments / (assignments.sum(dim=-1, keepdim=True) + self.eps)  # 归一化

        if return_probs:
            hard_assign = probs.argmax(dim=-1)  # [B]
            return hard_assign, probs
        else:
            hard_assign = probs.argmax(dim=-1)  # [B]
            return hard_assign, None

    def _gaussian_log_prob(self, z: torch.Tensor) -> torch.Tensor:
        """计算每个样本属于每个高斯分量的对数概率密度。"""
        B, D = z.shape
        K = self.n_clusters

        z_expanded = z.unsqueeze(1)          # [B, 1, D]
        centers_expanded = self.cluster_centers.unsqueeze(0)  # [1, K, D]

        diff = z_expanded - centers_expanded  # [B, K, D]

        # 对角协方差
        cov_diag = self.cluster_cov.diagonal(dim1=1, dim2=2) + self.eps  # [K, D]
        log_det_cov = cov_diag.log().sum(dim=-1)  # [K]
        cov_inv_diff = diff / cov_diag.unsqueeze(0)  # [B, K, D]
        mahal = (cov_inv_diff * diff).sum(dim=-1)  # [B, K]

        log_probs = -0.5 * (mahal + log_det_cov + D * math.log(2 * math.pi))

        # 混合权重
        log_probs = log_probs + self.cluster_weights.log().unsqueeze(0)

        return log_probs

    def compute_clustering_loss(
        self,
        z: torch.Tensor,
        alpha: float = 1.0,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        计算聚类变分目标（参考 PhytoCluster）。

        L_cluster = α · KL(q(k|z) || uniform) + entropy(q(k|z))

        即：鼓励聚类分布均匀，同时保持一定熵（不所有样本聚到同一簇）。

        参数
        ----
        z : [B, D]
        alpha : float  KL 项权重

        返回
        ----
        loss : 聚类损失
        assignment_probs : [B, K]  各样本属于各簇的概率
        """
        z_proj = self.projection(z)
        log_probs = self._gaussian_log_prob(z_proj)  # [B, K]
        probs = F.softmax(log_probs, dim=-1)         # [B, K]

        # 熵 H(q(k|z)) = - Σ_k q_k log q_k
        entropy = -(probs * probs.clamp(min=1e-8).log()).sum(dim=-1).mean()

        # 均匀先验的 KL
        uniform = torch.ones_like(probs) / self.n_clusters
        kl = (probs * (probs.clamp(min=1e-8).log() - uniform.clamp(min=1e-8).log())).sum(dim=-1).mean()

        loss = alpha * kl + entropy

        return loss, probs


# ---------------------------------------------------------------------------
# 2. Contrastive Cluster Head（SCAN 风格）
# ---------------------------------------------------------------------------

class ContrastiveClusterHead(Module):
    """
    对比聚类头（基于 SCAN / DCC 思想）。

    核心思想：
      - 在隐空间中，同一聚类的样本应该靠近，不同聚类的样本应该远离
      - 使用 NT-Xent（Normalized Temperature-scaled Cross Entropy）损失
      - 同时最大化簇内一致性和簇间差异性

    损失函数：
      L_contrast = - (1/N) Σ_i log( exp(sim(z_i, z_j)/τ) / Σ_k exp(sim(z_i, z_k)/τ) )

    其中 j 是与 i 同一聚类的正样本，τ 是温度参数。
    """

    def __init__(
        self,
        embed_dim: int,
        n_clusters: int,
        hidden_dim: int = 256,
        temperature: float = 0.1,
        n_heads: int = 1,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.n_clusters = n_clusters
        self.temperature = temperature
        self.n_heads = n_heads

        # 投影头（将 z 投影到归一化超球面）
        self.projection = nn.Sequential(
            nn.Linear(embed_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
        )
        self.register_buffer("prototype_vectors", torch.zeros(n_clusters, hidden_dim))

    def initialize(self, z: torch.Tensor, method: str = "kmeans"):
        """用 K-Means 初始化聚类原型向量。"""
        with torch.no_grad():
            z_proj = self.projection(z).detach()
            z_proj = F.normalize(z_proj, dim=-1)

            kmeans = KMeans(n_clusters=self.n_clusters, n_init=20, random_state=42)
            labels = kmeans.fit_predict(z_proj.cpu().numpy())

            for k in range(self.n_clusters):
                mask = torch.from_numpy(labels == k).bool()
                if mask.sum() > 0:
                    self.prototype_vectors[k] = z_proj[mask].mean(dim=0)

            self.prototype_vectors = F.normalize(self.prototype_vectors, dim=-1)

    def forward(self, z: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        前向传播：计算聚类分配概率和投影嵌入。

        返回
        ----
        assignments : [B]  硬分配标签
        z_proj : [B, hidden_dim]  投影后的嵌入
        """
        z_proj = self.projection(z)
        z_proj = F.normalize(z_proj, dim=-1)

        # 计算与聚类原型的相似度
        sim = z_proj @ self.prototype_vectors.T  # [B, K]
        assignments = sim.argmax(dim=-1)           # [B]

        return assignments, z_proj

    def compute_contrastive_loss(
        self,
        z: torch.Tensor,
        labels: torch.Tensor,
        cell_type_labels: Optional[torch.Tensor] = None,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        计算对比聚类损失。

        参数
        ----
        z : [B, D]
        labels : [B]  聚类分配（来自 forward）
        cell_type_labels : [B]  真实细胞类型（可选，用于辅助监督）

        返回
        ----
        loss : NT-Xent 损失
        z_proj : [B, hidden_dim]
        """
        B = z.size(0)
        z_proj = self.projection(z)
        z_proj = F.normalize(z_proj, dim=-1)

        # ---- Step 1: 与聚类原型计算相似度 ----
        sim_proto = z_proj @ self.prototype_vectors.T  # [B, K]
        sim_proto = sim_proto / self.temperature

        # ---- Step 2: 同一聚类的样本作为正样本 ----
        device = z.device

        # 构建正样本掩码
        labels_expanded_i = labels.unsqueeze(1)  # [B, 1]
        labels_expanded_j = labels.unsqueeze(0)   # [1, B]
        positive_mask = (labels_expanded_i == labels_expanded_j).float().to(device)  # [B, B]
        positive_mask.fill_diagonal_(0)  # 去掉自身

        # ---- Step 3: NT-Xent 损失 ----
        logits = sim_proto  # [B, K] 但我们需要 [B, B]
        # 使用原型相似度作为代理
        # 实际上用 z_i 与 z_j 的相似度

        # 重新计算：z_i · z_j
        sim_matrix = z_proj @ z_proj.T / self.temperature  # [B, B]

        # 数值稳定性
        sim_matrix = sim_matrix - sim_matrix.max(dim=-1, keepdim=True)[0]

        exp_sim = torch.exp(sim_matrix)  # [B, B]

        # 分母：所有相似度之和（除了自身）
        denom = exp_sim.sum(dim=-1, keepdim=True) - 1.0  # 去掉自身

        # 分子：正样本的相似度
        pos_sim = (exp_sim * positive_mask).sum(dim=-1)  # [B]
        pos_count = positive_mask.sum(dim=-1).clamp(min=1.0)  # [B]

        # NT-Xent
        loss = -torch.log(pos_sim / denom.squeeze(-1) + 1e-8).mean()

        return loss, z_proj


# ---------------------------------------------------------------------------
# 3. DEC-style Cluster Head（Deep Embedded Clustering）
# ---------------------------------------------------------------------------

class DECClusterHead(Module):
    """
    DEC 风格聚类头（基于 Xie et al. 2016）。

    核心思想：
      - 使用 Student's t 分布作为软分配核
      - q_{ik} = (1 + ||z_i - μ_k||²/α)^(-(α+1)/2) / Σ_j (1 + ||z_i - μ_j||²/α)^(-(α+1)/2)

    损失函数：
      L_KL = KL(P || Q) = Σ_i Σ_k p_{ik} log(p_{ik}/q_{ik})

    其中 P 是目标分布（用 q 的平方归一化得到）。
    """

    def __init__(
        self,
        embed_dim: int,
        n_clusters: int,
        alpha: float = 1.0,
    ):
        super().__init__()
        self.embed_dim = embed_dim
        self.n_clusters = n_clusters
        self.alpha = alpha

        # 可学习的聚类中心
        self.cluster_centers = nn.Parameter(torch.randn(n_clusters, embed_dim))

    def initialize(self, z: torch.Tensor, method: str = "kmeans"):
        """用 K-Means 初始化聚类中心。"""
        with torch.no_grad():
            kmeans = KMeans(n_clusters=self.n_clusters, n_init=20, random_state=42)
            labels = kmeans.fit_predict(z.cpu().numpy())
            self.cluster_centers.data.copy_(torch.from_numpy(kmeans.cluster_centers_).float())
        return labels

    def forward(self, z: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """
        前向传播：计算软分配概率。

        返回
        ----
        q : [B, K]  学生分布 q_{ik}
        p : [B, K]  目标分布 p_{ik}（用于 KL 损失）
        """
        # 计算 z 与聚类中心的距离
        dist = torch.cdist(z, self.cluster_centers, p=2) ** 2  # [B, K]

        # q_{ik} = (1 + d/α)^(-(α+1)/2)
        q = (1 + dist / self.alpha).pow(-(self.alpha + 1) / 2)
        q = q / (q.sum(dim=-1, keepdim=True) + 1e-8)  # 归一化

        # 目标分布 P：p_{ik} = q_{ik}² / Σ_k q_{ik}²
        p = q ** 2
        p = p / (p.sum(dim=-1, keepdim=True) + 1e-8)

        hard_assign = q.argmax(dim=-1)
        return hard_assign, q

    def compute_kl_loss(self, z: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        """计算 DEC 的 KL 散度损失。"""
        _, q = self.forward(z)
        _, p = self.forward(z)  # 重新计算 p（确保一致）

        # 重新计算 p（避免循环）
        dist = torch.cdist(z, self.cluster_centers, p=2) ** 2
        q_tmp = (1 + dist / self.alpha).pow(-(self.alpha + 1) / 2)
        q_tmp = q_tmp / (q_tmp.sum(dim=-1, keepdim=True) + 1e-8)
        p_tmp = q_tmp ** 2
        p_tmp = p_tmp / (p_tmp.sum(dim=-1, keepdim=True) + 1e-8)

        kl = p_tmp * (p_tmp.clamp(min=1e-8).log() - q_tmp.clamp(min=1e-8).log())
        loss = kl.sum(dim=-1).mean()

        return loss, q_tmp


# ---------------------------------------------------------------------------
# 工厂函数
# ---------------------------------------------------------------------------

class ClusterHeadFactory:
    """根据策略名称返回对应的聚类头。"""

    _HEAD_MAP = {
        "gmm": GMMClusterHead,
        "contrastive": ContrastiveClusterHead,
        "dec": DECClusterHead,
    }

    @classmethod
    def create(
        cls,
        strategy: Literal["gmm", "contrastive", "dec"],
        embed_dim: int,
        n_clusters: int,
        hidden_dim: int = 256,
        n_heads: int = 1,
        temperature: float = 0.1,
        alpha: float = 1.0,
        use_diagonal_cov: bool = True,
        eps: float = 1e-8,
        topk_k: int = 50,
        dropout: float = 0.1,
    ) -> Module:
        if strategy not in cls._HEAD_MAP:
            raise ValueError(f"Unknown strategy '{strategy}'. Available: {list(cls._HEAD_MAP.keys())}")

        if strategy == "gmm":
            return cls._HEAD_MAP[strategy](
                embed_dim=embed_dim,
                n_clusters=n_clusters,
                hidden_dim=hidden_dim,
                use_diagonal_cov=use_diagonal_cov,
                eps=eps,
            )
        elif strategy == "contrastive":
            return cls._HEAD_MAP[strategy](
                embed_dim=embed_dim,
                n_clusters=n_clusters,
                hidden_dim=hidden_dim,
                temperature=temperature,
                n_heads=n_heads,
            )
        elif strategy == "dec":
            return cls._HEAD_MAP[strategy](
                embed_dim=embed_dim,
                n_clusters=n_clusters,
                alpha=alpha,
            )


if __name__ == "__main__":
    # 简单测试
    import torch
    from sklearn.datasets import make_blobs

    B, D, K = 200, 64, 5
    device = "cuda" if torch.cuda.is_available() else "cpu"

    # 模拟数据
    X, _ = make_blobs(n_samples=B, n_features=D, centers=K, random_state=42)
    z = torch.from_numpy(X).float().to(device)

    for strategy in ["gmm", "contrastive", "dec"]:
        head = ClusterHeadFactory.create(
            strategy=strategy,
            embed_dim=D,
            n_clusters=K,
            hidden_dim=128,
        ).to(device)

        head.initialize(z, method="kmeans")

        if strategy in ["gmm", "dec"]:
            assign, probs = head(z[:32])
            print(f"{strategy}: assignments={assign.shape}, probs={probs.shape if probs is not None else None}")

            if strategy == "gmm":
                loss, _ = head.compute_clustering_loss(z[:32])
            elif strategy == "dec":
                loss, _ = head.compute_kl_loss(z[:32])
            print(f"  loss={loss:.4f}")

        elif strategy == "contrastive":
            assign, z_proj = head(z[:32])
            loss, _ = head.compute_contrastive_loss(z[:32], assign)
            print(f"{strategy}: loss={loss:.4f}")

    print("cluster_head.py test passed.")
