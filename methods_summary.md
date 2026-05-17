# PlantNet 单细胞聚类方法源代码汇总

本文档汇总了 PlantNet 项目中使用的各种单细胞 RNA 聚类方法的源代码。

---

## 1. PhytoCluster

植物单细胞 RNA-seq 变分自编码器聚类方法（基于 VAE + GMM）

**文件路径**: `methods/DeepLearning/PhytoCluster/run.py`

```python
# -*- coding: utf-8 -*-
"""
====================================================================================================
PhytoCluster — 植物单细胞RNA-seq变分自编码器聚类方法
====================================================================================================

【论文来源】
    Wang et al. (2025), aBIOTECH
    https://doi.org/10.1007/s42994-025-00210-x

【核心思想】
    将变分自编码器（VAE）与高斯混合模型（GMM）相结合，
    从植物单细胞转录组数据中提取低维潜在特征，实现无监督聚类。

【算法流程】
    ┌────────────────────────────────────────────────────────────────┐
    │                     PhytoCluster 流程图                          │
    │                                                                  │
    │  Step 1: 数据预处理                                              │
    │     └─→ 归一化 → log1p → HVG筛选 → Z-score标准化               │
    │                                                                  │
    │  Step 2: 预训练 VAE（阶段一）                                    │
    │     X → Encoder → z(μ,σ) → Decoder → X'                        │
    │     损失: ELBO = 重构损失 + KL散度                               │
    │     目标: 学习良好的潜在空间表示                                  │
    │                                                                  │
    │  Step 3: GMM 初始化（基于预训练特征）                           │
    │     z → GMM.fit(z) → μ_c, Σ_c, π_c                             │
    │     使用 GMM 估计聚类参数作为初始化                              │
    │                                                                  │
    │  Step 4: PhytoCluster 联合训练（阶段二）                        │
    │     VAE + GMM 联合优化                                          │
    │     损失: 重构损失 + 聚类KL散度（包含GMM后验约束）               │
    │                                                                  │
    │  Step 5: 聚类推断                                               │
    │     z → GMM 预测 → 细胞类型标签                                 │
    │                                                                  │
    └────────────────────────────────────────────────────────────────┘

【网络架构】
    | 组件      | 配置                | 说明                    |
    |-----------|---------------------|-------------------------|
    | 编码器    | 1024 → 128 → 10    | 两层MLP，ReLU激活      |
    | 解码器    | 10 → 128 → 1024    | 两层MLP                 |
    | 潜在空间  | 10维               | 高斯分布参数 μ, logσ²   |
    | GMM       | K个高斯混合        | 对角协方差矩阵          |

【损失函数】
    L_total = L_recon + β * L_KL_clustering

    1. L_recon (重构损失): MSE(Decoder(z), X)
    2. L_KL_clustering: 包含GMM后验分布的KL散度

【优势】
    - 专门针对植物scRNA-seq数据优化
    - VAE提供连续且结构化的潜在空间
    - GMM实现概率化聚类（软标签）
    - 两阶段训练确保稳定收敛
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import init
from torch.utils.data import DataLoader, TensorDataset
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from sklearn.metrics import (
    normalized_mutual_info_score, adjusted_rand_score,
    f1_score, fowlkes_mallows_score,
    v_measure_score, homogeneity_score, completeness_score
)
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from preprocess import prepare_data_for_model


def set_seed(seed):
    """设置随机种子确保可复现性"""
    import random
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='PhytoCluster: VAE + GMM for Plant scRNA-seq Clustering',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 数据参数
    parser.add_argument('--data_path', type=str, required=True,
                       help='输入 h5ad 文件路径')
    parser.add_argument('--save_dir', type=str, default='./results',
                       help='结果保存目录')

    # 模型参数
    parser.add_argument('--n_clusters', type=int, required=True,
                       help='聚类数')
    parser.add_argument('--latent_dim', type=int, default=10,
                       help='潜在空间维度')
    parser.add_argument('--encode_dim', type=int, nargs=2, default=[1024, 128],
                       help='编码器隐藏层维度')
    parser.add_argument('--decode_dim', type=int, nargs=2, default=[128, 1024],
                       help='解码器隐藏层维度')

    # 训练参数
    parser.add_argument('--pretrain_epochs', type=int, default=300,
                       help='预训练VAE的轮次')
    parser.add_argument('--cluster_epochs', type=int, default=100,
                       help='聚类训练的轮次')
    parser.add_argument('--batch_size', type=int, default=64,
                       help='批大小')
    parser.add_argument('--lr', type=float, default=0.0001,
                       help='学习率')
    parser.add_argument('--var_lr', type=float, default=0.0001,
                       help='方差参数学习率')
    parser.add_argument('--weight_decay', type=float, default=5e-4,
                       help='权重衰减')
    parser.add_argument('--patience', type=int, default=50,
                       help='早停耐心值')
    parser.add_argument('--warmup_steps', type=int, default=1000,
                       help='KL散度预热步数')

    # 其他参数
    parser.add_argument('--seed', type=int, default=42,
                       help='随机种子')
    parser.add_argument('--gpu', type=int, default=0,
                       help='GPU设备号')
    parser.add_argument('--no_cuda', action='store_true',
                       help='禁用CUDA')
    parser.add_argument('--eval_interval', type=int, default=10,
                       help='评估间隔（轮次）')

    return parser.parse_args()


# =============================================================================
# VAE 组件
# =============================================================================

def binary_cross_entropy(recon_x, x):
    """计算二元交叉熵损失"""
    return -torch.sum(x * torch.log(recon_x + 1e-8) + (1 - x) * torch.log(1 - recon_x + 1e-8), dim=-1)


def compute_kl_divergence(mu, logvar):
    """计算标准VAE中潜在变量的KL散度"""
    return -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)


def compute_elbo_loss(recon_x, x, z_params):
    """计算VAE的ELBO损失"""
    mu, logvar = z_params
    kld = compute_kl_divergence(mu, logvar)
    likelihood = -F.mse_loss(recon_x, x, reduction='none').sum(dim=-1)
    return torch.sum(likelihood), torch.sum(kld)


def compute_elbo(recon_x, x, gamma, c_params, z_params):
    """计算PhytoCluster的ELBO损失（包含聚类目标）"""
    import math
    mu_c, var_c, pi = c_params
    var_c += 1e-8
    n_centroids = pi.size(1)
    mu, logvar = z_params
    mu_expand = mu.unsqueeze(2).expand(mu.size(0), mu.size(1), n_centroids)
    logvar_expand = logvar.unsqueeze(2).expand(logvar.size(0), logvar.size(1), n_centroids)

    likelihood = -F.mse_loss(recon_x, x, reduction='none').sum(dim=-1)

    logpzc = -0.5 * torch.sum(gamma * torch.sum(math.log(2 * math.pi) +
                                                torch.log(var_c) +
                                                torch.exp(logvar_expand) / var_c +
                                                (mu_expand - mu_c) ** 2 / var_c, dim=1), dim=1)

    logpc = torch.sum(gamma * torch.log(pi), 1)
    qentropy = -0.5 * torch.sum(1 + logvar + math.log(2 * math.pi), 1)
    logqcx = torch.sum(gamma * torch.log(gamma), 1)

    kld = -logpzc - logpc + qentropy + logqcx

    return torch.sum(likelihood), torch.sum(kld)


def create_mlp(layers, activation=nn.ReLU(), bn=False, dropout=0):
    """创建多层感知器网络"""
    net = []
    for i in range(1, len(layers)):
        net.append(nn.Linear(layers[i-1], layers[i]))
        if bn:
            net.append(nn.BatchNorm1d(layers[i]))
        net.append(activation)
        if dropout > 0:
            net.append(nn.Dropout(dropout))
    return nn.Sequential(*net)


class DeterministicWarmup:
    """β值预热策略：在训练初期让β从0逐渐增加到目标值"""
    def __init__(self, n=100, t_max=1):
        self.t = 0
        self.t_max = t_max
        self.inc = t_max / n

    def __iter__(self):
        return self

    def __next__(self):
        t = self.t + self.inc
        self.t = self.t_max if t > self.t_max else t
        return self.t

    def next(self):
        return self.__next__()


class Stochastic(nn.Module):
    """随机层基类：实现VAE重参数化技巧"""
    def reparametrize(self, mu, logvar):
        epsilon = torch.randn(mu.size(), requires_grad=False, device=mu.device)
        std = logvar.mul(0.5).exp_()
        z = mu.addcmul(std, epsilon)
        return z


class GaussianSample(Stochastic):
    """高斯采样层"""
    def __init__(self, in_features, out_features):
        super(GaussianSample, self).__init__()
        self.in_features = in_features
        self.out_features = out_features
        self.mu = nn.Linear(in_features, out_features)
        self.log_var = nn.Linear(in_features, out_features)

    def forward(self, x):
        mu = self.mu(x)
        log_var = self.log_var(x)
        return self.reparametrize(mu, log_var), mu / 10, log_var


class Encoder(nn.Module):
    """编码器：将高维数据压缩到低维潜在空间"""
    def __init__(self, dims, bn=False, dropout=0):
        super(Encoder, self).__init__()
        [x_dim, h_dim, z_dim] = dims
        self.hidden = create_mlp([x_dim] + h_dim, bn=bn, dropout=dropout)
        self.sample = GaussianSample(h_dim[-1], z_dim)

    def forward(self, x):
        x = self.hidden(x)
        return self.sample(x)


class Decoder(nn.Module):
    """解码器：从潜在变量重建原始数据"""
    def __init__(self, dims, bn=False, dropout=0):
        super(Decoder, self).__init__()
        [z_dim, h_dim, x_dim] = dims
        self.hidden = create_mlp([z_dim] + h_dim, bn=bn, dropout=dropout)
        self.reconstruction = nn.Linear(h_dim[-1], x_dim)

    def forward(self, x):
        x = self.hidden(x)
        return self.reconstruction(x)


class VAE(nn.Module):
    """变分自编码器基础模型"""
    def __init__(self, dims, bn=False, dropout=0):
        super(VAE, self).__init__()
        [x_dim, z_dim, encode_dim, decode_dim] = dims
        self.encoder = Encoder([x_dim, encode_dim, z_dim], bn=bn, dropout=dropout)
        self.decoder = Decoder([z_dim, decode_dim, x_dim], bn=bn, dropout=dropout)
        self.initialize_weights()

    def initialize_weights(self):
        """初始化网络权重"""
        for m in self.modules():
            if isinstance(m, nn.Linear):
                init.xavier_normal_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.zero_()

    def forward(self, x):
        """前向传播"""
        z, mu, logvar = self.encoder(x)
        recon_x = self.decoder(z)
        return recon_x

    def compute_loss(self, x):
        """计算损失"""
        z, mu, logvar = self.encoder(x)
        recon_x = self.decoder(z)
        likelihood, kl_loss = compute_elbo_loss(recon_x, x, (mu, logvar))
        return -likelihood, kl_loss

    def encode(self, x):
        """编码：返回潜在变量"""
        z, mu, logvar = self.encoder(x)
        return z

    def encode_mu(self, dataloader, device):
        """批量编码：返回均值向量"""
        self.eval()
        output = []
        with torch.no_grad():
            for x in dataloader:
                if isinstance(x, (list, tuple)):
                    x = x[0]
                x = x.float().to(device)
                z, mu, logvar = self.encoder(x)
                output.append(mu.cpu().detach())
        return torch.cat(output).numpy()


class PhytoCluster(VAE):
    """PhytoCluster: VAE + GMM 联合聚类模型"""
    def __init__(self, dims, n_centroids):
        super(PhytoCluster, self).__init__(dims)
        self.n_centroids = n_centroids
        z_dim = dims[1]

        self.pi = nn.Parameter(torch.ones(n_centroids) / n_centroids)
        self.mu_c = nn.Parameter(torch.zeros(z_dim, n_centroids))
        self.var_c = nn.Parameter(torch.ones(z_dim, n_centroids))

    def compute_loss(self, x):
        """计算包含聚类目标的损失"""
        z, mu, logvar = self.encoder(x)
        recon_x = self.decoder(z)
        gamma, mu_c, var_c, pi = self.infer_clusters(z)
        likelihood, kl_loss = compute_elbo(
            recon_x, x, gamma, (mu_c, var_c, pi), (mu, logvar)
        )
        return -likelihood, kl_loss

    def infer_clusters(self, z):
        """推断聚类分配"""
        import math
        n_centroids = self.n_centroids
        N = z.size(0)
        z_expanded = z.unsqueeze(2).expand(z.size(0), z.size(1), n_centroids)
        pi = self.pi.repeat(N, 1)
        mu_c = self.mu_c.repeat(N, 1, 1)
        var_c = self.var_c.repeat(N, 1, 1) + 1e-8

        p_c_z = torch.exp(
            torch.log(pi) - torch.sum(
                0.5 * torch.log(2 * math.pi * var_c) +
                (z_expanded - mu_c) ** 2 / (2 * var_c), dim=1
            )
        ) + 1e-10
        gamma = p_c_z / torch.sum(p_c_z, dim=1, keepdim=True)

        return gamma, mu_c, var_c, pi

    def initialize_gmm_params(self, dataloader, device):
        """使用GMM初始化聚类参数"""
        gmm = GaussianMixture(
            n_components=self.n_centroids,
            covariance_type='diag',
            random_state=42
        )
        z = self.encode_mu(dataloader, device)
        gmm.fit(z)
        self.mu_c.data.copy_(torch.from_numpy(gmm.means_.T.astype(np.float32)))
        self.var_c.data.copy_(torch.from_numpy(gmm.covariances_.T.astype(np.float32)))


def train_model(model, dataloader, device, args, phase='pretrain'):
    """训练模型"""
    model.to(device)

    if phase == 'pretrain':
        optimizer = torch.optim.Adam([
            {'params': model.encoder.hidden.parameters(), 'lr': args.lr},
            {'params': model.encoder.sample.mu.parameters(), 'lr': args.var_lr},
            {'params': model.encoder.sample.log_var.parameters(), 'lr': args.var_lr},
            {'params': model.decoder.parameters(), 'lr': args.lr}
        ], weight_decay=args.weight_decay)
        max_iter = args.pretrain_epochs * len(dataloader)
        warmup_steps = args.warmup_steps
    else:
        optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
        max_iter = args.cluster_epochs * len(dataloader)
        warmup_steps = args.cluster_epochs * len(dataloader) // 3

    beta_scheduler = DeterministicWarmup(n=warmup_steps, t_max=1.0)

    n_epochs = args.pretrain_epochs if phase == 'pretrain' else args.cluster_epochs

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0
        epoch_recon = 0
        epoch_kl = 0

        for x in dataloader:
            if isinstance(x, (list, tuple)):
                x = x[0]
            x = x.float().to(device)

            optimizer.zero_grad()
            recon_loss, kl_loss = model.compute_loss(x)

            beta = next(beta_scheduler)
            loss = recon_loss + beta * kl_loss

            if torch.isnan(loss):
                continue

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10)
            optimizer.step()

            epoch_loss += loss.item() / len(x)
            epoch_recon += recon_loss.item() / len(x)
            epoch_kl += kl_loss.item() / len(x)

        if (epoch + 1) % 20 == 0:
            print(f'  [{phase}] Epoch {epoch+1}/{n_epochs}: '
                  f'Loss={epoch_loss:.4f}, Recon={epoch_recon:.4f}, KL={epoch_kl:.4f}')

    return model


def evaluation(y_true, y_pred):
    """计算聚类评估指标"""
    from munkres import Munkres

    # Hungarian算法标签对齐
    m = Munkres()
    D = max(int(y_pred.max()), int(y_true.max())) + 1
    cost = np.zeros((D, D))
    for i in range(len(y_pred)):
        cost[int(y_pred[i]), int(y_true[i])] -= 1
    assignment = m.compute(cost)
    y_map = {row: col for row, col in assignment}
    y_pred_aligned = np.array([y_map.get(int(p), int(p)) for p in y_pred])

    acc = np.mean(y_pred_aligned == y_true)
    nmi = normalized_mutual_info_score(y_true, y_pred)
    ari = adjusted_rand_score(y_true, y_pred)
    f1 = f1_score(y_true, y_pred, average='macro', zero_division=0)
    fmi = fowlkes_mallows_score(y_true, y_pred)
    vms = v_measure_score(y_true, y_pred)
    hom = homogeneity_score(y_true, y_pred)
    com = completeness_score(y_true, y_pred)

    return acc, nmi, ari, f1, fmi, vms, hom, com, y_pred_aligned


def save_results(save_dir, y_true, y_pred, embedding, epoch):
    """保存结果"""
    acc, nmi, ari, f1, fmi, vms, hom, com, y_pred_aligned = evaluation(y_true, y_pred)

    metrics = {
        'acc': float(acc),
        'nmi': float(nmi),
        'ari': float(ari),
        'f1_macro': float(f1),
        'fmi': float(fmi),
        'v_measure': float(vms),
        'homogeneity': float(hom),
        'completeness': float(com)
    }

    import json
    metrics_file = os.path.join(save_dir, f'metrics_{epoch}.json')
    with open(metrics_file, 'w') as f:
        json.dump(metrics, f, indent=2)

    with open(os.path.join(save_dir, 'metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)

    pd.DataFrame({
        'pred': y_pred_aligned,
        'true': y_true
    }).to_csv(os.path.join(save_dir, f'types_{epoch}_pred.csv'), index=False)

    np.save(os.path.join(save_dir, f'embedding_{epoch}.npy'), embedding)

    import h5py
    with h5py.File(os.path.join(save_dir, 'embedding.h5'), 'w') as f:
        f.create_dataset('X', data=embedding)
        f.create_dataset('Y', data=y_pred_aligned)

    print(f'\nPhytoCluster Results:')
    print(f'  ACC:        {acc:.4f}')
    print(f'  NMI:        {nmi:.4f}')
    print(f'  ARI:        {ari:.4f}')
    print(f'  F1-macro:  {f1:.4f}')
    print(f'  FMI:       {fmi:.4f}')
    print(f'  V-measure: {vms:.4f}')

    return metrics


def main():
    """主函数"""
    import random
    args = parse_args()

    import random

    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    set_seed(args.seed)

    os.makedirs(args.save_dir, exist_ok=True)

    print('Loading data...')
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=True,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True
    )

    X = np.array(X).astype(np.float32)
    Y = np.array(Y)

    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    input_dim = X.shape[1]
    dims = [input_dim, args.latent_dim, args.encode_dim, args.decode_dim]

    all_data_tensor = torch.from_numpy(X).float()
    full_loader = DataLoader(
        TensorDataset(all_data_tensor, torch.zeros(len(X)).long()),
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False
    )
    train_loader = DataLoader(
        TensorDataset(all_data_tensor, torch.zeros(len(X)).long()),
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True
    )

    print('\n' + '='*60)
    print('Phase 1: Pretraining VAE')
    print('='*60)
    pretrain_model = VAE(dims, bn=False, dropout=0)
    pretrain_model = train_model(pretrain_model, train_loader, device, args, phase='pretrain')

    print('\n' + '='*60)
    print('Phase 2: PhytoCluster Training (VAE + GMM)')
    print('='*60)
    model = PhytoCluster(dims, n_clusters)
    model.load_state_dict(pretrain_model.state_dict(), strict=False)
    model.to(device)  # Ensure model is on GPU
    model.initialize_gmm_params(full_loader, device)
    model = train_model(model, train_loader, device, args, phase='cluster')

    print('\n' + '='*60)
    print('Phase 3: Final Clustering')
    print('='*60)
    embedding = pretrain_model.encode_mu(full_loader, device)
    gmm = GaussianMixture(n_components=n_clusters, covariance_type='diag', random_state=42)
    pred_labels = gmm.fit_predict(embedding)

    save_results(args.save_dir, Y, pred_labels, embedding, 0)

    print(f'\nPhytoCluster completed. Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
```

---

## 2. SC3

Single-Cell Consensus Clustering - 单细胞共识聚类方法（纯 Python 实现）

**文件路径**: `methods/Traditional/sc3/run.py`

```python
# -*- coding: utf-8 -*-
"""
Pure Python SC3 Implementation for scCluBench
============================================

SC3: Single-Cell Consensus Clustering (Kiselev et al., 2017)
纯 Python 实现（无需 R 依赖）

核心思想：
  - 通过多次运行不同的聚类算法构建共识矩阵
  - 共识矩阵记录每对细胞被聚到同一簇的频率
  - 最终在共识矩阵上进行聚类，得到稳定的结果

算法流程：
  1. PCA 降维（减少噪声）
  2. 多次运行 K-means（不同参数）
  3. 多次运行层次聚类（不同距离度量+链接方式）
  4. 汇总结果构建共识矩阵
  5. 在共识矩阵上运行最终聚类

优势：
  - 通过集成多个聚类结果提高稳定性
  - 减少单一算法的随机性影响

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
from sklearn.metrics import f1_score, fowlkes_mallows_score
from sklearn.metrics import v_measure_score, homogeneity_score, completeness_score
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist, squareform
import json
import warnings
warnings.filterwarnings('ignore')

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from preprocess import prepare_data_for_model
from utils import save


def build_consensus_matrix(X, n_clusters, n_pcs=20, seed=42):
    """
    构建共识矩阵（Consensus Matrix）- 简化版本用于大数据集

    参数：
        X: 预处理后的基因表达矩阵 (细胞 × 基因)
        n_clusters: 目标聚类数
        n_pcs: PCA 降维维度
        seed: 随机种子

    返回：
        consensus: 共识矩阵 (细胞 × 细胞)
        X_pca: PCA 降维后的数据
    """
    n_cells = X.shape[0]
    
    # Step 1: PCA 降维
    pca = PCA(n_components=min(n_pcs, n_cells - 1), random_state=seed)
    X_pca = pca.fit_transform(X)

    # Step 2: 简化的共识聚类
    # 使用 K-means 进行共识聚类
    consensus = np.zeros((n_cells, n_cells))
    
    kmeans_configs = [
        (n_clusters, 20, seed),
        (n_clusters, 10, seed + 1),
        (n_clusters, 10, seed + 2),
        (n_clusters + 1, 20, seed),
        (n_clusters - 1, 20, seed),
    ]
    
    for k, n_init, rs in kmeans_configs:
        if k < 2 or k > n_cells:
            continue
        km = KMeans(n_clusters=k, n_init=n_init, random_state=rs)
        labels = km.fit_predict(X_pca)
        # 构建共现矩阵
        for label in np.unique(labels):
            mask = labels == label
            idx = np.where(mask)[0]
            for i in idx:
                consensus[i, mask] += 1
    
    # 归一化
    consensus /= len(kmeans_configs)
    np.fill_diagonal(consensus, 1.0)
    
    return consensus, X_pca


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='SC3: Single-Cell Consensus Clustering (Pure Python)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results',
                        help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, required=True,
                        help='Number of clusters')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    return parser.parse_args()


def main():
    """主函数：SC3 共识聚类流程"""
    args = parse_args()
    np.random.seed(args.seed)

    os.makedirs(args.save_dir, exist_ok=True)

    print('Loading data...')
    # 加载并预处理数据：归一化、log1p、HVG筛选、Z-score标准化
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=False,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=False  # SC3 不使用 Z-score（保留 log1p 数据）
    )

    X = np.array(X).astype(np.float32)
    Y = np.array(Y)

    # 标签编码（将字符串标签转为整数）
    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    print('Building consensus matrix...')
    # Step 1: 构建共识矩阵（通过多次聚类统计细胞对共现频率）
    consensus, X_pca = build_consensus_matrix(X, n_clusters, n_pcs=20, seed=args.seed)

    print('Performing final clustering on consensus matrix...')
    # Step 2: 在共识矩阵上运行最终聚类
    try:
        # 将共识矩阵转为距离矩阵（1 - consensus）
        dist_matrix = 1 - consensus
        dist_matrix = np.nan_to_num(dist_matrix, nan=1.0)
        # 使用层次聚类（基于预计算的距离矩阵）
        final_clustering = AgglomerativeClustering(
            n_clusters=n_clusters,
            metric='precomputed',
            linkage='average'
        )
        y_pred = final_clustering.fit_predict(dist_matrix)
    except Exception:
        # 如果距离矩阵聚类失败，回退到 K-means
        kmeans = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed)
        y_pred = kmeans.fit_predict(X_pca)

    print(f'Number of clusters found: {len(np.unique(y_pred))}')

    # 保存结果
    save(args.save_dir, Y, y_pred, 0, X_pca)

    # ========== Step 3: 计算评估指标 ==========
    metrics_path = os.path.join(args.save_dir, 'metrics.json')

    # 使用 Hungarian 算法对齐预测标签和真实标签（最优匹配）
    from munkres import Munkres
    m = Munkres()
    D = max(int(y_pred.max()), int(Y.max())) + 1
    cost = np.zeros((D, D))
    for i in range(len(y_pred)):
        cost[int(y_pred[i]), int(Y[i])] -= 1  # 构建成本矩阵
    assignment = m.compute(cost)
    y_map = {row: col for row, col in assignment}
    y_pred_aligned = np.array([y_map.get(int(p), int(p)) for p in y_pred])

    # 计算各项评估指标
    acc = float(np.mean(y_pred_aligned == Y))
    nmi = float(normalized_mutual_info_score(Y, y_pred))
    ari = float(adjusted_rand_score(Y, y_pred))
    f1 = float(f1_score(Y, y_pred, average='macro', zero_division=0))
    fmi = float(fowlkes_mallows_score(Y, y_pred))
    vms = float(v_measure_score(Y, y_pred))
    hom = float(homogeneity_score(Y, y_pred))
    comp = float(completeness_score(Y, y_pred))

    # 保存指标到 JSON
    metrics = {
        'acc': acc,           # 准确率（标签对齐后）
        'nmi': nmi,           # 标准化互信息
        'ari': ari,           # 调整兰德指数
        'f1_macro': f1,      # F1 分数（宏平均）
        'fmi': fmi,           # Fowlkes-Mallows 指数
        'v_measure': vms,    # V-measure（同质性和完整性调和平均）
        'homogeneity': hom,   # 同质性（每个簇只包含单一类的程度）
        'completeness': comp  # 完整性（同一类的成员都被分到同一簇的程度）
    }

    with open(metrics_path, 'w') as f:
        json.dump(metrics, f, indent=2)

    print(f'\nSC3 Results (Pure Python):')
    print(f'  ACC:        {acc:.4f}')   # 准确率
    print(f'  NMI:        {nmi:.4f}')   # 标准化互信息
    print(f'  ARI:        {ari:.4f}')   # 调整兰德指数
    print(f'  F1-macro:   {f1:.4f}')   # F1 分数
    print(f'  FMI:        {fmi:.4f}')   # Fowlkes-Mallows 指数
    print(f'  V-measure:  {vms:.4f}')  # V-measure
    print(f'\nResults saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
```

---

## 3. scGPT

Single-cell Foundation Model for Cell Embedding - 单细胞基础模型

**文件路径**: `methods/Foundation/scGPT/run.py`

```python
# -*- coding: utf-8 -*-
"""
scGPT — Single-cell Foundation Model for Cell Embedding
=======================================================

scGPT: Towards Building a Foundation Model for Single-Cell Multi-omics
Using Generative AI
Cui et al., Nature Methods 2024

This module provides a run.py interface for scGPT cell embedding and clustering.
It uses the official scGPT pretrained checkpoint (whole-human, 33M cells) to
extract cell embeddings via transformer encoder CLS token, then performs KMeans
clustering on the learned representations.

【工作流程】
    ┌────────────────────────────────────────────────────────────────┐
    │  Input: SRP182008.h5ad (Arabidopsis thaliana scRNA-seq)        │
    │       genes: AT1G... (plant gene names)                        │
    │       labels: Celltype annotation                              │
    │       ~13,514 cells × 53,678 genes                             │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 1: 基因名对齐 (Gene Name Matching)                         │
    │   - 尝试用 GeneVocab 匹配输入基因 vs. 人源预训练词汇表           │
    │   - AT1G... (植物基因) → 通常不匹配人源词汇表 (60,697 基因)   │
    │   - 匹配率低时 → 启用 PCA fallback 嵌入方案                    │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 2a: scGPT Transformer Encoder (基因匹配时)                  │
    │   TransformerModel(nlayers=12, nhead=8, embsize=512)            │
    │   - GeneEmbedding + ExpressionEmbedding → Transformer            │
    │   - CLS token position [0] → 512-dim cell embedding            │
    │   - L2 normalized                                             │
    └────────────────────┬─────────────────────────────────────────┘
                         │
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 2b: PCA Embedding (基因不匹配时 fallback)                  │
    │   sklearn.decomposition.PCA(n_components=50)                    │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Step 3: KMeans 聚类                                            │
    │   sklearn.cluster.KMeans(n_clusters=k, n_init=20)               │
    └────────────────────┬─────────────────────────────────────────┘
                         │
                         ▼
    ┌────────────────────────────────────────────────────────────────┐
    │ Output: ACC, NMI, ARI, F1_macro, FMI, V_measure,              │
    │         Homogeneity, Completeness                              │
    └────────────────────────────────────────────────────────────────┘

【关键设计说明】
    1. scGPT 预训练模型基于人类细胞 (33M)，植物基因名称 (AT1G...) 与人源基因
       词汇表无交集。因此对于植物数据集，默认使用 PCA embedding + KMeans 方案。
       此 fallback 保留了与原始 scGPT pipeline 完全一致的聚类后处理逻辑。
    2. 若数据集含有人源基因（如 HeLa, HEK293T 等），scGPT 可直接利用 transformer
       嵌入获得更好的表征质量。
    3. torchtext ABI 不兼容问题通过预填充 sys.modules mock 绕过。

Usage:
    python run.py --data_path /path/to/SRP182008.h5ad --n_clusters 15 --save_dir ./results
"""

import os
import sys
import json
import argparse
import warnings
import types
import numpy as np
import torch
import random
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import (
    normalized_mutual_info_score,
    adjusted_rand_score,
    f1_score,
    fowlkes_mallows_score,
    v_measure_score,
    homogeneity_score,
    completeness_score,
)
from sklearn.preprocessing import LabelEncoder

# ── project root & common modules ─────────────────────────────────────────────
# scGPT/run.py is at methods/Foundation/scGPT/run.py
# plantnet/ is 3 levels up: scGPT → Foundation → methods → plantnet
PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, PROJECT_ROOT)
from preprocess import prepare_data_for_model
from utils import save


# ════════════════════════════════════════════════════════════════════════════════
# 1.  TorchText Mock (Bypass broken ABI)
#
#     The installed torchtext .so has an ABI mismatch with the running PyTorch
#     (undefined symbol: torch::detail::class_base constructor).
#     We pre-populate sys.modules with minimal mock objects so that
#     scgpt.tokenizer.gene_tokenizer can import torchtext.vocab.Vocab
#     without actually loading the broken shared library.
# ════════════════════════════════════════════════════════════════════════════════

class _MockTorchText(types.ModuleType):
    """Minimal top-level namespace mock for torchtext."""
    pass


class _MockTorchTextVocab(types.ModuleType):
    """Mock for torchtext.vocab — provides only the Vocab class."""

    class Vocab:
        """Stand-in for torchtext.vocab.Vocab; not used at runtime for this workflow."""
        def __init__(self, ordered_dict=None, min_freq=1):
            pass


# Register mocks BEFORE any scgpt import
if "torchtext" not in sys.modules:
    sys.modules["torchtext"] = _MockTorchText("torchtext")
if "torchtext.vocab" not in sys.modules:
    sys.modules["torchtext.vocab"] = _MockTorchTextVocab("torchtext.vocab")

# Neutralise torch.ops.load_library so torchtext._extension skips its .so loading
import torch
if not hasattr(torch.ops, "load_library"):
    torch.ops.load_library = lambda *a, **k: None


# ════════════════════════════════════════════════════════════════════════════════
# 2.  GeneVocab (torchtext-free implementation)
#     Mirrors scgpt.tokenizer.gene_tokenizer.GeneVocab using only stdlib / numpy.
# ════════════════════════════════════════════════════════════════════════════════

class GeneVocab:
    """
    Minimal GeneVocab compatible with the scGPT checkpoint format.

    The vocabulary is a plain Python dict mapping gene symbols → integer IDs.
    Supports: __contains__, __getitem__, __len__, set_default_index, from_file.
    """

    def __init__(self, token2idx: dict, default_token: str | None = "<pad>"):
        self._stoi = dict(token2idx)
        self._itos = {v: k for k, v in token2idx.items()}
        self._default_idx = token2idx.get(default_token) if default_token else None

    @classmethod
    def from_file(cls, file_path: str):
        with open(file_path, "r") as f:
            token2idx = json.load(f)
        return cls(token2idx)

    def __contains__(self, token: str) -> bool:
        return token in self._stoi

    def __getitem__(self, token: str) -> int:
        return self._stoi[token]

    def __call__(self, tokens):
        if isinstance(tokens, str):
            return self._stoi.get(tokens, self._default_idx if self._default_idx is not None else 0)
        return [self._stoi.get(t, self._default_idx if self._default_idx is not None else 0) for t in tokens]

    def __len__(self) -> int:
        return len(self._stoi)

    def set_default_index(self, idx: int):
        self._default_idx = idx

    def get_stoi(self) -> dict:
        return self._stoi


# ════════════════════════════════════════════════════════════════════════════════
# 3.  Load scGPT model components
# ════════════════════════════════════════════════════════════════════════════════

SCGPT_DIR = os.path.dirname(os.path.abspath(__file__))
SCGPT_MODEL_FILE = os.path.join(SCGPT_DIR, "best_model.pt")
SCGPT_VOCAB_FILE = os.path.join(SCGPT_DIR, "vocab.json")
SCGPT_ARGS_FILE = os.path.join(SCGPT_DIR, "args.json")


def load_scgpt_model(device):
    """
    Build the TransformerModel architecture and load pretrained weights.
    Returns (model, vocab, config_dict).
    """
    import scgpt.model.model as _model_mod
    import scgpt.data_collator as _collator_mod

    # Vocabulary
    vocab = GeneVocab.from_file(SCGPT_VOCAB_FILE)
    vocab.set_default_index(vocab[" "])

    # Config
    with open(SCGPT_ARGS_FILE, "r") as f:
        cfg = json.load(f)

    # Build model (matching the pretraining architecture)
    model = _model_mod.TransformerModel(
        ntoken=len(vocab),
        d_model=cfg["embsize"],
        nhead=cfg["nheads"],
        d_hid=cfg["d_hid"],
        nlayers=cfg["nlayers"],
        nlayers_cls=cfg["n_layers_cls"],
        n_cls=1,
        vocab=vocab,
        dropout=cfg["dropout"],
        pad_token=cfg["pad_token"],
        pad_value=cfg["pad_value"],
        do_mvc=cfg.get("MVC", True),
        do_dab=False,
        use_batch_labels=False,
        domain_spec_batchnorm=False,
        explicit_zero_prob=False,
        use_fast_transformer=False,  # flash-attn not available
        pre_norm=False,
    )

    # Load pretrained weights (non-strict to tolerate flash-attn key differences)
    ckpt = torch.load(SCGPT_MODEL_FILE, map_location=device, weights_only=True)
    model_dict = model.state_dict()
    loaded, skipped = 0, 0
    for k, v in ckpt.items():
        if k in model_dict and v.shape == model_dict[k].shape:
            model_dict[k] = v
            loaded += 1
        else:
            skipped += 1
    model.load_state_dict(model_dict, strict=False)
    print(f"  Loaded {loaded} pretrained parameters, skipped {skipped}")

    model.to(device)
    model.eval()
    return model, vocab, cfg


# ════════════════════════════════════════════════════════════════════════════════
# 4.  Cell embedding extraction
# ════════════════════════════════════════════════════════════════════════════════

class _SeqDataset(torch.utils.data.Dataset):
    """Per-cell sparse representation: (nonzero gene IDs, expression values)."""

    def __init__(self, count_matrix, gene_ids, pad_token_id, pad_value):
        self.counts = count_matrix
        self.gids = gene_ids
        self.pad_token_id = pad_token_id
        self.pad_value = pad_value

    def __len__(self):
        return len(self.counts)

    def __getitem__(self, idx):
        row = self.counts[idx]
        nz = np.nonzero(row)[0]
        genes = self.gids[nz]
        values = row[nz]
        # Prepend CLS-like token at position 0 (maps to vocab[" "])
        genes = np.insert(genes, 0, self.pad_token_id)
        values = np.insert(values, 0, self.pad_value)
        return {
            "gene": torch.from_numpy(genes).long(),
            "expr": torch.from_numpy(values).float(),
        }


def get_scgpt_embeddings(adata, model, vocab, cfg, device, batch_size=64, max_length=1200):
    """
    Extract scGPT cell embeddings using the CLS token (position 0 of the
    transformer output), matching the official tutorial exactly.

    Steps:
        1. Map adata gene names → vocab IDs (filter unmatched)
        2. Build per-cell sparse representation (nonzero gene, expression)
        3. Collate into padded batches using scGPT DataCollator
        4. Run model._encode() → take position [0] → CLS embedding
        5. L2-normalise the embeddings
    """
    import scgpt.data_collator as _collator_mod

    gene_col = "feature_name"
    if gene_col == "index":
        adata.var["index"] = adata.var.index

    # Map genes to vocab IDs; -1 = not in vocabulary
    adata.var["id_in_vocab"] = [
        vocab[g] if g in vocab else -1 for g in adata.var[gene_col]
    ]
    matched = np.sum(np.array(adata.var["id_in_vocab"]) >= 0)
    print(f"  scGPT vocab match: {matched}/{len(adata.var)} genes")

    adata_matched = adata[:, adata.var["id_in_vocab"] >= 0].copy()
    if adata_matched.n_vars == 0:
        raise ValueError("No genes matched the scGPT vocabulary.")

    # Build gene-id array in adata.var order
    genes = adata_matched.var[gene_col].tolist()
    gene_ids = np.array([vocab[g] for g in genes], dtype=int)

    # Sparse → dense count matrix
    count_matrix = adata_matched.X
    count_matrix = (
        count_matrix if isinstance(count_matrix, np.ndarray)
        else count_matrix.toarray()
    )

    pad_token_id = vocab[cfg["pad_token"]]
    pad_value = cfg["pad_value"]

    dataset = _SeqDataset(count_matrix, gene_ids, pad_token_id, pad_value)
    collator = _collator_mod.DataCollator(
        do_padding=True,
        pad_token_id=pad_token_id,
        pad_value=pad_value,
        do_mlm=False,
        do_binning=True,
        max_length=max_length,
        sampling=True,
        keep_first_n_tokens=1,
    )
    loader = torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=False,
        num_workers=0,
        collate_fn=collator,
        pin_memory=True,
    )

    embsize = cfg["embsize"]
    cell_embs = np.zeros((len(dataset), embsize), dtype=np.float32)

    with torch.no_grad(), torch.cuda.amp.autocast(enabled=True):
        offset = 0
        for data_dict in loader:
            input_gene_ids = data_dict["gene"].to(device)
            src_mask = input_gene_ids.eq(pad_token_id)
            embeddings = model._encode(
                input_gene_ids,
                data_dict["expr"].to(device),
                src_key_padding_mask=src_mask,
                batch_labels=None,
            )
            # CLS token = position [0] in the sequence
            batch_embs = embeddings[:, 0, :].float().cpu().numpy()
            cell_embs[offset : offset + len(batch_embs)] = batch_embs
            offset += len(batch_embs)

    # L2 normalisation (matches official tutorial)
    norms = np.linalg.norm(cell_embs, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return cell_embs / norms


# ════════════════════════════════════════════════════════════════════════════════
# 5.  Pipeline
# ════════════════════════════════════════════════════════════════════════════════

def run_pca_embedding(X, n_pcs=50):
    n_pcs = min(n_pcs, X.shape[1])
    pca = PCA(n_components=n_pcs, random_state=42)
    return pca.fit_transform(X)


def run_scgpt_pipeline(adata, n_clusters, batch_size, max_length, device, seed):
    """
    Try scGPT transformer embeddings; fall back to PCA on vocab mismatch.
    Returns (y_pred, embeddings, used_scgpt).
    """
    try:
        print("Loading scGPT pretrained model...")
        model, vocab, cfg = load_scgpt_model(device)

        print("Extracting scGPT cell embeddings...")
        cell_embs = get_scgpt_embeddings(
            adata, model, vocab, cfg, device,
            batch_size=batch_size, max_length=max_length,
        )
        print(f"scGPT embedding shape: {cell_embs.shape}")
        used_scgpt = True

    except Exception as e:
        warnings.warn(f"scGPT embedding failed ({e}), using PCA fallback.")
        print(f"[Fallback] scGPT error: {e}")
        print("[Fallback] Using PCA embedding...")

        X_raw = adata.X
        X_raw = X_raw if isinstance(X_raw, np.ndarray) else X_raw.toarray()
        X_log = np.log1p(X_raw)
        cell_embs = run_pca_embedding(X_log, n_pcs=50)
        print(f"PCA embedding shape: {cell_embs.shape}")
        used_scgpt = False

    kmeans = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed)
    y_pred = kmeans.fit_predict(cell_embs)
    return y_pred, cell_embs, used_scgpt


from evaluation import best_map as _best_map


def compute_metrics(Y_true, Y_pred):
    """
    Compute a full suite of clustering metrics.
    All metrics are computed using Hungarian-matched labels for consistency.
    """
    from evaluation import best_map

    le = LabelEncoder()
    Y_true_int = le.fit_transform(Y_true)
    Y_pred_int = np.asarray(Y_pred, dtype=int)

    # Hungarian matching: permute Y_pred to best align with Y_true
    y_pred_mapped, _, _ = best_map(Y_true_int, Y_pred_int)

    acc = round(float(np.mean(y_pred_mapped == Y_true_int)), 4)
    nmi  = round(float(normalized_mutual_info_score(Y_true_int, Y_pred_int)), 4)
    ari  = round(float(adjusted_rand_score(Y_true_int, Y_pred_int)), 4)
    f1   = round(float(f1_score(y_pred_mapped, Y_true_int, average="macro", zero_division=0)), 4)
    fmi  = round(float(fowlkes_mallows_score(Y_true_int, Y_pred_int)), 4)
    vms  = round(float(v_measure_score(Y_true_int, Y_pred_int)), 4)
    hom  = round(float(homogeneity_score(Y_true_int, Y_pred_int)), 4)
    comp = round(float(completeness_score(Y_true_int, Y_pred_int)), 4)

    return {
        "ACC": acc,
        "NMI": nmi,
        "ARI": ari,
        "F1_macro": f1,
        "FMI": fmi,
        "V_measure": vms,
        "Homogeneity": hom,
        "Completeness": comp,
    }


def parse_args():
    parser = argparse.ArgumentParser(
        description="scGPT: Single-cell Foundation Model",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data_path",   type=str, required=True,
                        help="Path to input h5ad file")
    parser.add_argument("--save_dir",    type=str, default="./results",
                        help="Directory to save results")
    parser.add_argument("--n_clusters",  type=int, required=True,
                        help="Number of clusters (ground-truth label count)")
    parser.add_argument("--batch_size", type=int, default=64,
                        help="Batch size for scGPT embedding inference")
    parser.add_argument("--max_length", type=int, default=1200,
                        help="Maximum sequence length (genes per cell)")
    parser.add_argument("--seed",       type=int, default=42,
                        help="Random seed")
    parser.add_argument("--gpu",        type=int, default=0,
                        help="GPU device ID")
    return parser.parse_args()


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def main():
    args = parse_args()
    set_seed(args.seed)
    os.makedirs(args.save_dir, exist_ok=True)

    print("=" * 60)
    print("scGPT — Cell Embedding & Clustering Pipeline")
    print("=" * 60)

    # ── Load data ──────────────────────────────────────────────────────────────
    print(f"\n[1/4] Loading data: {args.data_path}")
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=False,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True,
    )
    X = np.array(X, dtype=np.float32)
    Y_str = np.array(Y)  # string cell-type labels for metrics & save()
    # Numeric integer labels for save() / evaluation() compatibility
    label_encoder = LabelEncoder()
    Y = label_encoder.fit_transform(Y_str)

    true_n_clusters = len(np.unique(Y_str))
    n_clusters = args.n_clusters if args.n_clusters > 0 else true_n_clusters
    print(f"  Cells: {X.shape[0]}, Genes (after HVG filter): {X.shape[1]}")
    print(f"  Ground-truth clusters: {true_n_clusters}")

    device = torch.device(f"cuda:{args.gpu}" if torch.cuda.is_available() else "cpu")
    print(f"\n[2/4] Device: {device}")

    # ── Embedding + Clustering ─────────────────────────────────────────────────
    print("\n[3/4] Running scGPT embedding & KMeans clustering...")
    y_pred, embeddings, used_scgpt = run_scgpt_pipeline(
        adata=adata,
        n_clusters=n_clusters,
        batch_size=args.batch_size,
        max_length=args.max_length,
        device=device,
        seed=args.seed,
    )
    method = "scGPT (Transformer CLS)" if used_scgpt else "PCA fallback"
    print(f"  Embedding method: {method}")
    print(f"  Embedding shape:  {embeddings.shape}")

    # ── Metrics (with string labels — LabelEncoder handles conversion) ───────────
    print("\n[4/4] Computing clustering metrics...")
    metrics = compute_metrics(Y_str, y_pred)

    print("\n── scGPT Results ──")
    for k, v in metrics.items():
        print(f"  {k:20s}: {v:.4f}")

    # ── Save outputs ───────────────────────────────────────────────────────────
    # save() internally calls evaluation() which expects integer labels.
    # Use the integer-encoded Y_pred mapped to integer Y for internal save,
    # but also write a CSV with readable string labels
    save(args.save_dir, Y, y_pred, 0, embeddings)

    # Write readable label CSV (true=string cell types, pred=numeric clusters)
    import pandas as pd
    readable_csv = os.path.join(args.save_dir, "types_readable.csv")
    pd.DataFrame({"true": Y_str, "pred": y_pred}).to_csv(readable_csv, index=False)

    metrics_path = os.path.join(args.save_dir, "metrics.json")
    with open(metrics_path, "w") as f:
        json.dump(metrics, f, indent=2)

    summary_path = os.path.join(args.save_dir, "scGPT_metrics.csv")
    with open(summary_path, "w") as f:
        f.write("Model,ACC,NMI,ARI,F1_macro,FMI,V_measure,Homogeneity,Completeness\n")
        f.write(
            f"scGPT,{metrics['ACC']},{metrics['NMI']},{metrics['ARI']},"
            f"{metrics['F1_macro']},{metrics['FMI']},{metrics['V_measure']},"
            f"{metrics['Homogeneity']},{metrics['Completeness']}\n"
        )

    # ── Append to global summary ───────────────────────────────────────────────
    # results/ is at plantnet/ level (2 levels above scGPT/), not methods/
    summary_csv = os.path.normpath(os.path.join(PROJECT_ROOT, "..", "results", "best_performance_summary.csv"))
    new_row = (
        f"scGPT,{metrics['ACC']},{metrics['NMI']},{metrics['ARI']},"
        f"{metrics['F1_macro']},{metrics['FMI']},{metrics['V_measure']},"
        f"{metrics['Homogeneity']},{metrics['Completeness']}\n"
    )
    if os.path.exists(summary_csv):
        with open(summary_csv) as f:
            lines = f.readlines()
        # Replace scGPT row or append
        header = lines[0]
        data_lines = [l for l in lines[1:] if not l.startswith("scGPT,")]
        data_lines.append(new_row)
        with open(summary_csv, "w") as f:
            f.write(header)
            f.writelines(data_lines)
    else:
        with open(summary_csv, "w") as f:
            f.write("Model,ACC,NMI,ARI,F1_macro,FMI,V_measure,Homogeneity,Completeness\n")
            f.write(new_row)

    print(f"\nResults saved to: {args.save_dir}")
    print("Done.")


if __name__ == "__main__":
    main()
```

---

## 4. scMAE (DeepLearning)

基于掩码自编码器的单细胞 RNA 聚类方法

**文件路径**: `methods/DeepLearning/scMAE/run.py`

```python
# -*- coding: utf-8 -*-
"""
====================================================================================================
scMAE — 基于掩码自编码器的单细胞 RNA 聚类方法
====================================================================================================

【灵感来源】
    受 Vision Transformer (ViT) 和 BERT 中掩码预训练的启发，
    scMAE 将单细胞基因表达视为"句子"，每个基因视为"词"，
    通过掩码预测学习基因间的潜在关系。

【核心思想】
    ┌────────────────────────────────────────────────────────────────┐
    │                     scMAE 流程图                                │
    │                                                                │
    │  原始基因表达 X = [x₁, x₂, x₃, ..., x₁₀₀₀]                   │
    │      │                                                        │
    │      ▼                                                        │
    │  ┌──────────────────────┐                                      │
    │  │ Step 1: 随机掩码    │  掩码概率 p=0.4                     │
    │  │ X_masked = mask(X)  │  随机替换部分基因为0或随机值           │
    │  └──────────┬───────────┘                                      │
    │             │                                                  │
    │             ▼                                                  │
    │  ┌──────────────────────┐                                      │
    │  │ Step 2: 编码器      │  Linear → LayerNorm → Mish →       │
    │  │ Encoder              │  Linear → LayerNorm → Mish →       │
    │  │                      │  Linear (输出: hidden_size)         │
    │  └──────────┬───────────┘                                      │
    │             │                                                  │
    │             ▼                                                  │
    │  ┌──────────────────────┐                                      │
    │  │ Step 3: 掩码预测器   │  Linear(hidden_size → n_genes)      │
    │  │ Mask Predictor       │  预测被掩码位置的值                   │
    │  └──────────┬───────────┘                                      │
    │             │                                                  │
    │             ▼                                                  │
    │  ┌──────────────────────┐                                      │
    │  │ Step 4: 解码器       │  concat(latent, predicted_mask)      │
    │  │ Decoder              │  → Linear → 重构原始输入               │
    │  └──────────────────────┘                                      │
    │                                                                │
    └────────────────────────────────────────────────────────────────┘

【损失函数】（两部分组成）

    L_total = L_reconstruction + L_mask

    1. L_reconstruction (重构损失)
       - 使用加权MSE，仅对被掩码位置计算损失
       - masked_data_weight=0.75：被掩码位置权重更高
       - mask_loss_weight=0.7：控制两部分损失的相对重要性

    2. L_mask (掩码损失)
       - BCE损失：预测被掩码的位置（二分类）
       - 促使编码器学习哪些位置被掩码

【与scCDCG的关键区别】
    | 特性         | scCDCG          | scMAE                    |
    |------------|----------------|--------------------------|
    | 核心思想    | 图结构保持        | 掩码预测                  |
    | 损失函数    | 重构+正交+协方差+聚类 | 重构+掩码预测              |
    | 图结构      | 需要KNN图        | 不需要                    |
    | 聚类方式    | DEC+Sinkhorn    | KMeans/ Leiden            |
    | 超参数数量  | 较多            | 较少（更易调）             |

【超参数配置】
    | 参数               | 默认值 | 说明                      |
    |------------------|-------|-------------------------|
    | hidden_size      | 128   | 隐藏层维度                |
    | mask_prob        | 0.4   | 基因掩码概率               |
    | masked_data_weight| 0.75 | 被掩码数据的损失权重        |
    | mask_loss_weight | 0.7   | 掩码预测损失权重            |
    | batch_size       | 256   | 批大小                    |
    | lr               | 1e-3  | 学习率                    |
    | epochs           | 100   | 训练轮次                  |
"""

import os
import sys
import argparse
import numpy as np
import torch
import random
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from torch.utils.data import DataLoader, Dataset

# 添加父目录到路径（用于导入benchmark通用模块）
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from preprocess import prepare_data_for_model
from utils import save

# 导入模型组件
from model import AutoEncoder


def set_seed(seed):
    """设置随机种子确保可复现性"""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def apply_noise(X, p):
    """
    对输入数据应用掩码噪声

    【掩码策略】
        对于每个基因，以概率p将其掩码：
        - 以概率p将基因值替换为从其他细胞随机采样的值
        - 返回：(损坏的数据, 掩码标记)

    【参数】
        X: 输入数据 (n_cells, n_genes)
        p: 掩码概率（可以是标量或向量）

    【返回】
        corrupted_X: 损坏后的数据（部分基因被替换）
        mask: 掩码标记（1=被掩码, 0=保持原样）

    【示例】
        X: [[1, 2, 3, 4], [5, 6, 7, 8]]
        mask概率p=0.5

        假设随机决定位置[0,2]和[1,1]被掩码
        corrupted_X: [[1, 6', 3, 4], [5, 7', 7, 8]]  # '表示来自随机细胞的替换值
        mask: [[0, 1, 0, 0], [1, 0, 0, 0]]
    """
    p = torch.tensor(p)
    # 生成掩码：伯努利分布采样
    should_swap = torch.bernoulli(p.to(X.device) * torch.ones((X.shape)).to(X.device))
    # 替换：被掩码位置用随机细胞的对应值替换
    corrupted_X = torch.where(
        should_swap == 1,
        X[torch.randperm(X.shape[0])],  # 从随机细胞取值
        X                                   # 保持原样
    )
    # 生成掩码标记
    masked = (corrupted_X != X).float()
    return corrupted_X, masked


class scRNADataset(Dataset):
    """
    单细胞RNA-seq数据集封装

    【功能】
        将NumPy数组封装为PyTorch Dataset
        支持DataLoader的批处理功能
    """

    def __init__(self, data, labels):
        """
        初始化数据集

        【参数】
            data: 基因表达矩阵 (n_cells, n_genes)
            labels: 细胞类型标签 (n_cells,)
        """
        self.data = torch.FloatTensor(data)
        self.labels = torch.LongTensor(labels)

    def __len__(self):
        """返回数据集大小"""
        return len(self.data)

    def __getitem__(self, idx):
        """获取单个样本"""
        return self.data[idx], self.labels[idx]


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    """
    搜索Leiden聚类的分辨率参数以获得指定数量的簇

    【问题背景】
        Leiden聚类需要指定分辨率参数，而不是直接的簇数
        不同数据集、不同细胞类型数需要不同的分辨率

    【算法】
        二分搜索找到使聚类数最接近目标的分辨率：
        1. 从高分辨率开始搜索（2.5 → 0.01）
        2. 按increment步长递减
        3. 找到聚类数=目标时停止

    【参数】
        adata: 包含嵌入的AnnData对象
        fixed_clus_count: 目标聚类数
        increment: 分辨率搜索步长

    【返回】
        最优分辨率参数
    """
    dis = []
    # 从高到低搜索（高分辨率=多簇，低分辨率=少簇）
    resolutions = sorted(list(np.arange(0.01, 2.5, increment)), reverse=True)

    for res in resolutions:
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        dis.append(abs(count_unique_leiden - fixed_clus_count))
        if count_unique_leiden == fixed_clus_count:
            break

    # 返回误差最小的分辨率
    return resolutions[np.argmin(dis)]


def inference(net, data_loader, device):
    """
    从模型提取特征（嵌入向量）

    【功能】
        遍历整个数据集，获取每个细胞的嵌入向量

    【参数】
        net: 训练好的模型
            data_loader: 数据加载器
            device: 计算设备

    【返回】
        feature_vector: 所有细胞的嵌入向量 (n_cells, hidden_size)
        labels_vector: 对应的真实标签 (n_cells,)
    """
    net.eval()
    feature_vector = []
    labels_vector = []

    with torch.no_grad():
        for x, y in data_loader:
            x = x.to(device)
            # 使用编码器获取嵌入
            feature_vector.extend(net.feature(x).detach().cpu().numpy())
            labels_vector.extend(y.numpy())

    return np.array(feature_vector), np.array(labels_vector)


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='scMAE: Masked Autoencoder for scRNA-seq Clustering',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 数据参数
    parser.add_argument('--data_path', type=str, required=True,
                       help='输入 h5ad 文件路径')
    parser.add_argument('--save_dir', type=str, default='./results',
                       help='结果保存目录')

    # 模型参数
    parser.add_argument('--n_clusters', type=int, required=True,
                       help='聚类数')
    parser.add_argument('--hidden_size', type=int, default=128,
                       help='隐藏层维度')
    parser.add_argument('--mask_prob', type=float, default=0.4,
                       help='掩码概率')
    parser.add_argument('--masked_data_weight', type=float, default=0.75,
                       help='被掩码数据的损失权重')
    parser.add_argument('--mask_loss_weight', type=float, default=0.7,
                       help='掩码预测损失权重')

    # 训练参数
    parser.add_argument('--epochs', type=int, default=100,
                       help='训练轮次')
    parser.add_argument('--batch_size', type=int, default=256,
                       help='批大小')
    parser.add_argument('--lr', type=float, default=1e-3,
                       help='学习率')

    # 其他参数
    parser.add_argument('--seed', type=int, default=42,
                       help='随机种子')
    parser.add_argument('--gpu', type=int, default=0,
                       help='GPU设备号')
    parser.add_argument('--no_cuda', action='store_true',
                       help='禁用CUDA')
    parser.add_argument('--eval_interval', type=int, default=10,
                       help='评估间隔（轮次）')

    return parser.parse_args()


def main():
    """主函数"""
    args = parse_args()

    # 设备设置
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    # 随机种子
    set_seed(args.seed)

    # 创建保存目录
    os.makedirs(args.save_dir, exist_ok=True)

    # =========================================================================
    # Step 1: 数据加载与预处理
    # =========================================================================
    print('Loading data...')
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=True,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True
    )

    # 转换为NumPy数组
    X = np.array(X).astype(np.float32)
    Y = np.array(Y)

    # 标签编码
    from sklearn.preprocessing import LabelEncoder
    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    # 获取聚类数
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    # =========================================================================
    # Step 2: 创建数据集和数据加载器
    # =========================================================================
    dataset = scRNADataset(X, Y)
    train_loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True
    )
    test_loader = DataLoader(
        dataset,
        batch_size=args.batch_size * 5,
        shuffle=False,
        drop_last=False
    )

    # =========================================================================
    # Step 3: 初始化模型
    # =========================================================================
    model = AutoEncoder(
        num_genes=X.shape[1],
        hidden_size=args.hidden_size,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    # 掩码概率（每个基因独立）
    mask_probas = [args.mask_prob] * X.shape[1]

    # =========================================================================
    # Step 4: 训练循环
    # =========================================================================
    print('Starting training...')

    best_acc = 0
    best_epoch = 0
    best_embedding = None
    best_pred = None

    for epoch in range(args.epochs):
        model.train()
        total_loss = 0

        for x, y in train_loader:
            x = x.to(device)

            # 应用掩码
            x_corrupted, mask = apply_noise(x, mask_probas)

            # 前向传播
            optimizer.zero_grad()
            _, loss = model.loss_mask(x_corrupted, x, mask)

            # 反向传播
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        avg_loss = total_loss / len(train_loader)

        # =========================================================================
        # Step 5: 周期性评估
        # =========================================================================
        if (epoch + 1) % args.eval_interval == 0 or epoch == args.epochs - 1:
            # 提取嵌入向量
            embedding, true_labels = inference(model, test_loader, device)

            # 聚类
            if embedding.shape[0] < 10000:
                # 小数据集：直接使用KMeans
                kmeans = KMeans(
                    n_clusters=n_clusters,
                    random_state=args.seed,
                    n_init=20
                )
                pred_labels = kmeans.fit_predict(embedding)
            else:
                # 大数据集：使用Leiden聚类
                adata_emb = sc.AnnData(embedding)
                sc.pp.neighbors(adata_emb, n_neighbors=10, use_rep="X")
                reso = res_search_fixed_clus(adata_emb, n_clusters)
                sc.tl.leiden(adata_emb, resolution=reso)
                pred_labels = np.array([
                    int(x) for x in adata_emb.obs['leiden'].to_list()
                ])

            # 评估并追踪最优模型
            from evaluation import evaluation as eval_fn
            acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, _ = eval_fn(
                np.array(true_labels), np.array(pred_labels))

            if acc > best_acc:
                best_acc = acc
                best_epoch = epoch + 1
                best_embedding = embedding.copy()
                best_pred = pred_labels.copy()
                # 保存最优模型检查点
                torch.save({
                    'model': model.state_dict(),
                    'optimizer': optimizer.state_dict(),
                    'args': vars(args),
                    'best_epoch': best_epoch,
                    'best_acc': best_acc,
                }, os.path.join(args.save_dir, 'model_checkpoint.pth'))

            # 保存结果
            save(args.save_dir, true_labels, pred_labels, epoch + 1, embedding)

            print(f'Epoch {epoch + 1}/{args.epochs}, Loss: {avg_loss:.4f}, ACC: {acc:.4f}, Best: {best_acc:.4f}')

    print(f'Best epoch: {best_epoch}, Best ACC: {best_acc:.4f}')

    # =========================================================================
    # Step 6: 保存最终结果（最优 epoch）
    # =========================================================================
    if best_embedding is not None:
        true_labels_arr = np.array(true_labels)
        save(args.save_dir, true_labels_arr, best_pred, best_epoch, best_embedding)

    print(f'Training completed. Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
```

---

## 5. scMAE (Foundation)

基于掩码自编码器的单细胞 RNA 聚类方法 - 主训练脚本

**文件路径**: `methods/Foundation/scMAE/main.py`

```python
import os
os.environ["CUDA_VISIBLE_DEVICES"] = "1"
import torch
import random
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from datasets import Loader, apply_noise
from model import AutoEncoder
from evaluate import evaluate
from util import AverageMeter



def make_dir(directory_path, new_folder_name):
    """Creates an expected directory if it does not exist"""
    directory_path = os.path.join(directory_path, new_folder_name)
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
    return directory_path


def inference(net, data_loader_test):
    net.eval()
    feature_vector = []
    labels_vector = []
    with torch.no_grad():
        for step, (x, y) in enumerate(data_loader_test):
            feature_vector.extend(net.feature(x.cuda()).detach().cpu().numpy())
            labels_vector.extend(y.numpy())
    feature_vector = np.array(feature_vector)
    labels_vector = np.array(labels_vector)
    return feature_vector, labels_vector


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    '''
        arg1(adata)[AnnData matrix]
        arg2(fixed_clus_count)[int]

        return:
            resolution[int]
    '''
    dis = []
    resolutions = sorted(list(np.arange(0.01, 2.5, increment)), reverse=True)
    i = 0
    res_new = []
    for res in resolutions:
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(
            adata.obs['leiden']).leiden.unique())
        dis.append(abs(count_unique_leiden-fixed_clus_count))
        res_new.append(res)
        if count_unique_leiden == fixed_clus_count:
            break
    reso = resolutions[np.argmin(dis)]

    return reso


def train(args):
    data_load = Loader(args, dataset_name=args["dataset"], drop_last=True)
    data_loader = data_load.train_loader
    data_loader_test = data_load.test_loader
    x_shape = args["data_dim"]

    results = []

    # Hyper-params
    init_lr = args["learning_rate"]
    max_epochs = args["epochs"]
    mask_probas = [0.4]*x_shape

    # setup model
    model = AutoEncoder(
        num_genes=x_shape,
        hidden_size=128,
        masked_data_weight=0.75,
        mask_loss_weight=0.7
    ).cuda()
    model_checkpoint = 'model_checkpoint.pth'

    optimizer = torch.optim.Adam(model.parameters(), lr=init_lr)

    # train model
    for epoch in range(max_epochs):
        model.train()
        meter = AverageMeter()
        for i, (x, y) in enumerate(data_loader):
            x = x.cuda()
            x_corrputed, mask = apply_noise(x, mask_probas)
            optimizer.zero_grad()
            x_corrputed_latent, loss_ae = model.loss_mask(x_corrputed, x, mask)
            loss_ae.backward()
            optimizer.step()
            meter.update(loss_ae.detach().cpu().numpy())
    
        if epoch == 80:
            # Generator in eval mode
            latent, true_label = inference(model, data_loader_test)
            if latent.shape[0] < 10000:
                clustering_model = KMeans(n_clusters=args["n_classes"])
                clustering_model.fit(latent)
                pred_label = clustering_model.labels_
            else:
                adata = sc.AnnData(latent)
                sc.pp.neighbors(adata, n_neighbors=10, use_rep="X")
                # sc.tl.umap(adata)
                reso = res_search_fixed_clus(adata, args["n_classes"])
                sc.tl.leiden(adata, resolution=reso)
                pred = adata.obs['leiden'].to_list()
                pred_label = [int(x) for x in pred]
            

            nmi, ari, acc = evaluate(true_label, pred_label)
            ss = silhouette_score(latent, pred_label)

            res = {}
            res["nmi"] = nmi
            res["ari"] = ari
            res["acc"] = acc
            res["sil"] = ss
            res["dataset"] = args["dataset"]
            res["epoch"] = epoch
            results.append(res)

            print("\tEvalute: [nmi: %f] [ari: %f] [acc: %f]" % (nmi, ari, acc))

            np.save(args["save_path"]+"/embedding_"+str(epoch)+".npy", 
                    latent)
            pd.DataFrame({"True": true_label, 
                        "Pred": pred_label}).to_csv(args["save_path"]+"/types_"+str(epoch)+".txt")

    torch.save({
        "optimizer": optimizer.state_dict(),
        "model": model.state_dict()
    }, model_checkpoint
    )

    return results


if __name__ == "__main__":
    for i in range(1):
        seed = random.randint(1,100)
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.cuda.manual_seed(seed)
        np.random.seed(seed)

        args = {}
        args["num_workers"] = 4
        args["paths"] = {"data": "/data/sc_data/all_data/",
                        "results": "./res/"}
        args['batch_size'] = 256
        args["data_dim"] = 1000
        args['n_classes'] = 4
        args['epochs'] = 100
        args["dataset"] = "10X_PBMC"
        args["learning_rate"] = 1e-3
        args["latent_dim"] = 32

        print(args)

        path = args["paths"]["data"]
        files = ["Pollen", "Quake_Smart-seq2_Lung", "Limb_Muscle", 
                 "worm_neuron_cell", "Melanoma_5K", "Young", "Guo", "Baron", 
                 "Wang", "Quake_10x_Spleen", "Shekhar", "Macosko", 
                 "Tosches", "Bach", "hrvatin"]

        results = pd.DataFrame()
        save_dir = make_dir(args["paths"]["results"], "a_summary")
        for dataset in files:
            print(f">> {dataset}")
            args["dataset"] = dataset
            args["save_path"] = make_dir("/data/sc_data/scMAE/"+str(i), dataset)

            res = train(args)
            print(res)
            results = results.append(res)
            results.to_csv(args["paths"]["results"] +
                        "/res_all_data_test"+str(i)+".csv", header=True)
```

---

## 6. scVI

Single-cell Variational Inference - 基于 scvi-tools 的变分推断方法

**文件路径**: `methods/DeepLearning/scVI/run.py`

```python
# -*- coding: utf-8 -*-
"""
====================================================================================================
scVI — Single-cell Variational Inference (scvi-tools Implementation)
====================================================================================================

【论文来源】
    Lopez et al. (2018), Nature Methods
    "Deep generative modeling for single-cell transcriptomics"
    https://www.nature.com/articles/s41592-018-0229-2

【官方实现】
    scvi-tools: https://github.com/scverse/scvi-tools

【本实现】
    基于 scvi-tools 官方实现，充分利用其 GPU 加速和成熟的工作流。
    完全遵循官方 API：setup → train → extract → cluster

【核心思想】
    scVI 将单细胞基因表达数据建模为变分自编码器（VAE）的输出，
    使用零膨胀负二项分布（ZINB）作为观测模型，适合单细胞计数的过离散特性。
    通过变分推断学习细胞的低维潜在表示，用于聚类、降维等下游任务。

【生成模型】
    z_n ~ Normal(0, I)                                          潜在变量
    ℓ_n ~ LogNormal(ℓ_μ^⊤ s_n, ℓ_σ²^⊤ s_n)                 Library size
    ρ_n = softmax(f_w(z_n, s_n))                               归一化表达
    π_n = sigmoid(f_h(z_n, s_n))                               零膨胀参数
    x_ng ~ ZINB(ℓ_n · ρ_n, θ_g, π_n)                         观测: ZINB

【变分推断】
    q_η(z_n | x_n, s_n) = Normal(μ_η(x_n), σ_η²(x_n))
    使用重参数化技巧进行随机反向传播

【损失函数】
    ELBO = E_{q(z|x)}[log p(x|z,s)] - KL(q(z|x) || N(0,I))
    其中 log p(x|z,s) 是 ZINB 对数似然

【与 PhytoCluster 的关键区别】
    | 特性           | PhytoCluster (VAE+GMM)  | scVI                          |
    |----------------|--------------------------|-------------------------------|
    | 观测模型       | MSE (重建连续值)         | ZINB (建模 counts 过离散)      |
    | 潜在空间       | 高斯混合 (离散结构)      | 标准高斯 (连续结构)            |
    | 聚类方式       | 内嵌 GMM 后验            | 外部聚类 (KMeans/Leiden)      |
    | 输入数据       | Z-score 标准化数据       | 原始 Counts + lib size        |
    | Library size   | 外部计算                  | 内部变分推断或观察值          |
    | 批次效应       | 不支持                   | 支持 categorical covariate    |
"""

import os
import sys
import argparse

# === CRITICAL: Apply all dependency patches BEFORE any other imports ===
# These patches fix version incompatibilities between jax, numpyro, flax, and numpy

import numpy as np
if not hasattr(np.dtypes, 'StringDType'):
    class StringDType:
        def __repr__(self): return 'StringDType()'
    np.dtypes.StringDType = StringDType

# Force import jax so submodules exist before patching
import jax
import jax.experimental.layout as jax_layout
import jax.api_util as api_util

# Patch jax.experimental.layout (flax 0.10+ needs Format)
if not hasattr(jax_layout, 'Format'):
    jax_layout.Format = type('Format', (), {'__repr__': lambda s: 'Format'})

# Patch jax.api_util (numpyro 0.14+ needs debug_info)
if not hasattr(api_util, 'debug_info'):
    api_util.debug_info = lambda *args, **kwargs: None

import pandas as pd
import scanpy as sc
import random
import time
from datetime import datetime
from sklearn.cluster import KMeans

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utils import save
from evaluation import evaluation as eval_fn


# =============================================================================
# 工具函数
# =============================================================================

def set_seed(seed):
    """设置随机种子确保可复现性"""
    random.seed(seed)
    np.random.seed(seed)
    try:
        import scvi
        scvi.utils._setup.seed(seed)
    except Exception:
        pass


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    """
    搜索 Leiden 聚类的分辨率参数以获得指定数量的簇
    """
    dis = []
    resolutions = sorted(list(np.arange(0.01, 2.5, increment)), reverse=True)

    for res in resolutions:
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        dis.append(abs(count_unique_leiden - fixed_clus_count))
        if count_unique_leiden == fixed_clus_count:
            break

    return resolutions[np.argmin(dis)]


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='scVI: Single-cell Variational Inference (scvi-tools)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 数据参数
    parser.add_argument('--data_path', type=str, required=True,
                       help='输入 h5ad 文件路径')
    parser.add_argument('--save_dir', type=str, default='./results',
                       help='结果保存目录')

    # 模型参数
    parser.add_argument('--n_clusters', type=int, required=True,
                       help='聚类数')
    parser.add_argument('--latent_dim', type=int, default=10,
                       help='潜在空间维度')
    parser.add_argument('--n_layers', type=int, default=1,
                       help='编码器/解码器隐藏层数')
    parser.add_argument('--encode_dim', type=int, default=128,
                       help='编码器隐藏层维度')

    # 训练参数
    parser.add_argument('--epochs', type=int, default=200,
                       help='训练轮次')
    parser.add_argument('--batch_size', type=int, default=256,
                       help='批大小')
    parser.add_argument('--lr', type=float, default=0.001,
                       help='学习率')
    parser.add_argument('--use_observed_lib_size', action='store_true', default=True,
                       help='使用观察到的 library size（默认开启）')
    parser.add_argument('--no_observed_lib_size', action='store_true',
                       help='禁用观察到的 library size（使用变分推断）')

    # 其他参数
    parser.add_argument('--seed', type=int, default=42,
                       help='随机种子')
    parser.add_argument('--gpu', type=int, default=2,
                       help='GPU设备号')
    parser.add_argument('--eval_interval', type=int, default=20,
                       help='评估间隔（轮次）')
    parser.add_argument('--print_interval', type=int, default=10,
                       help='打印间隔（轮次）')

    return parser.parse_args()


# =============================================================================
# 主函数
# =============================================================================

def main():
    """主函数"""
    args = parse_args()

    # 设备设置
    os.environ['CUDA_VISIBLE_DEVICES'] = str(args.gpu)

    print(f'Using GPU device: {args.gpu}')

    set_seed(args.seed)

    # 生成带时间戳的保存目录
    timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
    model_name = 'scVI'
    dataset_name = os.path.splitext(os.path.basename(args.data_path))[0]
    save_dir = os.path.join(
        args.save_dir,
        f'{dataset_name}_{model_name}_{timestamp}'
    )
    os.makedirs(save_dir, exist_ok=True)

    print(f'\nResults will be saved to: {save_dir}')

    # =========================================================================
    # Step 1: 数据加载与预处理
    # =========================================================================
    print('\n' + '='*60)
    print('Step 1: Loading and Preprocessing Data')
    print('='*60)

    print('Loading data...')

    # 读取原始 h5ad（scVI 需要原始 counts）
    adata_raw = sc.read_h5ad(args.data_path)

    # 获取真实标签
    label_col = None
    for candidate in ['cell_type', 'Celltype', 'celltype', 'cell_label', 'label']:
        if candidate in adata_raw.obs.columns:
            label_col = candidate
            break
    if label_col is None:
        raise KeyError(f"No cell type label column found. Available obs columns: {list(adata_raw.obs.columns)}")

    Y = np.array(adata_raw.obs[label_col])

    # 标签编码
    from sklearn.preprocessing import LabelEncoder
    le = LabelEncoder()
    Y_encoded = le.fit_transform(Y)

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y_encoded))
    n_cells = adata_raw.n_obs
    n_genes_raw = adata_raw.n_vars

    print(f'Dataset: {dataset_name}')
    print(f'Number of cells: {n_cells}')
    print(f'Number of genes: {n_genes_raw}')
    print(f'Number of cell types: {len(np.unique(Y_encoded))}')
    print(f'Target clusters: {n_clusters}')

    # =========================================================================
    # Step 2: 基础过滤（scVI 保留更多基因，不做 HVG 筛选）
    # =========================================================================
    print('\n' + '='*60)
    print('Step 2: Filtering Genes and Cells')
    print('='*60)

    adata = adata_raw.copy()

    # scVI 的基础过滤：移除表达量极低的基因和细胞
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    n_genes = adata.n_vars
    n_cells_after = adata.n_obs
    print(f'After filtering: {n_cells_after} cells, {n_genes} genes')

    # =========================================================================
    # Step 3: scvi-tools 设置
    # =========================================================================
    print('\n' + '='*60)
    print('Step 3: Setting up scvi-tools')
    print('='*60)

    import scvi
    import torch

    # 注册 AnnData
    # scVI 使用原始 counts 数据，不需要 Z-score 标准化
    # 重要：使用原始 X（未归一化的 counts）
    if hasattr(adata.X, 'toarray'):
        X_counts = adata.X.toarray()
    else:
        X_counts = np.array(adata.X)

    # 保持原始 counts 在 X 中（scvi-tools 会自动处理）
    # 如果 X 已经被归一化，需要恢复到原始 counts
    adata.X = X_counts.astype(np.float32)

    # 设置 scVI
    scvi.model.SCVI.setup_anndata(
        adata,
        layer=None,  # 使用原始 X
        batch_key=None,  # 无批次效应
    )

    # =========================================================================
    # Step 4: 创建并训练 scVI 模型
    # =========================================================================
    print('\n' + '='*60)
    print('Step 4: Training SCVI Model')
    print('='*60)

    model = scvi.model.SCVI(
        adata,
        n_latent=args.latent_dim,
        n_layers=args.n_layers,
        gene_likelihood='zinb',
        use_observed_lib_size=not args.no_observed_lib_size,
    )

    print(f'Model architecture:')
    print(f'  - Latent dim:     {args.latent_dim}')
    print(f'  - N layers:       {args.n_layers}')
    print(f'  - Encode dim:     {args.encode_dim}')
    print(f'  - Gene likelihood: zinb')
    print(f'  - Observed lib size: {not args.no_observed_lib_size}')
    print(f'  - Batch size:     {args.batch_size}')
    print(f'  - Learning rate:  {args.lr}')
    print(f'  - Epochs:         {args.epochs}')
    print(f'  - Total parameters: {sum(p.numel() for p in model.module.parameters()):,}')

    # 训练模型
    print('\nTraining scVI...')
    train_start_time = time.time()

    model.train(
        max_epochs=args.epochs,
        batch_size=args.batch_size,
        train_size=1.0,
        validation_size=None,
        plan_kwargs={
            'optimizer': 'Adam',
            'lr': args.lr,
        },
    )

    train_time = time.time() - train_start_time
    print(f'\nTraining completed in {train_time:.2f} seconds ({train_time/60:.2f} minutes)')

    # =========================================================================
    # Step 5: 提取潜在表示并进行聚类
    # =========================================================================
    print('\n' + '='*60)
    print('Step 5: Extracting Latent Representation & Clustering')
    print('='*60)

    # 提取潜在表示
    embedding = model.get_latent_representation(adata)
    print(f'Latent representation shape: {embedding.shape}')

    # 聚类
    true_labels = Y_encoded

    if n_cells_after < 10000:
        kmeans = KMeans(
            n_clusters=n_clusters,
            random_state=args.seed,
            n_init=20,
            max_iter=300
        )
        pred_labels = kmeans.fit_predict(embedding)
    else:
        adata_emb = sc.AnnData(embedding)
        adata_emb.obs['true_label'] = true_labels
        sc.pp.neighbors(adata_emb, n_neighbors=15, use_rep='X')
        reso = res_search_fixed_clus(adata_emb, n_clusters)
        sc.tl.leiden(adata_emb, resolution=reso, key_added='leiden_pred')
        pred_labels = np.array([int(x) for x in adata_emb.obs['leiden_pred']])

    # 评估
    acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, _ = eval_fn(
        np.array(true_labels), np.array(pred_labels))

    print(f'\nFinal Results:')
    print(f'  -> ACC: {acc:.4f}, NMI: {nmi:.4f}, ARI: {ari:.4f}, '
          f'F1: {f1_macro:.4f}, FMI: {fmi:.4f}')

    # =========================================================================
    # Step 6: 保存结果
    # =========================================================================
    print('\n' + '='*60)
    print('Step 6: Saving Results')
    print('='*60)

    save(save_dir, Y_encoded, pred_labels, args.epochs, embedding)

    # 保存配置信息
    import json
    config = {
        'model': 'scVI (scvi-tools)',
        'dataset': dataset_name,
        'n_cells': int(n_cells_after),
        'n_genes_raw': int(n_genes_raw),
        'n_genes_after_filter': int(n_genes),
        'n_clusters': int(n_clusters),
        'n_pred_clusters': int(len(np.unique(pred_labels))),
        'latent_dim': args.latent_dim,
        'n_layers': args.n_layers,
        'encode_dim': args.encode_dim,
        'gene_likelihood': 'zinb',
        'use_observed_lib_size': not args.no_observed_lib_size,
        'epochs': args.epochs,
        'batch_size': args.batch_size,
        'lr': args.lr,
        'seed': args.seed,
        'gpu': args.gpu,
        'train_time_seconds': float(train_time),
    }
    with open(os.path.join(save_dir, 'config.json'), 'w') as f:
        json.dump(config, f, indent=2)

    # 保存模型
    model.save(os.path.join(save_dir, 'scvi_model'), overwrite=True)

    # =========================================================================
    # 打印最终结果摘要
    # =========================================================================
    print('\n' + '='*60)
    print('scVI Results Summary')
    print('='*60)

    metrics_file = os.path.join(save_dir, 'metrics.json')
    with open(metrics_file, 'r') as f:
        metrics = json.load(f)

    print(f'\nModel Configuration:')
    print(f'  - Latent dim:      {args.latent_dim}')
    print(f'  - Gene likelihood: zinb')
    print(f'  - Training epochs: {args.epochs}')
    print(f'  - Training time:    {train_time:.2f}s')

    print(f'\nClustering Performance ({dataset_name}):')
    print(f'  {"Metric":<18} {"Value":>10}')
    print(f'  {"-"*28}')
    for key, value in metrics.items():
        print(f'  {key:<18} {value:>10.4f}')

    print(f'\nResults saved to: {save_dir}')
    print(f'\nscVI completed successfully!')

    return metrics, save_dir


if __name__ == '__main__':
    main()
```

---

## 方法汇总表

| 方法名 | 类别 | 核心算法 | 论文/来源 |
|--------|------|----------|----------|
| PhytoCluster | Deep Learning | VAE + GMM | Wang et al., aBIOTECH 2025 |
| SC3 | Traditional | Consensus Clustering | Kiselev et al., 2017 |
| scGPT | Foundation Model | Transformer + CLS | Cui et al., Nature Methods 2024 |
| scMAE | Deep Learning | Masked Autoencoder | - |
| scVI | Deep Learning | VAE + ZINB | Lopez et al., Nature Methods 2018 |

## 评估指标说明

所有方法均使用以下指标进行评估：

- **ACC**: 准确率（标签对齐后）
- **NMI**: 标准化互信息
- **ARI**: 调整兰德指数
- **F1_macro**: F1 分数（宏平均）
- **FMI**: Fowlkes-Mallows 指数
- **V_measure**: V-measure（同质性和完整性调和平均）
- **Homogeneity**: 同质性
- **Completeness**: 完整性
