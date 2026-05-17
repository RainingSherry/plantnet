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
    │  Step 5: 聚类推断（使用训练后的联合模型）                       │
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
    - 联合优化后的模型直接用于最终聚类推断（不是预训练VAE+sklearn GMM）
"""

import os
import sys
import argparse
import math
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


def binary_cross_entropy_loss(recon_x, x):
    """计算二元交叉熵损失（用于损失函数中）"""
    return -torch.sum(x * torch.log(recon_x + 1e-8) + (1 - x) * torch.log(1 - recon_x + 1e-8), dim=-1)


def compute_kl_divergence(mu, logvar):
    """计算标准VAE中潜在变量的KL散度"""
    return -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)


def compute_elbo_loss(recon_x, x, z_params, binary=True):
    """计算VAE的ELBO损失"""
    mu, logvar = z_params
    kld = compute_kl_divergence(mu, logvar)
    if binary:
        likelihood = -binary_cross_entropy_loss(recon_x, x)
    else:
        likelihood = -F.mse_loss(recon_x, x, reduction='none').sum(dim=-1)
    return torch.sum(likelihood), torch.sum(kld)


def compute_elbo(recon_x, x, gamma, c_params, z_params, binary=True):
    """计算PhytoCluster的ELBO损失（包含聚类目标）"""
    mu_c, var_c, pi = c_params
    var_c += 1e-8
    n_centroids = pi.size(1)
    mu, logvar = z_params
    mu_expand = mu.unsqueeze(2).expand(mu.size(0), mu.size(1), n_centroids)
    logvar_expand = logvar.unsqueeze(2).expand(mu.size(0), mu.size(1), n_centroids)

    if binary:
        likelihood = -binary_cross_entropy_loss(recon_x, x)
    else:
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
        self.sample = GaussianSample(([x_dim] + h_dim)[-1], z_dim)

    def forward(self, x):
        x = self.hidden(x)
        return self.sample(x)


class Decoder(nn.Module):
    """解码器：从潜在变量重建原始数据"""
    def __init__(self, dims, bn=False, dropout=0, output_activation=None):
        super(Decoder, self).__init__()
        [z_dim, h_dim, x_dim] = dims
        self.hidden = create_mlp([z_dim, *h_dim], bn=bn, dropout=dropout)
        self.reconstruction = nn.Linear([z_dim, *h_dim][-1], x_dim)
        self.output_activation = output_activation

    def forward(self, x):
        x = self.hidden(x)
        if self.output_activation is not None:
            return self.output_activation(self.reconstruction(x))
        return self.reconstruction(x)


class VAE(nn.Module):
    """变分自编码器基础模型"""
    def __init__(self, dims, bn=False, dropout=0, binary=False):
        super(VAE, self).__init__()
        [x_dim, z_dim, encode_dim, decode_dim] = dims
        self.binary = binary
        if binary:
            decode_activation = nn.Sigmoid()
        else:
            decode_activation = None

        self.encoder = Encoder([x_dim, encode_dim, z_dim], bn=bn, dropout=dropout)
        self.decoder = Decoder([z_dim, decode_dim, x_dim], bn=bn, dropout=dropout,
                               output_activation=decode_activation)
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
        likelihood, kl_loss = compute_elbo_loss(recon_x, x, (mu, logvar), binary=self.binary)
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

    def encode_batch(self, dataloader, device='cpu', output_type='z', transforms=None):
        """
        批量编码接口，支持多种输出类型。

        Args:
            dataloader: 数据加载器
            device: 计算设备
            output_type: 输出类型，'z'(重参数化样本)、'mu'(均值)、'log_var'(对数方差)
            transforms: 变换函数

        Returns:
            numpy array of encoded representations
        """
        output = []
        for x in dataloader:
            if isinstance(x, (list, tuple)):
                x = x[0]
            x = x.view(x.size(0), -1).float().to(device)
            z, mu, logvar = self.encoder(x)

            if output_type == 'z':
                output.append(z.detach().cpu())
            elif output_type == 'x':
                recon_x = self.decoder(z)
                output.append(recon_x.detach().cpu().data)
            elif output_type == 'mu':
                output.append(mu.cpu().detach().data)
            elif output_type == 'log_var':
                output.append(logvar.cpu().detach().data)

        output = torch.cat(output).numpy()
        return output


class PhytoCluster(VAE):
    """PhytoCluster: VAE + GMM 联合聚类模型

    核心创新：通过 VAE 编码器提取潜在特征 z，再用可学习的 GMM 参数
    (π_k, μ_k, Σ_k) 约束聚类目标，实现端到端的联合优化。
    训练完成后，使用联合优化后的模型直接进行聚类推断，而非重新拟合 GMM。
    """
    def __init__(self, dims, n_centroids, bn=False, dropout=0, binary=False):
        super(PhytoCluster, self).__init__(dims, bn=bn, dropout=dropout, binary=binary)
        self.n_centroids = n_centroids
        z_dim = dims[1]

        # GMM 聚类参数（可学习）
        # π: 各聚类的混合系数（先验概率）
        self.pi = nn.Parameter(torch.ones(n_centroids) / n_centroids)
        # μ_c: 各聚类的均值
        self.mu_c = nn.Parameter(torch.zeros(z_dim, n_centroids))
        # var_c: 各聚类的方差（对角协方差矩阵）
        self.var_c = nn.Parameter(torch.ones(z_dim, n_centroids))

    def compute_loss(self, x):
        """计算包含聚类目标的损失（ELBO with GMM posterior）"""
        z, mu, logvar = self.encoder(x)
        recon_x = self.decoder(z)
        gamma, mu_c, var_c, pi = self.infer_clusters(z)
        likelihood, kl_loss = compute_elbo(
            recon_x, x, gamma, (mu_c, var_c, pi), (mu, logvar), binary=self.binary
        )
        return -likelihood, kl_loss

    def infer_clusters(self, z):
        """
        基于潜在变量 z 推断聚类分配（软分配）。

        计算后验 q(c|x) = p(c) * p(z|c) / p(z)，
        其中 p(z|c) = N(z | μ_c, Σ_c)，p(c) = π_c。

        Args:
            z: 潜在变量，shape [N, z_dim]

        Returns:
            gamma: 后验概率矩阵 q(c|x)，shape [N, n_centroids]
            mu_c, var_c, pi: GMM 参数
        """
        n_centroids = self.n_centroids
        N = z.size(0)
        z_expanded = z.unsqueeze(2).expand(z.size(0), z.size(1), n_centroids)
        pi = self.pi.repeat(N, 1)
        mu_c = self.mu_c.repeat(N, 1, 1)
        var_c = self.var_c.repeat(N, 1, 1) + 1e-8

        # 计算 p(c) * p(z|c)
        p_c_z = torch.exp(
            torch.log(pi) - torch.sum(
                0.5 * torch.log(2 * math.pi * var_c) + (z_expanded - mu_c) ** 2 / (2 * var_c),
                dim=1
            )
        ) + 1e-10
        # 归一化得到后验 q(c|x)
        gamma = p_c_z / torch.sum(p_c_z, dim=1, keepdim=True)

        return gamma, mu_c, var_c, pi

    def initialize_gmm_params(self, dataloader, device='cpu'):
        """
        使用预训练 VAE 特征初始化 GMM 参数。

        先用 VAE encoder 提取潜在表示，再用 sklearn GaussianMixture 拟合，
        将拟合得到的均值和方差作为 GMM 参数的初始值。
        """
        gmm = GaussianMixture(
            n_components=self.n_centroids,
            covariance_type='diag',
            random_state=42
        )
        z = self.encode_batch(dataloader, device)
        gmm.fit(z)
        self.mu_c.data.copy_(torch.from_numpy(gmm.means_.T.astype(np.float32)))
        self.var_c.data.copy_(torch.from_numpy(gmm.covariances_.T.astype(np.float32)))

    def get_gamma(self, z):
        """获取聚类软分配（封装 infer_clusters 的第一返回值）"""
        return self.infer_clusters(z)[0]

    def predict(self, dataloader, device):
        """
        使用训练后的联合模型进行聚类预测。

        关键修复：不再使用预训练 VAE + sklearn GMM，
        而是使用阶段二联合优化后的 PhytoCluster 模型：
          1. 通过 encoder 获取潜在均值 μ
          2. 使用训练后的 GMM 参数 (π, μ_c, var_c) 计算后验分配
          3. 返回硬标签（argmax of gamma）

        Args:
            dataloader: 数据加载器
            device: 计算设备

        Returns:
            pred_labels: 聚类标签，shape [N,]
            embedding: 潜在均值，shape [N, z_dim]
            gamma: 软分配概率，shape [N, n_centroids]
        """
        self.eval()
        embedding = self.encode_batch(dataloader, device, output_type='mu')
        z_tensor = torch.from_numpy(embedding).float().to(device)

        with torch.no_grad():
            gamma, _, _, _ = self.infer_clusters(z_tensor)

        pred_labels = gamma.argmax(dim=1).cpu().numpy()
        return pred_labels, embedding, gamma.cpu().numpy()


def train_model(model, dataloader, device, args, phase='pretrain'):
    """训练模型（支持预训练阶段和联合聚类阶段）"""
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


def save_results(save_dir, y_true, y_pred, embedding, gamma, epoch):
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
    np.save(os.path.join(save_dir, f'gamma_{epoch}.npy'), gamma)

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


def save_trained_model(save_dir, model):
    """保存训练后的 PhytoCluster 联合模型"""
    model_path = os.path.join(save_dir, 'phytocluster_model.pt')
    torch.save(model.state_dict(), model_path)
    print(f'  PhytoCluster model saved to: {model_path}')

    # 保存 GMM 参数（用于分析）
    gmm_params = {
        'pi': model.pi.detach().cpu().numpy(),
        'mu_c': model.mu_c.detach().cpu().numpy(),
        'var_c': model.var_c.detach().cpu().numpy(),
        'n_centroids': model.n_centroids,
    }
    np.save(os.path.join(save_dir, 'gmm_params.npy'), gmm_params)
    print(f'  GMM parameters saved to: {os.path.join(save_dir, "gmm_params.npy")}')


def main():
    """主函数"""
    args = parse_args()

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
    print('Phase 2: PhytoCluster Training (VAE + GMM Joint Optimization)')
    print('='*60)
    model = PhytoCluster(dims, n_clusters, bn=False, dropout=0)
    model.load_state_dict(pretrain_model.state_dict(), strict=False)
    model.to(device)
    # 使用预训练 VAE 的特征初始化 GMM 参数
    model.initialize_gmm_params(full_loader, device)
    # 联合优化训练
    model = train_model(model, train_loader, device, args, phase='cluster')

    # 保存训练后的联合模型
    save_trained_model(args.save_dir, model)

    print('\n' + '='*60)
    print('Phase 3: Final Clustering (Using Jointly Trained Model)')
    print('='*60)
    # 关键修复：使用训练后的 PhytoCluster 联合模型进行聚类推断
    # 不再使用 pretrain_model.encode_mu() + sklearn GaussianMixture.fit_predict()
    # 而是使用 model.encode_batch() + model.infer_clusters() 直接推断
    pred_labels, embedding, gamma = model.predict(full_loader, device)

    save_results(args.save_dir, Y, pred_labels, embedding, gamma, 0)

    print(f'\nPhytoCluster completed. Results saved to: {args.save_dir}')
    print('\nNote: Clustering was performed using the jointly optimized')
    print('PhytoCluster model (VAE + GMM), not the pretrained VAE alone.')


if __name__ == '__main__':
    main()
