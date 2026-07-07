# -*- coding: utf-8 -*-
"""
====================================================================================================
scVI — Single-cell Variational Inference (Self-Contained PyTorch Implementation)
====================================================================================================

【论文来源】
    Lopez et al. (2018), Nature Methods
    "Deep generative modeling for single-cell transcriptomics"
    https://www.nature.com/articles/s41592-018-0229-2

【本实现】
    基于 PyTorch 的自包含实现，不依赖 scvi-tools/jax/numpyro。
    完全遵循原论文的 ZINB-VAE 模型：Encoder → Z → Decoder → ZINB likelihood

【生成模型】
    z_n ~ Normal(0, I)                                          潜在变量
    ℓ_n = sum(x_n)                                              Library size (observed)
    ρ_n = softmax(decoder(z_n))                                 归一化表达比例
    θ_g   (可学习)                                               基因 dropout 参数
    π_n = sigmoid(h(z_n))                                      Dropout 概率
    x_ng ~ ZINB(ℓ_n · ρ_ng, θ_g, π_n)                        观测: ZINB

【损失函数】
    ELBO = E_{q(z|x)}[log p(x|z)] - KL(q(z|x) || N(0,I))
    ZINB log-likelihood for each gene per cell.
"""

import os
import sys
import argparse
import random
import time
from datetime import datetime

import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.preprocessing import LabelEncoder

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
_DEEPDIR = os.path.dirname(SCRIPT_DIR)  # methods/DeepLearning/
_ROOT = os.path.dirname(os.path.dirname(_DEEPDIR))  # project root
for _p in [_ROOT, _DEEPDIR, SCRIPT_DIR]:
    if _p not in sys.path:
        sys.path.insert(0, _p)

try:
    from utils import save as _save_fn
except ImportError:
    _save_fn = None

try:
    from evaluation import evaluation as eval_fn
except ImportError:
    def eval_fn(y_true, y_pred):
        from sklearn.metrics import (
            accuracy_score, normalized_mutual_info_score,
            adjusted_rand_score, f1_score, fowlkes_mallows_score,
            v_measure_score, homogeneity_score, completeness_score
        )
        return (
            accuracy_score(y_true, y_pred),
            normalized_mutual_info_score(y_true, y_pred),
            adjusted_rand_score(y_true, y_pred),
            f1_score(y_true, y_pred, average='macro'),
            fowlkes_mallows_score(y_true, y_pred),
            v_measure_score(y_true, y_pred),
            homogeneity_score(y_true, y_pred),
            completeness_score(y_true, y_pred),
            {}
        )


# =============================================================================
# PyTorch ZINB-VAE 模型
# =============================================================================

def log_zinb_positive(x, mu, theta, pi):
    """Log ZINB probability. All inputs are positive."""
    theta = torch.clamp(theta, min=1e-6)
    mu = torch.clamp(mu, min=1e-6)
    pi = torch.clamp(pi, min=1e-6, max=1.0 - 1e-6)

    case_zero = torch.logaddexp(
        torch.log(pi),
        torch.log(1 - pi) + torch.lgamma(x + theta) - torch.lgamma(theta) - torch.lgamma(x + 1)
        + theta * (torch.log(theta) - torch.log(theta + mu))
        + x * (torch.log(mu) - torch.log(theta + mu))
    )
    case_nonzero = torch.log(1 - pi) + torch.lgamma(x + theta) - torch.lgamma(theta) - torch.lgamma(x + 1) \
        + theta * (torch.log(theta) - torch.log(theta + mu)) \
        + x * (torch.log(mu) - torch.log(theta + mu))
    mask = (x > 0).float()
    return mask * case_nonzero + (1 - mask) * case_zero


class Encoder(nn.Module):
    """Encodes gene expression into latent space (mean + log_var)."""
    def __init__(self, n_genes: int, latent_dim: int, hidden_dim: int, n_layers: int):
        super().__init__()
        layers = []
        dims = [n_genes] + [hidden_dim] * n_layers
        for i in range(n_layers):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            layers.append(nn.BatchNorm1d(dims[i + 1]))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.1))
        self.net = nn.Sequential(*layers)
        self.fc_mean = nn.Linear(hidden_dim, latent_dim)
        self.fc_log_var = nn.Linear(hidden_dim, latent_dim)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        h = self.net(x)
        mean = self.fc_mean(h)
        log_var = self.fc_log_var(h)
        return mean, log_var


class Decoder(nn.Module):
    """Decodes latent variable to gene-wise ZINB parameters."""
    def __init__(self, n_genes: int, latent_dim: int, hidden_dim: int, n_layers: int):
        super().__init__()
        layers = []
        dims = [latent_dim] + [hidden_dim] * n_layers
        for i in range(n_layers):
            layers.append(nn.Linear(dims[i], dims[i + 1]))
            layers.append(nn.BatchNorm1d(dims[i + 1]))
            layers.append(nn.ReLU())
            layers.append(nn.Dropout(0.1))
        layers.append(nn.Linear(hidden_dim, n_genes))
        self.net = nn.Sequential(*layers)

    def forward(self, z: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        out = self.net(z)
        mu = F.softmax(out, dim=-1)
        return mu


class ZINBVAE(nn.Module):
    """Zero-Inflated Negative Binomial Variational Autoencoder."""
    def __init__(self, n_genes: int, latent_dim: int, hidden_dim: int, n_layers: int, encode_dim: int):
        super().__init__()
        self.latent_dim = latent_dim
        self.n_genes = n_genes
        self.encoder = Encoder(n_genes, latent_dim, encode_dim, n_layers)
        self.decoder = Decoder(n_genes, latent_dim, encode_dim, n_layers)
        self.dropout_net = nn.Sequential(
            nn.Linear(latent_dim, encode_dim // 2),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(encode_dim // 2, n_genes),
        )
        # Gene-wise NB theta (inverse dispersion), initialized with NBMLE-like heuristic
        self.theta = nn.Parameter(torch.ones(n_genes) * 10.0)

    def reparameterize(self, mean: torch.Tensor, log_var: torch.Tensor) -> torch.Tensor:
        std = torch.exp(0.5 * log_var)
        eps = torch.randn_like(std)
        return mean + eps * std

    def forward(self, x: torch.Tensor, return_latent: bool = False):
        mean, log_var = self.encoder(x)
        z = self.reparameterize(mean, log_var)
        mu = self.decoder(z)
        pi = torch.sigmoid(self.dropout_net(z))

        if return_latent:
            return z, mean
        return z, mu, pi, mean, log_var

    def get_latent(self, x: torch.Tensor) -> np.ndarray:
        self.eval()
        with torch.no_grad():
            x_t = torch.as_tensor(x, dtype=torch.float32)
            if x_t.ndim == 1:
                x_t = x_t.unsqueeze(0)
            mean, _ = self.encoder(x_t)
            return mean.cpu().numpy()


# =============================================================================
# 训练函数
# =============================================================================

def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def zinb_loss(x: torch.Tensor, mu: torch.Tensor, theta: torch.Tensor,
              pi: torch.Tensor, mean: torch.Tensor, log_var: torch.Tensor,
              kl_weight: float = 1.0) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    """Compute ZINB reconstruction loss + KL divergence."""
    # Library size normalization
    lib_size = x.sum(dim=1, keepdim=True) + 1e-6
    mu_scaled = mu * lib_size

    # ZINB reconstruction loss (per-gene, then sum)
    ll = log_zinb_positive(x, mu_scaled, theta.exp(), pi)
    recon = -ll.sum(dim=1).mean()

    # KL divergence: KL(N(mean,var) || N(0,I))
    kl = -0.5 * (1 + log_var - mean.pow(2) - log_var.exp()).sum(dim=1).mean()

    total = recon + kl_weight * kl
    return total, recon, kl


def train_model(model, data_loader, optimizer, device, kl_weight=1.0):
    model.train()
    total_loss, total_recon, total_kl = 0.0, 0.0, 0.0
    n_batches = 0
    for batch in data_loader:
        if isinstance(batch, (list, tuple)):
            x = batch[0].to(device)
        else:
            x = batch.to(device)
        optimizer.zero_grad()
        _, mu, pi, mean, log_var = model(x)
        loss, recon, kl = zinb_loss(x, mu, model.theta, pi, mean, log_var, kl_weight)
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
        optimizer.step()
        total_loss += loss.item()
        total_recon += recon.item()
        total_kl += kl.item()
        n_batches += 1
    return total_loss / n_batches, total_recon / n_batches, total_kl / n_batches


# =============================================================================
# 参数解析
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description='scVI: Self-contained PyTorch ZINB-VAE implementation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--data_path', type=str, required=True)
    parser.add_argument('--save_dir', type=str, required=True)
    parser.add_argument('--n_clusters', type=int, required=True)
    parser.add_argument('--n_top_genes', type=int, default=2000)
    parser.add_argument('--latent_dim', type=int, default=10)
    parser.add_argument('--n_layers', type=int, default=1)
    parser.add_argument('--encode_dim', type=int, default=128)
    parser.add_argument('--epochs', type=int, default=200)
    parser.add_argument('--batch_size', type=int, default=256)
    parser.add_argument('--lr', type=float, default=0.001)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--gpu', type=int, default=1)
    return parser.parse_args()


# =============================================================================
# 主函数
# =============================================================================

def main():
    args = parse_args()
    set_seed(args.seed)

    device = torch.device(f'cuda:{args.gpu}' if torch.cuda.is_available() else 'cpu')
    print(f'Device: {device}')

    # Step 1: 加载数据
    print('\n' + '=' * 60)
    print('Step 1: Loading and Preprocessing Data')
    print('=' * 60)

    adata_raw = sc.read_h5ad(args.data_path)
    label_col = None
    for candidate in ['resolved_label', 'cell_type', 'Celltype', 'celltype', 'cell_label', 'label', 'maintype', 'type']:
        if candidate in adata_raw.obs.columns:
            label_col = candidate
            break
    if label_col is None:
        raise KeyError(f"No label column found. Available: {list(adata_raw.obs.columns)}")

    Y_encoded = LabelEncoder().fit_transform(np.array(adata_raw.obs[label_col]))
    n_clusters = args.n_clusters
    print(f'Cells: {adata_raw.n_obs}, Genes: {adata_raw.n_vars}, Clusters: {n_clusters}')

    # Step 2: 数据过滤（与 scvi-tools 保持一致）
    print('\n' + '=' * 60)
    print('Step 2: Filtering Genes and Cells')
    print('=' * 60)

    adata = adata_raw.copy()
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()

    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)

    if adata.n_vars > args.n_top_genes:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=args.n_top_genes, subset=True)

    print(f'After filtering: {adata.n_obs} cells, {adata.n_vars} genes')

    # Step 3: 准备训练数据
    print('\n' + '=' * 60)
    print('Step 3: Preparing Training Data')
    print('=' * 60)

    if hasattr(adata.X, 'toarray'):
        X = adata.X.toarray()
    else:
        X = np.array(adata.X)
    X = X.astype(np.float32)
    X = np.maximum(X, 0.0)  # 确保非负

    # Normalize to total counts per cell (library size)
    lib_sizes = X.sum(axis=1, keepdims=True)
    lib_sizes[lib_sizes == 0] = 1.0
    X_norm = X / lib_sizes * 1e4  # CPM-like normalization
    X_input = X_norm.astype(np.float32)

    n_genes = X_input.shape[1]
    n_cells = X_input.shape[0]
    print(f'Training data: {n_cells} cells x {n_genes} genes')

    # Step 4: 创建模型
    print('\n' + '=' * 60)
    print('Step 4: Building ZINB-VAE Model')
    print('=' * 60)

    model = ZINBVAE(
        n_genes=n_genes,
        latent_dim=args.latent_dim,
        hidden_dim=args.encode_dim,
        n_layers=args.n_layers,
        encode_dim=args.encode_dim,
    ).to(device)

    total_params = sum(p.numel() for p in model.parameters())
    print(f'Model parameters: {total_params:,}')
    print(f'Architecture: Input({n_genes}) -> Encoder -> z({args.latent_dim}) -> Decoder -> ZINB')

    # Step 5: 训练
    print('\n' + '=' * 60)
    print('Step 5: Training ZINB-VAE')
    print('=' * 60)

    dataset = torch.utils.data.TensorDataset(
        torch.from_numpy(X_input)
    )
    loader = torch.utils.data.DataLoader(
        dataset, batch_size=args.batch_size, shuffle=True, drop_last=False
    )

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-5)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    best_loss = float('inf')
    best_state = None
    start_time = time.time()

    for epoch in range(1, args.epochs + 1):
        loss, recon, kl = train_model(model, loader, optimizer, device, kl_weight=1.0)
        scheduler.step()

        if epoch % 20 == 0 or epoch == args.epochs:
            print(f'Epoch {epoch:3d}/{args.epochs} | loss={loss:.4f} | recon={recon:.4f} | kl={kl:.4f}')

        if loss < best_loss:
            best_loss = loss
            best_state = {k: v.cpu().clone() for k, v in model.state_dict().items()}

    train_time = time.time() - start_time
    print(f'\nTraining completed in {train_time:.1f}s')

    # 加载最佳模型
    if best_state is not None:
        model.load_state_dict(best_state)
        model.to(device)

    # Step 6: 提取 latent embedding
    print('\n' + '=' * 60)
    print('Step 6: Extracting Latent Embedding')
    print('=' * 60)

    model.eval()
    with torch.no_grad():
        X_tensor = torch.from_numpy(X_input).to(device)
        embeddings = []
        batch_size_eval = 512
        for i in range(0, len(X_tensor), batch_size_eval):
            batch = X_tensor[i:i+batch_size_eval]
            emb = model.get_latent(batch)
            embeddings.append(emb)
        embedding = np.concatenate(embeddings, axis=0)

    print(f'Embedding shape: {embedding.shape}')

    # Step 7: 聚类
    print('\n' + '=' * 60)
    print('Step 7: Clustering')
    print('=' * 60)

    if n_cells < 10000:
        kmeans = KMeans(n_clusters=n_clusters, random_state=args.seed, n_init=20, max_iter=300)
        pred_labels = kmeans.fit_predict(embedding)
    else:
        adata_emb = sc.AnnData(embedding)
        sc.pp.neighbors(adata_emb, n_neighbors=15, use_rep='X')
        sc.tl.leiden(adata_emb, resolution=1.0, key_added='leiden')
        pred_labels = np.array([int(x) for x in adata_emb.obs['leiden']])

    # Step 8: 评估
    print('\n' + '=' * 60)
    print('Step 8: Evaluation')
    print('=' * 60)

    acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, _ = eval_fn(
        np.array(Y_encoded), np.array(pred_labels)
    )
    print(f'ACC: {acc:.4f} | NMI: {nmi:.4f} | ARI: {ari:.4f}')
    print(f'F1:  {f1_macro:.4f} | FMI: {fmi:.4f} | V-measure: {v_measure:.4f}')

    # Step 9: 保存结果
    print('\n' + '=' * 60)
    print('Step 9: Saving Results')
    print('=' * 60)

    os.makedirs(args.save_dir, exist_ok=True)

    if _save_fn is not None:
        _save_fn(args.save_dir, Y_encoded, pred_labels, args.epochs, embedding)
    np.save(os.path.join(args.save_dir, 'embedding_final.npy'), embedding)

    import json
    metrics = {
        'acc': float(acc), 'nmi': float(nmi), 'ari': float(ari),
        'f1_macro': float(f1_macro), 'fmi': float(fmi),
        'v_measure': float(v_measure), 'homogeneity': float(hom),
        'completeness': float(com),
        'train_time_seconds': float(train_time),
    }
    with open(os.path.join(args.save_dir, 'metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=2)

    config = {
        'model': 'scVI (PyTorch ZINB-VAE)',
        'n_cells': int(n_cells),
        'n_genes': int(n_genes),
        'n_clusters': int(n_clusters),
        'latent_dim': args.latent_dim,
        'n_layers': args.n_layers,
        'encode_dim': args.encode_dim,
        'epochs': args.epochs,
        'batch_size': args.batch_size,
        'lr': args.lr,
        'seed': args.seed,
        'gpu': args.gpu,
        'train_time_seconds': float(train_time),
    }
    with open(os.path.join(args.save_dir, 'config.json'), 'w') as f:
        json.dump(config, f, indent=2)

    torch.save(model.state_dict(), os.path.join(args.save_dir, 'scvi_model.pt'))

    print(f'\nscVI completed successfully! Results saved to: {args.save_dir}')
    return metrics, args.save_dir


if __name__ == '__main__':
    main()
