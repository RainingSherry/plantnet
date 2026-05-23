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
    parser.add_argument('--n_top_genes', type=int, default=2000,
                       help='Number of genes expected by unified benchmark input')
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
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts'].copy()

    # scVI 的基础过滤：移除表达量极低的基因和细胞
    sc.pp.filter_genes(adata, min_counts=3)
    sc.pp.filter_cells(adata, min_counts=3)
    if args.n_top_genes and adata.n_vars > args.n_top_genes:
        sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=args.n_top_genes, subset=True)

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
