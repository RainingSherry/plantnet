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
from pathlib import Path
import numpy as np
import torch
import random
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from torch.utils.data import DataLoader, Dataset

# 添加项目根目录到路径（用于导入benchmark通用模块）
CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json, sanitize_anndata_for_write

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
    parser.add_argument('--dataset_name', type=str, default=None,
                       help='数据集名称（用于输出记录）')
    parser.add_argument('--method_name', type=str, default='scMAE',
                       help='方法名称（用于输出记录）')
    parser.add_argument('--variant_name', type=str, default='scmae_original_self_only',
                       help='变体名称（用于输出记录）')
    parser.add_argument('--label_key', type=str, default='auto',
                       help='标签列名；auto 时使用 scMAE-family 统一候选')
    parser.add_argument('--input_mode', type=str, default='auto', choices=['auto', 'raw', 'log1p'],
                       help='输入矩阵模式')
    parser.add_argument('--n_top_genes', type=int, default=1000,
                       help='scMAE-family 统一 HVG 数量，源码默认 data_dim=1000')
    parser.add_argument('--target_sum', type=float, default=10000.0,
                       help='raw counts normalize_total 目标总量')
    parser.add_argument('--scale_input', type=family.str2bool, default=True,
                       help='是否对输入做 Z-score scale')

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
    parser.add_argument('--dropout', type=float, default=0.0,
                       help='Dropout率；源码默认0')

    # 训练参数
    parser.add_argument('--epochs', type=int, default=80,
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
                       help='兼容旧参数；公平对照时最终统一 kmeans_known_k 评估')
    parser.add_argument('--skip_eval', type=family.str2bool, default=False,
                       help='是否跳过 runner 内部 kmeans_known_k 输出')
    parser.add_argument('--no_save_h5ad', action='store_true',
                       help='不保存带 embedding 的 h5ad')

    return parser.parse_args()


def main():
    args = parse_args()
    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    print(f'Using device: {device}')

    print('Loading data...')
    bundle = family.load_scmae_dataset(
        file_path=args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        scale_input=args.scale_input,
        label_key=args.label_key,
        seed=args.seed,
    )
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    with open(save_dir / "selected_genes.txt", "w", encoding="utf-8") as handle:
        for gene in bundle.gene_names:
            handle.write(f"{gene}\n")

    data_np = bundle.data
    labels = bundle.labels
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(labels))
    dataset_name = args.dataset_name or Path(args.data_path).stem
    print(f'Number of cells: {data_np.shape[0]}, Number of genes: {data_np.shape[1]}')
    print(f'Number of clusters: {n_clusters}')
    print(f'Variant: {args.variant_name}; preprocessing: scMAE_family')

    dataset = family.IndexedExpressionDataset(data_np, labels)
    generator = torch.Generator()
    generator.manual_seed(args.seed)
    train_loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True,
        generator=generator,
    )
    test_loader = DataLoader(
        dataset,
        batch_size=max(args.batch_size * 4, 512),
        shuffle=False,
        drop_last=False,
    )

    model = AutoEncoder(
        num_genes=data_np.shape[1],
        hidden_size=args.hidden_size,
        dropout=args.dropout,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight,
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    history = {"loss": []}
    print('Starting training...')

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        total_loss = 0.0
        n_batches = 0

        for _, x_cpu, _ in train_loader:
            x = x_cpu.to(device)
            x_corrupted, mask = family.apply_scmae_noise(x, args.mask_prob)
            optimizer.zero_grad()
            _, loss = model.loss_mask(x_corrupted, x, mask)
            loss.backward()
            optimizer.step()

            total_loss += float(loss.detach().cpu())
            n_batches += 1

        avg_loss = total_loss / max(1, n_batches)
        history["loss"].append(avg_loss)
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f'Epoch {epoch:03d}/{args.epochs} loss={avg_loss:.4f}')

    embedding, labels_out = family.extract_embedding(model, test_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save(
        {
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "args": vars(args),
            "gene_names": bundle.gene_names.astype(str),
        },
        save_dir / "model_checkpoint.pth",
    )

    result = None
    eval_extra = {
        "variant": args.variant_name,
        "baseline": "original_scmae_self_branch_only",
        "mask_ratio": float(args.mask_prob),
        "preprocessing": "scMAE_family",
    }
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(
            output_dir=save_dir,
            dataset=dataset_name,
            method=args.method_name,
            seed=args.seed,
            embedding=embedding,
            labels=labels_out,
            n_clusters=n_clusters,
            extra=eval_extra,
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))

    if not args.no_save_h5ad:
        bundle.adata.obsm["X_scmae"] = embedding
        bundle.adata.uns["scmae"] = {
            "method": args.method_name,
            "variant": args.variant_name,
            "mask_ratio": float(args.mask_prob),
            "preprocessing": "scMAE_family",
        }
        sanitize_anndata_for_write(bundle.adata)
        bundle.adata.write_h5ad(save_dir / "adata_scmae.h5ad", compression="gzip")

    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "variant": args.variant_name,
        "seed": int(args.seed),
        "n_cells": int(data_np.shape[0]),
        "n_genes": int(data_np.shape[1]),
        "n_clusters": int(n_clusters),
        "embedding_path": str((save_dir / "embedding_final.npy").resolve()),
        "fixed_metrics": result["fixed"] if result is not None else {},
        "note": "Original scMAE AutoEncoder and self masked branch; preprocessing/evaluation are shared with NeighborMix_scMAE.",
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(f'Training completed. Results saved to: {save_dir}')


if __name__ == '__main__':
    main()
