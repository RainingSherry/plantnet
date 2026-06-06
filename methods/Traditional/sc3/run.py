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
from scipy.optimize import linear_sum_assignment
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
    save(args.save_dir, Y, y_pred, 0, X_pca, args=vars(args))

    # ========== Step 3: 计算评估指标 ==========
    metrics_path = os.path.join(args.save_dir, 'metrics.json')

    # Use Hungarian algorithm (scipy) to align predicted labels with ground truth
    le = LabelEncoder()
    y_enc = le.fit_transform(y_pred)
    gt_enc = le.fit_transform(Y)
    D = max(int(y_enc.max()), int(gt_enc.max())) + 1
    cost = np.zeros((D, D), dtype=np.float64)
    for i in range(len(y_enc)):
        cost[int(y_enc[i]), int(gt_enc[i])] -= 1
    rows, cols = linear_sum_assignment(cost)
    y_map = {col: row for row, col in zip(rows, cols)}
    y_pred_aligned = np.array([y_map.get(int(p), int(p)) for p in y_enc])

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
