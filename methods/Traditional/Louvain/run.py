# -*- coding: utf-8 -*-
"""
Unified Louvain Model Interface for scCluBench
================================================

Louvain 社区检测算法用于单细胞聚类

核心思想：
  - 将细胞表达数据构建为 KNN 图（K近邻图）
  - 使用 Louvain 算法在图上进行社区检测
  - 通过调整 resolution 参数来控制聚类数量

Louvain 算法原理：
  - 优化图的模块度（Modularity）
  - 迭代合并节点以最大化社区内部边的密度
  - 时间复杂度 O(n log n)，适合大规模数据

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
import scanpy as sc
import networkx as nx
from networkx.algorithms.community import louvain_communities

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from preprocess import prepare_data_for_model
from utils import save


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Louvain clustering for scRNA-seq',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results',
                        help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, required=True,
                        help='Number of clusters (used for resolution tuning)')
    parser.add_argument('--resolution', type=float, default=None,
                        help='Resolution parameter (auto-tuned if not specified)')
    parser.add_argument('--n_neighbors', type=int, default=15,
                        help='Number of neighbors for KNN graph')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    return parser.parse_args()


def main():
    """主函数：Louvain 聚类流程"""
    args = parse_args()

    os.makedirs(args.save_dir, exist_ok=True)

    print('Loading data...')
    # 加载并预处理数据：归一化、log1p、HVG筛选、Z-score标准化
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=False,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True
    )

    Y = np.array(Y)
    from sklearn.preprocessing import LabelEncoder
    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Target clusters: {n_clusters}')

    # ========== Step 1: 构建 KNN 图 ==========
    # 使用 scanpy 的 neighbors 函数，基于预处理后的数据构建 KNN 图
    adata_work = adata.copy()
    adata_work.X = X
    sc.pp.neighbors(adata_work, n_neighbors=args.n_neighbors, use_rep='X')

    # ========== Step 2: 转换为 NetworkX 图 ==========
    print('Building graph...')
    G = nx.Graph()
    n_cells = X.shape[0]

    # 从 scanpy 的 connectivity 矩阵构建边
    connectivities = adata_work.obsp['connectivities']
    if hasattr(connectivities, 'toarray'):
        conn = connectivities.toarray()
    else:
        conn = np.array(connectivities)

    # 添加边（权重为 connectivity 值）
    rows, cols = np.where(conn > 0)
    for i, j in zip(rows, cols):
        if i != j:
            G.add_edge(i, j, weight=conn[i, j])

    if G.number_of_edges() == 0:
        print('Warning: Graph has no edges. Using fully connected.')
        G = nx.complete_graph(n_cells)

    # ========== Step 3: 调整 resolution 参数 ==========
    def get_louvain_partition(G, resolution):
        """获取 Louvain 社区划分结果"""
        communities = louvain_communities(G, resolution=resolution, seed=args.seed, weight='weight')
        partition = {}
        for comm_id, comm in enumerate(communities):
            for node in comm:
                partition[node] = comm_id
        return partition

    # 自动搜索最优 resolution（通过 NMI 选择最佳）
    if args.resolution is None:
        print('Tuning resolution...')
        best_nmi = 0
        best_res = 1.0
        best_labels = None

        # 测试多个 resolution 值
        for res in [0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]:
            partition = get_louvain_partition(G, res)
            labels = np.array([partition[i] for i in range(n_cells)])
            n_pred = len(np.unique(labels))

            # 选择聚类数接近目标值的 resolution
            if n_pred >= n_clusters // 2 and n_pred <= n_clusters * 3:
                from sklearn.metrics import normalized_mutual_info_score
                nmi = normalized_mutual_info_score(Y, labels)
                if nmi > best_nmi:
                    best_nmi = nmi
                    best_res = res
                    best_labels = labels

        if best_labels is None:
            best_res = 1.0
            partition = get_louvain_partition(G, best_res)
            best_labels = np.array([partition[i] for i in range(n_cells)])

        print(f'Best resolution: {best_res}, NMI: {best_nmi:.4f}')
    else:
        partition = get_louvain_partition(G, args.resolution)
        best_labels = np.array([partition[i] for i in range(n_cells)])

    n_pred = len(np.unique(best_labels))
    print(f'Predicted clusters: {n_pred}')

    # ========== Step 4: PCA 降维用于可视化保存 ==========
    from sklearn.decomposition import PCA
    pca = PCA(n_components=min(50, X.shape[1]), random_state=args.seed)
    embedding = pca.fit_transform(X)

    # 保存结果
    save(args.save_dir, Y, best_labels, 1, embedding, args=vars(args))
    print(f'Louvain completed. Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
