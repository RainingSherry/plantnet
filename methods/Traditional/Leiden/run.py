# -*- coding: utf-8 -*-
"""
Unified Leiden Model Interface for scCluBench
===============================================

Leiden 社区检测算法用于单细胞聚类

核心思想：
  - 将细胞表达数据构建为 KNN 图（K近邻图）
  - 使用 Leiden 算法在图上进行社区检测
  - 通过调整 resolution 参数来控制聚类数量

Leiden 算法原理：
  - 同样优化图的模块度（Modularity），与 Louvain 目标一致
  - 在节点合并策略上进行了改进：Louvain 可能产生"倒置"（dangling）节点，
    而 Leiden 通过在每次迭代后对社区进行精细化（refinement）来避免此问题
  - Leiden 保证产生的社区是连通的，而 Louvain 不保证
  - 时间复杂度 O(n log n)，与 Louvain 相当，但在稀疏图上通常更快

关键区别（Louvain vs Leiden）：
  - Louvain：贪心合并后直接进入下一轮，可能在模块度局部最优处停滞
  - Leiden：每轮合并后对社区做精细化（move → refinement → aggregation），
    得到更高质量的社区划分（更高模块度、更明确的连通性）
  - Leiden 的 refine 步骤会将节点在社区间迁移，以优化局部模块度

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
import scanpy as sc
import networkx as nx
import leidenalg as la
import igraph as ig

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from preprocess import prepare_data_for_model
from utils import save


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Leiden clustering for scRNA-seq',
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
    """主函数：Leiden 聚类流程"""
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

    # ========== Step 2: 转换为 igraph 图 ==========
    # Leidenalg 需要 igraph 对象作为输入，而非 networkx
    print('Building graph...')
    n_cells = X.shape[0]

    # 从 scanpy 的 connectivity 矩阵提取边
    connectivities = adata_work.obsp['connectivities']
    if hasattr(connectivities, 'toarray'):
        conn = connectivities.toarray()
    else:
        conn = np.array(connectivities)

    # 构建 igraph 的边列表（跳过自环）
    edges = []
    weights = []
    rows, cols = np.where(conn > 0)
    for i, j in zip(rows, cols):
        if i != j:
            edges.append((i, j))
            weights.append(conn[i, j])

    if len(edges) == 0:
        print('Warning: Graph has no edges. Using fully connected.')
        g = ig.Graph.Full(n_cells)
    else:
        #directed=False 表示无向图
        g = ig.Graph(n_cells, edges=edges, edge_attrs={'weight': weights}, directed=False)

    # ========== Step 3: 调整 resolution 参数 ==========
    def get_leiden_partition(g, resolution):
        """
        获取 Leiden 社区划分结果

        相比 Louvain，Leiden 多出一个 refinement 阶段：
          1. 局部移动（move）：将节点移入相邻社区以提升模块度
          2. 精细化（refinement）：对社区进行精细化调整
          3. 聚合（aggregation）：将 refinement 后的社区聚合为超级节点
        这三个阶段循环迭代，直到模块度不再提升。

        这里使用 RBConfigurationVertexPartition 作为目标函数，
        支持 resolution 参数控制社区粒度，且在模块度优化上比默认 CPM 更灵活。
        """
        partition = la.find_partition(
            g,
            partition_type=la.RBConfigurationVertexPartition,
            weights='weight',
            resolution_parameter=resolution,
            seed=args.seed
        )
        # partition 保存了每个节点的社区 ID
        labels = np.array(partition.membership)
        return labels

    # 自动搜索最优 resolution（通过 NMI 选择最佳）
    if args.resolution is None:
        print('Tuning resolution...')
        best_nmi = 0
        best_res = 1.0
        best_labels = None

        # 测试多个 resolution 值，范围与 Louvain 一致
        for res in [0.3, 0.5, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5, 2.0]:
            labels = get_leiden_partition(g, res)
            n_pred = len(np.unique(labels))

            # 选择聚类数接近目标值的 resolution
            if n_pred >= n_clusters // 2 and n_pred <= n_clusters * 3:
                from sklearn.metrics import normalized_mutual_info_score
                nmi = normalized_mutual_info_score(Y, labels)
                if nmi > best_nmi:
                    best_nmi = nmi
                    best_res = res
                    best_labels = labels

        # 若未找到合适 resolution，使用默认参数
        if best_labels is None:
            best_res = 1.0
            best_labels = get_leiden_partition(g, best_res)

        print(f'Best resolution: {best_res}, NMI: {best_nmi:.4f}')
    else:
        best_labels = get_leiden_partition(g, args.resolution)

    n_pred = len(np.unique(best_labels))
    print(f'Predicted clusters: {n_pred}')

    # ========== Step 4: PCA 降维用于可视化保存 ==========
    from sklearn.decomposition import PCA
    pca = PCA(n_components=min(50, X.shape[1]), random_state=args.seed)
    embedding = pca.fit_transform(X)

    # 保存结果
    save(args.save_dir, Y, best_labels, 1, embedding)
    print(f'Leiden completed. Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
