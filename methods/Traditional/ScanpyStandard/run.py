# -*- coding: utf-8 -*-
"""
====================================================================================================
Scanpy Standard Pipeline — 单细胞聚类基准方法
====================================================================================================

【目的】
    作为 benchmark 中的"生物学标准基线"方法，
    完整遵循 Scanpy 推荐的单细胞 RNA-seq 标准分析流程，
    用于与深度学习/传统聚类方法进行聚类性能对比。

【Scanpy 标准流程】
    Step 1: 读取 h5ad 数据
    Step 2: 质控（QC）— 过滤低质量细胞和基因
    Step 3: 归一化（normalize_total）
    Step 4: Log1p 变换（log(1+x)）
    Step 5: 高度可变基因（HVG）筛选
    Step 6: Z-score 标准化（可选，视下游任务而定）
    Step 7: PCA 降维
    Step 8: 构建邻域图（KNN）
    Step 9: UMAP 可视化降维
    Step 10: Leiden 社区检测聚类
    Step 11: 评估（vs ground truth）

【与其他 run.py 的核心区别】
    | 特性           | Benchmark 模型（scMAE/scVI/Leiden） | Scanpy 标准流程              |
    |---------------|----------------------------------|-----------------------------|
    | 预处理         | 固定 pipeline（normalize_per_cell） | Scanpy 官方推荐（normalize_total）|
    | HVG 参数       | 固定 n_top_genes=1000            | 可调参数                      |
    | Z-score       | 默认启用                         | 默认不启用（更接近原始生物学意义）    |
    | PCA           | 无                              | 有（step 7）                 |
    | UMAP          | 无                              | 有（step 9）                 |
    | 邻域图构建      | Leiden 用 igraph；scMAE 用 scanpy | scanpy pp.neighbors()       |
    | 聚类           | 各模型特定                      | Leiden（scanpy 内置）         |

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
import scanpy as sc

# NumPy 2.0 compatibility patch: np.string_ was removed
if not hasattr(np, 'string_'):
    np.string_ = np.bytes_

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from utils import save
from evaluation import evaluation as eval_fn


def parse_args():
    """解析命令行参数"""
    parser = argparse.ArgumentParser(
        description='Scanpy Standard Pipeline for scRNA-seq Clustering',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results',
                        help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, required=True,
                        help='Number of clusters (used for resolution tuning)')
    parser.add_argument('--min_genes', type=int, default=200,
                        help='Minimum genes per cell for QC filtering')
    parser.add_argument('--min_cells', type=int, default=3,
                        help='Minimum cells per gene for QC filtering')
    parser.add_argument('--target_sum', type=float, default=1e4,
                        help='Target total counts per cell after normalization')
    parser.add_argument('--n_top_genes', type=int, default=2000,
                        help='Number of highly variable genes to retain')
    parser.add_argument('--n_neighbors', type=int, default=15,
                        help='Number of neighbors for KNN graph')
    parser.add_argument('--n_pcs', type=int, default=50,
                        help='Number of principal components')
    parser.add_argument('--resolution', type=float, default=0.5,
                        help='Leiden clustering resolution parameter')
    parser.add_argument('--scale', action='store_true',
                        help='Enable Z-score scaling (default: disabled)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--figdir', type=str, default=None,
                        help='Figure output directory (default: <save_dir>/figures)')
    return parser.parse_args()


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    """
    搜索 Leiden 聚类的分辨率参数以获得指定数量的簇

    【算法】
        从高分辨率向低分辨率二分搜索，找到使聚类数最接近目标的 resolution。
    """
    import pandas as pd
    dis = []
    resolutions = sorted(list(np.arange(0.01, 2.5, increment)), reverse=True)

    for res in resolutions:
        sc.tl.leiden(adata, random_state=42, resolution=res, key_added='scanpy_leiden')
        count_unique_leiden = len(pd.DataFrame(adata.obs['scanpy_leiden']).scanpy_leiden.unique())
        dis.append(abs(count_unique_leiden - fixed_clus_count))
        if count_unique_leiden == fixed_clus_count:
            break

    return resolutions[np.argmin(dis)]


def main():
    """主函数：Scanpy 标准流程"""
    args = parse_args()

    os.makedirs(args.save_dir, exist_ok=True)

    # 设置图片保存路径
    if args.figdir:
        fig_dir = args.figdir
    else:
        fig_dir = os.path.join(args.save_dir, 'figures')
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    sc.settings.set_figure_params(dpi=80, facecolor='white')

    # =========================================================================
    # Step 1: 读取 h5ad 数据
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 1: Loading Data')
    print('=' * 60)

    adata = sc.read_h5ad(args.data_path)
    print(f'Raw data: {adata.n_obs} cells × {adata.n_vars} genes')
    print(f'obs columns: {list(adata.obs.columns)}')

    # 获取真实标签
    label_col = None
    for candidate in ['cell_type', 'Celltype', 'celltype', 'cell_label', 'label']:
        if candidate in adata.obs.columns:
            label_col = candidate
            break
    if label_col is None:
        raise KeyError(
            f"No cell type label column found. Available obs columns: {list(adata.obs.columns)}"
        )

    from sklearn.preprocessing import LabelEncoder
    Y_true = np.array(adata.obs[label_col])
    le = LabelEncoder()
    Y = le.fit_transform(Y_true)
    n_clusters_gt = len(np.unique(Y))
    n_clusters = args.n_clusters if args.n_clusters > 0 else n_clusters_gt
    print(f'Ground truth: {n_clusters_gt} cell types, target clusters: {n_clusters}')

    # =========================================================================
    # Step 2: 质控（QC）— 过滤低质量细胞和基因
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 2: Quality Control')
    print('=' * 60)

    # 计算 QC 指标（n_genes, n_counts）
    sc.pp.filter_cells(adata, min_genes=1)
    sc.pp.filter_genes(adata, min_cells=1)

    # 标记线粒体基因（如果基因名以 MT- 开头）
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True, percent_top=None)

    # 过滤低质量细胞和基因
    print(f'Before filtering: {adata.n_obs} cells, {adata.n_vars} genes')
    sc.pp.filter_cells(adata, min_genes=args.min_genes)
    sc.pp.filter_genes(adata, min_cells=args.min_cells)
    print(f'After filtering: {adata.n_obs} cells, {adata.n_vars} genes')

    # 保存原始数据到 adata.raw（用于后续基因表达量回溯）
    adata.raw = adata.copy()

    # =========================================================================
    # Step 3: 归一化（Per-cell normalization）
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 3: Normalization')
    print('=' * 60)

    sc.pp.normalize_total(adata, target_sum=args.target_sum)
    print(f'Normalized to target_sum={args.target_sum:.0e} per cell')

    # =========================================================================
    # Step 4: Log1p 变换
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 4: Log1p Transformation')
    print('=' * 60)

    sc.pp.log1p(adata)
    print('Applied log1p: log(1 + x)')

    # =========================================================================
    # Step 5: 高度可变基因（HVG）筛选
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 5: Highly Variable Genes Selection')
    print('=' * 60)

    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=args.n_top_genes,
        flavor='seurat'
    )
    n_hvg = adata.var['highly_variable'].sum()
    print(f'Selected {n_hvg} highly variable genes out of {adata.n_vars}')

    # 可视化 HVG
    sc.pl.highly_variable_genes(adata, save='_hvg.pdf')

    # 子集化，只保留 HVG
    adata = adata[:, adata.var['highly_variable']].copy()
    print(f'After HVG filtering: {adata.n_obs} cells × {adata.n_vars} genes')

    # =========================================================================
    # Step 6: Z-score 标准化（可选）
    # =========================================================================
    if args.scale:
        print('\n' + '=' * 60)
        print('Step 6: Z-score Scaling (optional)')
        print('=' * 60)
        sc.pp.scale(adata, max_value=10)
        print('Applied Z-score scaling with max_value=10')
    else:
        print('\n' + '=' * 60)
        print('Step 6: Skipping Z-score Scaling')
        print('=' * 60)
        print('Scale disabled — using normalized log1p data directly (recommended for clustering)')

    # =========================================================================
    # Step 7: PCA 降维
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 7: PCA Dimensionality Reduction')
    print('=' * 60)

    n_pcs_to_use = min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    sc.tl.pca(adata, svd_solver='arpack', n_comps=n_pcs_to_use)
    print(f'PCA computed with {n_pcs_to_use} components')
    print(f'Variance ratio (top 5 PCs): {adata.uns["pca"]["variance_ratio"][:5]}')

    # 可视化 PCA 方差
    sc.pl.pca_variance_ratio(adata, log=True, save='_pca.pdf')

    # =========================================================================
    # Step 8: 构建邻域图（KNN）
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 8: Building KNN Neighborhood Graph')
    print('=' * 60)

    n_pcs_neighbors = min(args.n_pcs, adata.n_vars - 1, adata.n_obs - 1)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, n_pcs=n_pcs_neighbors)
    print(f'KNN graph built: n_neighbors={args.n_neighbors}, n_pcs={n_pcs_neighbors}')

    # =========================================================================
    # Step 9: UMAP 可视化降维
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 9: UMAP Visualization')
    print('=' * 60)

    sc.tl.umap(adata)
    print('UMAP computed')

    # UMAP 可视化（用 ground truth 上色）
    sc.pl.umap(adata, color=[label_col], save='_groundtruth.pdf')
    print(f'UMAP saved with ground truth labels ({label_col})')

    # =========================================================================
    # Step 10: Leiden 社区检测聚类
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 10: Leiden Clustering')
    print('=' * 60)

    # 如果未指定 resolution，搜索以匹配目标聚类数
    if args.resolution is None or args.resolution <= 0:
        resolution = res_search_fixed_clus(adata, n_clusters)
        print(f'Auto-tuned resolution: {resolution:.4f} → target {n_clusters} clusters')
    else:
        resolution = args.resolution
        print(f'Using provided resolution: {resolution}')

    sc.tl.leiden(adata, resolution=resolution, random_state=42, key_added='scanpy_leiden')
    n_pred = adata.obs['scanpy_leiden'].nunique()
    print(f'Leiden clustering result: {n_pred} clusters predicted')

    # UMAP 可视化（用聚类结果上色）
    sc.pl.umap(adata, color='scanpy_leiden', legend_loc='on data', save='_leiden.pdf')

    # =========================================================================
    # Step 11: 评估
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 11: Evaluation')
    print('=' * 60)

    pred_labels = np.array([int(x) for x in adata.obs['scanpy_leiden'].to_list()])
    acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, _ = eval_fn(
        np.array(Y), pred_labels
    )

    print(f'\n{"=" * 50}')
    print(f'  Scanpy Standard Pipeline — Clustering Results')
    print(f'{"=" * 50}')
    print(f'  {"Metric":<18} {"Value":>10}')
    print(f'  {"-" * 30}')
    print(f'  {"Accuracy (ACC)":<18} {acc:>10.4f}')
    print(f'  {"NMI":<18} {nmi:>10.4f}')
    print(f'  {"ARI":<18} {ari:>10.4f}')
    print(f'  {"F1-macro":<18} {f1_macro:>10.4f}')
    print(f'  {"FMI":<18} {fmi:>10.4f}')
    print(f'  {"V-measure":<18} {v_measure:>10.4f}')
    print(f'  {"Homogeneity":<18} {hom:>10.4f}')
    print(f'  {"Completeness":<18} {com:>10.4f}')
    print(f'  {"-" * 30}')
    print(f'  Predicted clusters: {n_pred}')
    print(f'  Ground truth types: {n_clusters_gt}')
    print(f'  Resolution used: {resolution:.4f}')
    print(f'  HVG retained: {adata.n_vars}')
    print(f'  Cells after QC: {adata.n_obs}')
    print(f'  Scaled: {args.scale}')
    print(f'  {"=" * 50}\n')

    # =========================================================================
    # Step 12: 保存结果
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 12: Saving Results')
    print('=' * 60)

    # PCA embedding 作为降维表示（用于保存）
    embedding = adata.obsm['X_pca']

    # 保存聚类结果
    save(args.save_dir, Y, pred_labels, epoch=1, embedding=embedding)

    # 额外保存 UMAP embedding
    if 'X_umap' in adata.obsm:
        np.save(os.path.join(args.save_dir, 'umap_embedding.npy'), adata.obsm['X_umap'])

    # 保存配置信息
    import json
    config = {
        'model': 'Scanpy Standard Pipeline',
        'dataset': os.path.basename(args.data_path),
        'n_cells_raw': int(adata.raw.n_obs),
        'n_genes_raw': int(adata.raw.n_vars),
        'n_cells_after_qc': int(adata.n_obs),
        'n_genes_hvg': int(adata.n_vars),
        'n_clusters_gt': int(n_clusters_gt),
        'n_clusters_pred': int(n_pred),
        'n_top_genes': args.n_top_genes,
        'target_sum': args.target_sum,
        'n_neighbors': args.n_neighbors,
        'n_pcs': args.n_pcs,
        'resolution': float(resolution),
        'scale': args.scale,
        'min_genes': args.min_genes,
        'min_cells': args.min_cells,
        'seed': args.seed,
        'metrics': {
            'ACC': float(acc),
            'NMI': float(nmi),
            'ARI': float(ari),
            'F1_macro': float(f1_macro),
            'FMI': float(fmi),
            'V_measure': float(v_measure),
            'Homogeneity': float(hom),
            'Completeness': float(com),
        }
    }
    with open(os.path.join(args.save_dir, 'config.json'), 'w') as f:
        json.dump(config, f, indent=2)

    print(f'Results saved to: {args.save_dir}')
    print(f'Scanpy Standard Pipeline completed.')


if __name__ == '__main__':
    main()
