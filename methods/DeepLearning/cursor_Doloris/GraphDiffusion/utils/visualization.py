"""
utils/visualization.py
========================
可视化工具：UMAP/t-SNE、损失曲线、聚类结果、基因贡献热图。
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
import scanpy as sc
from typing import Optional, List


def plot_umap_sc(
    embeddings: np.ndarray,
    labels: np.ndarray,
    label_names: Optional[List[str]] = None,
    title: str = "UMAP",
    save_path: Optional[str] = None,
    figsize: tuple = (10, 8),
    point_size: float = 0.5,
    color_palette: str = "tab20",
):
    """
    使用 scanpy 绘制 UMAP 可视化。

    参数
    ----
    embeddings : [N, D]  嵌入矩阵
    labels : [N]  标签数组
    label_names : 标签名称列表
    title : str  图标题
    save_path : str  保存路径（.pdf / .png）
    figsize : tuple  图像大小
    point_size : float  点大小
    color_palette : str  调色板
    """
    try:
        adata = sc.AnnData(embeddings)
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X")
        sc.tl.umap(adata, random_state=42)

        label_key = "label"
        if label_names is not None:
            adata.obs[label_key] = [label_names[l] if l < len(label_names) else f"Type_{l}"
                                     for l in labels]
        else:
            adata.obs[label_key] = labels

        fig, ax = plt.subplots(figsize=figsize)
        sc.pl.umap(
            adata,
            color=label_key,
            ax=ax,
            show=False,
            title=title,
            size=point_size * 100,
            palette=color_palette,
        )
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"[UMAP] Saved to {save_path}")

        plt.close()
    except Exception as e:
        print(f"[UMAP] Failed: {e}")


def plot_loss_curves(
    loss_history: List[dict],
    save_path: Optional[str] = None,
    figsize: tuple = (14, 10),
):
    """
    绘制训练损失曲线。

    参数
    ----
    loss_history : List[dict]  每个 epoch 的损失字典
    save_path : str  保存路径
    figsize : tuple  图像大小
    """
    if not loss_history:
        return

    epochs = list(range(1, len(loss_history) + 1))

    loss_keys = list(loss_history[0].keys())
    n_plots = len(loss_keys)
    n_cols = 3
    n_rows = (n_plots + n_cols - 1) // n_cols

    fig, axes = plt.subplots(n_rows, n_cols, figsize=(figsize[0], n_rows * 4))
    axes = axes.flatten() if n_plots > 1 else [axes]

    for idx, key in enumerate(loss_keys):
        vals = [h.get(key, 0.0) for h in loss_history]
        axes[idx].plot(epochs, vals, label=key)
        axes[idx].set_title(f"{key} Loss")
        axes[idx].set_xlabel("Epoch")
        axes[idx].grid(True, alpha=0.3)
        axes[idx].legend()

    # 隐藏多余的子图
    for idx in range(n_plots, len(axes)):
        axes[idx].axis("off")

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"[Loss curves] Saved to {save_path}")

    plt.close()


def plot_gene_contribution_heatmap(
    gene_contributions: np.ndarray,  # [n_clusters, n_genes] or [n_cells, n_genes]
    gene_names: List[str],
    topk: int = 30,
    title: str = "Gene Contributions",
    save_path: Optional[str] = None,
    figsize: tuple = (14, 10),
):
    """
    绘制基因贡献热图。

    参数
    ----
    gene_contributions : [n_clusters, n_genes]  每个聚类对每个基因的贡献分数
    gene_names : 基因名列表
    topk : int  只显示 top-k 基因
    title : str  图标题
    save_path : str  保存路径
    figsize : tuple  图像大小
    """
    import seaborn as sns

    # 取每个聚类 top-k 基因
    n_clusters = gene_contributions.shape[0]

    # 取 top-k 基因
    top_gene_idx = []
    for k in range(n_clusters):
        top_k_idx = np.argsort(gene_contributions[k])[::-1][:topk]
        top_gene_idx.extend(top_k_idx)

    # 取并集（去重但保持顺序）
    top_gene_idx = list(dict.fromkeys(top_gene_idx))[:topk * n_clusters]
    sub_matrix = gene_contributions[:, top_gene_idx]
    sub_genes = [gene_names[i] for i in top_gene_idx]

    fig, ax = plt.subplots(figsize=figsize)
    sns.heatmap(
        sub_matrix,
        xticklabels=sub_genes,
        yticklabels=[f"Cluster_{k}" for k in range(n_clusters)],
        ax=ax,
        cmap="YlOrRd",
        cbar_kws={"label": "Contribution"},
    )
    ax.set_title(title)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right", fontsize=7)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"[Heatmap] Saved to {save_path}")

    plt.close()


def plot_cluster_comparison(
    embeddings: np.ndarray,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    true_names: Optional[List[str]] = None,
    pred_names: Optional[List[str]] = None,
    save_path: Optional[str] = None,
    figsize: tuple = (18, 5),
):
    """
    绘制聚类对比图：Ground Truth / Predicted / Leiden。
    """
    try:
        adata = sc.AnnData(embeddings)
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X")
        sc.tl.umap(adata, random_state=42)
        sc.tl.leiden(adata, resolution=0.5, random_state=42)

        label_true_key = "GroundTruth"
        label_pred_key = "Predicted"

        adata.obs[label_true_key] = [true_names[l] if true_names and l < len(true_names) else f"Type_{l}"
                                       for l in y_true]
        adata.obs[label_pred_key] = [f"Cluster_{l}" for l in y_pred]
        adata.obs["Leiden"] = adata.obs["leiden"].astype(str)

        fig, axes = plt.subplots(1, 3, figsize=figsize)
        sc.pl.umap(adata, color=label_true_key, ax=axes[0], show=False, title="Ground Truth", size=30)
        sc.pl.umap(adata, color=label_pred_key, ax=axes[1], show=False, title="Predicted", size=30)
        sc.pl.umap(adata, color="Leiden", ax=axes[2], show=False, title="Leiden", size=30)

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=150, bbox_inches="tight")
            print(f"[Comparison] Saved to {save_path}")

        plt.close()
    except Exception as e:
        print(f"[Comparison] Failed: {e}")


def plot_ablation_results(
    results: dict,  # {graph_type: {"acc":..., "nmi":..., "ari":...}}
    save_path: Optional[str] = None,
    figsize: tuple = (12, 6),
):
    """
    绘制消融实验结果对比图。
    """
    import seaborn as sns

    methods = list(results.keys())
    metrics = ["acc", "nmi", "ari", "f1_macro"]

    data_rows = []
    for method in methods:
        for metric in metrics:
            data_rows.append({
                "Method": method,
                "Metric": metric,
                "Value": results[method].get(metric, 0.0),
            })

    import pandas as pd
    df = pd.DataFrame(data_rows)

    fig, ax = plt.subplots(figsize=figsize)
    sns.barplot(data=df, x="Metric", y="Value", hue="Method", ax=ax, palette="Set2")
    ax.set_title("Ablation Study: Gene Graph Type")
    ax.set_ylim(0, 1)
    ax.grid(axis="y", alpha=0.3)

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=150, bbox_inches="tight")
        print(f"[Ablation] Saved to {save_path}")

    plt.close()
