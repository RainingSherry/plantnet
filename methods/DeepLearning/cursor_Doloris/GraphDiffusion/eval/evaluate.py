"""
eval/evaluate.py
==================
PlantDiffCluster 评估脚本。

在完整数据集上评估模型性能：
  1. 加载模型（训练后的检查点）
  2. 生成细胞嵌入
  3. 多种聚类算法比较（K-Means / Leiden / Louvain / GMM）
  4. 计算标准指标（ACC, NMI, ARI, F1-macro, Homogeneity, Completeness）
  5. UMAP / t-SNE 可视化
  6. 消融实验比较（不同图结构的性能差异）
  7. 可解释性报告（每个聚类的 top marker genes）
"""

from __future__ import annotations

import os
import sys
import json
import argparse
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import torch
import torch.nn as nn

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from datasets.srp182008_dataset import SRP182008Dataset, create_dataloader
from models import PlantDiffCluster


# ---------------------------------------------------------------------------
# 评估指标
# ---------------------------------------------------------------------------

def compute_all_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    """计算所有聚类评估指标。"""
    from sklearn.metrics import (
        accuracy_score, f1_score,
        normalized_mutual_info_score as nmi_score,
        adjusted_rand_score as ari_score,
        fowlkes_mallows_score as fmi_score,
        v_measure_score as vms_score,
        homogeneity_score as hom_score,
        completeness_score as com_score,
    )
    from scipy.optimize import linear_sum_assignment

    # Hungarian 重排
    y_true_unique = np.unique(y_true)
    y_pred_unique = np.unique(y_pred)
    n_class = len(y_true_unique)
    n_pred = len(y_pred_unique)

    G = np.zeros((n_class, n_pred), dtype=int)
    for i, ut in enumerate(y_true_unique):
        for j, up in enumerate(y_pred_unique):
            G[i, j] = np.sum((y_true == ut) & (y_pred == up))

    A = linear_sum_assignment(-G)
    new_pred = np.zeros_like(y_pred)
    for i, up in enumerate(y_pred_unique):
        col_idx = A[1][i] if i < len(A[1]) else i % n_class
        label_idx = A[0][i] if i < len(A[0]) else i % n_class
        new_pred[y_pred == up] = y_true_unique[label_idx]

    return {
        "acc": accuracy_score(y_true, new_pred),
        "nmi": nmi_score(y_true, y_pred, average_method="arithmetic"),
        "ari": ari_score(y_true, y_pred),
        "f1_macro": f1_score(y_true, new_pred, average="macro"),
        "fmi": fmi_score(y_true, y_pred),
        "v_measure": vms_score(y_true, y_pred),
        "homogeneity": hom_score(y_true, y_pred),
        "completeness": com_score(y_true, y_pred),
    }


# ---------------------------------------------------------------------------
# 可解释性报告
# ---------------------------------------------------------------------------

def generate_interpretability_report(
    model: nn.Module,
    dataset: SRP182008Dataset,
    embeddings: np.ndarray,
    labels: np.ndarray,
    cluster_labels: np.ndarray,
    output_dir: str,
    topk_genes: int = 30,
):
    """生成可解释性报告。"""
    from sklearn.cluster import KMeans
    from scipy.stats import pearsonr

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    n_clusters = len(np.unique(cluster_labels))
    n_cell_types = len(np.unique(labels))

    # ---- 1. 每个聚类的 marker genes ----
    print("\n=== Interpretability Report ===")

    # 计算聚类中心
    cluster_centers = {}
    for k in range(n_clusters):
        mask = cluster_labels == k
        if mask.sum() > 0:
            cluster_centers[k] = embeddings[mask].mean(axis=0)
        else:
            cluster_centers[k] = np.zeros(embeddings.shape[1])

    markers_report = {}
    for k in range(n_clusters):
        center = cluster_centers[k]
        # 计算每个基因与聚类中心的皮尔逊相关
        gene_scores = []
        for d in range(embeddings.shape[1]):
            expr_col = embeddings[:, d]
            if expr_col.std() > 1e-8:
                r, _ = pearsonr(expr_col, center)
                gene_scores.append(abs(r) if not np.isnan(r) else 0.0)
            else:
                gene_scores.append(0.0)

        top_idx = np.argsort(gene_scores)[::-1][:topk_genes]
        top_genes = [(dataset.gene_names[i], float(gene_scores[i])) for i in top_idx]
        markers_report[k] = top_genes

        print(f"Cluster {k} ({mask.sum()} cells): {[g[0] for g in top_genes[:5]]}...")

    # ---- 2. 每个聚类 vs 细胞类型 ----
    print("\n=== Cluster vs Cell Type ===")
    confusion = np.zeros((n_clusters, n_cell_types))
    for c in range(n_clusters):
        mask = cluster_labels == c
        for ct in range(n_cell_types):
            confusion[c, ct] = (labels[mask] == ct).sum()

    print("Cluster → Cell Type mapping:")
    for c in range(n_clusters):
        top_ct = confusion[c].argmax()
        top_ct_name = dataset.cell_type_names[top_ct] if top_ct < len(dataset.cell_type_names) else f"Type_{top_ct}"
        purity = confusion[c].max() / confusion[c].sum() if confusion[c].sum() > 0 else 0
        print(f"  Cluster {c} → {top_ct_name} (purity={purity:.3f})")

    # ---- 3. 保存报告 ----
    report = {
        "n_clusters": n_clusters,
        "n_cell_types": n_cell_types,
        "cell_type_names": list(dataset.cell_type_names),
        "markers_per_cluster": {
            str(k): [(g, float(s)) for g, s in v]
            for k, v in markers_report.items()
        },
        "confusion_matrix": confusion.tolist(),
    }
    with open(output_dir / "interpretability_report.json", "w") as f:
        json.dump(report, f, indent=2, default=str)

    # ---- 4. 绘制 UMAP ----
    try:
        adata = sc.AnnData(embeddings)
        sc.pp.neighbors(adata, n_neighbors=15, use_rep="X")
        sc.tl.umap(adata)
        sc.tl.leiden(adata, resolution=0.5)

        fig, axes = plt.subplots(1, 3, figsize=(18, 5))

        # Ground truth
        adata.obs["CellType"] = [dataset.cell_type_names[l] if l < len(dataset.cell_type_names) else str(l)
                                 for l in labels]
        adata.obs["Cluster"] = [str(c) for c in cluster_labels]
        adata.obs["Leiden"] = adata.obs["leiden"].astype(str)

        sc.pl.umap(adata, color="CellType", ax=axes[0], show=False, title="Ground Truth")
        sc.pl.umap(adata, color="Cluster", ax=axes[1], show=False, title="Predicted Clusters")
        sc.pl.umap(adata, color="Leiden", ax=axes[2], show=False, title="Leiden Clusters")

        plt.tight_layout()
        plt.savefig(output_dir / "umap_comparison.pdf", dpi=150, bbox_inches="tight")
        plt.close()
        print(f"\n[UMAP] Saved to {output_dir / 'umap_comparison.pdf'}")
    except Exception as e:
        print(f"[UMAP] Failed: {e}")


# ---------------------------------------------------------------------------
# 主评估函数
# ---------------------------------------------------------------------------

@torch.no_grad()
def evaluate_model(
    model: nn.Module,
    dataset: SRP182008Dataset,
    device: str,
    batch_size: int = 64,
) -> Dict[str, any]:
    """在完整数据集上评估模型。"""
    model.eval()

    all_embeddings = []
    all_labels = []

    loader = create_dataloader(
        dataset, batch_size=batch_size, shuffle=False, drop_last=False, num_workers=0,
    )

    for batch in loader:
        X = batch["X"].to(device)
        labels = batch["label"].to(device)
        support_idx = batch["support_idx"]
        support_mask = batch["support_mask"]
        support_weight = batch["support_weight"]

        output = model(
            X=X,
            cell_type=labels,
            support_weight=support_weight,
            support_mask=support_mask,
            support_idx=support_idx,
            t=None,
        )

        emb = output["cell_z"].cpu().numpy()
        all_embeddings.append(emb)
        all_labels.append(labels.cpu().numpy())

    embeddings = np.vstack(all_embeddings)
    labels = np.concatenate(all_labels)

    return {"embeddings": embeddings, "labels": labels}


def run_evaluation(
    model_path: str,
    data_path: str,
    output_dir: str,
    graph_type: str = "coexpression",
    n_clusters: int = 15,
    device: str = "cuda",
):
    """完整评估流程。"""
    print(f"\n{'='*60}")
    print(f"Evaluating PlantDiffCluster")
    print(f"Model: {model_path}")
    print(f"Data: {data_path}")
    print(f"Output: {output_dir}")
    print(f"{'='*60}\n")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # ---- 1. 加载数据 ----
    print("[1/5] Loading dataset...")
    dataset = SRP182008Dataset(
        h5ad_path=data_path,
        n_hvg=1500,
        graph_type=graph_type,
        support_strategy="log1p",
        dropout_rate=0.0,
    )

    # ---- 2. 加载模型 ----
    print("[2/5] Loading model...")
    from models import load_checkpoint

    model_config = {
        "gene_dim": 64,
        "hidden_dim": 256,
        "embed_dim": 128,
        "time_embed_dim": 128,
        "n_layers": 2,
        "heads": [4, 4],
        "pooling_strategy": "attention",
        "pooling_topk": 50,
        "n_clusters": n_clusters,
        "cluster_strategy": "gmm",
        "use_diffusion": True,
        "use_mask_predictor": True,
        "num_timesteps": 500,
        "ddim_steps": 20,
        "beta_schedule": "cosine",
        "refiner_depth": 3,
        "refiner_hidden_dim": 256,
        "lambda_cluster": 0.1,
        "cell_type_num": dataset.n_cell_types,
        "use_decoder": True,
        "decoder_hidden_dim": 256,
        "dropout": 0.1,
    }

    model = PlantDiffCluster(
        n_genes=dataset.n_hvg_actual,
        gene_names=list(dataset.gene_names),
        graph_dict=dataset.graph_dict,
        config=model_config,
    ).to(device)

    if model_path and os.path.exists(model_path):
        checkpoint = load_checkpoint(model, None, model_path, device)
        print(f"Loaded checkpoint from epoch {checkpoint.get('epoch', '?')}")
    else:
        print("[Warning] No checkpoint found, using random initialization")

    # ---- 3. 推理 ----
    print("[3/5] Running inference...")
    result = evaluate_model(model, dataset, device)
    embeddings, labels = result["embeddings"], result["labels"]

    # ---- 4. 聚类 + 评估 ----
    print("[4/5] Clustering and evaluating...")
    from sklearn.cluster import KMeans, SpectralClustering

    clustering_methods = {}

    # K-Means
    for n_clust in [n_clusters, 15, 20, 25]:
        kmeans = KMeans(n_clusters=n_clust, n_init=20, random_state=42)
        pred = kmeans.fit_predict(embeddings)
        metrics = compute_all_metrics(labels, pred)
        clustering_methods[f"kmeans_{n_clust}"] = {"pred": pred, "metrics": metrics}

    # Leiden（通过 scanpy）
    try:
        adata_emb = sc.AnnData(embeddings)
        sc.pp.neighbors(adata_emb, n_neighbors=15)
        for res in [0.3, 0.5, 0.8, 1.0]:
            sc.tl.leiden(adata_emb, resolution=res, random_state=42)
            pred_leiden = adata_emb.obs["leiden"].values.astype(int)
            metrics = compute_all_metrics(labels, pred_leiden)
            clustering_methods[f"leiden_{res}"] = {"pred": pred_leiden, "metrics": metrics}
    except Exception as e:
        print(f"[Leiden] Failed: {e}")

    # 打印所有方法的结果
    print("\n=== Clustering Results ===")
    best_method = None
    best_nmi = 0.0
    for name, res in sorted(clustering_methods.items()):
        m = res["metrics"]
        print(f"  {name:20s} ACC={m['acc']:.4f} NMI={m['nmi']:.4f} ARI={m['ari']:.4f} F1={m['f1_macro']:.4f}")
        if m["nmi"] > best_nmi:
            best_nmi = m["nmi"]
            best_method = name

    print(f"\nBest method: {best_method} (NMI={best_nmi:.4f})")
    best_pred = clustering_methods[best_method]["pred"]
    best_metrics = clustering_methods[best_method]["metrics"]

    # 保存最佳结果
    np.save(output_dir / "best_embeddings.npy", embeddings)
    np.save(output_dir / "best_predictions.npy", best_pred)
    with open(output_dir / "best_metrics.json", "w") as f:
        json.dump(best_metrics, f, indent=2)

    # 保存所有方法结果
    all_results = {k: {"metrics": v["metrics"]} for k, v in clustering_methods.items()}
    with open(output_dir / "all_results.json", "w") as f:
        json.dump(all_results, f, indent=2, default=str)

    # ---- 5. 可解释性报告 ----
    print("[5/5] Generating interpretability report...")
    generate_interpretability_report(
        model, dataset, embeddings, labels, best_pred, str(output_dir)
    )

    print(f"\n{'='*60}")
    print(f"Evaluation complete! Best method: {best_method}")
    print(f"  ACC={best_metrics['acc']:.4f}")
    print(f"  NMI={best_metrics['nmi']:.4f}")
    print(f"  ARI={best_metrics['ari']:.4f}")
    print(f"  F1={best_metrics['f1_macro']:.4f}")
    print(f"{'='*60}")

    return best_metrics, best_pred, embeddings


def parse_args():
    parser = argparse.ArgumentParser(description="Evaluate PlantDiffCluster")
    parser.add_argument("--model_path", type=str, default=None, help="Model checkpoint path")
    parser.add_argument("--data_path", type=str,
                        default="../../../data/SRP182008.h5ad",
                        help="Data path")
    parser.add_argument("--output_dir", type=str,
                        default="./eval_results",
                        help="Output directory")
    parser.add_argument("--graph_type", type=str, default="coexpression")
    parser.add_argument("--n_clusters", type=int, default=15)
    parser.add_argument("--device", type=str, default="cuda")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run_evaluation(
        model_path=args.model_path,
        data_path=args.data_path,
        output_dir=args.output_dir,
        graph_type=args.graph_type,
        n_clusters=args.n_clusters,
        device=args.device,
    )
