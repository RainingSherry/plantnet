import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from methods.evaluation import evaluation as benchmark_evaluation


# ---------- 聚类与评估 ----------
# 在嵌入空间中用KMeans聚类，然后与真实标签对比计算NMI/ARI等指标


def cluster_and_evaluate(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int = 42) -> dict:
    """在嵌入上执行KMeans聚类，计算与真实标签间的聚类评价指标.

    评估指标包括：
    - ACC（分类准确率，需标签对齐）
    - NMI（归一化互信息）
    - ARI（调整兰德指数）
    - F1-macro
    - FMI（Fowlkes-Mallows Index）
    - V-measure（同质性+完整性调和平均）
    - Silhouette Score（簇内紧密度与簇间分离度）
    """
    pred_labels = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    # benchmark_evaluation 内部做 Hungarian 标签对齐，返回准确率和对齐后的预测标签
    acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, y_pred_mapped = benchmark_evaluation(labels, pred_labels)
    if len(np.unique(pred_labels)) > 1 and embedding.shape[0] > len(np.unique(pred_labels)):
        try:
            silhouette = float(silhouette_score(embedding, pred_labels))
        except Exception:
            silhouette = float("nan")
    else:
        silhouette = float("nan")
    return {
        "pred_labels": pred_labels.astype(np.int64),
        "pred_labels_mapped": y_pred_mapped.astype(np.int64),
        "metrics": {
            "acc": float(acc),
            "nmi": float(nmi),
            "ari": float(ari),
            "f1_macro": float(f1_macro),
            "fmi": float(fmi),
            "v_measure": float(v_measure),
            "homogeneity": float(hom),
            "completeness": float(com),
            "silhouette": silhouette,
        },
    }


def evaluate_embedding_set(embeddings: dict, labels: np.ndarray, n_clusters: int, seed: int = 42) -> dict:
    """对多个嵌入字典中的每个嵌入分别聚类评估（用于对比不同方法的结果）."""
    return {
        name: cluster_and_evaluate(embedding, labels, n_clusters=n_clusters, seed=seed)
        for name, embedding in embeddings.items()
    }

