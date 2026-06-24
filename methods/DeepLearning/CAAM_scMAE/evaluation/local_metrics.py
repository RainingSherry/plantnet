from __future__ import annotations

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, f1_score, normalized_mutual_info_score


def best_map(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    true_values = np.unique(y_true)
    pred_values = np.unique(y_pred)
    n = max(len(true_values), len(pred_values))
    counts = np.zeros((n, n), dtype=np.int64)
    for i, t in enumerate(true_values):
        for j, p in enumerate(pred_values):
            counts[i, j] = int(np.sum((y_true == t) & (y_pred == p)))
    rows, cols = linear_sum_assignment(-counts)
    mapped = np.zeros_like(y_pred, dtype=np.int64)
    for row, col in zip(rows, cols):
        if row < len(true_values) and col < len(pred_values):
            mapped[y_pred == pred_values[col]] = true_values[row]
    return mapped


def clustering_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    mapped = best_map(y_true, y_pred.astype(np.int64))
    return {
        "acc": float(np.mean(mapped == y_true)),
        "nmi": float(normalized_mutual_info_score(y_true, y_pred)),
        "ari": float(adjusted_rand_score(y_true, y_pred)),
        "f1_macro": float(f1_score(y_true, mapped, average="macro", zero_division=0)),
        "n_pred_clusters": int(len(np.unique(y_pred))),
    }


def kmeans_known_k(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int) -> tuple[dict, np.ndarray]:
    pred = KMeans(n_clusters=int(n_clusters), n_init=20, random_state=int(seed)).fit_predict(embedding)
    metrics = clustering_metrics(labels, pred.astype(np.int64))
    metrics.update(
        {
            "uses_known_k": True,
            "oracle-K": True,
            "cluster_method": "kmeans_known_k",
        }
    )
    return metrics, pred.astype(np.int64)


def leiden_fixed(
    embedding: np.ndarray,
    labels: np.ndarray,
    *,
    resolution: float = 1.0,
    n_neighbors: int = 15,
    seed: int = 0,
) -> tuple[dict, np.ndarray]:
    """Run fixed-resolution Leiden on the learned embedding.

    This is the unknown-K protocol: it does not use the true number of clusters and
    does not perform label-selected/oracle resolution search.
    """
    import anndata as ad
    import scanpy as sc

    emb = np.asarray(embedding, dtype=np.float32)
    n = int(emb.shape[0])
    if n < 2:
        raise ValueError("leiden_fixed requires at least two cells.")
    adata = ad.AnnData(emb)
    adata.obsm["X_caam"] = emb
    sc.pp.neighbors(
        adata,
        n_neighbors=max(1, min(int(n_neighbors), n - 1)),
        use_rep="X_caam",
        random_state=int(seed),
    )
    sc.tl.leiden(adata, resolution=float(resolution), random_state=int(seed), key_added="leiden_fixed")
    categories = adata.obs["leiden_fixed"].astype("category").cat.codes.to_numpy().astype(np.int64)
    metrics = clustering_metrics(labels, categories)
    metrics.update(
        {
            "uses_known_k": False,
            "oracle-K": False,
            "cluster_method": "leiden_fixed",
            "resolution": float(resolution),
            "n_neighbors": int(n_neighbors),
        }
    )
    return metrics, categories
