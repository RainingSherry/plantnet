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


def kmeans_known_k(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int) -> tuple[dict, np.ndarray]:
    pred = KMeans(n_clusters=int(n_clusters), n_init=20, random_state=int(seed)).fit_predict(embedding)
    mapped = best_map(labels, pred.astype(np.int64))
    metrics = {
        "acc": float(np.mean(mapped == labels)),
        "nmi": float(normalized_mutual_info_score(labels, pred)),
        "ari": float(adjusted_rand_score(labels, pred)),
        "f1_macro": float(f1_score(labels, mapped, average="macro", zero_division=0)),
        "uses_known_k": True,
        "oracle-K": True,
        "cluster_method": "kmeans_known_k",
        "n_pred_clusters": int(len(np.unique(pred))),
    }
    return metrics, pred.astype(np.int64)

