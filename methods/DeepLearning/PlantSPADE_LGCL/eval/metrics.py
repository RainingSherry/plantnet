from __future__ import annotations

from typing import Optional, Tuple

import numpy as np
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    f1_score,
    fowlkes_mallows_score,
    homogeneity_score,
    normalized_mutual_info_score,
    silhouette_score,
    v_measure_score,
)


def best_map(y_true, y_pred) -> np.ndarray:
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    true_values = np.unique(y_true)
    pred_values = np.unique(y_pred)
    n = max(len(true_values), len(pred_values))
    overlap = np.zeros((n, n), dtype=np.int64)
    for i, true_label in enumerate(true_values):
        true_mask = y_true == true_label
        for j, pred_label in enumerate(pred_values):
            overlap[i, j] = int(np.sum(true_mask & (y_pred == pred_label)))

    rows, cols = linear_sum_assignment(-overlap)
    mapping = {}
    for row, col in zip(rows, cols):
        if row < len(true_values) and col < len(pred_values):
            mapping[pred_values[col]] = true_values[row]

    for pred_label in pred_values:
        if pred_label not in mapping:
            col = int(np.where(pred_values == pred_label)[0][0])
            best_true = true_values[int(np.argmax(overlap[: len(true_values), col]))]
            mapping[pred_label] = best_true

    mapped = np.empty(y_pred.shape, dtype=true_values.dtype)
    for pred_label, true_label in mapping.items():
        mapped[y_pred == pred_label] = true_label
    return mapped


def _safe_silhouette(
    embedding: Optional[np.ndarray],
    y_pred: np.ndarray,
    seed: int = 42,
    sample_size: Optional[int] = 3000,
) -> float:
    if embedding is None:
        return float("nan")
    embedding = np.asarray(embedding, dtype=np.float32)
    y_pred = np.asarray(y_pred)
    n_clusters = len(np.unique(y_pred))
    if embedding.shape[0] < 3 or n_clusters < 2 or n_clusters >= embedding.shape[0]:
        return float("nan")
    try:
        if sample_size is None or int(sample_size) <= 0:
            size = embedding.shape[0]
        else:
            size = min(int(sample_size), embedding.shape[0])
        kwargs = {"random_state": seed} if size < embedding.shape[0] else {}
        return float(silhouette_score(embedding, y_pred, sample_size=size, **kwargs))
    except Exception:
        return float("nan")


def compute_metrics(
    y_true,
    y_pred,
    embedding: Optional[np.ndarray] = None,
    seed: int = 42,
    silhouette_sample_size: Optional[int] = 3000,
) -> Tuple[dict, np.ndarray]:
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    mapped = best_map(y_true, y_pred)
    out = {
        "acc": float(np.mean(mapped == y_true)),
        "nmi": float(normalized_mutual_info_score(y_true, y_pred)),
        "ari": float(adjusted_rand_score(y_true, y_pred)),
        "f1_macro": float(f1_score(y_true, mapped, average="macro", zero_division=0)),
        "fmi": float(fowlkes_mallows_score(y_true, y_pred)),
        "v_measure": float(v_measure_score(y_true, y_pred)),
        "homogeneity": float(homogeneity_score(y_true, y_pred)),
        "completeness": float(completeness_score(y_true, y_pred)),
        "n_pred_clusters": int(len(np.unique(y_pred))),
        "silhouette": _safe_silhouette(embedding, y_pred, seed=seed, sample_size=silhouette_sample_size),
    }
    return out, mapped
