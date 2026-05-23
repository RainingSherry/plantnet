import json
import os
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
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
from sklearn.preprocessing import LabelEncoder


LABEL_CANDIDATES = [
    "cell_type",
    "Celltype",
    "celltype",
    "cell_label",
    "label",
    "Cluster",
    "cluster",
    "clusters",
    "Seurat_clusters",
]


def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return str(path)


def to_dense(x):
    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def find_label_key(adata):
    for key in LABEL_CANDIDATES:
        if key in adata.obs.columns:
            return key
    raise KeyError(f"No label column found. Available obs columns: {list(adata.obs.columns)}")


def labels_from_adata(adata):
    key = find_label_key(adata)
    labels = adata.obs[key].astype(str).to_numpy()
    encoded = LabelEncoder().fit_transform(labels).astype(np.int64)
    return encoded, key


def best_map(y_true, y_pred):
    true_values = np.unique(y_true)
    pred_values = np.unique(y_pred)
    n = max(len(true_values), len(pred_values))
    mat = np.zeros((n, n), dtype=np.int64)
    for i, true_label in enumerate(true_values):
        for j, pred_label in enumerate(pred_values):
            mat[i, j] = np.sum((y_true == true_label) & (y_pred == pred_label))
    rows, cols = linear_sum_assignment(-mat)
    mapped = np.zeros_like(y_pred, dtype=np.int64)
    for row, col in zip(rows, cols):
        if row < len(true_values) and col < len(pred_values):
            mapped[y_pred == pred_values[col]] = true_values[row]
    return mapped


def compute_metrics(y_true, y_pred, embedding=None):
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
    }
    out["silhouette"] = float("nan")
    return out, mapped


def leiden_labels(embedding, n_clusters, seed=42, n_neighbors=15):
    adata = sc.AnnData(np.asarray(embedding, dtype=np.float32))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    resolutions = np.asarray([1.0], dtype=float)
    best_res = resolutions[0]
    best_diff = float("inf")
    for res in resolutions:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            try:
                sc.tl.leiden(
                    adata,
                    random_state=seed,
                    resolution=float(res),
                    key_added="tmp",
                    flavor="igraph",
                    n_iterations=2,
                    directed=False,
                )
            except TypeError:
                sc.tl.leiden(adata, random_state=seed, resolution=float(res), key_added="tmp")
        count = adata.obs["tmp"].nunique()
        diff = abs(count - n_clusters)
        if diff < best_diff:
            best_diff = diff
            best_res = float(res)
        if count == n_clusters:
            break
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        try:
            sc.tl.leiden(
                adata,
                random_state=seed,
                resolution=best_res,
                key_added="leiden",
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
        except TypeError:
            sc.tl.leiden(adata, random_state=seed, resolution=best_res, key_added="leiden")
    return adata.obs["leiden"].astype(int).to_numpy(), best_res


def evaluate_embedding(embedding, labels, n_clusters=None, seed=42, n_neighbors=15):
    if n_clusters is None:
        n_clusters = len(np.unique(labels))
    results = {}
    km = KMeans(n_clusters=n_clusters, n_init=1, random_state=seed)
    pred_kmeans = km.fit_predict(embedding)
    results["kmeans"], mapped_kmeans = compute_metrics(labels, pred_kmeans, embedding)
    pred_leiden, res = leiden_labels(embedding, n_clusters, seed=seed, n_neighbors=n_neighbors)
    results["leiden"], mapped_leiden = compute_metrics(labels, pred_leiden, embedding)
    results["leiden"]["resolution"] = float(res)
    return results, {
        "kmeans": pred_kmeans.astype(np.int64),
        "kmeans_mapped": mapped_kmeans.astype(np.int64),
        "leiden": pred_leiden.astype(np.int64),
        "leiden_mapped": mapped_leiden.astype(np.int64),
    }


def save_json(obj, path):
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2, default=str)


def append_rows_csv(rows, path):
    df = pd.DataFrame(rows)
    if os.path.exists(path):
        old = pd.read_csv(path)
        df = pd.concat([old, df], ignore_index=True)
    df.to_csv(path, index=False)
