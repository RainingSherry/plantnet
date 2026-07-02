from __future__ import annotations

import numpy as np
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


def compute_local_reliability(target_expr, graph_space, k=15, pca_dim=50, floor=0.2, seed=42):
    """Per-cell-per-gene LOCAL reliability r_ig in [floor,1] from a DECOUPLED
    raw-data PCA-KNN graph (never the live embedding -> no circular collapse).

    r_ig high when gene g varies little among cell i's raw-KNN neighbors
    (trustworthy), low when it fluctuates (technical noise/dropout). Relative to
    per-gene median local variance (each gene judged vs its own baseline) so
    low-count rare markers are not globally flagged as noise; floor guarantees
    no gene/cell fully ignored.
    """
    n, g = target_expr.shape
    max_k = min(int(k), max(1, n - 1))
    dim = min(int(pca_dim), min(graph_space.shape) - 1)
    emb = PCA(n_components=max(2, dim), random_state=seed).fit_transform(graph_space.astype(np.float64)) if dim >= 2 else graph_space.astype(np.float64)
    nn = NearestNeighbors(n_neighbors=max_k + 1).fit(emb)
    idx = nn.kneighbors(emb, return_distance=False)[:, 1:max_k + 1]
    tgt = target_expr.astype(np.float32)
    local_var = np.empty((n, g), dtype=np.float32)
    for s in range(0, n, 2048):
        e = min(s + 2048, n)
        local_var[s:e] = tgt[idx[s:e]].var(axis=1)
    gene_med = np.median(local_var, axis=0, keepdims=True)
    rel = local_var / np.clip(gene_med, 1e-8, None)
    r = 1.0 / (1.0 + np.clip(rel - 1.0, 0.0, None))
    r = np.clip(floor + (1.0 - floor) * r, floor, 1.0).astype(np.float32)
    diag = {"reliability_mean": float(r.mean()), "reliability_min": float(r.min()),
            "reliability_p05": float(np.percentile(r, 5)), "reliability_p95": float(np.percentile(r, 95)),
            "neighbor_k": int(max_k)}
    return r, diag
