from __future__ import annotations

import numpy as np


def embedding_stats(z: np.ndarray) -> dict:
    z = np.asarray(z, dtype=np.float32)
    var = z.var(axis=0)
    centered = z - z.mean(axis=0, keepdims=True)
    try:
        s = np.linalg.svd(centered, compute_uv=False)
        p = s / (s.sum() + 1.0e-8)
        rank = float(np.exp(-(p * np.log(p + 1.0e-8)).sum()))
    except Exception:
        rank = 0.0
    norm = np.linalg.norm(z, axis=1, keepdims=True) + 1.0e-8
    zn = z / norm
    sample = zn[: min(1000, zn.shape[0])]
    cosine = float((sample @ sample.T).mean()) if sample.size else 0.0
    return {
        "per_dimension_variance_mean": float(var.mean()) if var.size else 0.0,
        "per_dimension_variance_min": float(var.min()) if var.size else 0.0,
        "effective_rank": rank,
        "mean_pairwise_cosine": cosine,
        "embedding_norm_mean": float(np.linalg.norm(z, axis=1).mean()) if z.size else 0.0,
    }

