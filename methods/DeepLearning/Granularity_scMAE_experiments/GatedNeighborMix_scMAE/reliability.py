from __future__ import annotations

import numpy as np


def compute_reliability(
    embedding: np.ndarray,
    q: np.ndarray,
    neighbor_indices: np.ndarray | None,
    k: int = 15,
    density_percentile: float = 20.0,
) -> tuple[np.ndarray, dict]:
    """Per-cell reliability field r_i in [0,1] = neighbor_agreement * local_density.

      agree_i    : fraction of KNN neighbors whose hard cluster assignment
                   (argmax q) matches this cell's own assignment. Low on
                   boundary cells whose neighbors disagree.
      density_i  : local density proxy from mean neighbor distance, robustly
                   min-max normalized. Rare/isolated cells sit in low-density
                   regions -> low density_i -> mixing suppressed.

    NOTE: absolute membership confidence max(q) is intentionally NOT multiplied
    in here. Early in training clusters are unformed and max(q)≈1/n_clusters for
    ALL cells, which would collapse the whole gate to ~0 and switch off both
    NeighborMix and DEC everywhere. The confidence dimension is already handled
    downstream by the DEC loss's `confidence_threshold` gate, so keeping it out
    of r_i avoids double-penalizing and avoids the early-training shutdown.
    Confidence is still returned as a diagnostic.

    A high r_i means "my neighborhood is trustworthy, smoothing is safe here".
    A low r_i means "I am rare/boundary, fall back to pure scMAE".

    Returns (r in [0,1] shape [n_cells], diagnostics dict).
    """
    n_cells = int(embedding.shape[0])
    q = np.asarray(q, dtype=np.float64)
    assign = q.argmax(axis=1)
    confidence = q.max(axis=1).astype(np.float64)

    if neighbor_indices is None or neighbor_indices.shape[1] == 0 or n_cells <= 2:
        agree = np.ones(n_cells, dtype=np.float64)
        density = np.ones(n_cells, dtype=np.float64)
    else:
        kk = min(int(k), neighbor_indices.shape[1])
        nb = neighbor_indices[:, :kk]
        agree = (assign[nb] == assign[:, None]).mean(axis=1).astype(np.float64)
        # local density from mean distance to the kk neighbors (smaller = denser).
        # Chunked to avoid a single (n_cells, kk, dim) broadcast that blows up
        # memory / OpenBLAS thread buffers on large datasets (e.g. 44k cells).
        mean_dist = np.empty(n_cells, dtype=np.float64)
        step = 4096
        for s in range(0, n_cells, step):
            e = min(s + step, n_cells)
            diff = embedding[s:e, None, :] - embedding[nb[s:e]]
            mean_dist[s:e] = np.linalg.norm(diff, axis=2).mean(axis=1)
        lo, hi = np.percentile(mean_dist, [5.0, 95.0])
        if hi - lo < 1e-8:
            density = np.ones(n_cells, dtype=np.float64)
        else:
            # invert: dense (small dist) -> high density score
            density = 1.0 - np.clip((mean_dist - lo) / (hi - lo), 0.0, 1.0)

    r = np.clip(agree * density, 0.0, 1.0).astype(np.float32)

    diag = {
        "reliability_mean": float(r.mean()),
        "reliability_min": float(r.min()),
        "reliability_max": float(r.max()),
        "core_fraction": float((r >= 0.5).mean()),
        "agree_mean": float(agree.mean()),
        "density_mean": float(density.mean()),
        "confidence_mean": float(confidence.mean()),
    }
    return r, diag

