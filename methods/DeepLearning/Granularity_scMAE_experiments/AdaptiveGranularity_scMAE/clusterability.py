from __future__ import annotations

import numpy as np


def compute_clusterability(
    embedding: np.ndarray,
    q: np.ndarray,
    neighbor_indices: np.ndarray | None,
    k: int = 15,
) -> tuple[np.ndarray, dict]:
    """Per-cell local clusterability c_i in [0,1] = neighbor_agreement * density.

      agree_i   : fraction of KNN neighbors whose hard assignment (argmax q)
                  matches this cell's own assignment. High = sits inside a
                  coherent cluster core; low = boundary between clusters.
      density_i : local density proxy (inverse mean neighbor distance, robustly
                  normalized). Low = rare / isolated cell.

    c_i high  -> cell is in a clean, dense cluster core -> safe to SHARPEN.
    c_i low   -> boundary or rare cell -> keep SOFT (fuzzy-rough boundary).

    Membership confidence is deliberately excluded (it collapses the field
    early in training when clusters are unformed; verified in the prior line).
    """
    n = int(embedding.shape[0])
    q = np.asarray(q, dtype=np.float64)
    assign = q.argmax(axis=1)
    confidence = q.max(axis=1).astype(np.float64)

    if neighbor_indices is None or neighbor_indices.shape[1] == 0 or n <= 2:
        agree = np.ones(n, dtype=np.float64)
        density = np.ones(n, dtype=np.float64)
    else:
        kk = min(int(k), neighbor_indices.shape[1])
        nb = neighbor_indices[:, :kk]
        agree = (assign[nb] == assign[:, None]).mean(axis=1).astype(np.float64)
        mean_dist = np.empty(n, dtype=np.float64)
        step = 4096
        for s in range(0, n, step):
            e = min(s + step, n)
            diff = embedding[s:e, None, :] - embedding[nb[s:e]]
            mean_dist[s:e] = np.linalg.norm(diff, axis=2).mean(axis=1)
        lo, hi = np.percentile(mean_dist, [5.0, 95.0])
        density = np.ones(n, dtype=np.float64) if hi - lo < 1e-8 else 1.0 - np.clip((mean_dist - lo) / (hi - lo), 0.0, 1.0)

    c = np.clip(agree * density, 0.0, 1.0).astype(np.float32)
    diag = {
        "clusterability_mean": float(c.mean()),
        "clusterability_min": float(c.min()),
        "clusterability_max": float(c.max()),
        "core_fraction": float((c >= 0.5).mean()),
        "agree_mean": float(agree.mean()),
        "density_mean": float(density.mean()),
        "confidence_mean": float(confidence.mean()),
    }
    return c, diag


def adaptive_target(q: np.ndarray, sharp_p: np.ndarray, c: np.ndarray) -> np.ndarray:
    """Per-cell interpolation between the sharp DEC target and the soft q.

        p_i^adaptive = c_i * sharp_p_i + (1 - c_i) * q_i

    Core cells follow the sharp DEC target; boundary/rare cells target their own
    q (KL(q||q)=0 -> no clustering pressure, left fuzzy).
    """
    c = np.asarray(c, dtype=np.float32)[:, None]
    p = c * np.asarray(sharp_p, dtype=np.float32) + (1.0 - c) * np.asarray(q, dtype=np.float32)
    return (p / np.clip(p.sum(axis=1, keepdims=True), 1e-8, None)).astype(np.float32)
