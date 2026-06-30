from __future__ import annotations

import numpy as np


def _mean_cosine(embedding: np.ndarray, max_cells: int = 2048) -> float:
    z = np.asarray(embedding, dtype=np.float32)
    if z.shape[0] > max_cells:
        rng = np.random.default_rng(0)
        z = z[rng.choice(z.shape[0], size=max_cells, replace=False)]
    denom = np.linalg.norm(z, axis=1, keepdims=True)
    z = z / np.maximum(denom, 1e-12)
    sim = z @ z.T
    if sim.shape[0] <= 1:
        return 1.0
    mask = ~np.eye(sim.shape[0], dtype=bool)
    return float(np.mean(sim[mask]))


def build_diagnostics(
    embedding: np.ndarray,
    pred_labels: np.ndarray | None,
    configured_mask_rate: float,
    effective_mask_rate: float,
    prototype_confidence_mean: float | None,
    neighbor_stats: dict | None,
    loss_components: dict,
) -> dict:
    z = np.nan_to_num(np.asarray(embedding, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    var = np.var(z, axis=0) if z.size else np.array([0.0], dtype=np.float32)
    if pred_labels is None:
        cluster_mass_max = 0.0
        cluster_mass_min = 0.0
        n_pred_clusters = 0
    else:
        _, counts = np.unique(pred_labels, return_counts=True)
        mass = counts.astype(np.float32) / max(1, int(np.sum(counts)))
        cluster_mass_max = float(np.max(mass)) if mass.size else 0.0
        cluster_mass_min = float(np.min(mass)) if mass.size else 0.0
        n_pred_clusters = int(len(counts))
    neighbor_stats = neighbor_stats or {}
    out = {
        "configured_mask_rate": float(configured_mask_rate),
        "effective_mask_rate": float(effective_mask_rate),
        "embedding_variance_mean": float(np.mean(var)),
        "embedding_variance_min": float(np.min(var)),
        "embedding_mean_cosine": _mean_cosine(z),
        "cluster_mass_max": cluster_mass_max,
        "cluster_mass_min": cluster_mass_min,
        "n_pred_clusters": n_pred_clusters,
        "prototype_confidence_mean": (
            float(prototype_confidence_mean) if prototype_confidence_mean is not None else float("nan")
        ),
        "neighbor_reliability_mean": float(neighbor_stats.get("neighbor_reliability_mean", float("nan"))),
        "neighbor_reliability_min": float(neighbor_stats.get("neighbor_reliability_min", float("nan"))),
        "neighbor_reliable_edge_fraction": float(neighbor_stats.get("neighbor_reliable_edge_fraction", float("nan"))),
        "loss_components": loss_components,
    }
    for key, value in neighbor_stats.items():
        if key.startswith("neighbor_") and key not in out:
            out[key] = float(value) if isinstance(value, (int, float, np.integer, np.floating)) else value
    return out
