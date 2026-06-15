from __future__ import annotations

import numpy as np
import torch

from methods.DeepLearning.CutAware_NeighborMix_scMAE.neighbor_graph import NeighborGraph


def make_neighbor_mixed_batch(
    data_np: np.ndarray,
    batch_indices: np.ndarray,
    batch_x: torch.Tensor,
    graph: NeighborGraph,
    edge_weights: np.ndarray,
    mix_neighbors: int,
    alpha: float,
) -> tuple[torch.Tensor, torch.Tensor, dict]:
    """Conservative NeighborMix expression branch for mix-plus-cut ablations."""
    if graph.indices.shape[1] == 0 or mix_neighbors <= 0:
        weight = torch.ones(batch_x.shape[0], dtype=batch_x.dtype, device=batch_x.device)
        return batch_x, weight, {"mean_mix_delta": 0.0, "mean_mix_neighbor_weight": 0.0}

    idx = np.asarray(batch_indices, dtype=np.int64)
    k = min(int(mix_neighbors), graph.indices.shape[1])
    nb = graph.indices[idx, :k]
    w = edge_weights[idx, :k].astype(np.float32)
    w = w / np.clip(w.sum(axis=1, keepdims=True), 1e-12, None)
    nb_expr = data_np[nb]
    mixed_neighbor = np.sum(nb_expr * w[:, :, None], axis=1).astype(np.float32)
    neighbor_t = torch.as_tensor(mixed_neighbor, dtype=batch_x.dtype, device=batch_x.device)
    alpha = float(np.clip(alpha, 0.0, 1.0))
    x_prime = alpha * batch_x + (1.0 - alpha) * neighbor_t
    delta = torch.linalg.vector_norm(x_prime - batch_x, dim=1)
    sample_weight = torch.ones(batch_x.shape[0], dtype=batch_x.dtype, device=batch_x.device)
    return x_prime, sample_weight, {
        "mean_mix_delta": float(delta.mean().detach().cpu()),
        "mean_mix_neighbor_weight": float(np.mean(np.max(w, axis=1))),
    }
