from __future__ import annotations

import numpy as np
import torch

from experimental_retired_models.CutAware_NeighborMix_scMAE.neighbor_graph import NeighborGraph


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


def make_gated_neighbor_mixed_batch(
    data_np: np.ndarray,
    batch_indices: np.ndarray,
    batch_x: torch.Tensor,
    graph: NeighborGraph,
    edge_weights: np.ndarray,
    model,
    anchor_latent: torch.Tensor,
    anchor_probs: torch.Tensor,
    mix_neighbors: int,
    alpha: float,
    gate_temperature: float,
    gate_min: float,
    cluster_temperature: float,
) -> tuple[torch.Tensor, torch.Tensor, dict, dict[str, torch.Tensor]]:
    """NeighborMix branch with learnable per-edge gates.

    Static graph weights, including cut-reweighted weights, remain the prior.
    The gate learns a second-stage reliability score from anchor/neighbor
    latents plus graph features. This mirrors recent attention/gating graph
    methods while keeping NeighborMix boundary-aware.
    """
    if graph.indices.shape[1] == 0 or mix_neighbors <= 0:
        weight = torch.ones(batch_x.shape[0], dtype=batch_x.dtype, device=batch_x.device)
        zero = torch.zeros((), dtype=batch_x.dtype, device=batch_x.device)
        return batch_x, weight, {"mean_mix_delta": 0.0, "mean_mix_neighbor_weight": 0.0, "mean_gate": 0.0}, {
            "gate_prior_loss": zero,
            "gate_entropy_loss": zero,
            "gate_cluster_loss": zero,
        }

    idx = np.asarray(batch_indices, dtype=np.int64)
    k = min(int(mix_neighbors), graph.indices.shape[1])
    nb = graph.indices[idx, :k]
    base_w = edge_weights[idx, :k].astype(np.float32)
    base_w = base_w / np.clip(base_w.sum(axis=1, keepdims=True), 1e-12, None)
    base_w_t = torch.as_tensor(base_w, dtype=batch_x.dtype, device=batch_x.device)

    nb_x = torch.as_tensor(data_np[nb.reshape(-1)], dtype=batch_x.dtype, device=batch_x.device)
    nb_latent = model.encoder(nb_x)
    src_latent = anchor_latent.repeat_interleave(k, dim=0)

    edge_feature_np = np.stack(
        [
            base_w.reshape(-1),
            graph.similarity[idx, :k].reshape(-1),
            graph.distance[idx, :k].reshape(-1),
            graph.mutual[idx, :k].astype(np.float32).reshape(-1),
            graph.snn[idx, :k].reshape(-1),
        ],
        axis=1,
    ).astype(np.float32)
    edge_features = torch.as_tensor(edge_feature_np, dtype=batch_x.dtype, device=batch_x.device)
    gate = model.edge_gate_scores(src_latent, nb_latent, edge_features, temperature=gate_temperature).view(batch_x.shape[0], k)
    gate_min = float(np.clip(gate_min, 0.0, 1.0))
    gated_raw = base_w_t * (gate_min + (1.0 - gate_min) * gate)
    gated_w = gated_raw / gated_raw.sum(dim=1, keepdim=True).clamp_min(1e-8)

    nb_expr = nb_x.view(batch_x.shape[0], k, batch_x.shape[1])
    mixed_neighbor = torch.sum(nb_expr * gated_w[:, :, None], dim=1)
    alpha = float(np.clip(alpha, 0.0, 1.0))
    x_prime = alpha * batch_x + (1.0 - alpha) * mixed_neighbor

    with torch.no_grad():
        prior = base_w_t / base_w_t.max(dim=1, keepdim=True).values.clamp_min(1e-8)
    gate_prior_loss = torch.mean((gate - prior) ** 2)
    gate_entropy_loss = -torch.mean(gate * torch.log(gate.clamp_min(1e-8)) + (1.0 - gate) * torch.log((1.0 - gate).clamp_min(1e-8)))
    q_dst = torch.softmax(model.cluster_logits(nb_latent) / max(float(cluster_temperature), 1e-6), dim=1).view(batch_x.shape[0], k, -1)
    same_prob = torch.sum(anchor_probs[:, None, :] * q_dst, dim=2)
    gate_cluster_loss = torch.mean(torch.sum(gated_w * (1.0 - same_prob), dim=1))

    delta = torch.linalg.vector_norm(x_prime - batch_x, dim=1)
    entropy = -torch.sum(gated_w * torch.log(gated_w.clamp_min(1e-8)), dim=1)
    sample_weight = torch.ones(batch_x.shape[0], dtype=batch_x.dtype, device=batch_x.device)
    return x_prime, sample_weight, {
        "mean_mix_delta": float(delta.mean().detach().cpu()),
        "mean_mix_neighbor_weight": float(gated_w.max(dim=1).values.mean().detach().cpu()),
        "mean_gate": float(gate.mean().detach().cpu()),
        "min_gate": float(gate.min().detach().cpu()),
        "max_gate": float(gate.max().detach().cpu()),
        "gate_effective_neighbors": float(torch.exp(entropy).mean().detach().cpu()),
    }, {
        "gate_prior_loss": gate_prior_loss,
        "gate_entropy_loss": gate_entropy_loss,
        "gate_cluster_loss": gate_cluster_loss,
    }
