from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize


@dataclass
class NeighborGraph:
    indices: np.ndarray
    probs: np.ndarray
    similarity: np.ndarray
    distance: np.ndarray
    embedding: np.ndarray
    mutual: np.ndarray
    snn: np.ndarray
    degree: np.ndarray
    profile: dict


def _safe_normalize(embedding: np.ndarray) -> np.ndarray:
    emb = np.nan_to_num(np.asarray(embedding, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    return normalize(emb, axis=1).astype(np.float32)


def build_pca_knn_graph(data_np: np.ndarray, k: int, pca_dim: int, tau: float, seed: int) -> NeighborGraph:
    data = np.asarray(data_np, dtype=np.float32)
    n_cells, n_genes = data.shape
    if k <= 0 or n_cells <= 1:
        return _empty_graph(n_cells)
    dim = max(1, min(int(pca_dim), n_genes, n_cells - 1))
    emb = PCA(n_components=dim, random_state=seed).fit_transform(data) if dim < min(data.shape) else data
    return build_embedding_knn_graph(emb, k=k, tau=tau, source="pca_input", pca_dim=dim)


def build_embedding_knn_graph(embedding: np.ndarray, k: int, tau: float, source: str = "embedding", pca_dim: int | None = None) -> NeighborGraph:
    emb = _safe_normalize(embedding)
    n_cells = int(emb.shape[0])
    if k <= 0 or n_cells <= 1:
        return _empty_graph(n_cells)

    k_eff = min(int(k), n_cells - 1)
    nn = NearestNeighbors(n_neighbors=k_eff + 1, metric="cosine")
    nn.fit(emb)
    distances, indices = nn.kneighbors(emb)
    indices = indices[:, 1 : k_eff + 1].astype(np.int64, copy=False)
    distances = distances[:, 1 : k_eff + 1].astype(np.float32, copy=False)
    similarity = (1.0 - distances).astype(np.float32)
    scaled = similarity / max(float(tau), 1e-8)
    scaled = scaled - scaled.max(axis=1, keepdims=True)
    exp_scaled = np.exp(scaled).astype(np.float32)
    probs = exp_scaled / np.clip(exp_scaled.sum(axis=1, keepdims=True), 1e-12, None)

    neighbor_sets = [set(row.tolist()) for row in indices]
    mutual = np.zeros_like(indices, dtype=bool)
    snn = np.zeros_like(similarity, dtype=np.float32)
    for i in range(n_cells):
        set_i = neighbor_sets[i]
        for pos, j in enumerate(indices[i]):
            mutual[i, pos] = i in neighbor_sets[j]
            union = set_i.union(neighbor_sets[j])
            snn[i, pos] = len(set_i.intersection(neighbor_sets[j])) / float(max(1, len(union)))

    degree = probs.sum(axis=1).astype(np.float32)
    profile = {
        "graph_source": source,
        "neighbor_k": int(k_eff),
        "tau": float(tau),
        "knn_pca_dim": int(pca_dim) if pca_dim is not None else None,
        "mean_neighbor_similarity": float(np.mean(similarity)),
        "mean_mutual_ratio": float(np.mean(mutual)),
        "mean_snn": float(np.mean(snn)),
        "mean_max_neighbor_prob": float(np.mean(np.max(probs, axis=1))),
    }
    return NeighborGraph(indices, probs.astype(np.float32), similarity, distances, emb, mutual, snn, degree, profile)


def compute_edge_weights(
    graph: NeighborGraph,
    mode: str,
    gamma_sim: float,
    gamma_mutual: float,
    gamma_snn: float,
    gamma_distance: float,
    prune_quantile: float = 0.0,
) -> tuple[np.ndarray, np.ndarray, dict]:
    if graph.indices.shape[1] == 0 or mode == "none":
        weights = graph.probs.copy()
        rel = np.ones_like(weights, dtype=np.float32)
        return rel, weights, summarize_edge_weights(weights, prune_quantile=0.0)

    rel = np.ones_like(graph.similarity, dtype=np.float32)
    if mode in {"sim", "sim_mutual", "sim_mutual_snn", "sim_mutual_snn_distance"}:
        rel *= np.exp(float(gamma_sim) * graph.similarity).astype(np.float32)
    if mode in {"sim_mutual", "sim_mutual_snn", "sim_mutual_snn_distance"}:
        rel *= 1.0 + float(gamma_mutual) * graph.mutual.astype(np.float32)
    if mode in {"sim_mutual_snn", "sim_mutual_snn_distance"}:
        rel *= 1.0 + float(gamma_snn) * graph.snn
    if mode == "sim_mutual_snn_distance":
        rel *= np.exp(-float(gamma_distance) * graph.distance).astype(np.float32)
    rel = np.clip(rel, 1e-6, 1e6).astype(np.float32)
    raw = graph.probs * rel

    pruned_fraction = 0.0
    if prune_quantile > 0.0 and raw.size:
        threshold = np.quantile(raw, float(prune_quantile))
        keep = raw >= threshold
        pruned_fraction = float(1.0 - np.mean(keep))
        raw = raw * keep.astype(np.float32)

    weights = raw / np.clip(raw.sum(axis=1, keepdims=True), 1e-12, None)
    return rel, weights.astype(np.float32), summarize_edge_weights(weights, prune_quantile=prune_quantile, pruned_fraction=pruned_fraction)


def apply_cluster_cut_reweight(
    graph: NeighborGraph,
    edge_weights: np.ndarray,
    n_clusters: int,
    cross_weight: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    """
    Cut-informed edge reweighting without backpropagating a cut loss.

    A provisional unsupervised graph partition is built from the graph embedding.
    Edges crossing this partition are treated as candidate cut edges and receive a
    smaller mixing weight. This adds the scCDCG-style "cut bad edges" idea to
    NeighborMix while avoiding unstable direct NCut/OT gradients.
    """
    weights = np.asarray(edge_weights, dtype=np.float32)
    if graph.indices.shape[1] == 0 or int(n_clusters) <= 1:
        labels = np.zeros(graph.indices.shape[0], dtype=np.int64)
        return labels, weights, {"cut_reweight_enabled": False}
    k = min(int(n_clusters), graph.embedding.shape[0])
    labels = KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(graph.embedding).astype(np.int64)
    same = labels[:, None] == labels[graph.indices]
    cross_weight = float(np.clip(cross_weight, 0.0, 1.0))
    factors = np.where(same, 1.0, cross_weight).astype(np.float32)
    raw = weights * factors
    row_sum = raw.sum(axis=1, keepdims=True)
    safe_raw = np.where(row_sum > 1e-12, raw / np.clip(row_sum, 1e-12, None), weights)
    cross_mask = ~same
    before_cross_mass = float(np.sum(weights * cross_mask) / max(1, weights.shape[0]))
    after_cross_mass = float(np.sum(safe_raw * cross_mask) / max(1, weights.shape[0]))
    summary = {
        "cut_reweight_enabled": True,
        "cut_partition_source": "kmeans_on_graph_embedding",
        "cut_cross_weight": cross_weight,
        "cut_partition_n_clusters": int(k),
        "fraction_candidate_cross_edges": float(np.mean(cross_mask)),
        "cross_edge_mass_before": before_cross_mass,
        "cross_edge_mass_after": after_cross_mass,
        "cross_edge_mass_reduction": float(before_cross_mass - after_cross_mass),
    }
    return labels, safe_raw.astype(np.float32), summary


def summarize_edge_weights(weights: np.ndarray, prune_quantile: float = 0.0, pruned_fraction: float = 0.0) -> dict:
    if weights.size == 0:
        return {
            "edge_weight_entropy": 0.0,
            "effective_neighbor_count": 0.0,
            "max_edge_weight_mean": 0.0,
            "max_edge_weight_p95": 0.0,
            "fraction_effective_neighbors_lt_2": 1.0,
            "prune_quantile": float(prune_quantile),
            "pruned_fraction": float(pruned_fraction),
        }
    entropy = -np.sum(weights * np.log(np.clip(weights, 1e-12, None)), axis=1)
    effective = np.exp(entropy)
    max_w = np.max(weights, axis=1)
    return {
        "edge_weight_entropy": float(np.mean(entropy)),
        "effective_neighbor_count": float(np.mean(effective)),
        "max_edge_weight_mean": float(np.mean(max_w)),
        "max_edge_weight_p95": float(np.percentile(max_w, 95)),
        "fraction_effective_neighbors_lt_2": float(np.mean(effective < 2.0)),
        "prune_quantile": float(prune_quantile),
        "pruned_fraction": float(pruned_fraction),
    }


def neighbor_tensors_for_batch(
    data_np: np.ndarray,
    batch_indices: np.ndarray,
    graph: NeighborGraph,
    edge_weights: np.ndarray,
    max_neighbors: int,
    device,
):
    import torch

    idx = np.asarray(batch_indices, dtype=np.int64)
    if graph.indices.shape[1] == 0 or max_neighbors <= 0:
        return None, None
    k = min(int(max_neighbors), graph.indices.shape[1])
    nb = graph.indices[idx, :k].reshape(-1)
    w = edge_weights[idx, :k].reshape(-1)
    src_rep = np.repeat(np.arange(idx.shape[0]), k)
    nb_x = torch.as_tensor(data_np[nb], dtype=torch.float32, device=device)
    src_rep_t = torch.as_tensor(src_rep, dtype=torch.long, device=device)
    weight_t = torch.as_tensor(w, dtype=torch.float32, device=device)
    return nb_x, src_rep_t, weight_t


def graph_cut_diagnostics(q: np.ndarray, graph: NeighborGraph, edge_weights: np.ndarray) -> dict:
    q = np.asarray(q, dtype=np.float32)
    if q.ndim != 2 or graph.indices.shape[1] == 0:
        return {
            "soft_cut": 0.0,
            "soft_association": 0.0,
            "ncut_surrogate": 0.0,
            "cluster_mass_min": 0.0,
            "cluster_mass_max": 0.0,
            "assignment_entropy": 0.0,
            "max_cluster_fraction": 0.0,
        }
    src = np.repeat(np.arange(graph.indices.shape[0]), graph.indices.shape[1])
    dst = graph.indices.reshape(-1)
    w = edge_weights.reshape(-1).astype(np.float32)
    same = np.sum(q[src] * q[dst], axis=1)
    soft_cut = float(np.sum(w * (1.0 - same)) / max(1, q.shape[0]))
    soft_assoc = float(np.sum(w * same) / max(1, q.shape[0]))
    degree = edge_weights.sum(axis=1).astype(np.float32)
    assoc = np.clip((degree[:, None] * q).sum(axis=0), 1e-8, None)
    cut_by_cluster = np.sum(w[:, None] * q[src] * (1.0 - q[dst]), axis=0)
    ncut = float(np.mean(cut_by_cluster / assoc))
    mean_q = q.mean(axis=0)
    pred = np.argmax(q, axis=1)
    return {
        "soft_cut": soft_cut,
        "soft_association": soft_assoc,
        "ncut_surrogate": ncut,
        "cluster_mass_min": float(np.min(mean_q)),
        "cluster_mass_max": float(np.max(mean_q)),
        "assignment_entropy": float(-np.mean(np.sum(q * np.log(np.clip(q, 1e-12, None)), axis=1))),
        "max_cluster_fraction": float(np.max(np.bincount(pred, minlength=q.shape[1])) / max(1, q.shape[0])),
    }


def _empty_graph(n_cells: int) -> NeighborGraph:
    empty_i = np.zeros((n_cells, 0), dtype=np.int64)
    empty_f = np.zeros((n_cells, 0), dtype=np.float32)
    return NeighborGraph(
        empty_i,
        empty_f,
        empty_f,
        empty_f,
        np.zeros((n_cells, 0), dtype=np.float32),
        empty_i.astype(bool),
        empty_f,
        np.zeros(n_cells, dtype=np.float32),
        {"neighbor_k": 0, "graph_source": "none"},
    )
