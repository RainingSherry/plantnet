from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors


@dataclass
class NeighborState:
    indices: np.ndarray
    reliability: np.ndarray
    similarity: np.ndarray
    shared_score: np.ndarray
    eligible: np.ndarray
    first_reliable: np.ndarray
    mean_reliable_count: float
    stats: dict


def _similarity_from_distance(distances: np.ndarray, metric: str) -> np.ndarray:
    if metric == "cosine":
        return np.clip(1.0 - distances.astype(np.float32), -1.0, 1.0)
    finite = distances[np.isfinite(distances)]
    scale = float(np.median(finite)) if finite.size else 1.0
    scale = max(scale, 1e-6)
    return np.exp(-distances.astype(np.float32) / scale)


def build_neighbor_state(
    embedding: np.ndarray,
    k: int,
    metric: str = "cosine",
    min_similarity: float = -1.0,
    min_shared_score: float = 0.0,
    score_threshold: float = 0.0,
    similarity_weight: float = 0.7,
    pseudo_filter: str = "none",
    pseudo_n_clusters: int | None = None,
    pseudo_seed: int = 0,
    pseudo_confidence_quantile: float = 0.0,
) -> NeighborState:
    z = np.nan_to_num(np.asarray(embedding, dtype=np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    n_cells = int(z.shape[0])
    pseudo_filter = str(pseudo_filter)
    pseudo_labels = None
    pseudo_confidence = None
    pseudo_confidence_threshold = 0.0
    if pseudo_filter != "none":
        if pseudo_n_clusters is None or int(pseudo_n_clusters) < 2:
            raise ValueError("pseudo_n_clusters must be at least 2 when NeighborMix pseudo filtering is enabled.")
        km = KMeans(n_clusters=int(pseudo_n_clusters), n_init=20, random_state=int(pseudo_seed))
        pseudo_labels = km.fit_predict(z).astype(np.int64)
        distances_to_centers = km.transform(z).astype(np.float32)
        if distances_to_centers.shape[1] > 1:
            nearest = np.partition(distances_to_centers, kth=1, axis=1)[:, :2]
            d1 = nearest[:, 0]
            d2 = np.maximum(nearest[:, 1], 1e-6)
            pseudo_confidence = np.clip(1.0 - d1 / d2, 0.0, 1.0).astype(np.float32)
        else:
            pseudo_confidence = np.ones(n_cells, dtype=np.float32)
        quantile = float(np.clip(pseudo_confidence_quantile, 0.0, 1.0))
        pseudo_confidence_threshold = float(np.quantile(pseudo_confidence, quantile))
    n_neighbors = max(2, min(int(k) + 1, n_cells))
    nn = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
    nn.fit(z)
    distances, neigh = nn.kneighbors(z, return_distance=True)
    rows = []
    dist_rows = []
    for i, (row, dist_row) in enumerate(zip(neigh, distances)):
        keep = row != i
        rows.append(row[keep][: int(k)])
        dist_rows.append(dist_row[keep][: int(k)])
    neigh = np.asarray(rows, dtype=np.int64)
    distances = np.asarray(dist_rows, dtype=np.float32)
    similarity = _similarity_from_distance(distances, metric=metric)
    neigh_sets = [set(row.tolist()) for row in neigh]
    reliability = np.zeros_like(neigh, dtype=np.float32)
    shared_score = np.zeros_like(neigh, dtype=np.float32)
    eligible = np.zeros_like(neigh, dtype=bool)
    first_reliable = np.arange(n_cells, dtype=np.int64)
    sim_weight = float(np.clip(similarity_weight, 0.0, 1.0))
    for i in range(n_cells):
        found = False
        for pos, j in enumerate(neigh[i]):
            if i in neigh_sets[int(j)]:
                shared = len(neigh_sets[i].intersection(neigh_sets[int(j)])) / max(1, int(k))
                shared_score[i, pos] = float(shared)
                sim_norm = 0.5 * (float(similarity[i, pos]) + 1.0) if metric == "cosine" else float(similarity[i, pos])
                score = sim_weight * sim_norm + (1.0 - sim_weight) * float(shared)
                reliability[i, pos] = float(score)
                passes = (
                    float(similarity[i, pos]) >= float(min_similarity)
                    and float(shared) >= float(min_shared_score)
                    and float(score) >= float(score_threshold)
                )
                if passes and pseudo_labels is not None:
                    same_cluster = pseudo_labels[i] == pseudo_labels[int(j)]
                    confident = (
                        float(pseudo_confidence[i]) >= pseudo_confidence_threshold
                        and float(pseudo_confidence[int(j)]) >= pseudo_confidence_threshold
                    )
                    if pseudo_filter == "same_cluster":
                        passes = bool(same_cluster)
                    elif pseudo_filter == "same_confident_cluster":
                        passes = bool(same_cluster and confident)
                    else:
                        raise ValueError(f"Unknown NeighborMix pseudo filter: {pseudo_filter}")
                eligible[i, pos] = passes
                if passes and not found:
                    first_reliable[i] = int(j)
                    found = True
    eligible_counts = eligible.sum(axis=1).astype(np.float32) if eligible.size else np.zeros(n_cells, dtype=np.float32)
    positive_scores = reliability[reliability > 0]
    eligible_scores = reliability[eligible]
    stats = {
        "neighbor_k": int(k),
        "neighbor_metric": metric,
        "neighbor_min_similarity": float(min_similarity),
        "neighbor_min_shared_score": float(min_shared_score),
        "neighbor_score_threshold": float(score_threshold),
        "neighbor_score_similarity_weight": float(sim_weight),
        "neighbor_pseudo_filter": pseudo_filter,
        "neighbor_pseudo_n_clusters": int(pseudo_n_clusters or 0),
        "neighbor_pseudo_confidence_quantile": float(np.clip(pseudo_confidence_quantile, 0.0, 1.0)),
        "neighbor_pseudo_confidence_threshold": float(pseudo_confidence_threshold),
        "neighbor_pseudo_confidence_mean": (
            float(np.mean(pseudo_confidence)) if pseudo_confidence is not None else 0.0
        ),
        "neighbor_reliability_mean": float(np.mean(eligible_scores)) if eligible_scores.size else 0.0,
        "neighbor_reliability_min": float(np.min(eligible_scores)) if eligible_scores.size else 0.0,
        "neighbor_reliability_all_edge_mean": float(np.mean(reliability)) if reliability.size else 0.0,
        "neighbor_mutual_edge_fraction": float(np.mean(reliability > 0)) if reliability.size else 0.0,
        "neighbor_reliable_edge_fraction": float(np.mean(eligible)) if eligible.size else 0.0,
        "neighbor_reliable_count_mean": float(np.mean(eligible_counts)) if eligible_counts.size else 0.0,
        "neighbor_reliable_count_min": float(np.min(eligible_counts)) if eligible_counts.size else 0.0,
        "neighbor_score_mean": float(np.mean(positive_scores)) if positive_scores.size else 0.0,
        "neighbor_score_p25": float(np.quantile(positive_scores, 0.25)) if positive_scores.size else 0.0,
        "neighbor_score_p50": float(np.quantile(positive_scores, 0.50)) if positive_scores.size else 0.0,
        "neighbor_score_p75": float(np.quantile(positive_scores, 0.75)) if positive_scores.size else 0.0,
        "neighbor_similarity_mean": float(np.mean(similarity[reliability > 0])) if positive_scores.size else 0.0,
        "neighbor_shared_score_mean": float(np.mean(shared_score[reliability > 0])) if positive_scores.size else 0.0,
        "cells_with_reliable_neighbor": float(np.mean(first_reliable != np.arange(n_cells))),
    }
    if pseudo_labels is not None:
        same_cluster_edges = pseudo_labels[np.arange(n_cells)[:, None]] == pseudo_labels[neigh]
        stats["neighbor_pseudo_same_cluster_edge_fraction"] = float(np.mean(same_cluster_edges[reliability > 0])) if positive_scores.size else 0.0
        stats["neighbor_pseudo_same_cluster_reliable_fraction"] = float(np.mean(same_cluster_edges[eligible])) if eligible_scores.size else 0.0
    else:
        stats["neighbor_pseudo_same_cluster_edge_fraction"] = 0.0
        stats["neighbor_pseudo_same_cluster_reliable_fraction"] = 0.0
    return NeighborState(
        indices=neigh,
        reliability=reliability,
        similarity=similarity,
        shared_score=shared_score,
        eligible=eligible,
        first_reliable=first_reliable,
        mean_reliable_count=stats["neighbor_reliable_count_mean"],
        stats=stats,
    )


def consensus_neighbor_state(states: list[NeighborState], min_hits: int = 2) -> NeighborState:
    if not states:
        raise ValueError("At least one NeighborState is required for consensus NeighborMix.")
    current = states[-1]
    window = len(states)
    min_hits = max(1, min(int(min_hits), window))
    if window == 1 or min_hits <= 1:
        current.stats["neighbor_consensus_window"] = int(window)
        current.stats["neighbor_consensus_min_hits"] = int(min_hits)
        current.stats["neighbor_consensus_hit_mean"] = 1.0
        return current

    n_cells, n_neighbors = current.indices.shape
    hit_count = np.zeros((n_cells, n_neighbors), dtype=np.int16)
    reliability_sum = np.zeros((n_cells, n_neighbors), dtype=np.float32)

    for state in states:
        for i in range(n_cells):
            keep = state.eligible[i]
            if not np.any(keep):
                continue
            old_ids = state.indices[i, keep]
            old_rel = state.reliability[i, keep]
            for pos, neighbor_id in enumerate(current.indices[i]):
                match = np.where(old_ids == neighbor_id)[0]
                if match.size:
                    hit_count[i, pos] += 1
                    reliability_sum[i, pos] += float(old_rel[int(match[0])])

    eligible = current.eligible & (hit_count >= min_hits)
    reliability = current.reliability.copy()
    reliability[eligible] = reliability_sum[eligible] / np.maximum(hit_count[eligible], 1)
    first_reliable = np.arange(n_cells, dtype=np.int64)
    for i in range(n_cells):
        positions = np.flatnonzero(eligible[i])
        if positions.size:
            first_reliable[i] = int(current.indices[i, int(positions[0])])

    eligible_counts = eligible.sum(axis=1).astype(np.float32) if eligible.size else np.zeros(n_cells, dtype=np.float32)
    eligible_scores = reliability[eligible]
    hit_scores = hit_count[eligible].astype(np.float32) if eligible_scores.size else np.array([], dtype=np.float32)
    stats = dict(current.stats)
    stats.update(
        {
            "neighbor_consensus_window": int(window),
            "neighbor_consensus_min_hits": int(min_hits),
            "neighbor_consensus_hit_mean": float(np.mean(hit_scores)) if hit_scores.size else 0.0,
            "neighbor_reliability_mean": float(np.mean(eligible_scores)) if eligible_scores.size else 0.0,
            "neighbor_reliability_min": float(np.min(eligible_scores)) if eligible_scores.size else 0.0,
            "neighbor_reliable_edge_fraction": float(np.mean(eligible)) if eligible.size else 0.0,
            "neighbor_reliable_count_mean": float(np.mean(eligible_counts)) if eligible_counts.size else 0.0,
            "neighbor_reliable_count_min": float(np.min(eligible_counts)) if eligible_counts.size else 0.0,
            "cells_with_reliable_neighbor": float(np.mean(first_reliable != np.arange(n_cells))),
        }
    )
    return NeighborState(
        indices=current.indices,
        reliability=reliability,
        similarity=current.similarity,
        shared_score=current.shared_score,
        eligible=eligible,
        first_reliable=first_reliable,
        mean_reliable_count=stats["neighbor_reliable_count_mean"],
        stats=stats,
    )


def adaptive_consensus_neighbor_state(
    states: list[NeighborState],
    loose_hits: int = 2,
    strict_hits: int | None = None,
    score_threshold: float = 0.84,
) -> NeighborState:
    if not states:
        raise ValueError("At least one NeighborState is required for adaptive consensus NeighborMix.")
    current = states[-1]
    window = len(states)
    strict_hits = window if strict_hits is None else int(strict_hits)
    loose_hits = max(1, min(int(loose_hits), window))
    strict_hits = max(loose_hits, min(strict_hits, window))
    score_threshold = float(score_threshold)
    if window == 1:
        current.stats["neighbor_adaptive_consensus_window"] = int(window)
        current.stats["neighbor_adaptive_loose_hits"] = int(loose_hits)
        current.stats["neighbor_adaptive_strict_hits"] = int(strict_hits)
        current.stats["neighbor_adaptive_score_threshold"] = float(score_threshold)
        current.stats["neighbor_adaptive_score_mean"] = (
            float(np.mean(current.reliability[current.eligible])) if np.any(current.eligible) else 0.0
        )
        current.stats["neighbor_adaptive_hit_mean"] = 1.0
        current.stats["neighbor_adaptive_core_edge_fraction"] = (
            float(np.mean(current.eligible)) if current.eligible.size else 0.0
        )
        current.stats["neighbor_adaptive_strict_edge_fraction"] = 0.0
        return current

    n_cells, n_neighbors = current.indices.shape
    hit_count = np.zeros((n_cells, n_neighbors), dtype=np.int16)
    reliability_sum = np.zeros((n_cells, n_neighbors), dtype=np.float32)

    for state in states:
        for i in range(n_cells):
            keep = state.eligible[i]
            if not np.any(keep):
                continue
            old_ids = state.indices[i, keep]
            old_rel = state.reliability[i, keep]
            for pos, neighbor_id in enumerate(current.indices[i]):
                match = np.where(old_ids == neighbor_id)[0]
                if match.size:
                    hit_count[i, pos] += 1
                    reliability_sum[i, pos] += float(old_rel[int(match[0])])

    reliability_avg = current.reliability.copy()
    hit_mask = hit_count > 0
    reliability_avg[hit_mask] = reliability_sum[hit_mask] / np.maximum(hit_count[hit_mask], 1)
    strict_edges = hit_count >= strict_hits
    core_edges = (hit_count >= loose_hits) & (reliability_avg >= score_threshold)
    eligible = current.eligible & (strict_edges | core_edges)

    hit_fraction = hit_count.astype(np.float32) / float(max(1, window))
    reliability = reliability_avg * np.clip(hit_fraction / max(loose_hits / float(window), 1e-6), 0.0, 1.0)
    reliability[~eligible] = current.reliability[~eligible]

    first_reliable = np.arange(n_cells, dtype=np.int64)
    for i in range(n_cells):
        positions = np.flatnonzero(eligible[i])
        if positions.size:
            first_reliable[i] = int(current.indices[i, int(positions[0])])

    eligible_counts = eligible.sum(axis=1).astype(np.float32) if eligible.size else np.zeros(n_cells, dtype=np.float32)
    eligible_scores = reliability[eligible]
    eligible_hits = hit_count[eligible].astype(np.float32) if eligible_scores.size else np.array([], dtype=np.float32)
    adaptive_scores = reliability_avg[eligible]
    stats = dict(current.stats)
    stats.update(
        {
            "neighbor_adaptive_consensus_window": int(window),
            "neighbor_adaptive_loose_hits": int(loose_hits),
            "neighbor_adaptive_strict_hits": int(strict_hits),
            "neighbor_adaptive_score_threshold": float(score_threshold),
            "neighbor_adaptive_score_mean": float(np.mean(adaptive_scores)) if adaptive_scores.size else 0.0,
            "neighbor_adaptive_hit_mean": float(np.mean(eligible_hits)) if eligible_hits.size else 0.0,
            "neighbor_adaptive_core_edge_fraction": float(np.mean(current.eligible & core_edges)) if eligible.size else 0.0,
            "neighbor_adaptive_strict_edge_fraction": float(np.mean(current.eligible & strict_edges)) if eligible.size else 0.0,
            "neighbor_reliability_mean": float(np.mean(eligible_scores)) if eligible_scores.size else 0.0,
            "neighbor_reliability_min": float(np.min(eligible_scores)) if eligible_scores.size else 0.0,
            "neighbor_reliable_edge_fraction": float(np.mean(eligible)) if eligible.size else 0.0,
            "neighbor_reliable_count_mean": float(np.mean(eligible_counts)) if eligible_counts.size else 0.0,
            "neighbor_reliable_count_min": float(np.min(eligible_counts)) if eligible_counts.size else 0.0,
            "cells_with_reliable_neighbor": float(np.mean(first_reliable != np.arange(n_cells))),
        }
    )
    return NeighborState(
        indices=current.indices,
        reliability=reliability,
        similarity=current.similarity,
        shared_score=current.shared_score,
        eligible=eligible,
        first_reliable=first_reliable,
        mean_reliable_count=stats["neighbor_reliable_count_mean"],
        stats=stats,
    )


def mix_batch(
    batch_indices: torch.Tensor,
    batch_x: torch.Tensor,
    full_data_cpu: torch.Tensor,
    state: NeighborState | None,
    alpha: float,
    mode: str = "first",
    soft_power: float = 1.0,
) -> tuple[torch.Tensor, float]:
    if state is None:
        return batch_x, 0.0
    idx_np = batch_indices.detach().cpu().numpy().astype(np.int64)
    if mode == "first":
        neighbor_idx = state.first_reliable[idx_np]
        has_neighbor = neighbor_idx != idx_np
        neighbor_x = full_data_cpu[neighbor_idx].to(batch_x.device)
        mixed = float(alpha) * batch_x + (1.0 - float(alpha)) * neighbor_x
        mixed = torch.where(torch.as_tensor(has_neighbor, device=batch_x.device).view(-1, 1), mixed, batch_x)
        return mixed, float(np.mean(has_neighbor)) if has_neighbor.size else 0.0
    if mode == "soft_first":
        neighbor_rows = state.indices[idx_np]
        reliability_rows = state.reliability[idx_np]
        eligible_rows = state.eligible[idx_np]
        has_neighbor = eligible_rows.any(axis=1)
        first_pos = np.argmax(eligible_rows, axis=1)
        row_ids = np.arange(len(idx_np))
        neighbor_idx = neighbor_rows[row_ids, first_pos]
        neighbor_idx = np.where(has_neighbor, neighbor_idx, idx_np)
        edge_weight = reliability_rows[row_ids, first_pos]
        edge_weight = np.where(has_neighbor, edge_weight, 0.0).astype(np.float32)
        edge_weight = np.power(np.clip(edge_weight, 0.0, 1.0), max(float(soft_power), 1e-6))
        neighbor_x = full_data_cpu[neighbor_idx].to(batch_x.device)
        beta = torch.as_tensor((1.0 - float(alpha)) * edge_weight, device=batch_x.device, dtype=batch_x.dtype).view(-1, 1)
        mixed = (1.0 - beta) * batch_x + beta * neighbor_x
        return mixed, float(np.mean(has_neighbor)) if has_neighbor.size else 0.0
    if mode not in {"mean", "weighted_mean"}:
        raise ValueError(f"Unknown NeighborMix mode: {mode}")
    neighbor_rows = state.indices[idx_np]
    reliability_rows = state.reliability[idx_np]
    eligible_rows = state.eligible[idx_np]
    mixed_neighbors = []
    has_neighbor = []
    for row, rel, elig in zip(neighbor_rows, reliability_rows, eligible_rows):
        keep = elig
        if np.any(keep):
            selected = full_data_cpu[row[keep]]
            if mode == "weighted_mean":
                weights = torch.as_tensor(rel[keep], dtype=selected.dtype).clamp_min(1e-6)
                weights = weights / weights.sum().clamp_min(1e-6)
                mixed_neighbors.append(torch.sum(selected * weights.view(-1, 1), dim=0))
            else:
                mixed_neighbors.append(selected.mean(dim=0))
            has_neighbor.append(True)
        else:
            mixed_neighbors.append(full_data_cpu[int(row[0])] * 0.0)
            has_neighbor.append(False)
    neighbor_x = torch.stack(mixed_neighbors, dim=0).to(batch_x.device)
    has_neighbor_t = torch.as_tensor(has_neighbor, device=batch_x.device, dtype=torch.bool).view(-1, 1)
    mixed = float(alpha) * batch_x + (1.0 - float(alpha)) * neighbor_x
    mixed = torch.where(has_neighbor_t, mixed, batch_x)
    return mixed, float(np.mean(has_neighbor)) if has_neighbor else 0.0
