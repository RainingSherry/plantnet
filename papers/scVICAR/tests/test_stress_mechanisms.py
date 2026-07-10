from __future__ import annotations

import numpy as np
import pytest
import torch

from experimental_retired_models.RG_NeighborMix_scMAE.mixing import make_pseudo_batch
from experimental_retired_models.RG_NeighborMix_scMAE.neighbor_graph import NeighborGraph, build_pca_knn_graph, inject_cross_label_edges


GRAPH_ARRAY_FIELDS = ("indices", "probs", "similarity", "distance", "embedding", "mutual", "snn")


def test_zero_bad_edge_ratio_is_numerically_identical() -> None:
    data = np.random.default_rng(17).normal(size=(10, 6)).astype(np.float32)
    labels = np.repeat(np.arange(2), 5)
    clean = build_pca_knn_graph(data, k=3, pca_dim=4, tau=0.2, seed=42)
    explicit_zero = build_pca_knn_graph(
        data,
        k=3,
        pca_dim=4,
        tau=0.2,
        seed=42,
        labels=labels,
        stress_bad_edge_ratio=0.0,
    )

    for field in GRAPH_ARRAY_FIELDS:
        np.testing.assert_array_equal(getattr(explicit_zero, field), getattr(clean, field))
    assert explicit_zero.profile == clean.profile
    assert explicit_zero.profile["label_leakage_diagnostic"] is False
    assert inject_cross_label_edges(clean, labels, ratio=0.0, tau=0.2, seed=42) is clean


def test_full_bad_edge_stress_makes_every_edge_cross_label_and_recomputes_statistics() -> None:
    data = np.random.default_rng(23).normal(size=(8, 5)).astype(np.float32)
    labels = np.asarray([0, 0, 0, 0, 1, 1, 1, 1])
    graph = build_pca_knn_graph(
        data,
        k=3,
        pca_dim=4,
        tau=0.3,
        seed=7,
        labels=labels,
        stress_bad_edge_ratio=1.0,
    )

    assert np.all(labels[graph.indices] != labels[:, None])
    expected_similarity = np.einsum("ij,ikj->ik", graph.embedding, graph.embedding[graph.indices])
    np.testing.assert_allclose(graph.similarity, expected_similarity, rtol=0.0, atol=1e-7)
    np.testing.assert_allclose(graph.distance, 1.0 - expected_similarity, rtol=0.0, atol=1e-7)
    scaled = expected_similarity / 0.3
    scaled -= scaled.max(axis=1, keepdims=True)
    expected_probs = np.exp(scaled) / np.exp(scaled).sum(axis=1, keepdims=True)
    np.testing.assert_allclose(graph.probs, expected_probs, rtol=1e-6, atol=1e-7)

    neighbor_sets = [set(row.tolist()) for row in graph.indices]
    for cell, row in enumerate(graph.indices):
        for position, neighbor in enumerate(row):
            assert graph.mutual[cell, position] == (cell in neighbor_sets[neighbor])
            union = neighbor_sets[cell].union(neighbor_sets[neighbor])
            expected_snn = len(neighbor_sets[cell].intersection(neighbor_sets[neighbor])) / max(1, len(union))
            assert graph.snn[cell, position] == pytest.approx(expected_snn)
    assert graph.profile["stress_bad_edge_ratio"] == 1.0
    assert graph.profile["stress_bad_edge_ratio_realized"] == 1.0
    assert graph.profile["stress_cross_label_edge_fraction"] == 1.0
    assert graph.profile["label_leakage_diagnostic"] is True


@pytest.mark.parametrize("ratio", [0.25, 0.5, 0.75])
def test_requested_bad_edge_fraction_is_realized_globally(ratio: float) -> None:
    data = np.random.default_rng(31).normal(size=(40, 8)).astype(np.float32)
    labels = np.repeat(np.arange(4), 10)
    clean = build_pca_knn_graph(data, k=5, pca_dim=6, tau=0.2, seed=13)
    stressed = inject_cross_label_edges(clean, labels, ratio=ratio, tau=0.2, seed=13)
    changed = stressed.indices != clean.indices
    assert changed.mean() == pytest.approx(ratio, abs=1.0 / changed.size)
    assert stressed.profile["stress_bad_edge_ratio_realized"] == pytest.approx(ratio, abs=1.0 / changed.size)
    assert np.all(labels[stressed.indices[changed]] != np.repeat(labels[:, None], 5, axis=1)[changed])


def estimator_fixture() -> tuple[np.ndarray, NeighborGraph]:
    data = np.asarray([[0.0], [1.0], [3.0]], dtype=np.float32)
    indices = np.asarray([[1, 2], [0, 2], [0, 1]], dtype=np.int64)
    probs = np.asarray([[0.25, 0.75], [0.5, 0.5], [0.5, 0.5]], dtype=np.float32)
    graph = NeighborGraph(
        indices=indices,
        probs=probs,
        similarity=np.zeros((3, 2), dtype=np.float32),
        distance=np.ones((3, 2), dtype=np.float32),
        embedding=np.zeros((3, 1), dtype=np.float32),
        mutual=np.ones((3, 2), dtype=bool),
        snn=np.zeros((3, 2), dtype=np.float32),
        profile={},
    )
    return data, graph


@pytest.mark.parametrize(
    ("estimator", "expected_mean"),
    [
        ("current", 2.6875),
        ("uniform_sample", 2.5),
        ("full", 2.5),
    ],
)
def test_neighbor_estimators_have_expected_mean_and_output_shape(estimator: str, expected_mean: float) -> None:
    data, graph = estimator_fixture()
    repeats = 10_000
    batch_indices = np.zeros(repeats, dtype=np.int64)
    batch_x = torch.zeros((repeats, 1), dtype=torch.float32)
    mixed, sample_weight, _ = make_pseudo_batch(
        data_np=data,
        batch_indices=batch_indices,
        batch_x=batch_x,
        mix_mode="fixed",
        graph=graph,
        edge_weights=graph.probs,
        node_gate=np.ones(3, dtype=np.float32),
        mix_neighbors=2,
        rng=np.random.default_rng(101),
        neighbor_estimator=estimator,
    )

    assert mixed.shape == batch_x.shape
    assert sample_weight.shape == (repeats,)
    assert float(mixed.mean()) == pytest.approx(expected_mean, abs=0.04)
    np.testing.assert_array_equal(sample_weight.numpy(), np.ones(repeats, dtype=np.float32))


def test_current_estimator_matches_legacy_weighting_exactly() -> None:
    data, graph = estimator_fixture()
    seed = 919
    reference_rng = np.random.default_rng(seed)
    choices = reference_rng.choice(2, size=2, replace=True, p=graph.probs[0])
    picked = graph.probs[0, choices]
    expected = np.sum(data[graph.indices[0, choices]] * (picked / picked.sum())[:, None], axis=0)

    mixed, _, _ = make_pseudo_batch(
        data_np=data,
        batch_indices=np.asarray([0]),
        batch_x=torch.zeros((1, 1), dtype=torch.float32),
        mix_mode="fixed",
        graph=graph,
        edge_weights=graph.probs,
        node_gate=np.ones(3, dtype=np.float32),
        mix_neighbors=2,
        rng=np.random.default_rng(seed),
        neighbor_estimator="current",
    )

    np.testing.assert_array_equal(mixed.numpy()[0], expected.astype(np.float32))
