from __future__ import annotations

import json

import anndata as ad
import numpy as np
import pytest

from papers.scVICAR.code.run_baseline import BASELINES, canonicalize_baseline


def test_baseline_matrix_has_all_dataset_method_seed_keys():
    matrix = json.loads(
        open("papers/scVICAR/experiments/baselines_v1/planned_matrix.json", encoding="utf-8").read()
    )
    assert len(matrix) == 108
    assert len({(row["dataset"], row["method"], row["seed"]) for row in matrix}) == 108
    assert {row["method"] for row in matrix} == set(BASELINES)


def test_canonicalize_uses_frozen_labels_and_checks_order(tmp_path):
    data_path = tmp_path / "canonical.h5ad"
    canonical = ad.AnnData(
        X=np.ones((10, 3), dtype=np.float32),
        obs={"resolved_label": ["b"] * 5 + ["a"] * 5},
    )
    canonical.obs_names = [f"cell_{index}" for index in range(10)]
    canonical.write_h5ad(data_path)
    raw = tmp_path / "raw"
    output = tmp_path / "out"
    raw.mkdir(); output.mkdir()
    np.save(raw / "embedding_final.npy", np.arange(20, dtype=np.float32).reshape(10, 2))
    np.save(raw / "labels.npy", np.array([7] * 5 + [3] * 5))
    np.save(raw / "cell_ids.npy", np.asarray(canonical.obs_names, dtype=str))
    result = canonicalize_baseline(raw, output, data_path, "Mouse_Pancreas_1", "pca_kmeans", 42)
    assert result["cell_set_complete"] is True
    saved = np.load(output / "clusters.npz")
    assert saved["labels"].tolist() == [1] * 5 + [0] * 5

    np.save(raw / "labels.npy", np.array([7, 3] * 5))
    with pytest.raises(ValueError, match="cell order"):
        canonicalize_baseline(raw, output, data_path, "Mouse_Pancreas_1", "pca_kmeans", 42)

    np.save(raw / "labels.npy", np.array([7] * 5 + [3] * 5))
    np.save(raw / "cell_ids.npy", np.asarray(list(reversed(canonical.obs_names)), dtype=str))
    with pytest.raises(ValueError, match="cell IDs"):
        canonicalize_baseline(raw, output, data_path, "Mouse_Pancreas_1", "pca_kmeans", 42)
