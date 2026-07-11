from __future__ import annotations

import numpy as np
import pandas as pd

from methods.DeepLearning import scMAE_family as family
from papers.scVICAR.code.config import DATASETS, SEEDS, VARIANTS
from papers.scVICAR.code.generate_manuscript_assets import write_full_label_sensitivity
from papers.scVICAR.code.secondary_aggregate import (
    CONFIRMATORY_METRICS,
    contrast_table,
    summarize,
    verify_arrays_and_metrics,
)


def test_secondary_arrays_are_independently_recomputed(tmp_path) -> None:
    labels = np.asarray([0, 0, 1, 1, 2, 2], dtype=np.int64)
    predicted = np.asarray([1, 1, 0, 0, 2, 2], dtype=np.int32)
    metrics, mapped = family.compute_kmeans_metrics(labels, predicted)
    metrics.update(
        cluster_method="leiden_fixed_resolution_1p0",
        uses_known_k=False,
        resolution=1.0,
        n_neighbors=15,
    )
    np.savez_compressed(tmp_path / "leiden_clusters.npz", labels=labels, predicted=predicted, mapped=mapped)
    recomputed, n_pred = verify_arrays_and_metrics(tmp_path, metrics, 6, 3, "test-run")
    assert n_pred == 3
    assert recomputed["ari"] == 1.0


def test_secondary_summaries_use_dataset_as_independent_unit() -> None:
    rows = []
    for dataset_index, dataset in enumerate(DATASETS):
        for variant_index, variant in enumerate(VARIANTS):
            for seed in SEEDS:
                base = 0.2 + 0.01 * dataset_index + 0.001 * variant_index + seed * 1e-8
                rows.append({
                    "dataset": dataset,
                    "variant": variant,
                    "model_seed": seed,
                    **{metric: base for metric in (
                        "ari", "nmi", "acc", "f1_macro", "fmi", "v_measure", "homogeneity", "completeness"
                    )},
                })
    run_level = pd.DataFrame(rows)
    dataset_level, overall = summarize(run_level)
    contrasts = contrast_table(run_level, dataset_level)
    assert len(dataset_level) == 36
    assert set(dataset_level["n_seeds"]) == {3}
    assert len(overall) == 6
    assert set(overall["n_datasets"]) == {6}
    assert len(contrasts) == len(CONFIRMATORY_METRICS) * 3
    assert set(contrasts["n_datasets"]) == {6}
    assert (contrasts[["wins", "ties", "losses"]].sum(axis=1) == 6).all()


def test_full_label_sensitivity_asset_requires_complete_three_variant_summary(tmp_path) -> None:
    variants = ["nomix", "fixed", "topology_full"]
    overall = pd.DataFrame([
        {
            "variant": variant, "n_datasets": 6,
            "ari_mean": 0.5 + index * 0.01, "ari_sd": 0.1,
            "nmi_mean": 0.6 + index * 0.01, "nmi_sd": 0.1,
            "f1_macro_mean": 0.4 + index * 0.01, "f1_macro_sd": 0.1,
        }
        for index, variant in enumerate(variants)
    ])
    contrasts = pd.DataFrame([
        {
            "metric": "ari", "contrast": contrast, "mean_delta": 0.01,
            "ci95_low": -0.01, "ci95_high": 0.03,
            "wins": 4, "ties": 0, "losses": 2, "holm_p": 0.5,
        }
        for contrast in ("F_vs_NoMix", "T_vs_NoMix", "T_vs_F")
    ])
    manuscript = tmp_path / "manuscript"; tables = tmp_path / "tables"
    manuscript.mkdir(); (manuscript / "generated").mkdir(); tables.mkdir()
    write_full_label_sensitivity(overall, contrasts, manuscript, tables)
    table = (tables / "full_label_sensitivity.tex").read_text()
    prose = (manuscript / "generated/full_label_sensitivity_results.tex").read_text()
    assert table.count("scVICAR-") == 2
    assert "wins/ties/losses 4/0/2" in prose
