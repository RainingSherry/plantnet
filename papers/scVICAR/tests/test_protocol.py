from __future__ import annotations

import json

import numpy as np

from papers.scVICAR.code.config import COMMON_MODEL_ARGS, DATASETS, SEEDS, VARIANTS, resolved_run_config
from papers.scVICAR.code.downstream import annotate_clusters, marker_sets, overlap_matrix
from papers.scVICAR.code.run_full_label_sensitivity import SENSITIVITY_VARIANTS


def test_matrix_has_108_unique_tasks() -> None:
    tasks = {(dataset, variant, seed) for dataset in DATASETS for variant in VARIANTS for seed in SEEDS}
    assert len(tasks) == 108


def test_full_label_sensitivity_is_a_separate_54_run_main_model_matrix() -> None:
    assert SENSITIVITY_VARIANTS == ("nomix", "fixed", "topology_full")
    tasks = {
        (dataset, variant, seed)
        for dataset in DATASETS for variant in SENSITIVITY_VARIANTS for seed in SEEDS
    }
    assert len(tasks) == 54


def test_matched_backbone_is_identical() -> None:
    allowed_variant_fields = {
        "method_name", "variant_name", "mix_mode", "gate_mode", "edge_reliability_mode",
        "neighbor_k", "mix_neighbors", "pseudo_weight",
    }
    configs = [resolved_run_config("Mouse_Pancreas_1", variant, 42) for variant in VARIANTS]
    for key, value in COMMON_MODEL_ARGS.items():
        if key in allowed_variant_fields:
            continue
        assert all(config[key] == value for config in configs)
    assert configs[2]["gate_max"] == 0.1


def test_labels_are_not_variant_configuration() -> None:
    serialized = json.dumps(VARIANTS).lower()
    assert "label" not in serialized
    assert "cell_type" not in serialized


def test_marker_overlap_and_deterministic_annotation() -> None:
    import pandas as pd

    reference = {
        "alpha": pd.DataFrame({"names": ["GCG", "TTR", "LOXL4"], "logfoldchanges": [3.0, 2.0, 1.0]}),
        "beta": pd.DataFrame({"names": ["INS", "IAPP", "PCSK1"], "logfoldchanges": [4.0, 3.0, 2.0]}),
    }
    predicted = {
        "0": pd.DataFrame({"names": ["INS", "IAPP", "PCSK1"], "logfoldchanges": [5.0, 4.0, 3.0]}),
        "1": pd.DataFrame({"names": ["GCG", "TTR", "OTHER"], "logfoldchanges": [5.0, 4.0, 0.1]}),
    }
    matrix = overlap_matrix(reference, predicted)
    assert matrix.loc["beta", "0"] == 3
    assert annotate_clusters(reference, predicted) == {"0": "beta", "1": "alpha"}
