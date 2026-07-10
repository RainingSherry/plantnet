import os

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


PROJECT_ROOT = Path(__file__).resolve().parents[3]
PAPER_ROOT = Path(__file__).resolve().parents[1]
PYTHON_BIN = Path("/data/luolie/conda/envs/scclubench-sccdcg-h100/bin/python")
MODEL_RUNNER = PROJECT_ROOT / "experimental_retired_models/RG_NeighborMix_scMAE/run.py"

REMOTE_DATA_ROOT = "<SCVICAR_DATA_ROOT>"
REMOTE_RESULT_ROOT = "<SCVICAR_RESULT_ROOT>"
REMOTE_HOST = os.getenv("SCVICAR_REMOTE_HOST", "<SCVICAR_REMOTE_HOST>")
REMOTE_PORT = 16335
REMOTE_USER = os.getenv("SCVICAR_REMOTE_USER", "<SCVICAR_REMOTE_USER>")
PROTOCOL_VERSION = "protocol_v1"

SEEDS = (42, 2024, 3407)
SPLIT_SEEDS = (11, 23, 37, 53, 71)


@dataclass(frozen=True)
class DatasetSpec:
    name: str
    source_path: str
    source_label_key: str
    expected_clusters: int
    exclude_labels: tuple[str, ...] = ("unknown", "nan", "unassigned")
    min_class_size: int = 10


DATASETS: dict[str, DatasetSpec] = {
    "Blood_BoneMarrow": DatasetSpec(
        "Blood_BoneMarrow", "data/其他/Blood_BoneMarrow.h5ad", "cell_type", 30
    ),
    "Human_Pancreas_1": DatasetSpec(
        "Human_Pancreas_1", "data/其他/Human_Pancreas_1.h5ad", "cell_type", 6
    ),
    "Human_Pancreas_3": DatasetSpec(
        "Human_Pancreas_3", "data/其他/Human_Pancreas_3.h5ad", "cell_type", 13
    ),
    "Mouse_Pancreas_1": DatasetSpec(
        "Mouse_Pancreas_1", "data/其他/Mouse_Pancreas_1.h5ad", "cell_type", 10
    ),
    "PRJNA895163": DatasetSpec(
        "PRJNA895163", "data/processed_benchmark/PRJNA895163.h5ad", "resolved_label", 12
    ),
    "TabulaSapiens_Pancreas": DatasetSpec(
        "TabulaSapiens_Pancreas", "data/其他/TabulaSapiens_Pancreas.h5ad", "cell_type", 16
    ),
}


# These defaults are deliberately shared across every variant. They match the
# untuned RG runner defaults and make the fixed variant exactly alpha=0.9.
COMMON_MODEL_ARGS: dict[str, Any] = {
    "label_key": "resolved_label",
    "input_mode": "auto",
    "n_top_genes": 1000,
    "target_sum": 10000.0,
    "scale_input": True,
    "hidden_size": 128,
    "dropout": 0.0,
    "masked_data_weight": 0.75,
    "mask_loss_weight": 0.7,
    "epochs": 80,
    "batch_size": 256,
    "lr": 1e-3,
    "mask_ratio": 0.4,
    "neighbor_k": 5,
    "mix_neighbors": 4,
    "neighbor_estimator": "current",
    "stress_bad_edge_ratio": 0.0,
    "tau": 0.2,
    "knn_pca_dim": 50,
    "pseudo_weight": 0.3,
    "gate_min": 0.0,
    "gate_max": 0.1,
    "gamma_sim": 1.0,
    "gamma_mutual": 1.0,
    "gamma_snn": 1.0,
    "gamma_distance": 1.0,
    "beta_mutual": 1.0,
    "beta_snn": 1.0,
    "beta_perturb": 2.0,
    "beta_uncertainty": 1.0,
    "contrast_weight": 0.0,
    "contrast_projection_dim": 0,
    "contrast_partition_mode": "none",
}


VARIANTS: dict[str, dict[str, Any]] = {
    "nomix": {
        "method_name": "scVICAR-NoMix",
        "variant_name": "mb_scmae_nomix",
        "mix_mode": "none",
        "gate_mode": "none",
        "edge_reliability_mode": "none",
        "neighbor_k": 0,
        "mix_neighbors": 0,
        "pseudo_weight": 0.0,
    },
    "random_mix": {
        "method_name": "scVICAR-Random",
        "variant_name": "mb_random_control",
        "mix_mode": "random",
        "gate_mode": "constant",
        "edge_reliability_mode": "none",
    },
    "fixed": {
        "method_name": "scVICAR-F",
        "variant_name": "mb_neighbormix_fixed",
        "mix_mode": "fixed",
        "gate_mode": "constant",
        "edge_reliability_mode": "none",
    },
    "topology_edge_only": {
        "method_name": "scVICAR-T-Edge",
        "variant_name": "mb_rg_edge_only",
        "mix_mode": "reliability",
        "gate_mode": "constant",
        "edge_reliability_mode": "sim_mutual_snn_distance",
    },
    "topology_gate_only": {
        "method_name": "scVICAR-T-Gate",
        "variant_name": "mb_rg_gate_only",
        "mix_mode": "reliability",
        "gate_mode": "topology",
        "edge_reliability_mode": "none",
    },
    "topology_full": {
        "method_name": "scVICAR-T",
        "variant_name": "mb_rg_full",
        "mix_mode": "reliability",
        "gate_mode": "topology",
        "edge_reliability_mode": "sim_mutual_snn_distance",
    },
}


def canonical_json(payload: Any) -> str:
    return json.dumps(payload, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def sha256_payload(payload: Any, length: int = 16) -> str:
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()[:length]


def resolved_run_config(dataset: str, variant: str, seed: int, epochs: int | None = None) -> dict[str, Any]:
    if dataset not in DATASETS:
        raise KeyError(f"Unknown dataset {dataset!r}; choose from {sorted(DATASETS)}")
    if variant not in VARIANTS:
        raise KeyError(f"Unknown variant {variant!r}; choose from {sorted(VARIANTS)}")
    if seed not in SEEDS:
        raise ValueError(f"Seed {seed} is outside the preregistered set {SEEDS}")
    spec = DATASETS[dataset]
    cfg = dict(COMMON_MODEL_ARGS)
    cfg.update(VARIANTS[variant])
    cfg.update(
        {
            "protocol_version": PROTOCOL_VERSION,
            "dataset_name": dataset,
            "n_clusters": spec.expected_clusters,
            "seed": int(seed),
        }
    )
    if epochs is not None:
        cfg["epochs"] = int(epochs)
        cfg["execution_mode"] = "smoke" if int(epochs) < COMMON_MODEL_ARGS["epochs"] else "formal"
    else:
        cfg["execution_mode"] = "formal"
    return cfg


def protocol_snapshot() -> dict[str, Any]:
    return {
        "protocol_version": PROTOCOL_VERSION,
        "seeds": list(SEEDS),
        "split_seeds": list(SPLIT_SEEDS),
        "datasets": {key: asdict(value) for key, value in DATASETS.items()},
        "common_model_args": COMMON_MODEL_ARGS,
        "variants": VARIANTS,
    }
