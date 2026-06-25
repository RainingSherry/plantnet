from __future__ import annotations

import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import scanpy as sc
import scipy.sparse as sp
from sklearn.preprocessing import LabelEncoder

from .dataset import EvaluationBundle


LABEL_CANDIDATES = [
    "resolved_label",
    "maintype",
    "cell_type",
    "Celltype",
    "celltype",
    "label",
    "labels",
    "cell_label",
    "Cluster",
    "cluster",
    "clusters",
    "Seurat_clusters",
    "true_label",
]


@dataclass
class CAAMDataBundle:
    x: np.ndarray
    gene_names: np.ndarray
    selected_gene_indices: np.ndarray
    batch_code: np.ndarray
    library_size: np.ndarray
    zero_ratio: np.ndarray
    evaluation: EvaluationBundle
    profile: dict[str, Any]
    preprocess_config: dict[str, Any]
    adata: sc.AnnData


def _dense_float32(matrix) -> np.ndarray:
    arr = matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)
    arr = np.asarray(arr, dtype=np.float32)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    arr[arr < 0.0] = 0.0
    return np.ascontiguousarray(arr, dtype=np.float32)


def _looks_like_raw_counts(matrix) -> bool:
    sample = matrix[: min(matrix.shape[0], 256)]
    values = sample.data if sp.issparse(sample) else np.asarray(sample).ravel()
    values = np.asarray(values, dtype=np.float32)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return False
    values = values[: min(values.size, 100000)]
    return bool(np.all(values >= 0.0) and np.allclose(values, np.round(values), atol=1e-4))


def _select_matrix(adata: sc.AnnData, input_mode: str):
    if input_mode in {"auto", "raw"} and "counts" in adata.layers:
        return adata.layers["counts"], "layers[counts]", "raw"
    if input_mode in {"auto", "raw"} and adata.raw is not None:
        return adata.raw.X, "adata.raw.X", "raw"
    inferred = "raw" if _looks_like_raw_counts(adata.X) else "log1p"
    if input_mode == "raw" and inferred != "raw":
        raise ValueError("--input_mode raw requested but no raw-looking source is available.")
    if input_mode == "log1p":
        inferred = "log1p"
    return adata.X, "adata.X", inferred


def _resolve_labels(adata: sc.AnnData) -> EvaluationBundle:
    for key in LABEL_CANDIDATES:
        if key in adata.obs.columns:
            raw = adata.obs[key].astype(str).to_numpy()
            enc = LabelEncoder()
            return EvaluationBundle(
                labels=enc.fit_transform(raw).astype(np.int64),
                label_names=enc.classes_.astype(str),
                label_key=key,
            )
    labels = np.zeros(adata.n_obs, dtype=np.int64)
    return EvaluationBundle(labels=labels, label_names=np.asarray(["unknown"]), label_key=None)


def _batch_code(adata: sc.AnnData) -> np.ndarray:
    for key in ("batch", "batch_id", "Batch", "sample", "sample_id", "donor"):
        if key in adata.obs.columns:
            raw = adata.obs[key].astype(str).to_numpy()
            return LabelEncoder().fit_transform(raw).astype(np.int64)
    return np.zeros(adata.n_obs, dtype=np.int64)


def _var_names(adata: sc.AnnData) -> np.ndarray:
    names = np.asarray(adata.var_names).astype(str)
    for col in ("gene_name", "features", "gene_symbols", "symbol", "_index"):
        if col in adata.var.columns:
            candidate = adata.var[col].astype(str).to_numpy()
            if len(candidate) == len(names) and len(np.unique(candidate)) == len(candidate):
                return candidate
    return names


def load_caam_data(
    data_path: str,
    *,
    input_mode: str,
    target_sum: float,
    n_top_genes: int,
    scale_input: bool,
    benchmark_mode: bool,
    seed: int,
) -> CAAMDataBundle:
    adata = sc.read_h5ad(data_path)
    work = adata.copy()
    original_index_key = "_caam_original_gene_index"
    work.var[original_index_key] = np.arange(work.n_vars, dtype=np.int64)
    source, source_name, inferred_mode = _select_matrix(work, input_mode)
    work.X = source.copy() if sp.issparse(source) else np.asarray(source).copy()

    if benchmark_mode:
        inferred_mode = "log1p"
        scale_input = False
    else:
        if inferred_mode == "raw" or (_looks_like_raw_counts(work.X) and input_mode in {"auto", "raw"}):
            sc.pp.normalize_total(work, target_sum=target_sum)
            sc.pp.log1p(work)
            inferred_mode = "raw->log1p"

    if n_top_genes and n_top_genes > 0:
        if work.n_vars > n_top_genes:
            sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=n_top_genes, subset=False)
            score_key = "dispersions_norm" if "dispersions_norm" in work.var.columns else "dispersions"
            scores = work.var[score_key].to_numpy(dtype=np.float64)
            scores = np.nan_to_num(scores, nan=-np.inf, posinf=np.inf, neginf=-np.inf)
            selected = np.argsort(-scores, kind="mergesort")[: int(n_top_genes)]
            keep = np.zeros(work.n_vars, dtype=bool)
            keep[selected] = True
            work = work[:, keep].copy()
            feature_space_source = "hvg"
        else:
            feature_space_source = "all_genes_below_hvg_target"
    else:
        feature_space_source = "full_gene_stress" if benchmark_mode else "all_genes"

    if scale_input:
        warnings.warn(
            "WARNING: scale_input=True changes the semantic meaning of zero values.",
            RuntimeWarning,
        )
        sc.pp.scale(work)

    x = _dense_float32(work.X)
    gene_names = _var_names(work)
    selected_gene_indices = work.var[original_index_key].to_numpy(dtype=np.int64)
    batch_code = _batch_code(work)
    library_size = x.sum(axis=1).astype(np.float32)
    zero_ratio = (x <= 0.0).mean(axis=1).astype(np.float32)
    evaluation = _resolve_labels(work)
    label_counts: dict[str, int] = {}
    if evaluation.label_key is not None:
        counts = work.obs[evaluation.label_key].astype(str).value_counts(dropna=False).sort_index()
        label_counts = {str(k): int(v) for k, v in counts.items()}
    profile = {
        "dataset_name": Path(data_path).stem,
        "n_cells": int(work.n_obs),
        "n_genes": int(work.n_vars),
        "n_genes_original": int(adata.n_vars),
        "source": source_name,
        "input_mode": inferred_mode,
        "benchmark_mode": bool(benchmark_mode),
        "scale_input": bool(scale_input),
        "feature_space_source": feature_space_source,
        "actual_n_genes_after_selection": int(work.n_vars),
        "label_key": evaluation.label_key,
        "n_cell_types": int(len(evaluation.label_names)),
        "cell_type_counts": label_counts,
    }
    preprocess_config = {
        "data_path": str(data_path),
        "input_mode": input_mode,
        "resolved_input_mode": inferred_mode,
        "target_sum": float(target_sum),
        "n_top_genes": int(n_top_genes),
        "scale_input": bool(scale_input),
        "benchmark_mode": bool(benchmark_mode),
        "feature_space_source": feature_space_source,
        "actual_n_genes_after_selection": int(work.n_vars),
        "seed": int(seed),
    }
    return CAAMDataBundle(
        x=x,
        gene_names=gene_names,
        selected_gene_indices=selected_gene_indices,
        batch_code=batch_code,
        library_size=library_size,
        zero_ratio=zero_ratio,
        evaluation=evaluation,
        profile=profile,
        preprocess_config=preprocess_config,
        adata=work,
    )
