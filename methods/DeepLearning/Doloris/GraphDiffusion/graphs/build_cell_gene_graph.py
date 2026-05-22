from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import numpy as np
import scanpy as sc
import scipy.sparse as sp
from sklearn.feature_extraction.text import TfidfTransformer

from .build_gene_graph import GeneGraphBundle, load_adata, select_hvgs


@dataclass
class CellSupportBundle:
    x_dense: np.ndarray
    support_indices: List[np.ndarray]
    support_weights: List[np.ndarray]
    labels: np.ndarray
    label_names: List[str]
    gene_names: List[str]
    size_factors: np.ndarray
    cell_names: List[str]


def _read_obs_dataset(ds) -> np.ndarray:
    if hasattr(ds, "keys"):
        if "categories" in ds.keys() and "codes" in ds.keys():
            categories = _read_obs_dataset(ds["categories"])
            codes = np.asarray(ds["codes"])
            return np.array([categories[int(code)] for code in codes])
        return {key: _read_obs_dataset(ds[key]) for key in ds.keys()}
    if hasattr(ds, "dtype") and ds.dtype.kind in {"S", "O", "U"}:
        arr = []
        for i in range(ds.len()):
            val = ds[i]
            if isinstance(val, bytes):
                val = val.decode("utf-8")
            elif isinstance(val, np.bytes_):
                val = val.decode("utf-8") if val else ""
            arr.append(val)
        return np.array(arr)
    return np.array(ds)


def _dense_array(x):
    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def _infer_label_col(adata: sc.AnnData, label_col: Optional[str]) -> str:
    if label_col and label_col in adata.obs.columns:
        return label_col
    for candidate in ["cell_type", "Celltype", "celltype", "celltype_after", "cell_label", "label", "Seurat_clusters"]:
        if candidate in adata.obs.columns:
            return candidate
    raise KeyError(f"No label column found. Available columns: {list(adata.obs.columns)}")


def _encode_labels(raw_labels: np.ndarray) -> Tuple[np.ndarray, List[str]]:
    unique = list(np.unique(raw_labels))
    mapping = {label: idx for idx, label in enumerate(unique)}
    encoded = np.array([mapping[x] for x in raw_labels], dtype=np.int64)
    return encoded, [str(x) for x in unique]


def compute_support_weights(
    matrix: np.ndarray,
    mode: str = "log1p_count",
) -> np.ndarray:
    matrix = np.asarray(matrix, dtype=np.float32)
    if mode == "log1p_count":
        weights = np.log1p(np.clip(matrix, a_min=0.0, a_max=None))
    elif mode == "normalized_count":
        denom = matrix.sum(axis=1, keepdims=True)
        denom[denom == 0] = 1.0
        weights = matrix / denom
    elif mode == "rank_weight":
        weights = np.zeros_like(matrix, dtype=np.float32)
        for row_idx in range(matrix.shape[0]):
            nz = np.where(matrix[row_idx] > 0)[0]
            if nz.size == 0:
                continue
            ordered = nz[np.argsort(matrix[row_idx, nz])[::-1]]
            rank_values = 1.0 / np.log2(np.arange(2, ordered.size + 2))
            weights[row_idx, ordered] = rank_values.astype(np.float32)
    elif mode == "tfidf":
        transformer = TfidfTransformer(norm="l1", use_idf=True, smooth_idf=True, sublinear_tf=True)
        weights = transformer.fit_transform(matrix).toarray().astype(np.float32)
    else:
        raise ValueError(f"Unsupported support weight mode: {mode}")
    return weights.astype(np.float32)


def build_cell_support_bundle(
    data_path: str,
    gene_graph: GeneGraphBundle,
    label_col: Optional[str] = "Celltype",
    support_weight_mode: str = "log1p_count",
    target_sum: float = 1e4,
) -> CellSupportBundle:
    adata = load_adata(data_path)
    label_col = _infer_label_col(adata, label_col)

    raw = adata.raw.X.copy() if adata.raw is not None else adata.X.copy()
    work = adata.copy()
    work.X = raw
    work = work[:, gene_graph.gene_names].copy()

    counts = _dense_array(work.X).astype(np.float32)
    size_factors = counts.sum(axis=1).astype(np.float32)
    size_factors = size_factors / np.median(size_factors[size_factors > 0])

    sc.pp.normalize_total(work, target_sum=target_sum)
    sc.pp.log1p(work)
    x_dense = _dense_array(work.X).astype(np.float32)
    weights = compute_support_weights(counts, mode=support_weight_mode)

    support_indices: List[np.ndarray] = []
    support_weights: List[np.ndarray] = []
    for row_idx in range(counts.shape[0]):
        nz = np.where(counts[row_idx] > 0)[0]
        if nz.size == 0:
            nz = np.array([int(np.argmax(counts[row_idx]))], dtype=np.int64)
        support_indices.append(nz.astype(np.int64))
        row_w = weights[row_idx, nz]
        if row_w.sum() <= 0:
            row_w = np.ones_like(row_w, dtype=np.float32)
        row_w = row_w / row_w.sum()
        support_weights.append(row_w.astype(np.float32))

    labels, label_names = _encode_labels(np.asarray(work.obs[label_col]))
    return CellSupportBundle(
        x_dense=x_dense,
        support_indices=support_indices,
        support_weights=support_weights,
        labels=labels,
        label_names=label_names,
        gene_names=list(work.var_names),
        size_factors=size_factors,
        cell_names=work.obs_names.tolist(),
    )
