from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import torch
from torch.utils.data import Dataset


@dataclass
class APABundle:
    x: np.ndarray
    labels: np.ndarray | None
    label_key: str | None
    gene_stats: np.ndarray
    prototypes: np.ndarray
    selected_gene_indices: np.ndarray
    gene_names: np.ndarray
    preprocess_config: dict[str, Any]


class APAExpressionDataset(Dataset):
    """Training dataset. It intentionally returns no labels."""

    def __init__(self, x: np.ndarray) -> None:
        self.x = torch.as_tensor(np.asarray(x, dtype=np.float32))

    def __len__(self) -> int:
        return int(self.x.shape[0])

    def __getitem__(self, index: int) -> dict[str, torch.Tensor]:
        return {"index": torch.tensor(index, dtype=torch.long), "x": self.x[index]}


def assert_no_training_labels(batch: dict[str, Any]) -> None:
    forbidden = {"label", "labels", "y", "y_true", "cell_type", "celltype", "n_clusters"}
    leaked = sorted(forbidden.intersection(batch.keys()))
    if leaked:
        raise RuntimeError(f"Training batch leaked label-like fields: {leaked}")


def _dense(x: Any) -> np.ndarray:
    import scipy.sparse as sp

    if sp.issparse(x):
        return x.toarray()
    return np.asarray(x)


def _looks_like_counts(x: np.ndarray) -> bool:
    if x.size == 0:
        return False
    sample = x if x.size <= 1_000_000 else x.reshape(-1)[:1_000_000]
    return bool(np.nanmin(sample) >= 0 and np.allclose(sample, np.round(sample), atol=1.0e-6))


LABEL_CANDIDATES = (
    "cell_type",
    "celltype",
    "CellType",
    "cell type",
    "cell_types",
    "celltype_l1",
    "celltype.l1",
    "celltype_l2",
    "celltype.l2",
    "cell_type1",
    "Cell.type",
    "Cell Types",
    "cell_types",
    "cell_ontology_class",
    "cell_ontology",
    "annotation",
    "annotations",
    "Annotation",
    "manual_annotation",
    "assigned_cell_type",
    "cell_label",
    "cell_labels",
    "label",
    "labels",
    "Label",
    "class",
    "Class",
    "type",
    "y",
    "Y",
)


def _label_array(adata: Any) -> tuple[np.ndarray | None, str | None]:
    for key in LABEL_CANDIDATES:
        if key in adata.obs:
            values = np.asarray(adata.obs[key].astype("category").cat.codes, dtype=np.int64)
            return values, key
    return None, None


def _raw_count_source(adata: Any) -> tuple[np.ndarray | None, np.ndarray | None, str | None]:
    if "counts" in adata.layers:
        return _dense(adata.layers["counts"]).astype(np.float32), np.asarray(adata.var_names, dtype=str), 'layers["counts"]'
    if adata.raw is not None:
        return _dense(adata.raw.X).astype(np.float32), np.asarray(adata.raw.var_names, dtype=str), "adata.raw.X"
    x = _dense(adata.X).astype(np.float32)
    if _looks_like_counts(x):
        return x, np.asarray(adata.var_names, dtype=str), "adata.X"
    return None, None, None


def _hvg_indices(x: np.ndarray, n_top_genes: int) -> np.ndarray:
    n_genes = int(x.shape[1])
    if n_top_genes <= 0 or n_top_genes >= n_genes:
        return np.arange(n_genes, dtype=np.int64)
    var = np.var(x, axis=0)
    return np.argsort(-var, kind="mergesort")[: int(n_top_genes)].astype(np.int64)


def compute_gene_stats(x: np.ndarray) -> np.ndarray:
    mean = x.mean(axis=0)
    var = x.var(axis=0)
    zero_rate = (np.abs(x) <= 1.0e-8).mean(axis=0)
    ranks = np.empty(x.shape[1], dtype=np.float32)
    ranks[np.argsort(-var, kind="mergesort")] = np.arange(x.shape[1], dtype=np.float32)
    denom = max(1, x.shape[1] - 1)
    ranks = ranks / float(denom)
    stats = np.stack([mean, var, zero_rate, ranks], axis=1).astype(np.float32)
    stats[:, :2] = (stats[:, :2] - stats[:, :2].mean(axis=0, keepdims=True)) / (
        stats[:, :2].std(axis=0, keepdims=True) + 1.0e-6
    )
    return stats


def build_prototypes(x: np.ndarray, n_prototypes: int, pca_dim: int, seed: int) -> np.ndarray:
    from sklearn.cluster import KMeans
    from sklearn.decomposition import PCA

    n_cells = int(x.shape[0])
    k = max(1, min(int(n_prototypes), n_cells))
    dim = max(1, min(int(pca_dim), int(x.shape[1]), max(1, n_cells - 1)))
    if n_cells <= 1:
        reduced = np.zeros((n_cells, dim), dtype=np.float32)
    else:
        reduced = PCA(n_components=dim, random_state=seed).fit_transform(x).astype(np.float32)
    if k == 1:
        centers = reduced.mean(axis=0, keepdims=True)
    else:
        centers = KMeans(n_clusters=k, random_state=seed, n_init=10).fit(reduced).cluster_centers_
    return np.asarray(centers, dtype=np.float32)


def load_apa_data(
    data_path: str,
    *,
    input_mode: str,
    target_sum: float,
    n_top_genes: int,
    scale_input: bool,
    n_prototypes: int,
    pca_dim: int,
    seed: int,
    require_labels: bool = True,
) -> APABundle:
    import anndata as ad

    adata = ad.read_h5ad(data_path)
    x_adata = _dense(adata.X).astype(np.float32)
    labels, label_key = _label_array(adata)
    if labels is None and require_labels:
        raise ValueError(
            "No label column found in adata.obs. Tried: "
            + ", ".join(LABEL_CANDIDATES)
            + ". Use --skip_eval true only when labels are intentionally unavailable."
        )
    count_x, count_gene_names, count_source = _raw_count_source(adata)
    mode = input_mode
    if mode == "auto":
        mode = "raw" if count_x is not None else "log1p"
    if mode == "raw":
        if count_x is None:
            raise ValueError("input_mode='raw' requested, but no raw count source was found in layers['counts'], adata.raw.X, or adata.X.")
        x = count_x.copy()
        source_gene_names = count_gene_names
        totals = x.sum(axis=1, keepdims=True)
        totals[totals <= 0] = 1.0
        x = np.log1p(x / totals * float(target_sum)).astype(np.float32)
        matrix_source = count_source
    elif mode == "log1p":
        x = x_adata.copy()
        source_gene_names = np.asarray(adata.var_names, dtype=str)
        matrix_source = "adata.X"
    else:
        raise ValueError(f"Unsupported input_mode={input_mode!r}")
    selected = _hvg_indices(x, int(n_top_genes))
    x = x[:, selected].astype(np.float32, copy=False)
    if scale_input:
        x = ((x - x.mean(axis=0, keepdims=True)) / (x.std(axis=0, keepdims=True) + 1.0e-6)).astype(np.float32)
    gene_names = np.asarray(source_gene_names, dtype=str)[selected]
    gene_stats = compute_gene_stats(x)
    prototypes = build_prototypes(x, int(n_prototypes), int(pca_dim), int(seed))
    return APABundle(
        x=x,
        labels=labels,
        label_key=label_key,
        gene_stats=gene_stats,
        prototypes=prototypes,
        selected_gene_indices=selected,
        gene_names=gene_names,
        preprocess_config={
            "input_mode_requested": input_mode,
            "input_mode_resolved": mode,
            "raw_count_source": matrix_source if mode == "raw" else None,
            "label_key": label_key,
            "labels_available": labels is not None,
            "target_sum": float(target_sum),
            "n_top_genes": int(n_top_genes),
            "actual_n_genes": int(x.shape[1]),
            "scale_input": bool(scale_input),
            "feature_space_source": "hvg_variance",
        },
    )


class ScMAEShuffleCorruption:
    """Gene-wise scMAE-style replacement values; the generator never creates V."""

    def __init__(self, x_full: torch.Tensor, *, seed: int, atol: float, rtol: float) -> None:
        self.x_full = x_full
        self.n_cells, self.n_genes = [int(v) for v in x_full.shape]
        self.atol = float(atol)
        self.rtol = float(rtol)
        rng = np.random.default_rng(int(seed) + 7919)
        perms = np.empty((self.n_cells, self.n_genes), dtype=np.int64)
        for gene_id in range(self.n_genes):
            perms[:, gene_id] = rng.permutation(self.n_cells)
        self.permutations = torch.as_tensor(perms, dtype=torch.long, device=x_full.device)

    def sample(self, batch_indices: torch.Tensor) -> dict[str, torch.Tensor]:
        idx = batch_indices.long().to(self.x_full.device)
        donor = self.permutations[idx]
        gene_ids = torch.arange(self.n_genes, device=self.x_full.device).view(1, self.n_genes).expand_as(donor)
        replacement = self.x_full[donor, gene_ids]
        original = self.x_full[idx]
        effective = ~torch.isclose(replacement, original, atol=self.atol, rtol=self.rtol)
        return {"replacement": replacement, "effective": effective.float(), "donor_indices": donor}


def save_text_lines(path: Path, values: np.ndarray) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(map(str, values)) + "\n")
