"""Data loading and preprocessing for maskdiffusion.

This module follows the unified preprocessing philosophy from the design doc:
  1. Load raw counts
  2. Construct support from raw counts BEFORE normalization
  3. Normalize + log1p for expression values
  4. Select HVG and scale to [0, 1]
  5. Preserve support = (raw_counts > 0) throughout

This means support is always constructed from the raw, unnormalized data,
while expression values are normalized. The two are never mixed.
"""

import warnings
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch


# ── Sparse matrix helpers ──────────────────────────────────────────────────────


def _to_dense(X) -> np.ndarray:
    """Convert sparse or dense matrix to dense numpy array."""
    if sp.issparse(X):
        return X.toarray().astype(np.float32)
    if hasattr(X, "toarray"):
        return X.toarray().astype(np.float32)
    return np.asarray(X, dtype=np.float32)


def _is_likely_raw(X: np.ndarray, sample_n: int = 512) -> bool:
    """Heuristic: raw counts should be close to integer values."""
    flat = X[:sample_n].flatten()
    return np.all(flat >= 0) and np.allclose(flat, np.round(flat), atol=1e-4)


# ── Main dataclass ────────────────────────────────────────────────────────────


@dataclass
class SCDatasetBundle:
    """Container for a preprocessed single-cell dataset.

    Holds three parallel arrays:
      - values:  normalized gene expression, shape (n_cells, n_genes), range [0, 1]
      - support: binary support from RAW counts, shape (n_cells, n_genes), in {0, 1}
      - labels:  integer cell-type labels, shape (n_cells,)

    And metadata: adata, gene_names, input_mode (raw/log1p).
    """

    adata: sc.AnnData
    values: np.ndarray          # (n_cells, n_genes), normalized, [0, 1]
    support: np.ndarray         # (n_cells, n_genes), binary from raw counts
    labels: np.ndarray         # (n_cells,), integer-encoded cell types
    gene_names: np.ndarray     # (n_genes,)
    x_max: float               # max value used for [0,1] scaling
    input_mode: str            # "raw" or "log1p"
    label_key: str             # obs column name for labels


# ── Label extraction ─────────────────────────────────────────────────────────


def _detect_label_key(adata: sc.AnnData) -> str:
    """Find the most likely cell-type label column in adata.obs.

    Priority: cell_type > label > Cluster > cluster > leiden > any
    """
    priority = [
        "cell_type", "celltype", "cell_type_ontology_term_id",
        "cell_label", "label",
        "Cluster", "cluster",
        "leiden",
        "celltype_parsed",
        "cell_type.l1",
        "Celltype",  # SRP182008 / PlantNet format
    ]
    for key in priority:
        if key in adata.obs.columns:
            return key
    # Fall back to any categorical column
    cats = [c for c in adata.obs.columns if adata.obs[c].dtype.name == "category"]
    if cats:
        return cats[0]
    # Fall back to any column with few unique values
    n_uniq = adata.obs.nunique()
    candidates = n_uniq[n_uniq < adata.n_obs * 0.5].index.tolist()
    if candidates:
        return candidates[0]
    raise ValueError("Could not detect a cell-type label column in adata.obs. "
                     "Please pass --label_key explicitly.")


def _extract_labels(adata: sc.AnnData, label_key: Optional[str], n_clusters: int) -> Tuple[np.ndarray, str, int]:
    """Extract integer-encoded labels from adata.obs.

    Returns: (labels, label_key, n_clusters)
    """
    if label_key is None:
        label_key = _detect_label_key(adata)

    col = adata.obs[label_key]

    # Convert to string for uniform handling (works for both categorical and string cols)
    col_str = col.astype(str).values
    unique_vals = sorted(set(v for v in col_str if v != "nan" and v != "None"))
    label_map = {cat: i for i, cat in enumerate(unique_vals)}
    labels = np.array([label_map.get(v, -1) for v in col_str], dtype=np.int64)

    n_clusters_detected = len(unique_vals)
    if n_clusters > 0:
        assert n_clusters == n_clusters_detected, (
            f"Label column '{label_key}' has {n_clusters_detected} categories "
            f"but --n_clusters={n_clusters} was specified."
        )

    return labels, label_key, n_clusters_detected


# ── Main loading function ──────────────────────────────────────────────────────


def load_sc_dataset(
    file_path: str,
    n_top_genes: int = 2000,
    input_mode: str = "auto",
    normalize_total: float = 1e4,
    min_genes: int = 200,
    min_cells: int = 3,
    label_key: Optional[str] = None,
    n_clusters: int = 0,
    target_sum: float = 1e4,
) -> SCDatasetBundle:
    """Load and preprocess a single-cell dataset for maskdiffusion.

    Processing steps (matching design doc unification):
      1. Load h5ad
      2. Filter low-coverage cells and genes
      3. Construct SUPPORT from raw counts (before normalization!)
      4. Normalize total counts + log1p
      5. Select top N highly variable genes (seurat flavor)
      6. Scale expression values to [0, 1] using global x_max
      7. Extract labels

    Args:
        file_path: Path to .h5ad file.
        n_top_genes: Number of highly variable genes to select.
        input_mode: "auto" to detect raw/log1p from data, "raw" or "log1p" to force.
        normalize_total: Total count to normalize to (per cell).
        min_genes: Minimum genes per cell to keep.
        min_cells: Minimum cells per gene to keep.
        label_key: obs column name for labels (auto-detected if None).
        n_clusters: Number of expected clusters (for validation).
        target_sum: Synonym for normalize_total.

    Returns:
        SCDatasetBundle with .values, .support, .labels, etc.
    """
    # ── 1. Load ────────────────────────────────────────────────────────────────
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        adata = sc.read_h5ad(file_path, backed=False)

    # Handle legacy h5ad formats
    if sp.issparse(adata.X) and not hasattr(adata.X, "toarray"):
        adata.X = sp.csc_matrix(adata.X)

    # ── 2. Light filtering (before any normalization) ────────────────────────
    adata.raw = adata.copy()  # Keep raw layer

    n_obs0, n_var0 = adata.n_obs, adata.n_vars
    if min_genes and min_genes > 0:
        sc.pp.filter_cells(adata, min_genes=min_genes)
    if min_cells and min_cells > 0:
        sc.pp.filter_genes(adata, min_cells=min_cells)

    if adata.n_obs < 50:
        raise ValueError(f"Too few cells ({adata.n_obs}) after filtering. "
                         "Try lowering --min_genes or --min_cells.")
    if adata.n_vars < 100:
        raise ValueError(f"Too few genes ({adata.n_vars}) after filtering. "
                         "Try lowering --min_cells.")

    # ── 3. Detect input mode (raw vs log1p) from adata.X ────────────────────
    # We detect from the MAIN layer (adata.X), not raw.X, because:
    #   - If adata.X is already log1p normalized, we should NOT re-normalize
    #   - raw.X is only used for support construction
    if input_mode == "auto":
        X_main = _to_dense(adata.X[: min(256, adata.n_obs)])
        detected = "raw" if _is_likely_raw(X_main) else "log1p"
        print(f"[load_sc_dataset] Detected input mode from adata.X: '{detected}'")
        input_mode = detected

    # ── 4. Normalize + log1p for expression values ───────────────────────────
    # Work on a copy for the values pipeline
    work = adata.copy()
    if input_mode == "raw":
        sc.pp.normalize_total(work, target_sum=target_sum)
        sc.pp.log1p(work)
    # If already log1p, skip normalization

    # ── 5. Select highly variable genes ───────────────────────────────────────
    n_hvg = min(n_top_genes, work.n_vars)
    try:
        if input_mode == "raw":
            sc.pp.highly_variable_genes(
                work,
                flavor="seurat_v3",
                n_top_genes=n_hvg,
                subset=True,
            )
        else:
            # Already log1p: use standard seurat flavor
            sc.pp.highly_variable_genes(
                work,
                flavor="seurat",
                n_top_genes=n_hvg,
                subset=True,
            )
    except Exception as e:
        print(f"  Warning: HVG selection failed ({e}), falling back to variance-based selection")
        from sklearn.feature_selection import VarianceThreshold
        vt = VarianceThreshold(threshold=0.0)
        X_dense = _to_dense(work.X)
        if X_dense.shape[1] > n_hvg:
            vt.fit(X_dense)
            var_order = vt.variances_.argsort()[::-1][:n_hvg]
            work = work[:, var_order]

    if work.n_vars < 100:
        raise ValueError(
            f"Only {work.n_vars} genes remain after HVG selection. "
            "Try reducing --n_top_genes."
        )

    # ── 6. Construct SUPPORT from raw counts, SUBSETTED to HVG genes ──────────
    # Key unification: support is always from raw counts, matched to the selected HVG set.
    # We need to subset the raw matrix to the same gene indices as work.
    if hasattr(work.var, "_varlier"):  # Seurat-style HVG
        hvg_names = work.var_names.tolist()
    else:
        hvg_names = work.var_names.tolist()

    if "support" in work.layers:
        support = _to_dense(work.layers["support"]).astype(np.float32)
    elif "counts" in work.layers:
        support = (_to_dense(work.layers["counts"]) > 0).astype(np.float32)
    elif adata.raw is not None and adata.raw.X is not None:
        # Get raw counts for the full gene set
        X_raw_full = _to_dense(adata.raw.X)
        # Find the indices of the HVG genes in the original var_names
        full_var_names = adata.raw.var_names.tolist() if hasattr(adata.raw, "var_names") and adata.raw.var_names is not None else None
        if full_var_names:
            # Map HVG gene names to indices in the raw layer
            raw_gene_to_idx = {g: i for i, g in enumerate(full_var_names)}
            hvg_indices = [raw_gene_to_idx[g] for g in hvg_names]
            X_raw_subset = X_raw_full[:, hvg_indices]
        else:
            # Fall back: assume raw.var_names matches adata.var_names ordering
            X_raw_subset = X_raw_full[:, :work.n_vars]
    else:
        # No raw layer: use the main X matrix
        X_raw_subset = _to_dense(work.X)
        support = (X_raw_subset > 0).astype(np.float32)  # (n_cells, n_genes_HVG)
    if "support" not in work.layers and "counts" not in work.layers:
        support = (X_raw_subset > 0).astype(np.float32)  # (n_cells, n_genes_HVG)
    print(f"[load_sc_dataset] Support sparsity: {1.0 - support.mean():.3f} "
          f"(support density: {support.mean():.3f})")

    # ── 7. Extract expression values ─────────────────────────────────────────
    values = _to_dense(work.X)  # shape: (n_cells, n_genes_selected)

    # Scale to [0, 1] by global x_max
    x_max = values.max()
    if x_max > 0:
        values = values / x_max
    else:
        values = values.astype(np.float32)

    # Clip any stray negative values (can happen with some log1p data)
    values = np.clip(values, 0.0, 1.0).astype(np.float32)

    # ── 8. Extract labels ────────────────────────────────────────────────────
    labels, label_key, n_clusters_detected = _extract_labels(
        adata, label_key, n_clusters
    )
    n_clusters = n_clusters or n_clusters_detected

    gene_names = work.var_names.to_numpy() if hasattr(work.var_names, "to_numpy") else np.array(work.var_names)

    print(f"[load_sc_dataset] Final: {values.shape[0]} cells × {values.shape[1]} genes, "
          f"{n_clusters} clusters, input_mode={input_mode}")

    return SCDatasetBundle(
        adata=work,
        values=values,
        support=support,
        labels=labels,
        gene_names=gene_names,
        x_max=float(x_max),
        input_mode=input_mode,
        label_key=label_key,
    )


# ── PyTorch DataLoader helper ─────────────────────────────────────────────────


def build_dataloader(
    values: np.ndarray,
    support: Optional[np.ndarray] = None,
    labels: Optional[np.ndarray] = None,
    batch_size: int = 256,
    shuffle: bool = True,
    drop_last: bool = False,
):
    """Build a PyTorch DataLoader from numpy arrays.

    Each batch is a dict: {'values': tensor, 'support': tensor, 'labels': tensor}.
    """
    from torch import tensor

    data = {"values": tensor(values, dtype=torch.float32)}
    if support is not None:
        data["support"] = tensor(support, dtype=torch.float32)
    if labels is not None:
        data["labels"] = tensor(labels, dtype=torch.long)

    class SCDataset:
        def __init__(self, data_dict):
            self.data = data_dict
            self.n = len(data_dict["values"])

        def __len__(self):
            return self.n

        def __getitem__(self, idx):
            return {k: v[idx] for k, v in self.data.items()}

    dataset = SCDataset(data)
    return torch.utils.data.DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=drop_last,
        num_workers=0,
    )
