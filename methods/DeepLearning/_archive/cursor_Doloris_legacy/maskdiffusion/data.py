# -*- coding: utf-8 -*-
"""
Data loading and preprocessing utilities for maskdiffusion.
"""

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import torch
from torch.utils.data import Dataset, DataLoader


def _read_legacy_h5ad(data_path: str) -> sc.AnnData:
    """
    Read legacy h5ad format that is incompatible with current anndata version.

    This handles older h5ad files that use np.string_ instead of np.bytes_.

    Args:
        data_path: Path to h5ad file

    Returns:
        AnnData object
    """
    with h5py.File(data_path, 'r') as f:
        # Read X matrix
        X = f['X'][:]

        # Read obs
        obs_data = f['obs']
        obs_dict = {}
        for name in obs_data.dtype.names:
            if obs_data.dtype[name].kind in ('U', 'S'):
                # String fields
                obs_dict[name] = [s.decode() if isinstance(s, bytes) else s
                                  for s in obs_data[name]]
            else:
                obs_dict[name] = obs_data[name][:]
        obs = pd.DataFrame(obs_dict)
        obs.index = obs['index'].astype(str)
        obs = obs.drop('index', axis=1)

        # Read var
        var_data = f['var']
        var_dict = {}
        for name in var_data.dtype.names:
            if var_data.dtype[name].kind in ('U', 'S'):
                # String fields
                var_dict[name] = [s.decode() if isinstance(s, bytes) else s
                                  for s in var_data[name]]
            else:
                var_dict[name] = var_data[name][:]
        var = pd.DataFrame(var_dict)
        var.index = var['index'].astype(str)
        var = var.drop('index', axis=1)

        # Create AnnData
        adata = sc.AnnData(X=X, obs=obs, var=var)

        return adata


class ScRNADataset(Dataset):
    """
    PyTorch Dataset for single-cell RNA-seq data.
    """

    def __init__(self, X: np.ndarray, y: np.ndarray = None):
        """
        Args:
            X: Gene expression matrix (n_cells, n_genes)
            y: Cell type labels (n_cells,) or None
        """
        self.X = torch.from_numpy(X.astype(np.float32))
        self.y = torch.from_numpy(y.astype(np.int64)) if y is not None else None

    def __len__(self):
        return len(self.X)

    def __getitem__(self, idx):
        if self.y is not None:
            return self.X[idx], self.y[idx]
        return self.X[idx]


def load_and_preprocess(
    data_path: str,
    n_top_genes: int = 1000,
    log_transform: bool = True,
    normalize: bool = True,
    scale: bool = False,  # Default: no scaling to preserve sparsity
    preserve_sparsity: bool = True,  # New: preserve sparsity structure
) -> tuple:
    """
    Load and preprocess h5ad data for maskdiffusion.

    Key optimization: Preserve sparsity structure for better mask prediction.
    Instead of z-score scaling (which destroys sparsity), we clip to [0, 1].

    Args:
        data_path: Path to h5ad file
        n_top_genes: Number of highly variable genes to keep
        log_transform: Whether to apply log1p transformation
        normalize: Whether to normalize per cell
        scale: Whether to scale genes (NOT recommended - destroys sparsity)
        preserve_sparsity: Whether to preserve sparsity (recommended)

    Returns:
        tuple: (X, Y, adata) where X is expression matrix, Y is labels, adata is AnnData
    """
    # Try loading with anndata first
    try:
        adata = sc.read_h5ad(data_path)
    except (AttributeError, ValueError) as e:
        # Fallback: read legacy h5ad format directly
        print(f"Warning: Failed to read with anndata ({e}), trying legacy format...")
        adata = _read_legacy_h5ad(data_path)

    # Raw backup
    adata.raw = adata.copy()

    # Smart filtering: only filter if data has many genes/cells
    # Some preprocessed datasets have already been filtered
    if adata.n_vars > 500 and adata.n_obs > 100:
        sc.pp.filter_cells(adata, min_genes=50)  # Lower threshold for prefiltered data
        sc.pp.filter_genes(adata, min_cells=3)
        print(f"Filtered data: {adata.n_obs} cells, {adata.n_vars} genes")
    else:
        print(f"Using prefiltered data: {adata.n_obs} cells, {adata.n_vars} genes")

    # Normalize per cell
    if normalize:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

    # Log transform
    if log_transform:
        sc.pp.log1p(adata)

    # Select highly variable genes
    if n_top_genes > 0:
        # Ensure n_top_genes doesn't exceed available genes
        n_genes_available = adata.n_vars
        n_top_genes_actual = min(n_top_genes, n_genes_available)

        if n_top_genes_actual < n_genes_available:
            print(f"Selecting top {n_top_genes_actual} genes from {n_genes_available} available")

        if n_top_genes_actual > 1:
            try:
                sc.pp.highly_variable_genes(
                    adata,
                    n_top_genes=n_top_genes_actual,
                    flavor='seurat_v3' if hasattr(sc.pp, 'high') else 'seurat'
                )
                adata = adata[:, adata.var.highly_variable]
            except ValueError as e:
                # If HVG selection fails (e.g., too few genes), use all genes
                print(f"Warning: HVG selection failed ({e}), using all genes")
                n_top_genes_actual = adata.n_vars
        else:
            # If only 1 gene available, can't do HVG selection
            print(f"Warning: Only {n_genes_available} genes available, using all genes")

    # Get expression matrix BEFORE scaling
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = X.astype(np.float32)

    # CRITICAL: Preserve sparsity by clipping to [0, 1] instead of z-score scaling
    if preserve_sparsity:
        # This preserves the sparsity structure: zeros stay zeros, 
        # expressed values are normalized to [0, 1]
        X = np.clip(X, 0, None)  # Remove negative values first
        # Normalize each gene to [0, 1] range for better model training
        gene_max = X.max(axis=0, keepdims=True)
        gene_max = np.where(gene_max > 0, gene_max, 1.0)  # Avoid division by zero
        X = X / gene_max
        # Any value that was 0 stays 0, expressed values are in [0, 1]
        adata.X = X
    elif scale:
        # Fallback to original scaling
        sc.pp.scale(adata, max_value=10)
        X = adata.X
        if hasattr(X, 'toarray'):
            X = X.toarray()
        X = X.astype(np.float32)

    # Final safety check for NaN/inf
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    # Get labels
    label_col = None
    label_candidates = [
        'cell_type', 'Celltype', 'celltype', 'cell_label', 'label', 'Cluster',
        'cell_type', 'clusters', 'Cluster', 'cluster', 'celltype',
        'cell_cluster', 'celltype', 'cell_label',
        'paul15_clusters',  # Special case for paul15 dataset
    ]
    for col in label_candidates:
        if col in adata.obs.columns:
            label_col = col
            break

    if label_col is None:
        Y = np.zeros(len(adata), dtype=np.int64)
    else:
        labels, indices = np.unique(adata.obs[label_col].values, return_inverse=True)
        Y = indices.astype(np.int64)
        print(f"Detected {len(labels)} clusters from column '{label_col}'")

    return X, Y, adata


def create_dataloader(
    X: np.ndarray,
    y: np.ndarray = None,
    batch_size: int = 128,
    shuffle: bool = True,
    drop_last: bool = True,
) -> DataLoader:
    """
    Create PyTorch DataLoader from expression matrix.

    Args:
        X: Gene expression matrix
        y: Cell type labels
        batch_size: Batch size
        shuffle: Whether to shuffle
        drop_last: Whether to drop last incomplete batch

    Returns:
        DataLoader
    """
    dataset = ScRNADataset(X, y)
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        drop_last=drop_last,
    )


def get_expression_mask(X: np.ndarray, threshold: float = 0.0) -> np.ndarray:
    """
    Get binary mask for expressed genes.

    Args:
        X: Expression matrix
        threshold: Expression threshold for "expressed"

    Returns:
        Binary mask where 1 = expressed, 0 = zero
    """
    return (X > threshold).astype(np.float32)


def split_train_val(
    X: np.ndarray,
    y: np.ndarray = None,
    val_ratio: float = 0.1,
    seed: int = 42,
) -> tuple:
    """
    Split data into train and validation sets.

    Args:
        X: Expression matrix
        y: Labels
        val_ratio: Validation ratio
        seed: Random seed

    Returns:
        tuple of (X_train, X_val, y_train, y_val)
    """
    np.random.seed(seed)
    n_samples = len(X)
    indices = np.random.permutation(n_samples)
    n_val = int(n_samples * val_ratio)

    train_idx = indices[n_val:]
    val_idx = indices[:n_val]

    X_train = X[train_idx]
    X_val = X[val_idx]

    if y is not None:
        y_train = y[train_idx]
        y_val = y[val_idx]
        return X_train, X_val, y_train, y_val

    return X_train, X_val
