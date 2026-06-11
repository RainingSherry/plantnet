#!/usr/bin/env python3
"""
prepare_dataset.py
=================
Convert raw .h5 / .h5ad files to a canonical, annotated .h5ad for benchmarking.

Accepts:
  - .h5  : scMAE-style HDF5 (X / Y / cell_names / gene_names)
  - .h5ad: Already AnnData, re-annotated only

Outputs a canonical .h5ad with:
  - adata.X           : expression matrix (float)
  - adata.obs         : metadata including resolved_label_key
  - adata.var_names   : gene names
  - adata.obs_names   : cell names
  - adata.uns["source_format"]   : "h5" or "h5ad"
  - adata.uns["source_file"]     : original path
  - adata.uns["resolved_label_key"]
  - adata.uns["matrix_key"]
  - adata.uns["n_clusters"]

Usage:
    # Auto-detect everything
    python scripts/prepare_dataset.py \
        --input_path data/scMAE/Quake_Smart-seq2_Lung.h5 \
        --dataset_name Quake_Smart-seq2_Lung \
        --output_dir data/processed

    # With explicit keys
    python scripts/prepare_dataset.py \
        --input_path data/scMAE/Wang.h5 \
        --dataset_name Wang \
        --output_dir data/processed \
        --matrix_key X \
        --label_key Y \
        --n_clusters 14

    # From .h5ad (re-annotate only)
    python scripts/prepare_dataset.py \
        --input_path data/SRP182008.h5ad \
        --dataset_name SRP182008 \
        --output_dir data/processed
"""

import argparse
import json
import os
import sys

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp


CANDIDATE_MATRIX_KEYS = ["X", "data", "exprs", "matrix", "counts", "raw", "value"]
CANDIDATE_LABEL_KEYS = [
    "Y", "label", "labels", "cell_type", "celltype", "cell_type_label",
    "cell_label", "cluster", "clusters", "type", "group", "group_id", "Group",
]
CANDIDATE_CELL_KEYS = ["cell_names", "obs_names", "barcode", "barcodes", "cell_id", "cells"]
CANDIDATE_GENE_KEYS = ["gene_names", "var_names", "gene", "genes", "gene_id"]


# ─── Schema detection ──────────────────────────────────────────────────────────

def auto_detect_matrix(f):
    """Return (key, shape) of the expression matrix."""
    for key in CANDIDATE_MATRIX_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset):
                return key, ds.shape
            if isinstance(ds, h5py.Group):
                if "data" in ds and "indices" in ds:
                    return key, ds["shape"][...]
                if "data" in ds and "indptr" in ds:
                    return key, ds["shape"][...]
    # Fallback: find largest 2-D array
    candidates = []
    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset) and len(obj.shape) >= 2:
            candidates.append((name, obj.shape))
    f.visititems(visitor)
    candidates.sort(key=lambda x: -np.prod(x[1]))
    return (candidates[0] if candidates else (None, None))


def auto_detect_labels(f):
    """Return (key, shape) of the label array."""
    for key in CANDIDATE_LABEL_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset) and len(ds.shape) == 1:
                return key, ds.shape
    # Fallback: find largest 1-D array
    candidates = []
    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset) and len(obj.shape) == 1:
            candidates.append((name, obj.shape))
    f.visititems(visitor)
    candidates.sort(key=lambda x: -x[1][0])
    return (candidates[0] if candidates else (None, None))


def auto_detect_names(f, n_cells, n_genes):
    """Return (cell_key, gene_key) that match n_cells / n_genes."""
    cell_key = None
    for key in CANDIDATE_CELL_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset) and ds.shape[0] == n_cells:
                cell_key = key
                break
    gene_key = None
    for key in CANDIDATE_GENE_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset) and ds.shape[0] == n_genes:
                gene_key = key
                break
    return cell_key, gene_key


def infer_n_clusters(labels):
    """Count unique labels. Labels can be int, float, or string."""
    try:
        unique = pd.unique(labels)
        return len(unique)
    except Exception:
        return None


# ─── .h5 → AnnData conversion ────────────────────────────────────────────────

def h5_to_anndata(
    h5_path,
    matrix_key=None,
    label_key=None,
    cell_key=None,
    gene_key=None,
    label_is_integer=False,
    verbose=True,
):
    """
    Load a scMAE-style .h5 file and return an AnnData object.
    """
    if verbose:
        print(f"Reading {h5_path} ...")

    with h5py.File(h5_path, "r") as f:
        # ── Matrix ────────────────────────────────────────────────
        if matrix_key:
            mat_ds = f[matrix_key]
        else:
            detected, shape = auto_detect_matrix(f)
            matrix_key = detected
            if matrix_key is None:
                raise ValueError(
                    f"Cannot detect expression matrix. "
                    f"Available keys: {sorted(f.keys())}. "
                    f"Please provide --matrix_key manually."
                )
            mat_ds = f[matrix_key]

        if isinstance(mat_ds, h5py.Group):
            # Sparse CSR/CSC
            if "data" in mat_ds and "indptr" in mat_ds:
                mat = sp.csr_matrix(
                    (mat_ds["data"][...], mat_ds["indices"][...], mat_ds["indptr"][...]),
                    shape=mat_ds["shape"][...]
                )
            elif "data" in mat_ds and "indices" in mat_ds:
                mat = sp.csr_matrix(
                    (mat_ds["data"][...], mat_ds["indices"][...], mat_ds["indptr"][...]),
                    shape=mat_ds["shape"][...]
                )
            else:
                raise ValueError(f"Group {matrix_key} does not look like a sparse matrix.")
            n_cells, n_genes = mat.shape
            X = mat
        else:
            # Dense
            X = mat_ds[...]
            n_cells, n_genes = X.shape

        X = np.asarray(X, dtype=np.float32)

        # ── Labels ────────────────────────────────────────────────
        if label_key:
            label_ds = f[label_key]
        else:
            detected, shape = auto_detect_labels(f)
            label_key = detected
            if label_key is None:
                raise ValueError(
                    f"Cannot detect label field. "
                    f"Available keys: {sorted(f.keys())}. "
                    f"Please provide --label_key manually."
                )
            label_ds = f[label_key]

        raw_labels = label_ds[...]
        if hasattr(raw_labels, "tolist"):
            raw_labels = raw_labels.tolist()

        # Convert to string labels for consistent handling
        str_labels = np.array(raw_labels).astype(str)

        # Detect integer label mode
        label_is_integer = label_is_integer or np.issubdtype(label_ds.dtype, np.integer)

        n_clusters = infer_n_clusters(str_labels)
        if verbose:
            print(f"  Matrix: {matrix_key} → shape=({n_cells}, {n_genes})")
            print(f"  Labels: {label_key} → dtype={label_ds.dtype}, n_clusters={n_clusters}")

        # ── Cell names ────────────────────────────────────────────
        if cell_key:
            cell_names_ds = f[cell_key]
        else:
            cell_key, _ = auto_detect_names(f, n_cells, n_genes)
            cell_names_ds = f[cell_key] if cell_key else None

        if cell_names_ds is not None:
            raw_names = cell_names_ds[...]
            if hasattr(raw_names, "tolist"):
                raw_names = raw_names.tolist()
            obs_names = np.asarray(raw_names).astype(str)
            if len(obs_names) != n_cells:
                obs_names = np.arange(n_cells).astype(str)
        else:
            obs_names = np.arange(n_cells).astype(str)

        # ── Gene names ────────────────────────────────────────────
        if gene_key:
            gene_names_ds = f[gene_key]
        else:
            _, gene_key = auto_detect_names(f, n_cells, n_genes)
            gene_names_ds = f[gene_key] if gene_key else None

        if gene_names_ds is not None:
            raw_gene_names = gene_names_ds[...]
            if hasattr(raw_gene_names, "tolist"):
                raw_gene_names = raw_gene_names.tolist()
            var_names = np.asarray(raw_gene_names).astype(str)
            if len(var_names) != n_genes:
                var_names = np.arange(n_genes).astype(str)
        else:
            var_names = np.arange(n_genes).astype(str)

    # ── Build AnnData ────────────────────────────────────────────────────────
    obs = pd.DataFrame(index=pd.Index(obs_names, name="cell"))
    obs["resolved_label"] = str_labels
    obs["_label_is_integer"] = label_is_integer

    # Compute n_counts from the expression matrix (required by normalize_sc)
    if sp.issparse(X):
        counts_per_cell = np.asarray(X.sum(axis=1)).flatten()
    else:
        counts_per_cell = np.asarray(X).sum(axis=1)
    obs["n_counts"] = counts_per_cell.astype(np.float32)

    var = pd.DataFrame(index=pd.Index(var_names, name="gene"))

    # Detect if data is already normalized (float + non-integer values).
    # Check both dtype and value range: raw counts are typically int or have
    # large max values; normalized/log1p data is float in [0, ~20].
    X_max = float(np.max(X)) if np.prod(X.shape) > 0 else 0.0
    is_float_not_counts = not np.issubdtype(X.dtype, np.integer) and X_max < 100

    if sp.issparse(X):
        X_sparse = X
    else:
        X_sparse = sp.csr_matrix(X)

    adata = sc.AnnData(X=X_sparse, obs=obs, var=var)
    adata.uns["source_format"] = "h5"
    adata.uns["source_file"] = str(h5_path)
    adata.uns["matrix_key"] = matrix_key
    adata.uns["resolved_label_key"] = "resolved_label"
    adata.uns["n_clusters"] = n_clusters
    adata.uns["label_key"] = label_key
    adata.uns["label_is_integer"] = label_is_integer
    adata.uns["cell_names_key"] = cell_key
    adata.uns["gene_names_key"] = gene_key
    adata.uns["is_pre_normalized"] = is_float_not_counts
    if is_float_not_counts:
        adata.uns["pre_normalized_note"] = (
            f"Data is already normalized (float, max={X_max:.2f}). "
            "normalize_sc with logtrans_input=True will be skipped."
        )

    return adata, n_clusters


# ─── .h5ad re-annotation ─────────────────────────────────────────────────────

def h5ad_annotate(h5ad_path, label_key="auto", verbose=True):
    """
    Load an existing .h5ad and annotate its uns with metadata.
    If label_key is "auto", tries common column names.
    """
    adata = sc.read_h5ad(h5ad_path)

    if label_key == "auto":
        found = None
        for candidate in ["cell_type", "Celltype", "celltype", "cell_label",
                          "label", "cluster", "Group", "cell_type"]:
            if candidate in adata.obs.columns:
                found = candidate
                break
        if found is None:
            raise ValueError(
                f"Cannot auto-detect label column. "
                f"Available obs columns: {list(adata.obs.columns)}. "
                f"Please provide --label_key manually."
            )
        label_key = found

    # Ensure resolved_label exists
    if "resolved_label" not in adata.obs:
        adata.obs["resolved_label"] = adata.obs[label_key].astype(str)

    n_clusters = infer_n_clusters(adata.obs["resolved_label"])

    adata.uns["source_format"] = "h5ad"
    adata.uns["source_file"] = str(h5ad_path)
    adata.uns["resolved_label_key"] = "resolved_label"
    adata.uns["n_clusters"] = n_clusters
    adata.uns["label_key"] = label_key

    if verbose:
        print(f"  Loaded: shape=({adata.n_obs}, {adata.n_vars})")
        print(f"  Label column: {label_key} → n_clusters={n_clusters}")

    return adata, n_clusters


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Prepare a canonical .h5ad from raw .h5 or existing .h5ad",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Auto-detect everything from .h5
  python scripts/prepare_dataset.py \\
      --input_path data/scMAE/Pollen.h5 \\
      --dataset_name Pollen \\
      --output_dir data/processed

  # Auto-detect everything from .h5ad
  python scripts/prepare_dataset.py \\
      --input_path data/SRP182008.h5ad \\
      --dataset_name SRP182008 \\
      --output_dir data/processed

  # Manual key override
  python scripts/prepare_dataset.py \\
      --input_path data/scMAE/Wang.h5 \\
      --dataset_name Wang \\
      --output_dir data/processed \\
      --matrix_key X --label_key Y --n_clusters 14
""",
    )
    parser.add_argument("--input_path", type=str, required=True,
                        help="Path to input .h5 or .h5ad file")
    parser.add_argument("--dataset_name", type=str, default=None,
                        help="Dataset name (default: inferred from input filename)")
    parser.add_argument("--output_dir", type=str, default="data/processed",
                        help="Output directory (default: data/processed)")
    parser.add_argument("--label_key", type=str, default="auto",
                        help="Label column key in .h5 or obs column in .h5ad. "
                             "'auto' = try common names (default: auto)")
    parser.add_argument("--matrix_key", type=str, default=None,
                        help="Matrix key in .h5 file. 'auto' = auto-detect (default: auto)")
    parser.add_argument("--n_clusters", type=str, default="auto",
                        help="Number of clusters. 'auto' = infer from labels (default: auto)")
    parser.add_argument("--force", action="store_true",
                        help="Overwrite existing output .h5ad")
    parser.add_argument("--verbose", "-v", action="store_true")
    args = parser.parse_args()

    input_path = os.path.abspath(args.input_path)
    if not os.path.exists(input_path):
        print(f"ERROR: File not found: {input_path}")
        sys.exit(1)

    dataset_name = args.dataset_name or os.path.splitext(os.path.basename(input_path))[0]
    output_dir = os.path.abspath(args.output_dir)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"{dataset_name}.h5ad")

    if os.path.exists(output_path) and not args.force:
        print(f"Output already exists: {output_path}")
        print("Use --force to overwrite.")
        sys.exit(0)

    ext = os.path.splitext(input_path)[1].lower()

    if ext == ".h5":
        print(f"\n[{dataset_name}] Converting .h5 → .h5ad")
        print(f"  Input:  {input_path}")

        adata, n_clusters = h5_to_anndata(
            input_path,
            matrix_key=args.matrix_key,
            label_key=args.label_key if args.label_key != "auto" else None,
            verbose=args.verbose,
        )

        # Override n_clusters if user provided it
        if args.n_clusters != "auto":
            n_clusters = int(args.n_clusters)
            adata.uns["n_clusters"] = n_clusters
        elif adata.uns.get("n_clusters"):
            n_clusters = adata.uns["n_clusters"]

    elif ext == ".h5ad":
        print(f"\n[{dataset_name}] Re-annotating existing .h5ad")
        print(f"  Input:  {input_path}")

        adata, n_clusters = h5ad_annotate(
            input_path,
            label_key=args.label_key,
            verbose=args.verbose,
        )

        if args.n_clusters != "auto":
            n_clusters = int(args.n_clusters)
            adata.uns["n_clusters"] = n_clusters
        elif adata.uns.get("n_clusters"):
            n_clusters = adata.uns["n_clusters"]

    else:
        print(f"ERROR: Unsupported file extension: {ext} (expected .h5 or .h5ad)")
        sys.exit(1)

    print(f"  Output: {output_path}")
    print(f"  Shape:  ({adata.n_obs}, {adata.n_vars})")
    print(f"  n_clusters: {n_clusters}")
    print(f"  Label col: {adata.uns.get('resolved_label_key')}")

    adata.write_h5ad(output_path)
    print(f"\nSaved: {output_path}")

    # Write a companion .json sidecar with resolved metadata
    input_stat = os.stat(input_path)
    meta = {
        "dataset_name": dataset_name,
        "input_path": str(input_path),
        "input_size": int(input_stat.st_size),
        "input_mtime_ns": int(input_stat.st_mtime_ns),
        "output_path": str(output_path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_clusters": n_clusters,
        "source_format": adata.uns.get("source_format"),
        "source_file": adata.uns.get("source_file"),
        "matrix_key": adata.uns.get("matrix_key"),
        "label_key": adata.uns.get("label_key"),
        "resolved_label_key": adata.uns.get("resolved_label_key"),
        "label_is_integer": adata.uns.get("label_is_integer", False),
        "conversion_label_key": args.label_key,
        "conversion_matrix_key": args.matrix_key or "auto",
        "conversion_n_clusters": args.n_clusters,
    }
    meta_path = output_path.replace(".h5ad", ".meta.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=2)
    print(f"Metadata: {meta_path}")


if __name__ == "__main__":
    main()
