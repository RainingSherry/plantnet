#!/usr/bin/env python3
"""
inspect_h5_dataset.py
====================
Inspect the schema of a raw .h5 file used by scMAE-style pipelines.

Usage:
    python scripts/inspect_h5_dataset.py --input_path data/scMAE/Pollen.h5
    python scripts/inspect_h5_dataset.py --input_path data/scMAE/Quake_Smart-seq2_Lung.h5 --verbose
"""

import argparse
import h5py
import numpy as np


CANDIDATE_MATRIX_KEYS = ["X", "data", "exprs", "matrix", "counts", "raw", "value"]
CANDIDATE_LABEL_KEYS = ["Y", "label", "labels", "cell_type", "celltype", "cell_type_label",
                          "cell_label", "cluster", "clusters", "type", "group", "group_id"]
CANDIDATE_CELL_KEYS = ["cell_names", "obs_names", "barcode", "barcodes", "cell_id", "cells"]
CANDIDATE_GENE_KEYS = ["gene_names", "var_names", "gene", "genes", "gene_id"]


def _print_keys(prefix, obj, verbose, indent=0):
    """Recursively print all keys with shapes and dtypes."""
    if isinstance(obj, h5py.Dataset):
        shape = obj.shape
        dtype = obj.dtype
        try:
            if shape[0] < 5:
                sample = obj[...].tolist()
            else:
                sample = obj[:3].tolist()
        except Exception:
            sample = None
        line = f"{'  ' * indent}{prefix}: shape={shape}, dtype={dtype}"
        if sample is not None and verbose:
            line += f", sample={str(sample)[:80]}"
        print(line)
    elif isinstance(obj, h5py.Group):
        print(f"{'  ' * indent}{prefix}/ (group, keys={list(obj.keys())})")
        for k in sorted(obj.keys()):
            _print_keys(k, obj[k], verbose, indent + 1)


def auto_detect_matrix(f):
    """Find the expression matrix dataset."""
    # 1. Check known candidates directly at root
    for key in CANDIDATE_MATRIX_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset):
                return key, ds.shape
            # Group containing a sparse matrix
            if isinstance(ds, h5py.Group):
                if "data" in ds and "indices" in ds:
                    return key, ds["shape"][...]
                if "data" in ds and "indptr" in ds:
                    return key, ds["shape"][...]
                return key, None
    # 2. Search all datasets and find the largest 2-D array
    candidates = []
    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset) and len(obj.shape) >= 2:
            candidates.append((name, obj.shape, obj.dtype))
    f.visititems(visitor)
    candidates.sort(key=lambda x: -np.prod(x[1]) if x[1] else 0)
    if candidates:
        return candidates[0][0], candidates[0][1]
    return None, None


def auto_detect_labels(f):
    """Find the label array (1-D, length = n_cells)."""
    # 1. Check known candidates
    for key in CANDIDATE_LABEL_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset) and len(ds.shape) == 1:
                return key, ds.shape
    # 2. Search all 1-D arrays
    candidates = []
    def visitor(name, obj):
        if isinstance(obj, h5py.Dataset) and len(obj.shape) == 1:
            candidates.append((name, obj.shape, obj.dtype))
    f.visititems(visitor)
    # Prefer 1-D arrays whose length matches a detected matrix row count
    candidates.sort(key=lambda x: x[1][0], reverse=True)
    if candidates:
        return candidates[0][0], candidates[0][1]
    return None, None


def auto_detect_names(f, n_cells, n_genes):
    """Find cell and gene name arrays."""
    cell_key, gene_key = None, None
    for key in CANDIDATE_CELL_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset) and ds.shape[0] == n_cells:
                cell_key = key
                break
    for key in CANDIDATE_GENE_KEYS:
        if key in f:
            ds = f[key]
            if isinstance(ds, h5py.Dataset) and ds.shape[0] == n_genes:
                gene_key = key
                break
    return cell_key, gene_key


def infer_n_clusters(label_ds, h5_file):
    """Read labels and count unique values."""
    try:
        labels = label_ds[...]
        if hasattr(labels, 'tolist'):
            labels = labels.tolist()
        if hasattr(labels, 'astype'):
            labels = labels.astype(str)
        else:
            labels = np.array(labels).astype(str)
        unique = np.unique(labels)
        return len(unique), list(unique[:10])
    except Exception as e:
        return None, str(e)


def main():
    parser = argparse.ArgumentParser(
        description="Inspect .h5 dataset schema",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/inspect_h5_dataset.py --input_path data/scMAE/Pollen.h5
  python scripts/inspect_h5_dataset.py --input_path data/scMAE/Quake_Smart-seq2_Lung.h5 --verbose
""",
    )
    parser.add_argument("--input_path", type=str, required=True,
                        help="Path to .h5 file")
    parser.add_argument("--verbose", "-v", action="store_true",
                        help="Show sample values for each dataset")
    args = parser.parse_args()

    import os
    if not os.path.exists(args.input_path):
        print(f"ERROR: File not found: {args.input_path}")
        return

    print(f"{'=' * 60}")
    print(f"File: {args.input_path}")
    print(f"{'=' * 60}")

    with h5py.File(args.input_path, "r") as f:
        # 1. Show all keys
        print(f"\nTop-level keys: {sorted(f.keys())}")
        print("\n--- Full schema ---")
        for k in sorted(f.keys()):
            _print_keys(k, f[k], args.verbose)

        # 2. Auto-detect matrix
        matrix_key, matrix_shape = auto_detect_matrix(f)
        print(f"\n--- Auto-detection ---")
        print(f"Detected matrix key: {matrix_key!r}  shape={matrix_shape}")

        # 3. Auto-detect labels
        label_key, label_shape = auto_detect_labels(f)
        print(f"Detected label key:  {label_key!r}  shape={label_shape}")

        # 4. Auto-detect names
        n_cells = matrix_shape[0] if matrix_shape else None
        n_genes = matrix_shape[1] if matrix_shape and len(matrix_shape) == 2 else None
        if n_cells and n_genes:
            cell_key, gene_key = auto_detect_names(f, n_cells, n_genes)
            print(f"Detected cell names key: {cell_key!r}")
            print(f"Detected gene names key: {gene_key!r}")
        else:
            cell_key, gene_key = None, None

        # 5. Infer n_clusters
        if label_key:
            n_clusters, sample_labels = infer_n_clusters(f[label_key], args.input_path)
            print(f"Inferred n_clusters: {n_clusters}")
            if sample_labels:
                print(f"  Sample label values: {sample_labels}")

        # 6. Summary
        print(f"\n--- Summary ---")
        print(f"  Matrix:   {matrix_key}  (shape={matrix_shape})")
        print(f"  Labels:   {label_key}  (shape={label_shape})")
        print(f"  Cells:    {cell_key}")
        print(f"  Genes:    {gene_key}")
        if n_clusters:
            print(f"  Clusters: {n_clusters}")
        print()


if __name__ == "__main__":
    main()
