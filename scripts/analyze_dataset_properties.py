#!/usr/bin/env python3
"""
analyze_dataset_properties.py
==============================
Compute structural properties of each dataset that may explain why NeighborMix
succeeds on some datasets but not others.

Computed metrics (all based on raw input before any preprocessing):

  1. Basic stats
     - n_cells, n_genes, n_clusters
     - sparsity = zero_fraction

  2. Gene detection
     - median_detected_genes, mean_detected_genes (genes with > 0 count per cell)

  3. Label distribution
     - class_imbalance_ratio = largest_class / smallest_class
     - label_entropy = H(labels) / H_max  (normalized to [0, 1])

  4. KNN graph purity (neighborhood reliability)
     - knn_purity_k{K} = mean fraction of k-NN that share the same label
       Computed on PCA(50) of raw expression matrix.

  5. Cross-type edge ratio
     - cross_type_edge_ratio_k{K} = fraction of k-NN edges that connect
       cells of different types. High ratio → unreliable neighborhood.

  6. Silhouette score (global separation)
     - silhouette_pca = silhouette score on PCA(50) embeddings
     - silhouette_hvg_pca = silhouette score on PCA(50) of HVG-subset

Usage:
    python scripts/analyze_dataset_properties.py \
        --data_paths \
            data/processed/Quake_Smart-seq2_Lung.h5ad \
            data/processed/Wang.h5ad \
            data/processed/Pollen.h5ad \
            data/hrvatin_geo_maintype_counts.h5ad \
            data/SRP182008.h5ad \
            data/SRP171040.h5ad \
            data/SRP235541.h5ad \
        --label_key auto \
        --out_dir results/analysis
"""

import argparse
import csv
import json
import os
import sys
import warnings
from pathlib import Path
from typing import Optional

import numpy as np
import scanpy as sc
from scipy.sparse import issparse
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)


# ─── Metric computation ────────────────────────────────────────────────────────

def _get_label_col(adata, label_key) -> str:
    """Resolve label column, trying candidates if label_key='auto'."""
    if label_key != "auto":
        return label_key
    candidates = [
        "resolved_label", "cell_type", "Celltype", "celltype",
        "cell_label", "label", "Group", "cluster",
    ]
    for c in candidates:
        if c in adata.obs.columns:
            return c
    raise ValueError(
        f"Cannot auto-detect label column. "
        f"Available: {list(adata.obs.columns)}. "
        f"Use --label_key to specify."
    )


def compute_zero_fraction(X) -> float:
    """Fraction of zero entries in the expression matrix."""
    if issparse(X):
        X = X.toarray()
    return float((X == 0).sum()) / X.size


def compute_gene_detection(X):
    """Per-cell and summary gene detection stats."""
    if issparse(X):
        X = X.toarray()
    genes_per_cell = (X > 0).sum(axis=1)
    return {
        "median_detected_genes": float(np.median(genes_per_cell)),
        "mean_detected_genes": float(np.mean(genes_per_cell)),
    }


def compute_class_imbalance_ratio(labels) -> float:
    """Ratio of largest class size to smallest class size."""
    _, counts = np.unique(labels, return_counts=True)
    return float(counts.max() / counts.min()) if len(counts) > 1 else 1.0


def compute_label_entropy(labels) -> float:
    """Normalized Shannon entropy of label distribution. 0 = pure, 1 = uniform."""
    _, counts = np.unique(labels, return_counts=True)
    probs = counts / counts.sum()
    H = -np.sum(probs * np.log(probs + 1e-12))
    H_max = np.log(len(counts) + 1e-12)
    return float(H / H_max) if H_max > 0 else 0.0


def compute_knn_metrics(X, labels, ks=(5, 10, 20), n_pcs=50, random_state=42):
    """
    Compute KNN purity and cross-type edge ratio at multiple k values.
    Uses PCA(50) of raw X for distance computation.
    """
    # PCA on raw expression
    if issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = X
    X_dense = X_dense.astype(np.float32)

    # Subsample if too large (for efficiency)
    max_for_pca = 10000
    if X_dense.shape[0] > max_for_pca:
        idx = np.random.RandomState(random_state).choice(
            X_dense.shape[0], max_for_pca, replace=False
        )
        X_sub = X_dense[idx]
        labels_sub = labels[idx]
    else:
        X_sub = X_dense
        labels_sub = labels

    # PCA
    n_comp = min(n_pcs, X_sub.shape[0] - 1, X_sub.shape[1])
    pca = PCA(n_components=n_comp, random_state=random_state)
    X_pca = pca.fit_transform(X_sub)

    results = {}
    for k in ks:
        nn = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=1)
        nn.fit(X_pca)
        dists, indices = nn.kneighbors(X_pca)

        # Exclude self (index 0)
        neighbor_labels = labels_sub[indices[:, 1:]]
        same_label = (neighbor_labels == labels_sub[:, np.newaxis])
        purity = same_label.mean(axis=1).mean()
        cross_ratio = 1.0 - purity

        results[f"knn_purity_k{k}"] = round(float(purity), 4)
        results[f"cross_type_edge_ratio_k{k}"] = round(float(cross_ratio), 4)

    return results


def compute_silhouette(X, labels, n_pcs=50, random_state=42):
    """
    Compute silhouette scores on PCA embeddings (full matrix and HVG-subset).
    """
    if issparse(X):
        X_dense = X.toarray()
    else:
        X_dense = X
    X_dense = X_dense.astype(np.float32)

    # Encode labels
    le = LabelEncoder()
    y = le.fit_transform(labels)

    # Subsample for silhouette if too large
    max_for_sil = 5000
    if X_dense.shape[0] > max_for_sil:
        idx = np.random.RandomState(random_state).choice(
            X_dense.shape[0], max_for_sil, replace=False
        )
        X_sub = X_dense[idx]
        y_sub = y[idx]
    else:
        X_sub = X_dense
        y_sub = y

    n_comp = min(n_pcs, X_sub.shape[0] - 2, X_sub.shape[1])
    pca = PCA(n_components=n_comp, random_state=random_state)
    X_pca = pca.fit_transform(X_sub)

    try:
        sil_full = silhouette_score(X_pca, y_sub)
    except Exception:
        sil_full = np.nan

    # HVG-based silhouette
    try:
        adata_tmp = sc.AnnData(X=X_sub)
        sc.pp.highly_variable_genes(
            adata_tmp, n_top_genes=min(1000, X_sub.shape[1]),
            flavor="seurat", inplace=True
        )
        X_hvg = X_sub[:, adata_tmp.var["highly_variable"]]
        pca_hvg = PCA(n_components=min(50, X_hvg.shape[1], X_hvg.shape[0] - 2),
                      random_state=random_state)
        X_hvg_pca = pca_hvg.fit_transform(X_hvg)
        sil_hvg = silhouette_score(X_hvg_pca, y_sub)
    except Exception:
        sil_hvg = np.nan

    return {
        "silhouette_pca": round(float(sil_full), 4) if not np.isnan(sil_full) else None,
        "silhouette_hvg_pca": round(float(sil_hvg), 4) if not np.isnan(sil_hvg) else None,
    }


def compute_properties(adata, label_key="auto", random_state=42) -> dict:
    """Compute all properties for a single dataset."""
    label_col = _get_label_col(adata, label_key)
    labels = adata.obs[label_col].astype(str).values
    X = adata.X

    props = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_clusters": int(len(np.unique(labels))),
        "zero_fraction": round(compute_zero_fraction(X), 4),
    }

    # Gene detection
    detection = compute_gene_detection(X)
    props.update(detection)

    # Class imbalance
    props["class_imbalance_ratio"] = round(compute_class_imbalance_ratio(labels), 4)
    props["label_entropy"] = round(compute_label_entropy(labels), 4)

    # KNN metrics
    knn = compute_knn_metrics(X, labels, ks=(5, 10, 20), random_state=random_state)
    props.update(knn)

    # Silhouette
    sil = compute_silhouette(X, labels, random_state=random_state)
    props.update(sil)

    props["label_key"] = label_col
    return props


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Analyze structural properties of datasets",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/analyze_dataset_properties.py \\
      --data_paths data/processed/Quake_Smart-seq2_Lung.h5ad \\
                  data/processed/Wang.h5ad \\
                  data/processed/Pollen.h5ad \\
      --label_key auto \\
      --out_dir results/analysis
""",
    )
    parser.add_argument("--data_paths", type=str, nargs="+", required=True,
                        help="Paths to .h5ad files")
    parser.add_argument("--label_key", type=str, default="auto",
                        help="Label column key (default: auto)")
    parser.add_argument("--out_dir", type=str, default="results/analysis",
                        help="Output directory")
    parser.add_argument("--seeds", type=int, nargs="+", default=[42, 2024, 3407],
                        help="Seeds for KNN subsampling (default: 42 2024 3407)")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    all_rows = []
    out_csv = Path(args.out_dir) / "dataset_property_summary.csv"
    out_json = Path(args.out_dir) / "dataset_property_summary.json"

    # Column order for the CSV
    FIELD_ORDER = [
        "dataset",
        "n_cells", "n_genes", "n_clusters",
        "zero_fraction",
        "median_detected_genes", "mean_detected_genes",
        "class_imbalance_ratio", "label_entropy",
        "knn_purity_k5", "knn_purity_k10", "knn_purity_k20",
        "cross_type_edge_ratio_k5", "cross_type_edge_ratio_k10",
        "cross_type_edge_ratio_k20",
        "silhouette_pca", "silhouette_hvg_pca",
        "label_key",
    ]

    for data_path in args.data_paths:
        if not os.path.exists(data_path):
            print(f"WARNING: File not found: {data_path} — skipping")
            continue

        dataset_name = os.path.splitext(os.path.basename(data_path))[0]
        print(f"\n[{dataset_name}] Loading {data_path} ...")

        adata = sc.read_h5ad(data_path)

        # Compute properties using the first seed; report consistency note
        seed = args.seeds[0]
        np.random.seed(seed)
        props = compute_properties(adata, label_key=args.label_key, random_state=seed)
        props["dataset"] = dataset_name

        print(f"  n_cells={props['n_cells']}, n_genes={props['n_genes']}, "
              f"n_clusters={props['n_clusters']}")
        print(f"  zero_fraction={props['zero_fraction']:.4f}")
        print(f"  knn_purity_k10={props['knn_purity_k10']:.4f}")
        print(f"  cross_type_edge_ratio_k10={props['cross_type_edge_ratio_k10']:.4f}")
        print(f"  silhouette_pca={props['silhouette_pca']}")
        print(f"  class_imbalance_ratio={props['class_imbalance_ratio']:.2f}")
        print(f"  label_entropy={props['label_entropy']:.4f}")

        all_rows.append(props)

    if not all_rows:
        print("ERROR: No datasets processed successfully.")
        sys.exit(1)

    # Write CSV
    with open(out_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=FIELD_ORDER, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(all_rows)
    print(f"\nSaved: {out_csv}")

    # Write JSON (machine-readable, all fields)
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(all_rows, f, indent=2)
    print(f"Saved: {out_json}")

    # Print summary table
    print(f"\n{'Dataset':<35} {'cells':>8} {'genes':>8} {'k':>3} "
          f"{'zero%':>6} {'purity10':>8} {'cross10':>8} {'sil':>6}")
    print("-" * 95)
    for row in all_rows:
        print(
            f"{row['dataset']:<35} "
            f"{row['n_cells']:>8} "
            f"{row['n_genes']:>8} "
            f"{row['n_clusters']:>3} "
            f"{row['zero_fraction']*100:>5.1f}% "
            f"{row['knn_purity_k10']:>8.4f} "
            f"{row['cross_type_edge_ratio_k10']:>8.4f} "
            f"{str(row.get('silhouette_pca','')):>6}"
        )


if __name__ == "__main__":
    main()
