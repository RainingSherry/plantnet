#!/usr/bin/env python3
"""
Compute metrics.json for models that have embeddings/predictions but missing metrics.
Handles:
  1. cursor_Doloris_maskdiffusion (has results.json)
  2. cursor_Doloris_GraphDiffusion (has config.json + loss_history.json, embeddings need regeneration)
  3. Doloris_DiffusionBridge (has bridge_embedding.npy + bridge_kmeans_labels.npy)
  4. cursor_Doloris_DiffusionBridge (has bridge_embedding.npy + bridge_kmeans_labels.npy)
"""

import os
import sys
import json
import h5py
import numpy as np
from pathlib import Path
from sklearn.metrics import (
    normalized_mutual_info_score as nmi_score,
    adjusted_rand_score as ari_score,
    fowlkes_mallows_score as fmi_score,
    v_measure_score as vms_score,
    homogeneity_score as hom_score,
    completeness_score as com_score,
    f1_score,
    accuracy_score,
)
from scipy.optimize import linear_sum_assignment


def get_label_column(h5ad_path):
    """Get ground truth label column from h5ad file."""
    with h5py.File(h5ad_path, 'r') as f:
        obs_keys = list(f['obs'].keys())
        if 'Celltype' in obs_keys:
            labels_raw = f['obs/Celltype'][:]
            if labels_raw.dtype.kind == 'O':
                return np.array([s.decode() if isinstance(s, bytes) else s for s in labels_raw])
            return np.array(labels_raw).astype(str)
        elif 'cell_type' in obs_keys:
            codes = f['obs/cell_type/codes'][:]
            cats = f['obs/cell_type/categories'][:]
            cats_str = np.array([s.decode() if isinstance(s, bytes) else s for s in cats])
            return cats_str[codes]
    return None


def encode_labels(labels):
    """Encode string labels to integers."""
    unique = np.unique(labels)
    label_map = {v: i for i, v in enumerate(unique)}
    return np.array([label_map[v] for v in labels]), len(unique)


def compute_metrics(y_true, y_pred):
    """Compute all clustering metrics using Hungarian label matching."""
    y_true_enc, _ = encode_labels(y_true)
    y_pred_enc = y_pred.astype(int)

    # Build confusion matrix for Hungarian matching
    true_unique = np.unique(y_true_enc)
    pred_unique = np.unique(y_pred_enc)
    n_class = len(true_unique)
    n_pred = len(pred_unique)

    G = np.zeros((n_class, n_pred), dtype=int)
    for i, ut in enumerate(true_unique):
        for j, up in enumerate(pred_unique):
            G[i, j] = np.sum((y_true_enc == ut) & (y_pred_enc == up))

    # Hungarian assignment
    A = linear_sum_assignment(-G)
    new_pred = np.zeros_like(y_pred_enc)
    for i, up in enumerate(pred_unique):
        col_idx = A[1][i] if i < len(A[1]) else i % n_class
        label_idx = A[0][i] if i < len(A[0]) else i % n_class
        new_pred[y_pred_enc == up] = true_unique[label_idx]

    return {
        "acc": float(accuracy_score(y_true_enc, new_pred)),
        "nmi": float(nmi_score(y_true_enc, y_pred_enc, average_method="arithmetic")),
        "ari": float(ari_score(y_true_enc, y_pred_enc)),
        "f1_macro": float(f1_score(y_true_enc, new_pred, average='macro', zero_division=0)),
        "fmi": float(fmi_score(y_true_enc, y_pred_enc)),
        "v_measure": float(vms_score(y_true_enc, y_pred_enc)),
        "homogeneity": float(hom_score(y_true_enc, y_pred_enc)),
        "completeness": float(com_score(y_true_enc, y_pred_enc)),
    }


def compute_maskdiffusion_metrics(results_path, data_path, output_path):
    """Handle cursor_Doloris_maskdiffusion - has results.json with embeddings/predictions."""
    print(f"\n=== cursor_Doloris_maskdiffusion @ {results_path} ===")
    results_path = Path(results_path)
    data_path = Path(data_path)
    output_path = Path(output_path) / "metrics.json"

    if output_path.exists():
        print(f"  metrics.json already exists, skipping")
        return

    # Load results
    with open(results_path / "results.json") as f:
        results = json.load(f)

    # Get predictions
    if (results_path / "pred_labels.npy").exists():
        y_pred = np.load(results_path / "pred_labels.npy")
    elif "pred_labels" in results:
        y_pred = np.array(results["pred_labels"])
    else:
        print(f"  ERROR: No pred_labels found")
        return

    # Get true labels
    y_true = get_label_column(data_path)
    if y_true is None:
        print(f"  ERROR: Could not find label column")
        return

    # Ensure matching length
    if len(y_true) != len(y_pred):
        print(f"  WARNING: Length mismatch {len(y_true)} vs {len(y_pred)}, truncating to min")
        n = min(len(y_true), len(y_pred))
        y_true = y_true[:n]
        y_pred = y_pred[:n]

    metrics = compute_metrics(y_true, y_pred)

    # Add metadata
    metrics["source"] = "cursor_Doloris_maskdiffusion"
    if "best_nmi" in results:
        metrics["best_nmi_from_training"] = results["best_nmi"]
    if "n_clusters" in results:
        metrics["n_clusters"] = results["n_clusters"]

    with open(output_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"  Saved: {output_path}")
    print(f"  NMI={metrics['nmi']:.4f} ARI={metrics['ari']:.4f} ACC={metrics['acc']:.4f}")


def compute_graphdiffusion_metrics(results_path, data_path, output_path):
    """Handle cursor_Doloris_GraphDiffusion - has config + loss_history, embeddings missing."""
    print(f"\n=== cursor_Doloris_GraphDiffusion @ {results_path} ===")
    results_path = Path(results_path)
    data_path = Path(data_path)
    output_path = Path(output_path) / "metrics.json"

    if output_path.exists():
        print(f"  metrics.json already exists, skipping")
        return

    # Load config
    with open(results_path / "config.json") as f:
        config = json.load(f)

    n_clusters = config.get("n_clusters", 15)

    # Get true labels
    y_true = get_label_column(data_path)
    if y_true is None:
        print(f"  ERROR: Could not find label column")
        return

    # Try to get predictions from kmeans_n_clust files
    # The eval script saves predictions as "best_predictions.npy" or similar
    # Let's check for existing prediction files
    kmeans_files = sorted(results_path.glob("kmeans_*.npy"))
    pred_files = sorted(results_path.glob("*pred*.npy"))

    if pred_files:
        # Use best predictions found
        for pf in pred_files:
            if "best" in pf.name.lower() or "kmeans" in pf.name.lower():
                y_pred = np.load(pf)
                break
        else:
            y_pred = np.load(pred_files[0])
    elif kmeans_files:
        # Use kmeans with n_clusters
        target_file = results_path / f"kmeans_{n_clusters}.npy"
        if target_file.exists():
            y_pred = np.load(target_file)
        else:
            y_pred = np.load(kmeans_files[0])
    else:
        # Need to regenerate from embeddings - check if embeddings exist
        emb_files = sorted(results_path.glob("embedding*.npy"))
        if not emb_files:
            print(f"  ERROR: No embeddings or predictions found. Need to re-run training.")
            print(f"  Available files: {list(results_path.glob('*'))}")
            return
        print(f"  Embeddings found: {len(emb_files)} files. Using last one.")
        emb = np.load(emb_files[-1])
        from sklearn.cluster import KMeans
        km = KMeans(n_clusters=n_clusters, n_init=20, random_state=42)
        y_pred = km.fit_predict(emb)

    # Ensure matching length
    if len(y_true) != len(y_pred):
        print(f"  WARNING: Length mismatch {len(y_true)} vs {len(y_pred)}, truncating to min")
        n = min(len(y_true), len(y_pred))
        y_true = y_true[:n]
        y_pred = y_pred[:n]

    metrics = compute_metrics(y_true, y_pred)

    # Add metadata
    metrics["source"] = "cursor_Doloris_GraphDiffusion"
    if (results_path / "loss_history.json").exists():
        with open(results_path / "loss_history.json") as f:
            lh = json.load(f)
        metrics["best_epoch"] = lh.get("best_epoch", "?")
        metrics["best_nmi_from_training"] = lh.get("best_nmi", "?")

    with open(output_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"  Saved: {output_path}")
    print(f"  NMI={metrics['nmi']:.4f} ARI={metrics['ari']:.4f} ACC={metrics['acc']:.4f}")


def compute_diffusionbridge_metrics(results_path, data_path, output_path):
    """Handle Doloris_DiffusionBridge and cursor_Doloris_DiffusionBridge."""
    print(f"\n=== DiffusionBridge @ {results_path} ===")
    results_path = Path(results_path)
    data_path = Path(data_path)
    output_path = Path(output_path) / "metrics.json"

    if output_path.exists():
        print(f"  metrics.json already exists, skipping")
        return

    # Get true labels
    y_true = get_label_column(data_path)
    if y_true is None:
        print(f"  ERROR: Could not find label column")
        return

    # Get kmeans predictions (or dec labels)
    kmeans_file = results_path / "bridge_kmeans_labels.npy"
    dec_file = results_path / "bridge_dec_labels.npy"

    if kmeans_file.exists():
        y_pred = np.load(kmeans_file)
        source = "kmeans"
    elif dec_file.exists():
        y_pred = np.load(dec_file)
        source = "dec_labels"
    else:
        print(f"  ERROR: No predictions found (bridge_kmeans_labels.npy or bridge_dec_labels.npy)")
        return

    # Ensure matching length
    if len(y_true) != len(y_pred):
        print(f"  WARNING: Length mismatch {len(y_true)} vs {len(y_pred)}, truncating to min")
        n = min(len(y_true), len(y_pred))
        y_true = y_true[:n]
        y_pred = y_pred[:n]

    metrics = compute_metrics(y_true, y_pred)
    metrics["source"] = f"DiffusionBridge_{source}"

    with open(output_path, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"  Saved: {output_path}")
    print(f"  NMI={metrics['nmi']:.4f} ARI={metrics['ari']:.4f} ACC={metrics['acc']:.4f}")


def main():
    base = Path("/home/luolie/biopipeline/dimension-reduction/plantnet")
    data_dir = base / "data"
    results_dir = base / "results"

    # Dataset configurations
    datasets = {
        "SRP182008": {
            "h5ad": data_dir / "SRP182008.h5ad",
            "label_col": "Celltype",
            "n_cells": 13514,
            "n_cell_types": 15,
        },
        "Mouse_Pancreas_1": {
            "h5ad": data_dir / "Mouse_Pancreas_1.h5ad",
            "label_col": "cell_type",
            "n_cells": 1886,
            "n_cell_types": 13,
        },
    }

    # Model configurations: (results_subdir, compute_func, extra_args)
    models = [
        # cursor_Doloris_maskdiffusion
        ("cursor_Doloris_maskdiffusion/SRP182008", "maskdiffusion", "cursor_Doloris_maskdiffusion"),
        ("cursor_Doloris_maskdiffusion/Mouse_Pancreas_1", "maskdiffusion", "cursor_Doloris_maskdiffusion"),
        # cursor_Doloris_GraphDiffusion
        ("cursor_Doloris_GraphDiffusion/SRP182008", "graphdiffusion", "cursor_Doloris_GraphDiffusion"),
        ("cursor_Doloris_GraphDiffusion/Mouse_Pancreas_1", "graphdiffusion", "cursor_Doloris_GraphDiffusion"),
        # Doloris_DiffusionBridge
        ("Doloris_DiffusionBridge/SRP182008", "diffusionbridge", "Doloris_DiffusionBridge"),
        ("Doloris_DiffusionBridge/Mouse_Pancreas_1", "diffusionbridge", "Doloris_DiffusionBridge"),
        # cursor_Doloris_DiffusionBridge
        ("cursor_Doloris_DiffusionBridge/SRP182008", "diffusionbridge", "cursor_Doloris_DiffusionBridge"),
        ("cursor_Doloris_DiffusionBridge/Mouse_Pancreas_1", "diffusionbridge", "cursor_Doloris_DiffusionBridge"),
    ]

    for model_path, model_type, model_name in models:
        result_path = results_dir / model_path
        dataset_name = model_path.split("/")[1]
        dataset_info = datasets.get(dataset_name)
        if dataset_info is None:
            print(f"\nUnknown dataset: {dataset_name}")
            continue

        if not result_path.exists():
            print(f"\nSkipping (not found): {result_path}")
            continue

        metrics_file = result_path / "metrics.json"
        if metrics_file.exists():
            print(f"\nSkipping (metrics.json exists): {model_name}/{dataset_name}")
            continue

        h5ad_path = dataset_info["h5ad"]
        if not h5ad_path.exists():
            print(f"\nERROR: Data file not found: {h5ad_path}")
            continue

        if model_type == "maskdiffusion":
            compute_maskdiffusion_metrics(result_path, h5ad_path, result_path)
        elif model_type == "graphdiffusion":
            compute_graphdiffusion_metrics(result_path, h5ad_path, result_path)
        elif model_type == "diffusionbridge":
            compute_diffusionbridge_metrics(result_path, h5ad_path, result_path)


if __name__ == "__main__":
    os.environ["TMPDIR"] = "/data/tmp"
    main()
