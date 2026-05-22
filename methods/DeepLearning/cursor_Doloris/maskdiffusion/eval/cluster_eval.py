# -*- coding: utf-8 -*-
"""
Clustering evaluation metrics.

Provides comprehensive clustering evaluation including:
- Standard metrics: ACC, NMI, ARI, F1-macro
- Advanced metrics: Silhouette Score
- Rare cell preservation analysis
"""

import numpy as np
import pandas as pd
from sklearn.metrics import (
    normalized_mutual_info_score,
    adjusted_rand_score,
    accuracy_score,
    fowlkes_mallows_score,
    v_measure_score,
    homogeneity_score,
    completeness_score,
    silhouette_score,
)
from sklearn.metrics.cluster import fowlkes_mallows_score as fmi


def compute_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    """
    Compute best-case accuracy by matching clusters to labels optimally.

    Uses Hungarian algorithm for optimal assignment.
    """
    from scipy.optimize import linear_sum_assignment

    # Get unique labels and predictions
    true_labels = np.unique(y_true)
    pred_labels = np.unique(y_pred)

    # Build cost matrix
    n_true = len(true_labels)
    n_pred = len(pred_labels)
    cost_matrix = np.zeros((n_true, n_pred))

    for i, tl in enumerate(true_labels):
        for j, pl in enumerate(pred_labels):
            # Count matches
            mask_true = y_true == tl
            mask_pred = y_pred == pl
            cost_matrix[i, j] = -np.sum(mask_true & mask_pred)  # Negative for max

    # Find optimal assignment
    row_ind, col_ind = linear_sum_assignment(cost_matrix)

    # Compute accuracy
    y_pred_mapped = np.zeros_like(y_pred)
    for i, j in zip(row_ind, col_ind):
        y_pred_mapped[y_pred == pred_labels[j]] = true_labels[i]

    acc = np.mean(y_pred_mapped == y_true)
    return acc


def evaluate_clustering(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    X: np.ndarray = None,
    verbose: bool = True,
) -> dict:
    """
    Compute comprehensive clustering evaluation metrics.

    Args:
        y_true: Ground truth labels
        y_pred: Predicted cluster labels
        X: Optional feature matrix for silhouette score
        verbose: Whether to print results

    Returns:
        dict with metrics
    """
    metrics = {}

    # Standard clustering metrics
    metrics['nmi'] = normalized_mutual_info_score(y_true, y_pred)
    metrics['ari'] = adjusted_rand_score(y_true, y_pred)
    metrics['v_measure'] = v_measure_score(y_true, y_pred)
    metrics['homogeneity'] = homogeneity_score(y_true, y_pred)
    metrics['completeness'] = completeness_score(y_true, y_pred)
    metrics['fmi'] = fmi(y_true, y_pred)

    # Best-case accuracy
    metrics['acc'] = compute_accuracy(y_true, y_pred)

    # F1-macro
    from sklearn.metrics import f1_score
    metrics['f1_macro'] = f1_score(y_true, y_pred, average='macro', zero_division=0)

    # Silhouette score (if features provided)
    if X is not None:
        try:
            metrics['silhouette'] = silhouette_score(X, y_pred)
        except Exception:
            metrics['silhouette'] = None
    else:
        metrics['silhouette'] = None

    if verbose:
        print('=' * 50)
        print('Clustering Evaluation Results')
        print('=' * 50)
        print(f"NMI:          {metrics['nmi']:.4f}")
        print(f"ARI:          {metrics['ari']:.4f}")
        print(f"V-Measure:    {metrics['v_measure']:.4f}")
        print(f"Homogeneity:  {metrics['homogeneity']:.4f}")
        print(f"Completeness: {metrics['completeness']:.4f}")
        print(f"FMI:          {metrics['fmi']:.4f}")
        print(f"ACC:          {metrics['acc']:.4f}")
        print(f"F1-macro:     {metrics['f1_macro']:.4f}")
        if metrics['silhouette'] is not None:
            print(f"Silhouette:   {metrics['silhouette']:.4f}")
        print('=' * 50)

    return metrics


def evaluate_rare_cell_preservation(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    min_cluster_size: int = 10,
    verbose: bool = True,
) -> dict:
    """
    Evaluate preservation of rare cell types.

    Args:
        y_true: Ground truth labels
        y_pred: Predicted cluster labels
        min_cluster_size: Minimum size to be considered non-rare
        verbose: Whether to print results

    Returns:
        dict with rare cell preservation metrics
    """
    metrics = {}

    # Find rare cell types (small populations in ground truth)
    unique_true, counts_true = np.unique(y_true, return_counts=True)
    rare_types = unique_true[counts_true < min_cluster_size]

    if len(rare_types) == 0:
        if verbose:
            print('No rare cell types found (threshold: {})'.format(min_cluster_size))
        return {'rare_preservation_rate': 1.0, 'n_rare_types': 0}

    # For each rare type, check if it's preserved in predictions
    preserved = 0
    for rt in rare_types:
        mask = y_true == rt
        pred_labels_rt = y_pred[mask]
        unique_pred, counts_pred = np.unique(pred_labels_rt, return_counts=True)

        # A rare type is preserved if:
        # 1. Most cells go to one cluster
        # 2. That cluster is not dominated by other types
        if len(unique_pred) > 0:
            dominant_cluster = unique_pred[np.argmax(counts_pred)]
            n_in_cluster = np.sum(y_pred == dominant_cluster)
            n_from_rare = np.sum(pred_labels_rt == dominant_cluster)

            # Check purity
            purity = n_from_rare / n_in_cluster
            if purity > 0.5 and n_from_rare >= min(min_cluster_size // 2, 3):
                preserved += 1

    metrics['n_rare_types'] = len(rare_types)
    metrics['rare_preserved'] = preserved
    metrics['rare_preservation_rate'] = preserved / len(rare_types)

    # Per-rare-type details
    metrics['rare_types'] = rare_types.tolist()
    metrics['rare_sizes'] = counts_true[counts_true < min_cluster_size].tolist()

    if verbose:
        print('=' * 50)
        print('Rare Cell Preservation Analysis')
        print('=' * 50)
        print(f"Rare types (< {min_cluster_size} cells): {len(rare_types)}")
        print(f"Preserved: {preserved}/{len(rare_types)}")
        print(f"Preservation rate: {metrics['rare_preservation_rate']:.4f}")
        print('=' * 50)

    return metrics


def compare_methods(
    results: dict,
    verbose: bool = True,
) -> pd.DataFrame:
    """
    Compare multiple clustering methods.

    Args:
        results: dict mapping method name to (y_true, y_pred) tuples
        verbose: Whether to print results

    Returns:
        DataFrame with comparison results
    """
    import pandas as pd

    rows = []
    for method, (y_true, y_pred) in results.items():
        metrics = evaluate_clustering(y_true, y_pred, verbose=False)
        row = {'method': method}
        row.update(metrics)
        rows.append(row)

    df = pd.DataFrame(rows)

    if verbose:
        print('=' * 80)
        print('Method Comparison')
        print('=' * 80)
        print(df.to_string(index=False))
        print('=' * 80)

    return df


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    """
    Search Leiden resolution to achieve target cluster count.

    Args:
        adata: AnnData object with embedding in obsm
        fixed_clus_count: Target number of clusters
        increment: Resolution increment step

    Returns:
        Optimal resolution parameter
    """
    import scanpy as sc
    import pandas as pd

    dis = []
    resolutions = sorted(list(np.arange(0.01, 2.5, increment)), reverse=True)

    for res in resolutions:
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        dis.append(abs(count_unique_leiden - fixed_clus_count))
        if count_unique_leiden == fixed_clus_count:
            break

    return resolutions[np.argmin(dis)]


def run_all_evaluations(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    X: np.ndarray = None,
    min_rare_size: int = 10,
    verbose: bool = True,
) -> dict:
    """
    Run all clustering evaluations.

    Args:
        y_true: Ground truth labels
        y_pred: Predicted cluster labels
        X: Optional feature matrix
        min_rare_size: Threshold for rare cell identification
        verbose: Whether to print results

    Returns:
        dict with all metrics
    """
    results = {}

    # Standard metrics
    results['clustering'] = evaluate_clustering(y_true, y_pred, X, verbose=verbose)

    # Rare cell preservation
    results['rare_cells'] = evaluate_rare_cell_preservation(
        y_true, y_pred, min_cluster_size=min_rare_size, verbose=verbose
    )

    return results
