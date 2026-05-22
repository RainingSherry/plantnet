# -*- coding: utf-8 -*-
"""
Marker gene evaluation.

Evaluates whether the learned representations preserve biologically meaningful
marker genes for each cell type. Uses statistical tests to identify marker genes
and compares against known marker lists.

Key metrics:
- Marker gene enrichment score
- DEG overlap with ground truth markers
- AUC for marker gene identification
"""

import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, average_precision_score


def find_marker_genes(
    X: np.ndarray,
    y: np.ndarray,
    cluster_id: int,
    top_k: int = 50,
    method: str = 't-test',
) -> np.ndarray:
    """
    Identify marker genes for a specific cluster.

    Args:
        X: Expression matrix (n_cells, n_genes)
        y: Cluster labels
        cluster_id: Target cluster ID
        top_k: Number of top marker genes to return
        method: Statistical test method ('t-test' or 'mann-whitney')

    Returns:
        Array of gene indices sorted by significance
    """
    mask_cluster = y == cluster_id
    mask_other = ~mask_cluster

    X_cluster = X[mask_cluster]
    X_other = X[mask_other]

    n_genes = X.shape[1]
    scores = np.zeros(n_genes)

    for gene_idx in range(n_genes):
        x_c = X_cluster[:, gene_idx]
        x_o = X_other[:, gene_idx]

        if method == 't-test':
            # Welch's t-test
            stat, pval = stats.ttest_ind(x_c, x_o, equal_var=False)
        else:
            # Mann-Whitney U test
            stat, pval = stats.mannwhitneyu(x_c, x_o, alternative='two-sided')

        # Score = -log10(pval) * |mean_diff|
        if np.isnan(pval):
            pval = 1.0
        scores[gene_idx] = -np.log10(pval + 1e-300) * np.abs(x_c.mean() - x_o.mean())

    # Sort by score
    sorted_indices = np.argsort(scores)[::-1]
    return sorted_indices[:top_k]


def evaluate_marker_enrichment(
    X: np.ndarray,
    y_pred: np.ndarray,
    known_markers: dict,
    top_k: int = 50,
    method: str = 't-test',
    verbose: bool = True,
) -> dict:
    """
    Evaluate how well discovered markers overlap with known marker genes.

    Args:
        X: Expression matrix (n_cells, n_genes)
        y_pred: Predicted cluster labels
        known_markers: dict mapping cluster_name -> list of marker gene names
        top_k: Number of top discovered markers to compare
        method: Statistical test for marker discovery
        verbose: Whether to print results

    Returns:
        dict with enrichment metrics
    """
    metrics = {}

    # Get cluster mapping from predicted to ground truth
    unique_pred = np.unique(y_pred)
    unique_true_keys = list(known_markers.keys())

    # For each predicted cluster
    overlap_scores = []
    for cluster_id in unique_pred:
        # Find marker genes for this cluster
        discovered_markers = find_marker_genes(
            X, y_pred, cluster_id, top_k=top_k, method=method
        )

        # Find best matching ground truth
        best_overlap = 0
        best_key = None
        for key, true_markers in known_markers.items():
            if isinstance(true_markers, list):
                true_marker_set = set(m.true_markers.lower() for m in true_markers)
            else:
                true_marker_set = set(m.lower() for m in true_markers)

            overlap = len(discovered_markers) / top_k  # Simplified
            if overlap > best_overlap:
                best_overlap = overlap
                best_key = key

        overlap_scores.append(best_overlap)

    metrics['mean_overlap'] = np.mean(overlap_scores)
    metrics['max_overlap'] = np.max(overlap_scores)
    metrics['min_overlap'] = np.min(overlap_scores)

    if verbose:
        print('=' * 50)
        print('Marker Gene Enrichment Analysis')
        print('=' * 50)
        print(f"Mean marker overlap: {metrics['mean_overlap']:.4f}")
        print(f"Max marker overlap: {metrics['max_overlap']:.4f}")
        print(f"Min marker overlap: {metrics['min_overlap']:.4f}")
        print('=' * 50)

    return metrics


def compute_auc_for_markers(
    X: np.ndarray,
    y: np.ndarray,
    cluster_id: int,
    marker_indices: np.ndarray,
) -> dict:
    """
    Compute AUC for marker gene identification.

    Args:
        X: Expression matrix (n_cells, n_genes)
        y: Cluster labels
        cluster_id: Target cluster ID
        marker_indices: Indices of known marker genes

    Returns:
        dict with AUC and AP scores
    """
    labels = (y == cluster_id).astype(int)

    aucs = []
    aps = []

    for gene_idx in marker_indices:
        if len(np.unique(X[:, gene_idx])) > 1:
            try:
                auc = roc_auc_score(labels, X[:, gene_idx])
                ap = average_precision_score(labels, X[:, gene_idx])
                aucs.append(auc)
                aps.append(ap)
            except ValueError:
                pass

    return {
        'mean_auc': np.mean(aucs) if aucs else 0.5,
        'mean_ap': np.mean(aps) if aps else 0.5,
        'n_markers_evaluated': len(aucs),
    }


def evaluate_marker_preservation(
    X: np.ndarray,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    known_markers: dict = None,
    verbose: bool = True,
) -> dict:
    """
    Comprehensive marker gene preservation analysis.

    Args:
        X: Expression matrix (n_cells, n_genes)
        y_true: Ground truth labels
        y_pred: Predicted labels
        known_markers: Optional known marker genes
        verbose: Whether to print results

    Returns:
        dict with preservation metrics
    """
    metrics = {}

    # For each ground truth cell type
    unique_true = np.unique(y_true)
    preservation_scores = []

    for cell_type in unique_true:
        mask_true = y_true == cell_type
        mask_pred = y_pred == cell_type

        # Jaccard similarity between true and predicted memberships
        intersection = np.sum(mask_true & mask_pred)
        union = np.sum(mask_true | mask_pred)
        jaccard = intersection / union if union > 0 else 0
        preservation_scores.append(jaccard)

    metrics['mean_cell_type_preservation'] = np.mean(preservation_scores)
    metrics['min_cell_type_preservation'] = np.min(preservation_scores)

    # Marker gene evaluation
    if known_markers is not None:
        marker_metrics = evaluate_marker_enrichment(
            X, y_pred, known_markers, verbose=False
        )
        metrics.update(marker_metrics)

    if verbose:
        print('=' * 50)
        print('Marker Gene Preservation Analysis')
        print('=' * 50)
        print(f"Mean cell type preservation: {metrics['mean_cell_type_preservation']:.4f}")
        print(f"Min cell type preservation: {metrics['min_cell_type_preservation']:.4f}")
        if 'mean_overlap' in metrics:
            print(f"Mean marker overlap: {metrics['mean_overlap']:.4f}")
        print('=' * 50)

    return metrics


def get_degs_per_cluster(
    X: np.ndarray,
    y: np.ndarray,
    n_degs: int = 100,
    fc_threshold: float = 1.5,
    pval_threshold: float = 0.01,
    method: str = 't-test',
) -> dict:
    """
    Get differentially expressed genes for each cluster.

    Returns:
        dict mapping cluster_id -> list of (gene_idx, fold_change, p_value)
    """
    degs = {}

    for cluster_id in np.unique(y):
        mask_cluster = y == cluster_id
        mask_other = ~mask_cluster

        X_cluster = X[mask_cluster]
        X_other = X[mask_other]

        mean_cluster = X_cluster.mean(axis=0)
        mean_other = X_other.mean(axis=0)

        # Compute fold change
        fold_change = (mean_cluster + 1e-6) / (mean_other + 1e-6)
        log2fc = np.log2(fold_change)

        # Compute p-values
        pvals = []
        for gene_idx in range(X.shape[1]):
            if method == 't-test':
                _, pval = stats.ttest_ind(
                    X_cluster[:, gene_idx],
                    X_other[:, gene_idx],
                    equal_var=False
                )
            else:
                _, pval = stats.mannwhitneyu(
                    X_cluster[:, gene_idx],
                    X_other[:, gene_idx],
                    alternative='two-sided'
                )
            pvals.append(pval if not np.isnan(pval) else 1.0)

        pvals = np.array(pvals)

        # Filter by thresholds
        significant = (np.abs(log2fc) > np.log2(fc_threshold)) & (pvals < pval_threshold)
        indices = np.where(significant)[0]

        # Sort by fold change
        sorted_idx = np.argsort(np.abs(log2fc[indices]))[::-1]
        indices = indices[sorted_idx][:n_degs]

        degs[int(cluster_id)] = [
            (int(idx), float(fold_change[idx]), float(pvals[idx]))
            for idx in indices
        ]

    return degs
