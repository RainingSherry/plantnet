# -*- coding: utf-8 -*-
"""
Sparsity evaluation.

Evaluates how well the model handles sparsity in single-cell data:
- Zero-inflation analysis
- Dropout recovery
- Sparsity-aware metrics
"""

import numpy as np
from scipy import stats


def compute_sparsity_stats(X: np.ndarray) -> dict:
    """
    Compute sparsity statistics for expression matrix.

    Args:
        X: Expression matrix (n_cells, n_genes)

    Returns:
        dict with sparsity statistics
    """
    # Total sparsity
    total_elements = X.shape[0] * X.shape[1]
    n_zeros = np.sum(X == 0)
    n_nonzeros = np.sum(X != 0)

    sparsity = n_zeros / total_elements

    # Per-cell sparsity
    cell_sparsity = np.mean(np.sum(X == 0, axis=1) / X.shape[1])

    # Per-gene sparsity
    gene_sparsity = np.mean(np.sum(X == 0, axis=0) / X.shape[0])

    # Sparsity quartiles
    sparsity_per_cell = np.sum(X == 0, axis=1) / X.shape[1]
    sparsity_q25 = np.percentile(sparsity_per_cell, 25)
    sparsity_q50 = np.percentile(sparsity_per_cell, 50)
    sparsity_q75 = np.percentile(sparsity_per_cell, 75)

    return {
        'total_sparsity': float(sparsity),
        'total_zeros': int(n_zeros),
        'total_nonzeros': int(n_nonzeros),
        'mean_cell_sparsity': float(cell_sparsity),
        'mean_gene_sparsity': float(gene_sparsity),
        'sparsity_q25': float(sparsity_q25),
        'sparsity_q50': float(sparsity_q50),
        'sparsity_q75': float(sparsity_q75),
    }


def evaluate_dropout_recovery(
    X_true: np.ndarray,
    X_pred: np.ndarray,
    dropout_rate: float = 0.3,
    verbose: bool = True,
) -> dict:
    """
    Evaluate dropout recovery performance.

    Simulates dropout and measures recovery.

    Args:
        X_true: Ground truth expression
        X_pred: Predicted expression
        dropout_rate: Rate of simulated dropout
        verbose: Whether to print results

    Returns:
        dict with recovery metrics
    """
    # Identify true zeros and dropouts
    true_zeros = X_true == 0
    true_nonzeros = X_true != 0

    # For predicted values at true zero positions
    if true_zeros.sum() > 0:
        pred_at_true_zeros = X_pred[true_zeros]
        recovery_at_zeros = np.mean(pred_at_true_zeros == 0)

        # MSE at true zero positions
        mse_at_zeros = np.mean(pred_at_true_zeros ** 2)
    else:
        recovery_at_zeros = 0
        mse_at_zeros = 0

    # For predicted values at true non-zero positions
    if true_nonzeros.sum() > 0:
        pred_at_nonzeros = X_pred[true_nonzeros]
        true_at_nonzeros = X_true[true_nonzeros]

        # Correlation
        corr, _ = stats.pearsonr(pred_at_nonzeros, true_at_nonzeros)

        # MSE
        mse_at_nonzeros = np.mean((pred_at_nonzeros - true_at_nonzeros) ** 2)

        # R-squared
        ss_res = np.sum((true_at_nonzeros - pred_at_nonzeros) ** 2)
        ss_tot = np.sum((true_at_nonzeros - true_at_nonzeros.mean()) ** 2)
        r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    else:
        corr = 0
        mse_at_nonzeros = 0
        r2 = 0

    metrics = {
        'recovery_at_zeros': float(recovery_at_zeros),
        'mse_at_zeros': float(mse_at_zeros),
        'pearson_corr_at_nonzeros': float(corr),
        'mse_at_nonzeros': float(mse_at_nonzeros),
        'r2_at_nonzeros': float(r2),
    }

    if verbose:
        print('=' * 50)
        print('Dropout Recovery Evaluation')
        print('=' * 50)
        print(f"Zero recovery rate:     {metrics['recovery_at_zeros']:.4f}")
        print(f"MSE at true zeros:      {metrics['mse_at_zeros']:.6f}")
        print(f"Pearson corr (nonzero): {metrics['pearson_corr_at_nonzeros']:.4f}")
        print(f"MSE at true nonzeros:  {metrics['mse_at_nonzeros']:.6f}")
        print(f"R2 at true nonzeros:   {metrics['r2_at_nonzeros']:.4f}")
        print('=' * 50)

    return metrics


def evaluate_sparsity_handling(
    X_original: np.ndarray,
    X_reconstructed: np.ndarray,
    mask: np.ndarray = None,
    verbose: bool = True,
) -> dict:
    """
    Comprehensive evaluation of sparsity handling.

    Args:
        X_original: Original expression matrix
        X_reconstructed: Reconstructed expression matrix
        mask: Optional mask indicating active genes (1=active, 0=inactive)
        verbose: Whether to print results

    Returns:
        dict with evaluation metrics
    """
    metrics = {}

    # Sparsity preservation
    orig_sparsity = np.mean(X_original == 0)
    recon_sparsity = np.mean(X_reconstructed == 0)
    sparsity_diff = abs(orig_sparsity - recon_sparsity)

    metrics['orig_sparsity'] = float(orig_sparsity)
    metrics['recon_sparsity'] = float(recon_sparsity)
    metrics['sparsity_diff'] = float(sparsity_diff)

    # Zero patterns
    orig_zeros = X_original == 0
    recon_zeros = X_reconstructed == 0

    # Agreement
    agreement = np.mean(orig_zeros == recon_zeros)
    metrics['zero_pattern_agreement'] = float(agreement)

    # False positive rate (predicting non-zero when true is zero)
    fp = np.sum(~orig_zeros & recon_zeros) / orig_zeros.sum()
    metrics['false_positive_rate'] = float(fp)

    # False negative rate (predicting zero when true is non-zero)
    fn = np.sum(orig_zeros & ~recon_zeros) / (~orig_zeros).sum()
    metrics['false_negative_rate'] = float(fn)

    # MSE on active genes
    if mask is not None:
        active_mask = mask > 0
        if active_mask.sum() > 0:
            mse_active = np.mean(
                ((X_reconstructed - X_original) ** 2)[active_mask]
            )
            metrics['mse_on_active_genes'] = float(mse_active)
        else:
            metrics['mse_on_active_genes'] = None
    else:
        active_mask = X_original != 0
        if active_mask.sum() > 0:
            mse_active = np.mean(
                ((X_reconstructed - X_original) ** 2)[active_mask]
            )
            metrics['mse_on_active_genes'] = float(mse_active)
        else:
            metrics['mse_on_active_genes'] = None

    if verbose:
        print('=' * 50)
        print('Sparsity Handling Evaluation')
        print('=' * 50)
        print(f"Original sparsity:        {metrics['orig_sparsity']:.4f}")
        print(f"Reconstruction sparsity:   {metrics['recon_sparsity']:.4f}")
        print(f"Sparsity difference:     {metrics['sparsity_diff']:.4f}")
        print(f"Zero pattern agreement:   {metrics['zero_pattern_agreement']:.4f}")
        print(f"False positive rate:     {metrics['false_positive_rate']:.4f}")
        print(f"False negative rate:     {metrics['false_negative_rate']:.4f}")
        if metrics['mse_on_active_genes'] is not None:
            print(f"MSE on active genes:    {metrics['mse_on_active_genes']:.6f}")
        print('=' * 50)

    return metrics


def analyze_zero_distribution(
    X: np.ndarray,
    y: np.ndarray = None,
    verbose: bool = True,
) -> dict:
    """
    Analyze the distribution of zeros across cell types.

    Args:
        X: Expression matrix
        y: Optional cluster labels
        verbose: Whether to print results

    Returns:
        dict with zero distribution analysis
    """
    metrics = {}

    # Overall sparsity
    metrics.update(compute_sparsity_stats(X))

    # Per-cluster sparsity
    if y is not None:
        unique_labels = np.unique(y)
        per_cluster_sparsity = {}

        for label in unique_labels:
            mask = y == label
            sparsity = np.mean(X[mask] == 0)
            per_cluster_sparsity[int(label)] = float(sparsity)

        metrics['per_cluster_sparsity'] = per_cluster_sparsity
        metrics['cluster_sparsity_std'] = float(np.std(list(per_cluster_sparsity.values())))

        if verbose:
            print('=' * 50)
            print('Zero Distribution Analysis')
            print('=' * 50)
            print(f"Overall sparsity: {metrics['total_sparsity']:.4f}")
            print(f"Per-cluster sparsity (std: {metrics['cluster_sparsity_std']:.4f}):")
            for label, sparsity in sorted(per_cluster_sparsity.items()):
                print(f"  Cluster {label}: {sparsity:.4f}")
            print('=' * 50)
    else:
        if verbose:
            print('=' * 50)
            print('Zero Distribution Analysis')
            print('=' * 50)
            print(f"Overall sparsity: {metrics['total_sparsity']:.4f}")
            print('=' * 50)

    return metrics


def compare_sparsity_patterns(
    X1: np.ndarray,
    X2: np.ndarray,
    verbose: bool = True,
) -> dict:
    """
    Compare sparsity patterns between two expression matrices.

    Useful for comparing original vs. denoised data.

    Args:
        X1: First expression matrix
        X2: Second expression matrix
        verbose: Whether to print results

    Returns:
        dict with comparison metrics
    """
    metrics = {}

    # Sparsity comparison
    sparsity1 = np.mean(X1 == 0)
    sparsity2 = np.mean(X2 == 0)
    metrics['sparsity1'] = float(sparsity1)
    metrics['sparsity2'] = float(sparsity2)
    metrics['sparsity_change'] = float(sparsity2 - sparsity1)

    # Correlation of non-zero patterns
    mask1 = (X1 != 0).astype(float)
    mask2 = (X2 != 0).astype(float)

    corr_mask, _ = stats.pearsonr(mask1.flatten(), mask2.flatten())
    metrics['nonzero_pattern_correlation'] = float(corr_mask)

    # Jaccard similarity of non-zero elements
    intersection = np.sum((X1 != 0) & (X2 != 0))
    union = np.sum((X1 != 0) | (X2 != 0))
    jaccard = intersection / union if union > 0 else 0
    metrics['nonzero_jaccard'] = float(jaccard)

    if verbose:
        print('=' * 50)
        print('Sparsity Pattern Comparison')
        print('=' * 50)
        print(f"Sparsity 1:       {metrics['sparsity1']:.4f}")
        print(f"Sparsity 2:       {metrics['sparsity2']:.4f}")
        print(f"Sparsity change:  {metrics['sparsity_change']:+.4f}")
        print(f"Pattern corr:     {metrics['nonzero_pattern_correlation']:.4f}")
        print(f"Nonzero Jaccard:  {metrics['nonzero_jaccard']:.4f}")
        print('=' * 50)

    return metrics
