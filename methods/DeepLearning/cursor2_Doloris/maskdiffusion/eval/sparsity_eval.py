"""Support sparsity evaluation: measure how well the mask network predicts gene activation."""

import numpy as np


def evaluate_support_predictions(
    true_support: np.ndarray,
    predicted_probs: np.ndarray,
    threshold: float = 0.5,
) -> dict:
    """Evaluate mask network performance at predicting expression support.

    Args:
        true_support: Ground truth binary support, shape (n_cells, n_genes),
                      values in {0, 1}. Must be constructed from raw counts.
        predicted_probs: Predicted activation probabilities, shape (n_cells, n_genes),
                         values in [0, 1].
        threshold: Threshold for binarizing predicted probabilities.

    Returns:
        dict with precision, recall, F1, AUROC, sparsity metrics.
    """
    from sklearn.metrics import roc_auc_score, precision_recall_curve, auc

    pred_binary = (predicted_probs >= threshold).astype(np.float32)

    # Basic TP/FP/FN
    tp = float((pred_binary * true_support).sum())
    fp = float((pred_binary * (1 - true_support)).sum())
    fn = float(((1 - pred_binary) * true_support).sum())
    tn = float(((1 - pred_binary) * (1 - true_support)).sum())

    precision = tp / (tp + fp + 1e-10)
    recall = tp / (tp + fn + 1e-10)
    f1 = 2 * precision * recall / (precision + recall + 1e-10)

    # AUROC
    try:
        auroc = roc_auc_score(
            true_support.flatten(),
            predicted_probs.flatten(),
        )
    except ValueError:
        auroc = 0.5

    # AUPR (Area Under Precision-Recall curve)
    try:
        precision_curve, recall_curve, _ = precision_recall_curve(
            true_support.flatten(), predicted_probs.flatten()
        )
        aupr = auc(recall_curve, precision_curve)
    except ValueError:
        aupr = 0.0

    # Sparsity metrics
    predicted_sparsity = 1.0 - pred_binary.mean()
    true_sparsity = 1.0 - true_support.mean()

    return {
        "support_precision": float(precision),
        "support_recall": float(recall),
        "support_f1": float(f1),
        "support_auroc": float(auroc),
        "support_aupr": float(aupr),
        "predicted_sparsity": float(predicted_sparsity),
        "true_sparsity": float(true_sparsity),
        "sparsity_error": float(abs(predicted_sparsity - true_sparsity)),
        "threshold": threshold,
    }


def compute_gene_importance_scores(
    true_support: np.ndarray,
    predicted_probs: np.ndarray,
) -> dict:
    """Compute per-gene importance scores aggregated across cells.

    Args:
        true_support: (n_cells, n_genes)
        predicted_probs: (n_cells, n_genes)

    Returns:
        dict with per-gene precision, recall, and "importance" (fraction expressed).
    """
    n_genes = true_support.shape[1]
    gene_precision = []
    gene_recall = []
    gene_importance = true_support.mean(axis=0)

    for g in range(n_genes):
        tp = ((predicted_probs[:, g] >= 0.5) * true_support[:, g]).sum()
        fp = ((predicted_probs[:, g] >= 0.5) * (1 - true_support[:, g])).sum()
        fn = ((predicted_probs[:, g] < 0.5) * true_support[:, g]).sum()
        prec = tp / (tp + fp + 1e-10)
        rec = tp / (tp + fn + 1e-10)
        gene_precision.append(prec)
        gene_recall.append(rec)

    return {
        "gene_precision": np.array(gene_precision),
        "gene_recall": np.array(gene_recall),
        "gene_importance": gene_importance,
    }
