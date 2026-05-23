import numpy as np


def evaluate_support_predictions(true_support, predicted_probs, threshold=0.5):
    pred_support = (predicted_probs >= threshold).astype(np.float32)
    true_support = true_support.astype(np.float32)
    tp = float((pred_support * true_support).sum())
    fp = float((pred_support * (1 - true_support)).sum())
    fn = float(((1 - pred_support) * true_support).sum())
    precision = tp / max(tp + fp, 1e-8)
    recall = tp / max(tp + fn, 1e-8)
    f1 = 2 * precision * recall / max(precision + recall, 1e-8)
    sparsity = float(1.0 - pred_support.mean())
    return {
        "support_precision": precision,
        "support_recall": recall,
        "support_f1": f1,
        "predicted_sparsity": sparsity,
    }
