from __future__ import annotations

import numpy as np


def rare_cell_summary(labels: np.ndarray, pred_labels: np.ndarray | None = None, mask_rate_per_cell: np.ndarray | None = None, rare_fraction: float = 0.05) -> dict:
    """Offline label-dependent rare-cell diagnostics. Never call from training."""
    labels = np.asarray(labels)
    values, counts = np.unique(labels, return_counts=True)
    threshold = max(1, int(np.ceil(float(rare_fraction) * labels.size)))
    rare_values = values[counts <= threshold]
    rare_mask = np.isin(labels, rare_values)
    out = {
        "status": "ok",
        "rare_fraction_threshold": float(rare_fraction),
        "rare_label_count": int(len(rare_values)),
        "rare_cell_count": int(rare_mask.sum()),
        "rare_labels": [str(v) for v in rare_values.tolist()],
    }
    if mask_rate_per_cell is not None:
        m = np.asarray(mask_rate_per_cell, dtype=np.float32)
        out["rare_cell_mask_rate"] = float(m[rare_mask].mean()) if rare_mask.any() else float("nan")
        out["nonrare_cell_mask_rate"] = float(m[~rare_mask].mean()) if (~rare_mask).any() else float("nan")
    if pred_labels is not None:
        pred = np.asarray(pred_labels)
        out["rare_cell_raw_match_rate"] = float((pred[rare_mask] == labels[rare_mask]).mean()) if rare_mask.any() else float("nan")
    return out
