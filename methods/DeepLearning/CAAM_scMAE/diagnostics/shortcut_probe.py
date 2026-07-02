from __future__ import annotations

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.multiclass import OneVsRestClassifier
from sklearn.preprocessing import OneHotEncoder


def _safe_metrics(y_true, score) -> dict:
    try:
        auroc = float(roc_auc_score(y_true, score))
    except Exception:
        auroc = float("nan")
    try:
        auprc = float(average_precision_score(y_true, score))
    except Exception:
        auprc = float("nan")
    return {"auroc": auroc, "auprc": auprc}


def run_shortcut_probe(x_tilde: np.ndarray, mask: np.ndarray, max_positions: int = 200000, seed: int = 0) -> dict:
    """Single-value probe: input is only corrupted value and gene id."""
    x_tilde = np.asarray(x_tilde, dtype=np.float32)
    mask = np.asarray(mask, dtype=np.int64)
    n, g = x_tilde.shape
    values = x_tilde.reshape(-1, 1)
    gene_ids = np.tile(np.arange(g, dtype=np.int64), n).reshape(-1, 1)
    y = mask.reshape(-1)
    rng = np.random.default_rng(seed)
    if y.size > max_positions:
        keep = rng.choice(y.size, size=max_positions, replace=False)
        values = values[keep]
        gene_ids = gene_ids[keep]
        y = y[keep]
    if len(np.unique(y)) < 2:
        return {"auroc": float("nan"), "auprc": float("nan"), "status": "single_class"}
    enc = OneHotEncoder(handle_unknown="ignore", sparse_output=False)
    gene_onehot = enc.fit_transform(gene_ids)
    features = np.concatenate([values, gene_onehot.astype(np.float32)], axis=1)
    x_train, x_test, y_train, y_test = train_test_split(features, y, test_size=0.3, random_state=seed, stratify=y)
    clf = LogisticRegression(max_iter=200, class_weight="balanced")
    clf.fit(x_train, y_train)
    score = clf.predict_proba(x_test)[:, 1]
    return {**_safe_metrics(y_test, score), "status": "ok", "probe": "single_value"}


def run_context_probe(x_tilde: np.ndarray, mask: np.ndarray, max_cells: int = 5000, seed: int = 0) -> dict:
    """Context probe: input is the whole corrupted cell, output is per-gene mask label."""
    x_tilde = np.asarray(x_tilde, dtype=np.float32)
    mask = np.asarray(mask, dtype=np.int64)
    rng = np.random.default_rng(seed)
    if x_tilde.shape[0] > max_cells:
        keep = rng.choice(x_tilde.shape[0], size=max_cells, replace=False)
        x_tilde = x_tilde[keep]
        mask = mask[keep]
    if np.unique(mask).size < 2:
        return {"auroc": float("nan"), "auprc": float("nan"), "status": "single_class"}
    x_train, x_test, y_train, y_test = train_test_split(x_tilde, mask, test_size=0.3, random_state=seed)
    clf = OneVsRestClassifier(LogisticRegression(max_iter=200, class_weight="balanced"))
    clf.fit(x_train, y_train)
    score = clf.predict_proba(x_test)
    return {**_safe_metrics(y_test.ravel(), np.asarray(score).ravel()), "status": "ok", "probe": "context"}
