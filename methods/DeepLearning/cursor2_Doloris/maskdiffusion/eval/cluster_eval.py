"""Clustering evaluation for scRNA-seq analysis.

Provides:
    - benchmark_evaluation(): core metrics with Hungarian label mapping
    - cluster_and_evaluate(): full clustering pipeline (KMeans + Leiden + all metrics)
"""

import os
from typing import Optional, Tuple

import numpy as np
import scanpy as sc
from sklearn.cluster import KMeans
from sklearn.metrics import (
    silhouette_score,
    normalized_mutual_info_score as NMI,
    adjusted_rand_score as ARI,
    f1_score as sklearn_f1_score,
    fowlkes_mallows_score as FMI,
    homogeneity_score as Homogeneity,
    completeness_score as Completeness,
    v_measure_score as VMeasure,
)


# ── Try to import PlantNet's unified evaluation ─────────────────────────────────

try:
    from methods.evaluation import evaluation as _plantnet_evaluation
    from methods.utils import save as _plantnet_save
    _HAS_PLANTNET = True
except ImportError:
    _plantnet_evaluation = None
    _plantnet_save = None
    _HAS_PLANTNET = False


# ── Hungarian label mapping ─────────────────────────────────────────────────────


def _hungarian_accuracy(y_true: np.ndarray, y_pred: np.ndarray) -> Tuple[float, np.ndarray]:
    """Compute clustering accuracy after optimal label assignment (Hungarian algorithm).

    This is CRITICAL for fair metric computation. Without label mapping,
    F1-macro and ACC can be arbitrarily bad even when the clustering
    structure is correct (just the integer labels are permuted).
    """
    from scipy.optimize import linear_sum_assignment

    n_clusters = max(y_true.max(), y_pred.max()) + 1
    cost = np.zeros((n_clusters, n_clusters), dtype=np.float64)

    for i in range(n_clusters):
        for j in range(n_clusters):
            cost[i, j] = -np.sum((y_true == i) & (y_pred == j))

    row_ind, col_ind = linear_sum_assignment(cost)
    mapping = np.full(n_clusters, -1, dtype=np.int64)
    mapping[col_ind] = row_ind

    y_pred_mapped = np.vectorize(lambda x: mapping[x] if x < n_clusters else x)(y_pred)
    accuracy = np.mean(y_pred_mapped == y_true)
    return accuracy, y_pred_mapped


# ── Core metric computation ────────────────────────────────────────────────────


def benchmark_evaluation(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    use_plantnet: bool = True,
) -> Tuple[float, float, float, float, float, float, float, float, np.ndarray]:
    """Compute all clustering metrics with Hungarian label mapping.

    All metrics (including F1-macro!) are computed on the MAPPED labels.
    This is critical for fair comparison with Claude/Doloris benchmarks.

    Args:
        y_true: Ground truth labels, shape (n_cells,).
        y_pred: Predicted cluster labels, shape (n_cells,).
        use_plantnet: If True and PlantNet evaluation is available, use it.

    Returns:
        (acc, nmi, ari, f1_macro, fmi, v_measure, homogeneity, completeness, y_pred_mapped)
    """
    if use_plantnet and _HAS_PLANTNET:
        acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, y_pred_mapped = _plantnet_evaluation(
            y_true, y_pred
        )
        return acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, y_pred_mapped

    # Standalone implementation
    acc, y_pred_mapped = _hungarian_accuracy(y_true, y_pred)
    nmi = NMI(y_true, y_pred_mapped)
    ari = ARI(y_true, y_pred_mapped)
    fmi = FMI(y_true, y_pred_mapped)
    v_measure = VMeasure(y_true, y_pred_mapped)
    homogeneity = Homogeneity(y_true, y_pred_mapped)
    completeness = Completeness(y_true, y_pred_mapped)

    # F1-macro on MAPPED labels (not raw y_pred!)
    f1_macro = sklearn_f1_score(y_true, y_pred_mapped, average="macro", zero_division=0)

    return acc, nmi, ari, f1_macro, fmi, v_measure, homogeneity, completeness, y_pred_mapped


# ── Full clustering pipeline ────────────────────────────────────────────────────


def cluster_and_evaluate(
    embeddings: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    methods: tuple = ("kmeans", "leiden"),
    resolution_range: tuple = (0.1, 2.5),
    n_resolutions: int = 60,
    metric: str = "nmi",
    return_all: bool = True,
    save_dir: str = None,
    save_key: str = "clustering",
) -> dict:
    """Run clustering (KMeans + Leiden) and evaluate against ground truth.

    This function runs both KMeans and Leiden clustering, picks the best
    result (by the specified metric), and returns all metrics for both.

    Args:
        embeddings: Cell embeddings, shape (n_cells, latent_dim).
        labels: Ground truth labels, shape (n_cells,).
        n_clusters: Number of ground-truth clusters.
        methods: Which clustering methods to try.
        resolution_range: Leiden resolution search range.
        n_resolutions: Number of resolutions to search.
        metric: Which metric to use for best-method selection.
        return_all: If True, return metrics for all methods.
        save_dir: Optional directory to save results.
        save_key: Prefix for saved files.

    Returns:
        dict with keys for each method's metrics + best_method.
    """
    results = {}

    # ── KMeans ───────────────────────────────────────────────────────────────
    if "kmeans" in methods:
        kmeans = KMeans(n_clusters=n_clusters, n_init=20, random_state=42)
        pred_kmeans = kmeans.fit_predict(embeddings)
        acc, nmi, ari, f1, fmi, vm, hom, com, pred_mapped = benchmark_evaluation(
            labels, pred_kmeans
        )
        sil = silhouette_score(embeddings, pred_mapped)
        results["kmeans"] = {
            "pred_labels": pred_kmeans,
            "acc": float(acc),
            "nmi": float(nmi),
            "ari": float(ari),
            "f1_macro": float(f1),
            "fmi": float(fmi),
            "v_measure": float(vm),
            "homogeneity": float(hom),
            "completeness": float(com),
            "silhouette": float(sil),
        }

    # ── Leiden ─────────────────────────────────────────────────────────────
    if "leiden" in methods:
        adata = sc.AnnData(embeddings.astype(np.float32))
        sc.pp.neighbors(adata, n_neighbors=10, use_rep="X", random_state=42)

        # Search resolution to match ground-truth cluster count
        # Use igraph flavor for 10-100x speedup
        resolutions = np.linspace(resolution_range[0], resolution_range[1], n_resolutions)
        best_res = resolutions[0]
        best_diff = float("inf")
        try:
            leiden_flavor = "igraph"
        except NameError:
            leiden_flavor = "leidenalg"

        for res in resolutions:
            sc.tl.leiden(
                adata,
                random_state=42,
                resolution=float(res),
                key_added="tmp_leiden",
                flavor=leiden_flavor,
                n_iterations=2,
            )
            n_found = adata.obs["tmp_leiden"].nunique()
            diff = abs(n_found - n_clusters)
            if diff < best_diff:
                best_diff = diff
                best_res = float(res)
            if n_found == n_clusters:
                break

        sc.tl.leiden(
            adata,
            random_state=42,
            resolution=best_res,
            key_added="leiden",
            flavor=leiden_flavor,
            n_iterations=2,
        )
        pred_leiden = adata.obs["leiden"].astype(int).to_numpy()

        acc, nmi, ari, f1, fmi, vm, hom, com, pred_mapped = benchmark_evaluation(
            labels, pred_leiden
        )
        sil = silhouette_score(embeddings, pred_mapped)
        results["leiden"] = {
            "pred_labels": pred_leiden,
            "acc": float(acc),
            "nmi": float(nmi),
            "ari": float(ari),
            "f1_macro": float(f1),
            "fmi": float(fmi),
            "v_measure": float(vm),
            "homogeneity": float(hom),
            "completeness": float(com),
            "silhouette": float(sil),
            "resolution": best_res,
        }

    # ── Select best method ─────────────────────────────────────────────────
    if not results:
        raise ValueError("No clustering methods ran successfully.")

    best_method = max(
        results.keys(),
        key=lambda m: results[m].get(metric, -1.0),
    )
    results["best_method"] = best_method
    results["best_metrics"] = results[best_method]

    # ── Save ───────────────────────────────────────────────────────────────
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        for method_name, res in results.items():
            if isinstance(res, dict) and "pred_labels" in res:
                np.save(
                    os.path.join(save_dir, f"{save_key}_{method_name}_pred_labels.npy"),
                    res["pred_labels"],
                )
        if _HAS_PLANTNET and _plantnet_save:
            best_pred = results[best_method]["pred_labels"]
            best_metrics = results[best_method]
            _plantnet_save(
                save_dir,
                labels,
                best_pred,
                epoch="final",
                embedding=embeddings,
            )

    return results


# ── Rare cell preservation ─────────────────────────────────────────────────────


def rare_cell_preservation(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    min_cluster_size: int = 5,
) -> dict:
    """Measure how well rare cell populations are preserved in clustering.

    Rare cell groups are defined as ground-truth clusters with fewer than
    min_cluster_size cells. For each such group, check if it's captured
    as a single cluster (purity) or split across multiple (split_ratio).
    """
    from collections import Counter

    true_counts = Counter(y_true)
    rare_true_labels = {k for k, v in true_counts.items() if v <= min_cluster_size}

    if not rare_true_labels:
        return {"rare_preservation_score": 1.0, "n_rare": 0}

    pred_counts = Counter(y_pred)

    # For each rare true cluster, find its dominant predicted cluster
    preservation_scores = []
    for true_label in rare_true_labels:
        mask = y_true == true_label
        sub_pred = y_pred[mask]
        pred_counts_sub = Counter(sub_pred)
        dominant_pred = pred_counts_sub.most_common(1)[0][0]
        purity = pred_counts_sub[dominant_pred] / len(sub_pred)
        preservation_scores.append(purity)

    avg_preservation = np.mean(preservation_scores)
    return {
        "rare_preservation_score": float(avg_preservation),
        "n_rare": len(rare_true_labels),
        "rare_preservation_details": {
            int(k): float(v) for k, v in zip(rare_true_labels, preservation_scores)
        },
    }
