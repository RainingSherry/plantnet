"""Marker gene enrichment evaluation.

Measures how well the predicted clusters recover the same marker genes
as the ground-truth cell types, using Wilcoxon rank-sum test.
"""

from typing import Optional, List
import numpy as np


def marker_gene_enrichment(
    adata,
    pred_key: str = "scspade_cluster",
    label_key: Optional[str] = None,
    top_n: int = 20,
) -> dict:
    """Evaluate marker gene recovery between ground-truth and predicted clusters.

    Uses scanpy's rank_genes_groups with Wilcoxon test for both ground-truth
    and predicted clusterings, then computes overlap of top-N markers.

    Args:
        adata: AnnData object with .X (normalized expression), .obs[pred_key],
               and .obs[label_key].
        pred_key: Column in adata.obs with predicted cluster labels.
        label_key: Column in adata.obs with ground-truth labels.
                   If None, auto-detected.
        top_n: Number of top marker genes to compare per cluster.

    Returns:
        dict with marker_overlap, n_shared_markers, etc.
    """
    import scanpy as sc

    if label_key is None:
        # Auto-detect ground-truth label key
        for candidate in ["cell_type", "label", "Cluster", "cluster", "celltype"]:
            if candidate in adata.obs.columns:
                label_key = candidate
                break

    if label_key is None or label_key not in adata.obs.columns:
        return {
            "marker_overlap": None,
            "n_shared_markers": None,
            "note": f"Could not find ground-truth label key",
        }

    # Make a copy to avoid modifying the original
    work = adata.copy()
    for key in [pred_key, label_key]:
        if key in work.obs.columns:
            work.obs[key] = work.obs[key].astype(str).astype("category")

    # Rank genes by Wilcoxon for ground truth
    try:
        sc.tl.rank_genes_groups(
            work,
            groupby=label_key,
            method="wilcoxon",
            use_raw=False,
            n_genes=work.n_vars,
        )
        true_marker_sets = {}
        for grp in work.uns["rank_genes_groups"]["names"]:
            true_marker_sets[grp] = set(
                work.uns["rank_genes_groups"]["names"][grp][:top_n].astype(str).tolist()
            )
    except Exception:
        true_marker_sets = {}

    # Rank genes by Wilcoxon for predictions
    try:
        sc.tl.rank_genes_groups(
            work,
            groupby=pred_key,
            method="wilcoxon",
            use_raw=False,
            n_genes=work.n_vars,
        )
        pred_marker_sets = {}
        for grp in work.uns["rank_genes_groups"]["names"]:
            pred_marker_sets[grp] = set(
                work.uns["rank_genes_groups"]["names"][grp][:top_n].astype(str).tolist()
            )
    except Exception:
        pred_marker_sets = {}

    if not true_marker_sets or not pred_marker_sets:
        return {
            "marker_overlap": None,
            "note": "Could not compute marker genes (too few cells per cluster?)",
        }

    # Compute pairwise Jaccard overlap between matching clusters
    # (clusters are matched by size similarity)
    from scipy.optimize import linear_sum_assignment

    common_labels = list(true_marker_sets.keys())
    common_preds = list(pred_marker_sets.keys())

    if not common_labels or not common_preds:
        return {"marker_overlap": 0.0, "note": "Empty cluster labels"}

    n_l, n_p = len(common_labels), len(common_preds)
    cost = np.zeros((n_l, n_p), dtype=np.float64)

    for i, label in enumerate(common_labels):
        for j, pred in enumerate(common_preds):
            overlap = len(true_marker_sets[label] & pred_marker_sets[pred])
            union = len(true_marker_sets[label] | pred_marker_sets[pred])
            cost[i, j] = -overlap / (union + 1e-10)

    row_ind, col_ind = linear_sum_assignment(cost)
    overlaps = [-cost[i, j] for i, j in zip(row_ind, col_ind)]
    jaccard_scores = [o / (2 * top_n - o + 1e-10) for o in overlaps]

    avg_jaccard = np.mean(jaccard_scores)
    max_jaccard = np.max(jaccard_scores)

    return {
        "marker_overlap": float(avg_jaccard),
        "max_marker_overlap": float(max_jaccard),
        "n_shared_markers": int(np.sum(jaccard_scores)),
        "per_cluster_jaccard": {
            str(common_labels[i]): float(jaccard_scores[row_ind.tolist().index(i)])
            for i in range(len(row_ind))
            if i in row_ind
        },
        "top_n": top_n,
    }


def gene_set_overlap_score(
    marker_genes_true: List[str],
    marker_genes_pred: List[str],
) -> dict:
    """Simple Jaccard overlap between two marker gene lists.

    Args:
        marker_genes_true: Ground-truth marker gene names.
        marker_genes_pred: Predicted marker gene names.
    Returns:
        dict with jaccard, overlap_count, etc.
    """
    set_true = set(marker_genes_true)
    set_pred = set(marker_genes_pred)
    overlap = set_true & set_pred
    union = set_true | set_pred
    jaccard = len(overlap) / (len(union) + 1e-10)
    return {
        "jaccard": float(jaccard),
        "overlap_count": len(overlap),
        "precision": len(overlap) / (len(set_pred) + 1e-10),
        "recall": len(overlap) / (len(set_true) + 1e-10),
    }
