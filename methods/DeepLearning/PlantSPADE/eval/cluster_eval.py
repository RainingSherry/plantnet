import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
import scanpy as sc

from methods.evaluation import evaluation as benchmark_evaluation


def _rare_cell_preservation(labels, pred_labels):
    unique, counts = np.unique(labels, return_counts=True)
    if len(unique) == 0:
        return float("nan")
    threshold = max(1, int(0.01 * len(labels)))
    rare_labels = unique[counts <= threshold]
    if len(rare_labels) == 0:
        rare_labels = np.array([unique[np.argmin(counts)]])

    scores = []
    for rare_label in rare_labels:
        mask = labels == rare_label
        rare_pred = pred_labels[mask]
        if rare_pred.size == 0:
            continue
        majority = np.bincount(rare_pred.astype(int)).argmax()
        scores.append(float((rare_pred == majority).mean()))
    return float(np.mean(scores)) if scores else float("nan")


def cluster_and_evaluate(embedding, labels, n_clusters, method="leiden"):
    if method == "kmeans":
        pred_labels = KMeans(n_clusters=n_clusters, n_init=20, random_state=0).fit_predict(embedding)
    else:
        adata = sc.AnnData(embedding)
        sc.pp.neighbors(adata, n_neighbors=10, use_rep="X")
        resolutions = np.linspace(0.1, 2.5, 60)[::-1]
        best_res = resolutions[0]
        best_diff = float("inf")
        for res in resolutions:
            sc.tl.leiden(adata, random_state=0, resolution=float(res))
            count = adata.obs["leiden"].nunique()
            diff = abs(count - n_clusters)
            if diff < best_diff:
                best_diff = diff
                best_res = float(res)
            if count == n_clusters:
                break
        sc.tl.leiden(adata, random_state=0, resolution=best_res)
        pred_labels = adata.obs["leiden"].astype(int).to_numpy()

    acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, y_pred_ = benchmark_evaluation(labels, pred_labels)
    sil = float(silhouette_score(embedding, pred_labels)) if len(np.unique(pred_labels)) > 1 else float("nan")
    rare_preservation = _rare_cell_preservation(labels, pred_labels)
    return {
        "pred_labels": pred_labels,
        "metrics": {
            "acc": float(acc),
            "nmi": float(nmi),
            "ari": float(ari),
            "f1_macro": float(f1_macro),
            "fmi": float(fmi),
            "v_measure": float(v_measure),
            "homogeneity": float(hom),
            "completeness": float(com),
            "silhouette": sil,
            "rare_cell_preservation": rare_preservation,
        },
        "y_pred_mapped": y_pred_,
    }
