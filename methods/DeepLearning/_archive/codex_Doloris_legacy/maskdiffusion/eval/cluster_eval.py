import numpy as np
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

from methods.evaluation import evaluation as benchmark_evaluation


def cluster_and_evaluate(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int = 42) -> dict:
    pred_labels = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    acc, nmi, ari, f1_macro, fmi, v_measure, hom, com, y_pred_mapped = benchmark_evaluation(labels, pred_labels)
    if len(np.unique(pred_labels)) > 1 and embedding.shape[0] > len(np.unique(pred_labels)):
        try:
            silhouette = float(silhouette_score(embedding, pred_labels))
        except Exception:
            silhouette = float("nan")
    else:
        silhouette = float("nan")
    return {
        "pred_labels": pred_labels.astype(np.int64),
        "pred_labels_mapped": y_pred_mapped.astype(np.int64),
        "metrics": {
            "acc": float(acc),
            "nmi": float(nmi),
            "ari": float(ari),
            "f1_macro": float(f1_macro),
            "fmi": float(fmi),
            "v_measure": float(v_measure),
            "homogeneity": float(hom),
            "completeness": float(com),
            "silhouette": silhouette,
        },
    }


def evaluate_embedding_set(embeddings: dict, labels: np.ndarray, n_clusters: int, seed: int = 42) -> dict:
    return {
        name: cluster_and_evaluate(embedding, labels, n_clusters=n_clusters, seed=seed)
        for name, embedding in embeddings.items()
    }

