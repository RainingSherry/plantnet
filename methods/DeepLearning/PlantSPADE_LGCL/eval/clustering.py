from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Dict, Iterable, Optional

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans

from .metrics import compute_metrics


def _format_resolution(resolution: float) -> str:
    return str(float(resolution)).replace(".", "p")


def _neighbors(embedding: np.ndarray, n_neighbors: int):
    adata = sc.AnnData(np.asarray(embedding, dtype=np.float32))
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    return adata


def _leiden(adata, resolution: float, seed: int, key: str) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        try:
            sc.tl.leiden(
                adata,
                random_state=seed,
                resolution=float(resolution),
                key_added=key,
                flavor="igraph",
                n_iterations=2,
                directed=False,
            )
        except TypeError:
            sc.tl.leiden(adata, random_state=seed, resolution=float(resolution), key_added=key)
    return adata.obs[key].astype(int).to_numpy()


def _louvain(adata, resolution: float, seed: int, key: str) -> Optional[np.ndarray]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FutureWarning)
        try:
            sc.tl.louvain(adata, random_state=seed, resolution=float(resolution), key_added=key)
            return adata.obs[key].astype(int).to_numpy()
        except Exception:
            pass
    try:
        import networkx as nx
        from networkx.algorithms.community import louvain_communities

        conn = adata.obsp["connectivities"].tocoo()
        graph = nx.Graph()
        graph.add_nodes_from(range(adata.n_obs))
        for row, col, weight in zip(conn.row, conn.col, conn.data):
            if row < col and weight > 0:
                graph.add_edge(int(row), int(col), weight=float(weight))
        communities = louvain_communities(graph, resolution=float(resolution), seed=seed, weight="weight")
        labels = np.empty(adata.n_obs, dtype=np.int64)
        for cluster_id, nodes in enumerate(communities):
            labels[list(nodes)] = cluster_id
        return labels
    except Exception:
        return None


def _empty_metric_row(error: str) -> dict:
    return {
        "acc": float("nan"),
        "nmi": float("nan"),
        "ari": float("nan"),
        "f1_macro": float("nan"),
        "fmi": float("nan"),
        "v_measure": float("nan"),
        "homogeneity": float("nan"),
        "completeness": float("nan"),
        "n_pred_clusters": 0,
        "silhouette": float("nan"),
        "error": error,
    }


def evaluate_embedding_protocol(
    embedding: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    seed: int,
    n_neighbors: int = 15,
    leiden_fixed_resolution: float = 1.0,
    louvain_fixed_resolution: Optional[float] = 1.0,
    leiden_sweep_resolutions: Iterable[float] = (0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
    include_louvain: bool = True,
    sweep_max_cells: Optional[int] = 10000,
) -> Dict:
    embedding = np.asarray(embedding, dtype=np.float32)
    labels = np.asarray(labels)
    fixed = {}
    oracle = {}
    sweep_rows = []
    preds = {}
    mapped_preds = {}

    pred_kmeans = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    vals, mapped = compute_metrics(labels, pred_kmeans, embedding=embedding, seed=seed)
    vals.update({"protocol": "fixed", "cluster_method": "kmeans_known_k", "uses_known_k": True})
    fixed["kmeans_known_k"] = vals
    preds["kmeans_known_k"] = pred_kmeans.astype(np.int64)
    mapped_preds["kmeans_known_k"] = mapped.astype(np.int64)

    adata = _neighbors(embedding, n_neighbors=n_neighbors)
    pred_leiden_fixed = _leiden(adata, leiden_fixed_resolution, seed=seed, key="leiden_fixed")
    vals, mapped = compute_metrics(labels, pred_leiden_fixed, embedding=embedding, seed=seed)
    vals.update(
        {
            "protocol": "fixed",
            "cluster_method": "leiden_fixed",
            "resolution": float(leiden_fixed_resolution),
            "uses_known_k": False,
        }
    )
    fixed["leiden_fixed"] = vals
    preds["leiden_fixed"] = pred_leiden_fixed.astype(np.int64)
    mapped_preds["leiden_fixed"] = mapped.astype(np.int64)

    if include_louvain and louvain_fixed_resolution is not None:
        pred_louvain = _louvain(adata, louvain_fixed_resolution, seed=seed, key="louvain_fixed")
        if pred_louvain is None:
            vals = _empty_metric_row("scanpy.tl.louvain unavailable")
            vals.update(
                {
                    "protocol": "fixed",
                    "cluster_method": "louvain_fixed",
                    "resolution": float(louvain_fixed_resolution),
                    "uses_known_k": False,
                }
            )
            fixed["louvain_fixed"] = vals
        else:
            vals, mapped = compute_metrics(labels, pred_louvain, embedding=embedding, seed=seed)
            vals.update(
                {
                    "protocol": "fixed",
                    "cluster_method": "louvain_fixed",
                    "resolution": float(louvain_fixed_resolution),
                    "uses_known_k": False,
                }
            )
            fixed["louvain_fixed"] = vals
            preds["louvain_fixed"] = pred_louvain.astype(np.int64)
            mapped_preds["louvain_fixed"] = mapped.astype(np.int64)

    best_key = None
    best_nmi = -np.inf
    best_pred = None
    best_mapped = None
    if sweep_max_cells is not None and embedding.shape[0] > sweep_max_cells:
        rng = np.random.default_rng(seed)
        subsample_idx = rng.choice(embedding.shape[0], sweep_max_cells, replace=False)
        sweep_embedding = embedding[subsample_idx]
        sweep_labels = labels[subsample_idx]
        sweep_adata = _neighbors(sweep_embedding, n_neighbors=n_neighbors)
    else:
        sweep_adata = adata
        sweep_labels = labels
    for resolution in leiden_sweep_resolutions:
        key = f"leiden_res_{_format_resolution(float(resolution))}"
        pred = _leiden(sweep_adata, float(resolution), seed=seed, key=key)
        vals, mapped = compute_metrics(sweep_labels, pred, embedding=sweep_embedding if sweep_max_cells else embedding, seed=seed)
        vals.update(
            {
                "protocol": "full_sweep",
                "cluster_method": key,
                "resolution": float(resolution),
                "uses_known_k": False,
            }
        )
        sweep_rows.append(vals)
        preds[key] = pred.astype(np.int64)
        mapped_preds[key] = mapped.astype(np.int64)
        if vals["nmi"] > best_nmi:
            best_nmi = vals["nmi"]
            best_key = key
            best_pred = pred
            best_mapped = mapped

    if best_key is not None:
        best_vals = dict([row for row in sweep_rows if row["cluster_method"] == best_key][0])
        best_vals.update(
            {
                "protocol": "oracle",
                "cluster_method": "leiden_oracle_best",
                "selected_from": best_key,
                "oracle_selection_metric": "nmi",
            }
        )
        oracle["leiden_oracle_best"] = best_vals
        preds["leiden_oracle_best"] = best_pred.astype(np.int64)
        mapped_preds["leiden_oracle_best"] = best_mapped.astype(np.int64)

    return {
        "fixed": fixed,
        "oracle": oracle,
        "sweep_rows": sweep_rows,
        "preds": preds,
        "mapped_preds": mapped_preds,
    }


def _rows_from_payload(dataset: str, method: str, seed: int, section: Dict[str, dict], extra: Optional[dict] = None):
    rows = []
    for cluster_method, vals in section.items():
        row = {"dataset": dataset, "method": method, "seed": int(seed), "cluster_method": cluster_method}
        if extra:
            row.update(extra)
        row.update(vals)
        rows.append(row)
    return rows


def write_evaluation_outputs(
    output_dir: str,
    dataset: str,
    method: str,
    seed: int,
    embedding: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    n_neighbors: int = 15,
    leiden_fixed_resolution: float = 1.0,
    louvain_fixed_resolution: Optional[float] = 1.0,
    leiden_sweep_resolutions: Iterable[float] = (0.2, 0.4, 0.6, 0.8, 1.0, 1.2),
    include_louvain: bool = True,
    sweep_max_cells: Optional[int] = 10000,
    prefix: str = "eval",
    extra: Optional[dict] = None,
) -> Dict:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    result = evaluate_embedding_protocol(
        embedding=embedding,
        labels=labels,
        n_clusters=n_clusters,
        seed=seed,
        n_neighbors=n_neighbors,
        leiden_fixed_resolution=leiden_fixed_resolution,
        louvain_fixed_resolution=louvain_fixed_resolution,
        leiden_sweep_resolutions=leiden_sweep_resolutions,
        include_louvain=include_louvain,
        sweep_max_cells=sweep_max_cells,
    )

    fixed_rows = _rows_from_payload(dataset, method, seed, result["fixed"], extra=extra)
    oracle_rows = _rows_from_payload(dataset, method, seed, result["oracle"], extra=extra)
    sweep_rows = []
    for vals in result["sweep_rows"]:
        row = {"dataset": dataset, "method": method, "seed": int(seed)}
        if extra:
            row.update(extra)
        row.update(vals)
        sweep_rows.append(row)

    pd.DataFrame(fixed_rows).to_csv(output / f"{prefix}_fixed.csv", index=False)
    pd.DataFrame(oracle_rows).to_csv(output / f"{prefix}_oracle.csv", index=False)
    pd.DataFrame(sweep_rows).to_csv(output / f"{prefix}_sweep.csv", index=False)

    payload = {
        "dataset": dataset,
        "method": method,
        "seed": int(seed),
        "n_clusters": int(n_clusters),
        "fixed": result["fixed"],
        "oracle": result["oracle"],
        "sweep": result["sweep_rows"],
    }
    with open(output / f"{prefix}_metrics.json", "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, default=str)

    for name, pred in result["preds"].items():
        np.save(output / f"{prefix}_{name}.npy", pred.astype(np.int64))
    for name, pred in result["mapped_preds"].items():
        np.save(output / f"{prefix}_{name}_mapped.npy", pred.astype(np.int64))

    return result
