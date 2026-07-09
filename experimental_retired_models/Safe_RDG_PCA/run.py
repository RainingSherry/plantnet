#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import math
import os
import random
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy.sparse.csgraph import connected_components
from scipy.sparse.linalg import eigsh
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import LabelEncoder
from sklearn.preprocessing import StandardScaler

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family  # noqa: E402
from methods.shared_utils import ensure_dir, save_json  # noqa: E402


VARIANTS = {
    "pca_kmeans",
    "pca_spectral_kmeans",
    "rdg_cell_only",
    "rdg_gene_only",
    "rdg_concat_kmeans",
    "rdg_always_on",
    "safe_rdg_heuristic",
    "neg_random_cell_graph",
    "neg_degree_shuffle_graph",
    "neg_shuffled_gene_cell_graph",
    "stage_a_all",
    "negative_controls_all",
}

STAGE_A_VARIANTS = [
    "pca_kmeans",
    "pca_spectral_kmeans",
    "rdg_cell_only",
    "rdg_gene_only",
    "rdg_concat_kmeans",
    "rdg_always_on",
    "safe_rdg_heuristic",
]

NEGATIVE_CONTROL_VARIANTS = [
    "neg_random_cell_graph",
    "neg_degree_shuffle_graph",
    "neg_shuffled_gene_cell_graph",
]


@dataclass
class Config:
    raw_pca_dim: int = 128
    cell_knn_k: int = 15
    cell_spectral_dim: int = 30
    final_knn_k: int = 15
    spectral_dim: int = 0
    kmeans_n_init: int = 50
    epsilon_jaccard: float = 1e-3
    mutual_boost_eta: float = 0.5
    gene_bootstrap_B: int = 20
    gene_bootstrap_cell_fraction: float = 0.8
    gene_knn_k: int = 10
    gene_edge_stability_threshold: float = 0.5
    min_module_size: int = 3
    max_module_size: int = 300
    gene_pca_dim: int = 20
    lambda_cell: float = 1.0
    lambda_gene: float = 1.0
    w0: float = 1.0
    wc: float = 1.0
    wg: float = 1.0
    heuristic_threshold: float = 0.45
    eps: float = 1e-8


def register_null_h5ad_reader() -> None:
    try:
        import h5py as _h5py
        from anndata._io.specs.registry import IOSpec, _REGISTRY

        def _read_null(*args, **kwargs):
            return None

        for typ in (_h5py.Dataset, _h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


def looks_like_scmae_xy_h5(path: Path) -> bool:
    if path.suffix.lower() not in {".h5", ".hdf5"}:
        return False
    try:
        with h5py.File(path, "r") as handle:
            return "X" in handle and "Y" in handle
    except Exception:
        return False


def materialize_scmae_xy_h5ad(path: Path, out_dir: Path, dataset_name: str) -> Path:
    out_path = out_dir / "_prepared_input" / f"{dataset_name}.h5ad"
    if out_path.exists():
        return out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "r") as handle:
        x = np.asarray(handle["X"])
        y_raw = np.asarray(handle["Y"])
    if y_raw.dtype.kind in {"S", "O", "U"}:
        labels_text = np.asarray([v.decode("utf-8") if isinstance(v, bytes) else str(v) for v in y_raw], dtype=object)
        label_is_integer = False
    else:
        labels_text = np.asarray(y_raw).astype(str)
        label_is_integer = bool(np.allclose(np.asarray(y_raw, dtype=float), np.round(np.asarray(y_raw, dtype=float)), equal_nan=False))
    obs = pd.DataFrame(
        {
            "resolved_label": labels_text,
            "_label_is_integer": label_is_integer,
            "n_counts": np.asarray(x, dtype=np.float64).sum(axis=1),
        },
        index=[f"cell_{i}" for i in range(x.shape[0])],
    )
    var = pd.DataFrame(index=[f"gene_{j}" for j in range(x.shape[1])])
    adata = sc.AnnData(X=np.asarray(x, dtype=np.float32), obs=obs, var=var)
    adata.uns["source_format"] = "scmae_xy_h5"
    adata.uns["source_file"] = str(path)
    adata.write_h5ad(out_path)
    meta = {
        "dataset_name": dataset_name,
        "input_path": str(path),
        "output_path": str(out_path),
        "source_format": "h5",
        "matrix_key": "X",
        "label_key": "Y",
        "resolved_label_key": "resolved_label",
        "n_cells": int(x.shape[0]),
        "n_genes": int(x.shape[1]),
        "n_clusters": int(len(LabelEncoder().fit(labels_text).classes_)),
    }
    save_json(meta, str(out_path.with_suffix(".meta.json")))
    return out_path


def load_dataset_compat(args: argparse.Namespace, save_dir: Path, dataset_name: str) -> family.DataBundle:
    data_path = Path(args.data_path)
    load_path = data_path
    if looks_like_scmae_xy_h5(data_path):
        load_path = materialize_scmae_xy_h5ad(data_path, save_dir, dataset_name)
    bundle = family.load_scmae_dataset(
        str(load_path),
        args.input_mode,
        args.n_top_genes,
        args.target_sum,
        args.scale_input,
        args.label_key,
        args.seed,
    )
    bundle.preprocess_config = {
        **bundle.preprocess_config,
        "original_data_path": str(data_path),
        "loaded_data_path": str(load_path),
        "compat_loader": "scmae_xy_h5_to_h5ad" if load_path != data_path else "scmae_family",
    }
    bundle.profile = {
        **bundle.profile,
        "original_data_path": str(data_path),
        "loaded_data_path": str(load_path),
        "compat_loader": bundle.preprocess_config["compat_loader"],
    }
    return bundle


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)


def zscore_columns(x: np.ndarray) -> np.ndarray:
    if x.size == 0:
        return x.astype(np.float32)
    out = StandardScaler(with_mean=True, with_std=True).fit_transform(np.asarray(x, dtype=np.float32))
    return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)


def compute_pca(x: np.ndarray, dim: int, seed: int) -> tuple[np.ndarray, list[float]]:
    n_comp = max(1, min(int(dim), x.shape[0] - 1, x.shape[1] - 1))
    pca = PCA(n_components=n_comp, random_state=seed, svd_solver="randomized")
    z = pca.fit_transform(zscore_columns(x)).astype(np.float32)
    return zscore_columns(z), [float(v) for v in pca.explained_variance_ratio_]


def effective_rank(z: np.ndarray, eps: float = 1e-8) -> float:
    if z.size == 0 or z.shape[1] == 0:
        return 0.0
    var = np.var(z, axis=0).astype(np.float64)
    total = float(var.sum())
    if total <= eps:
        return 0.0
    p = var / total
    return float(np.exp(-np.sum(p * np.log(p + eps))))


def compute_knn(data: np.ndarray, k: int, metric: str = "cosine") -> tuple[np.ndarray, np.ndarray]:
    n = int(data.shape[0])
    kk = max(1, min(int(k), n - 1))
    nn = NearestNeighbors(n_neighbors=kk + 1, metric=metric)
    nn.fit(data)
    distances, indices = nn.kneighbors(data, return_distance=True)
    return distances[:, 1:].astype(np.float32), indices[:, 1:].astype(np.int64)


def sym_max(a: sp.csr_matrix) -> sp.csr_matrix:
    return a.maximum(a.T).tocsr()


def build_soft_knn_graph(data: np.ndarray, k: int, metric: str = "cosine", eps: float = 1e-8) -> tuple[sp.csr_matrix, dict]:
    distances, indices = compute_knn(data, k, metric=metric)
    n, kk = indices.shape
    sigma = distances[:, -1] + eps
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    for i in range(n):
        for pos in range(kk):
            j = int(indices[i, pos])
            d = float(distances[i, pos])
            w = math.exp(-(d * d) / float(sigma[i] * sigma[j] + eps))
            rows.append(i)
            cols.append(j)
            vals.append(w)
    graph = sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float32)
    graph = sym_max(graph)
    graph.eliminate_zeros()
    return graph, {"indices": indices, "distances": distances}


def random_cell_graph_like(reference: sp.csr_matrix, seed: int) -> sp.csr_matrix:
    rng = np.random.default_rng(seed)
    ref = reference.tocsr()
    n = ref.shape[0]
    edge_count = max(1, ref.nnz // 2)
    max_edges = max(1, n * (n - 1) // 2)
    edge_count = min(edge_count, max_edges)
    rows: list[int] = []
    cols: list[int] = []
    seen: set[tuple[int, int]] = set()
    while len(seen) < edge_count:
        batch = max(1024, 2 * (edge_count - len(seen)))
        u = rng.integers(0, n, size=batch)
        v = rng.integers(0, n, size=batch)
        for a, b in zip(u, v):
            if a == b:
                continue
            i, j = (int(a), int(b)) if a < b else (int(b), int(a))
            if (i, j) not in seen:
                seen.add((i, j))
                rows.extend([i, j])
                cols.extend([j, i])
                if len(seen) >= edge_count:
                    break
    weights = np.asarray(ref.data, dtype=np.float32)
    if weights.size == 0:
        vals = np.ones(len(rows), dtype=np.float32)
    else:
        vals = rng.choice(weights, size=len(rows), replace=True).astype(np.float32)
    graph = sp.csr_matrix((vals, (rows, cols)), shape=ref.shape, dtype=np.float32)
    graph.eliminate_zeros()
    return graph


def degree_preserving_shuffle_graph(reference: sp.csr_matrix, seed: int) -> sp.csr_matrix:
    rng = np.random.default_rng(seed)
    coo = sp.triu(reference.tocsr(), k=1).tocoo()
    if coo.nnz == 0:
        return sp.csr_matrix(reference.shape, dtype=np.float32)
    endpoints = np.concatenate([coo.row, coo.col]).astype(np.int64)
    rng.shuffle(endpoints)
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    seen: set[tuple[int, int]] = set()
    weights = np.asarray(coo.data, dtype=np.float32).copy()
    rng.shuffle(weights)
    weight_pos = 0
    for pos in range(0, len(endpoints) - 1, 2):
        a = int(endpoints[pos])
        b = int(endpoints[pos + 1])
        if a == b:
            continue
        i, j = (a, b) if a < b else (b, a)
        if (i, j) in seen:
            continue
        seen.add((i, j))
        w = float(weights[weight_pos % len(weights)])
        weight_pos += 1
        rows.extend([i, j])
        cols.extend([j, i])
        vals.extend([w, w])
    graph = sp.csr_matrix((np.asarray(vals, dtype=np.float32), (rows, cols)), shape=reference.shape, dtype=np.float32)
    graph.eliminate_zeros()
    return graph


def shuffled_cell_graph_from_embedding(z: np.ndarray, k: int, seed: int, eps: float = 1e-8) -> sp.csr_matrix:
    if z.size == 0 or z.shape[1] == 0:
        return sp.csr_matrix((z.shape[0], z.shape[0]), dtype=np.float32)
    rng = np.random.default_rng(seed)
    shuffled = np.asarray(z, dtype=np.float32).copy()
    for col in range(shuffled.shape[1]):
        rng.shuffle(shuffled[:, col])
    graph, _ = build_soft_knn_graph(shuffled, k, eps=eps)
    return graph


def build_reliable_cell_graph(knn: dict, cfg: Config) -> tuple[sp.csr_matrix, dict]:
    indices = knn["indices"]
    distances = knn["distances"]
    n, kk = indices.shape
    sigma = distances[:, -1] + cfg.eps
    neighbor_sets = [set(map(int, indices[i])) for i in range(n)]
    directed = {(i, int(j)) for i in range(n) for j in indices[i]}
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    jaccards: list[float] = []
    mutual_count = 0
    for i in range(n):
        set_i = neighbor_sets[i]
        for pos in range(kk):
            j = int(indices[i, pos])
            set_j = neighbor_sets[j]
            inter = len(set_i.intersection(set_j))
            union = max(1, len(set_i.union(set_j)))
            jac = inter / union
            mutual = 1.0 if (j, i) in directed else 0.0
            if mutual:
                mutual_count += 1
            d = float(distances[i, pos])
            rbf = math.exp(-(d * d) / float(sigma[i] * sigma[j] + cfg.eps))
            boost = cfg.mutual_boost_eta + (1.0 - cfg.mutual_boost_eta) * mutual
            w = rbf * (jac + cfg.epsilon_jaccard) * boost
            rows.append(i)
            cols.append(j)
            vals.append(w)
            jaccards.append(jac)
    graph = sp.csr_matrix((vals, (rows, cols)), shape=(n, n), dtype=np.float32)
    graph = sym_max(graph)
    graph.eliminate_zeros()
    info = {
        "mutual_edge_ratio": float(mutual_count / max(1, n * kk)),
        "mean_shared_neighbor_jaccard": float(np.mean(jaccards) if jaccards else 0.0),
        "edge_weight_mean": float(np.mean(vals) if vals else 0.0),
        "edge_weight_p10": float(np.percentile(vals, 10) if vals else 0.0),
        "edge_weight_p90": float(np.percentile(vals, 90) if vals else 0.0),
    }
    return graph, info


def normalize_graph_symmetric(a: sp.csr_matrix, eps: float = 1e-8) -> sp.csr_matrix:
    a = a.astype(np.float32).tocsr()
    a = a + sp.eye(a.shape[0], dtype=np.float32, format="csr")
    degree = np.asarray(a.sum(axis=1)).ravel().astype(np.float64)
    inv = 1.0 / np.sqrt(degree + eps)
    d_inv = sp.diags(inv.astype(np.float32), format="csr")
    out = d_inv @ a @ d_inv
    out.eliminate_zeros()
    return out.tocsr()


def spectral_embedding_from_affinity(a_norm: sp.csr_matrix, dim: int, seed: int) -> np.ndarray:
    n = a_norm.shape[0]
    if n <= 2:
        return np.zeros((n, 1), dtype=np.float32)
    k = max(2, min(int(dim) + 1, n - 1))
    try:
        vals, vecs = eigsh(a_norm.astype(np.float64), k=k, which="LA", tol=1e-3, maxiter=2000)
        order = np.argsort(vals)[::-1]
        vecs = vecs[:, order]
        emb = vecs[:, 1 : min(k, int(dim) + 1)]
    except Exception:
        rng = np.random.default_rng(seed)
        emb = rng.normal(size=(n, max(1, min(int(dim), n - 1))))
    return zscore_columns(np.asarray(emb, dtype=np.float32))


def kmeans_labels(z: np.ndarray, n_clusters: int, seed: int, n_init: int) -> np.ndarray:
    return KMeans(n_clusters=int(n_clusters), n_init=int(n_init), random_state=int(seed)).fit_predict(z).astype(np.int64)


def graph_diagnostics(a: sp.csr_matrix) -> dict:
    a = a.tocsr()
    n = int(a.shape[0])
    undirected_edges = int(a.nnz // 2)
    degree = np.asarray(a.sum(axis=1)).ravel().astype(np.float64)
    n_comp, comp = connected_components(a, directed=False, connection="weak")
    counts = np.bincount(comp) if comp.size else np.array([0])
    sorted_degree = np.sort(degree)
    if sorted_degree.sum() > 0:
        idx = np.arange(1, len(sorted_degree) + 1)
        gini = float((2 * np.sum(idx * sorted_degree) / (len(sorted_degree) * sorted_degree.sum())) - (len(sorted_degree) + 1) / len(sorted_degree))
    else:
        gini = 0.0
    return {
        "n_nodes": n,
        "n_edges": undirected_edges,
        "density": float(a.nnz / max(1, n * (n - 1))),
        "n_connected_components": int(n_comp),
        "largest_component_ratio": float(counts.max() / max(1, n)),
        "degree_mean": float(degree.mean()) if degree.size else 0.0,
        "degree_std": float(degree.std()) if degree.size else 0.0,
        "degree_max": float(degree.max()) if degree.size else 0.0,
        "degree_gini": gini,
    }


def spectral_gap_proxy(a_norm: sp.csr_matrix) -> float:
    try:
        k = min(4, a_norm.shape[0] - 1)
        if k < 2:
            return 0.0
        vals = eigsh(a_norm.astype(np.float64), k=k, which="LA", return_eigenvectors=False, tol=1e-2, maxiter=1000)
        vals = np.sort(vals)[::-1]
        return float(max(0.0, vals[1] - vals[2] if len(vals) > 2 else vals[0] - vals[1]))
    except Exception:
        return 0.0


def build_bootstrap_stable_gene_graph(x: np.ndarray, cfg: Config, seed: int) -> tuple[sp.csr_matrix, dict]:
    rng = np.random.default_rng(seed)
    n_cells, n_genes = x.shape
    if n_genes < cfg.min_module_size:
        return sp.csr_matrix((n_genes, n_genes), dtype=np.float32), {"stable_edge_ratio": 0.0}
    sample_size = max(3, min(n_cells, int(round(cfg.gene_bootstrap_cell_fraction * n_cells))))
    edge_counts: dict[tuple[int, int], int] = {}
    edge_sims: dict[tuple[int, int], list[float]] = {}
    for _ in range(max(1, int(cfg.gene_bootstrap_B))):
        cells = rng.choice(n_cells, size=sample_size, replace=False)
        genes_as_samples = zscore_columns(x[cells, :].T)
        kk = max(1, min(cfg.gene_knn_k, n_genes - 1))
        nn = NearestNeighbors(n_neighbors=kk + 1, metric="cosine")
        nn.fit(genes_as_samples)
        distances, indices = nn.kneighbors(genes_as_samples, return_distance=True)
        seen: dict[tuple[int, int], float] = {}
        for u in range(n_genes):
            for pos in range(1, kk + 1):
                v = int(indices[u, pos])
                sim = max(0.0, 1.0 - float(distances[u, pos]))
                if sim <= 0:
                    continue
                key = (u, v) if u < v else (v, u)
                seen[key] = max(seen.get(key, 0.0), sim)
        for key, sim in seen.items():
            edge_counts[key] = edge_counts.get(key, 0) + 1
            edge_sims.setdefault(key, []).append(sim)
    rows: list[int] = []
    cols: list[int] = []
    vals: list[float] = []
    threshold = float(cfg.gene_edge_stability_threshold)
    for (u, v), count in edge_counts.items():
        stability = count / max(1, int(cfg.gene_bootstrap_B))
        if stability >= threshold:
            weight = stability * float(np.mean(edge_sims[(u, v)]))
            rows.extend([u, v])
            cols.extend([v, u])
            vals.extend([weight, weight])
    graph = sp.csr_matrix((vals, (rows, cols)), shape=(n_genes, n_genes), dtype=np.float32)
    graph.eliminate_zeros()
    possible = max(1, n_genes * min(cfg.gene_knn_k, n_genes - 1))
    info = {
        "stable_edge_ratio": float((len(vals) / 2) / possible),
        "candidate_edge_count": int(len(edge_counts)),
        "stable_edge_count": int(len(vals) / 2),
    }
    return graph, info


def detect_gene_modules(a_gene: sp.csr_matrix, cfg: Config, seed: int) -> tuple[list[np.ndarray], dict]:
    n = a_gene.shape[0]
    if n == 0 or a_gene.nnz == 0:
        return [], {"method": "none", "n_gene_modules": 0}
    n_comp, comp = connected_components(a_gene, directed=False, connection="weak")
    modules: list[np.ndarray] = []
    oversized: list[np.ndarray] = []
    for c in range(n_comp):
        idx = np.where(comp == c)[0]
        if cfg.min_module_size <= len(idx) <= cfg.max_module_size:
            modules.append(idx)
        elif len(idx) > cfg.max_module_size:
            oversized.append(idx)
    method = "connected_components"
    if oversized:
        try:
            import igraph as ig
            import leidenalg as la

            rows, cols = a_gene.nonzero()
            edges = [(int(i), int(j)) for i, j in zip(rows, cols) if i < j]
            weights = [float(a_gene[i, j]) for i, j in edges]
            graph = ig.Graph(n=n, edges=edges, directed=False)
            graph.es["weight"] = weights
            part = la.find_partition(
                graph,
                la.RBConfigurationVertexPartition,
                weights="weight",
                resolution_parameter=1.0,
                seed=seed,
            )
            modules = []
            labels = np.asarray(part.membership)
            for lab in np.unique(labels):
                idx = np.where(labels == lab)[0]
                if cfg.min_module_size <= len(idx) <= cfg.max_module_size:
                    modules.append(idx)
            method = "leiden_fixed"
        except Exception:
            for idx in oversized:
                for start in range(0, len(idx), cfg.max_module_size):
                    chunk = idx[start : start + cfg.max_module_size]
                    if len(chunk) >= cfg.min_module_size:
                        modules.append(chunk)
            method = "component_chunks"
    sizes = [len(m) for m in modules]
    if sizes:
        probs = np.asarray(sizes, dtype=np.float64) / max(1.0, float(np.sum(sizes)))
        entropy = float(-np.sum(probs * np.log(probs + cfg.eps)) / np.log(len(probs) + cfg.eps)) if len(probs) > 1 else 0.0
    else:
        entropy = 0.0
    return modules, {
        "method": method,
        "n_gene_modules": int(len(modules)),
        "module_size_min": int(min(sizes)) if sizes else 0,
        "module_size_median": float(np.median(sizes)) if sizes else 0.0,
        "module_size_max": int(max(sizes)) if sizes else 0,
        "module_size_entropy": entropy,
    }


def compute_module_eigengenes(x: np.ndarray, modules: list[np.ndarray], seed: int) -> np.ndarray:
    parts: list[np.ndarray] = []
    for module in modules:
        sub = np.asarray(x[:, module], dtype=np.float32)
        if sub.shape[1] < 3 or float(np.var(sub)) <= 1e-12:
            continue
        z, _ = compute_pca(sub, 1, seed)
        pc = z[:, 0]
        try:
            corr = np.corrcoef(pc, zscore_columns(sub).T)[0, 1:]
            if np.nanmean(corr) < 0:
                pc = -pc
        except Exception:
            pass
        parts.append(pc.reshape(-1, 1).astype(np.float32))
    if not parts:
        return np.zeros((x.shape[0], 0), dtype=np.float32)
    return zscore_columns(np.concatenate(parts, axis=1))


def cluster_diagnostics(labels: np.ndarray) -> dict:
    counts = np.bincount(labels.astype(np.int64))
    probs = counts / max(1, counts.sum())
    entropy = float(-np.sum(probs * np.log(probs + 1e-8)) / np.log(len(counts) + 1e-8)) if len(counts) > 1 else 0.0
    return {
        "n_pred_clusters": int(len(counts)),
        "cluster_size_min": int(counts.min()) if counts.size else 0,
        "cluster_size_median": float(np.median(counts)) if counts.size else 0.0,
        "cluster_size_max": int(counts.max()) if counts.size else 0,
        "cluster_size_entropy": entropy,
        "singleton_cluster_count": int(np.sum(counts == 1)),
    }


def q_scores(cell_info: dict, cell_diag: dict, cell_gap: float, gene_info: dict, gene_diag: dict, gene_erank: float) -> dict:
    mutual = float(cell_info.get("mutual_edge_ratio", 0.0))
    jacc = float(cell_info.get("mean_shared_neighbor_jaccard", 0.0))
    lcc = float(cell_diag.get("largest_component_ratio", 0.0))
    hub = float(cell_diag.get("degree_max", 0.0)) / max(float(cell_diag.get("degree_mean", 0.0)), 1e-8)
    hub_penalty = min(1.0, max(0.0, (hub - 5.0) / 20.0))
    jacc_scaled = jacc / (jacc + 0.10) if jacc > 0 else 0.0
    gap_scaled = cell_gap / (cell_gap + 0.05) if cell_gap > 0 else 0.0
    q_cell = np.clip(0.30 * mutual + 0.30 * jacc_scaled + 0.25 * lcc + 0.15 * gap_scaled - 0.20 * hub_penalty, 0.0, 1.0)

    stable = float(gene_info.get("stable_edge_ratio", 0.0))
    n_modules = float(gene_diag.get("n_gene_modules", 0.0))
    entropy = float(gene_diag.get("module_size_entropy", 0.0))
    module_score = min(1.0, n_modules / 20.0)
    erank_score = min(1.0, gene_erank / 10.0)
    q_gene = np.clip(0.35 * min(1.0, stable * 20.0) + 0.25 * module_score + 0.25 * entropy + 0.15 * erank_score, 0.0, 1.0)
    q_total = float(0.5 * q_cell + 0.5 * q_gene)
    return {"q_cell": float(q_cell), "q_gene": float(q_gene), "q_total": q_total}


def save_embedding_h5(path: Path, embedding: np.ndarray, labels: np.ndarray) -> None:
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=embedding.astype(np.float32))
        handle.create_dataset("labels", data=labels.astype(np.int64))


def write_outputs(
    save_dir: Path,
    dataset: str,
    method_name: str,
    variant: str,
    seed: int,
    embedding: np.ndarray,
    labels_true: np.ndarray,
    labels_pred: np.ndarray,
    n_clusters: int,
    diagnostics: dict,
    gate: dict,
    args: argparse.Namespace,
    preprocess_config: dict,
) -> dict:
    metrics, mapped = family.compute_kmeans_metrics(labels_true, labels_pred.astype(np.int64))
    metrics["variant"] = variant
    metrics["method_raw"] = variant
    metrics["runtime_seconds"] = float(diagnostics.get("runtime_seconds", 0.0))
    save_json(metrics, str(save_dir / "metrics.json"))
    row = {"dataset": dataset, "method": method_name, "seed": seed, "variant": variant, **metrics}
    pd.DataFrame([row]).to_csv(save_dir / "eval_fixed.csv", index=False)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_true.astype(np.int64))
    np.save(save_dir / "pred_labels.npy", labels_pred.astype(np.int64))
    np.save(save_dir / "pred_labels_mapped.npy", mapped.astype(np.int64))
    save_embedding_h5(save_dir / "embedding.h5", embedding, labels_true)
    save_json(diagnostics, str(save_dir / "diagnostics.json"))
    save_json(gate, str(save_dir / "gate_decision.json"))
    save_json(vars(args), str(save_dir / "args.json"))
    save_json(preprocess_config, str(save_dir / "preprocess_config.json"))
    return metrics


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default="Safe-RDG-PCA")
    parser.add_argument("--variant_name", default="safe_rdg_heuristic", choices=sorted(VARIANTS))
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--raw_pca_dim", type=int, default=128)
    parser.add_argument("--cell_knn_k", type=int, default=15)
    parser.add_argument("--cell_spectral_dim", type=int, default=30)
    parser.add_argument("--final_knn_k", type=int, default=15)
    parser.add_argument("--spectral_dim", type=int, default=0)
    parser.add_argument("--kmeans_n_init", type=int, default=50)
    parser.add_argument("--gene_bootstrap_B", type=int, default=20)
    parser.add_argument("--gene_bootstrap_cell_fraction", type=float, default=0.8)
    parser.add_argument("--gene_knn_k", type=int, default=10)
    parser.add_argument("--gene_edge_stability_threshold", type=float, default=0.5)
    parser.add_argument("--lambda_cell", type=float, default=1.0)
    parser.add_argument("--lambda_gene", type=float, default=1.0)
    parser.add_argument("--heuristic_threshold", type=float, default=0.45)
    parser.add_argument("--include_negative_controls", type=family.str2bool, default=False)
    return parser.parse_args()


def cfg_from_args(args: argparse.Namespace) -> Config:
    return Config(
        raw_pca_dim=args.raw_pca_dim,
        cell_knn_k=args.cell_knn_k,
        cell_spectral_dim=args.cell_spectral_dim,
        final_knn_k=args.final_knn_k,
        spectral_dim=args.spectral_dim,
        kmeans_n_init=args.kmeans_n_init,
        gene_bootstrap_B=args.gene_bootstrap_B,
        gene_bootstrap_cell_fraction=args.gene_bootstrap_cell_fraction,
        gene_knn_k=args.gene_knn_k,
        gene_edge_stability_threshold=args.gene_edge_stability_threshold,
        lambda_cell=args.lambda_cell,
        lambda_gene=args.lambda_gene,
        heuristic_threshold=args.heuristic_threshold,
    )


def main() -> int:
    args = parse_args()
    register_null_h5ad_reader()
    set_seed(args.seed)
    start = time.time()
    save_dir = Path(ensure_dir(args.save_dir))
    dataset_name = args.dataset_name or Path(args.data_path).stem
    cfg = cfg_from_args(args)

    bundle = load_dataset_compat(args, save_dir, dataset_name)
    x = np.asarray(bundle.data, dtype=np.float32)
    y = np.asarray(bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(y)))

    z_raw, raw_var = compute_pca(x, cfg.raw_pca_dim, args.seed)
    a_pca, pca_knn = build_soft_knn_graph(z_raw, cfg.cell_knn_k, eps=cfg.eps)
    a_pca_bar = normalize_graph_symmetric(a_pca, eps=cfg.eps)

    a_cell, cell_info = build_reliable_cell_graph(pca_knn, cfg)
    a_cell_bar = normalize_graph_symmetric(a_cell, eps=cfg.eps)
    z_cell = spectral_embedding_from_affinity(a_cell_bar, cfg.cell_spectral_dim, args.seed)

    a_gene_raw, gene_info = build_bootstrap_stable_gene_graph(x, cfg, args.seed)
    gene_modules, gene_module_info = detect_gene_modules(a_gene_raw, cfg, args.seed)
    module_matrix = compute_module_eigengenes(x, gene_modules, args.seed)
    if module_matrix.shape[1] > 0:
        z_gene, gene_var = compute_pca(module_matrix, min(cfg.gene_pca_dim, module_matrix.shape[1]), args.seed)
        a_gene_cell, _ = build_soft_knn_graph(z_gene, cfg.final_knn_k, eps=cfg.eps)
        a_gene_cell_bar = normalize_graph_symmetric(a_gene_cell, eps=cfg.eps)
    else:
        z_gene = np.zeros((x.shape[0], 0), dtype=np.float32)
        gene_var = []
        a_gene_cell = sp.csr_matrix((x.shape[0], x.shape[0]), dtype=np.float32)
        a_gene_cell_bar = normalize_graph_symmetric(a_gene_cell, eps=cfg.eps)

    cell_diag = graph_diagnostics(a_cell)
    gene_cell_diag = graph_diagnostics(a_gene_cell)
    final_pca_diag = graph_diagnostics(a_pca)
    cell_gap = spectral_gap_proxy(a_cell_bar)
    qs = q_scores(cell_info, cell_diag, cell_gap, gene_info, gene_module_info, effective_rank(z_gene))
    graph_enabled = bool(qs["q_total"] >= cfg.heuristic_threshold and (qs["q_cell"] >= 0.35 or qs["q_gene"] >= 0.35))

    z_dual_parts = [z_raw]
    if z_cell.shape[1] > 0:
        z_dual_parts.append(cfg.lambda_cell * z_cell)
    if z_gene.shape[1] > 0:
        z_dual_parts.append(cfg.lambda_gene * z_gene)
    z_dual = zscore_columns(np.concatenate(z_dual_parts, axis=1))

    a_final_from_dual, _ = build_soft_knn_graph(z_dual, cfg.final_knn_k, eps=cfg.eps)
    a_final_dual_bar = normalize_graph_symmetric(a_final_from_dual, eps=cfg.eps)
    a_safe = (cfg.w0 * a_pca_bar + cfg.wc * qs["q_cell"] * a_cell_bar + cfg.wg * qs["q_gene"] * a_gene_cell_bar).tocsr()
    a_always = (cfg.w0 * a_pca_bar + cfg.wc * a_cell_bar + cfg.wg * a_gene_cell_bar).tocsr()

    spectral_dim = int(cfg.spectral_dim) if cfg.spectral_dim and cfg.spectral_dim > 0 else max(n_clusters, 10)
    variant = args.variant_name
    gate = {"variant": variant, **qs, "heuristic_threshold": cfg.heuristic_threshold, "graph_enabled": graph_enabled}

    embedding_cache: dict[str, np.ndarray] = {}
    negative_graph_cache: dict[str, sp.csr_matrix] = {}

    def negative_graph(name: str) -> sp.csr_matrix:
        if name in negative_graph_cache:
            return negative_graph_cache[name]
        if name == "random_cell":
            graph = random_cell_graph_like(a_cell, args.seed + 1009)
        elif name == "degree_shuffle":
            graph = degree_preserving_shuffle_graph(a_cell, args.seed + 2027)
        elif name == "shuffled_gene_cell":
            graph = shuffled_cell_graph_from_embedding(z_gene, cfg.final_knn_k, args.seed + 4099, eps=cfg.eps)
        else:
            raise ValueError(name)
        negative_graph_cache[name] = graph
        negative_graph_cache[f"{name}_bar"] = normalize_graph_symmetric(graph, eps=cfg.eps)
        return graph

    def negative_graph_bar(name: str) -> sp.csr_matrix:
        if f"{name}_bar" not in negative_graph_cache:
            negative_graph(name)
        return negative_graph_cache[f"{name}_bar"]

    def select_variant(v: str) -> tuple[np.ndarray, np.ndarray, dict]:
        local_gate = {
            **gate,
            "variant": v,
            "graph_enabled": False,
            "uses_graph_clustering": False,
            "uses_rdg_graph": False,
        }
        if v == "pca_kmeans":
            emb = z_raw
        elif v == "pca_spectral_kmeans":
            local_gate["uses_graph_clustering"] = True
            emb = embedding_cache.setdefault("pca_spectral", spectral_embedding_from_affinity(a_pca_bar, spectral_dim, args.seed))
        elif v == "rdg_cell_only":
            local_gate.update({"graph_enabled": True, "uses_graph_clustering": True, "uses_rdg_graph": True})
            emb = embedding_cache.setdefault(
                "cell_only",
                spectral_embedding_from_affinity((a_pca_bar + a_cell_bar).tocsr(), spectral_dim, args.seed),
            )
        elif v == "rdg_gene_only":
            local_gate.update({"graph_enabled": True, "uses_graph_clustering": True, "uses_rdg_graph": True})
            emb = embedding_cache.setdefault(
                "gene_only",
                spectral_embedding_from_affinity((a_pca_bar + a_gene_cell_bar).tocsr(), spectral_dim, args.seed),
            )
        elif v == "rdg_concat_kmeans":
            emb = z_dual
        elif v == "rdg_always_on":
            local_gate.update({"graph_enabled": True, "uses_graph_clustering": True, "uses_rdg_graph": True})
            emb = embedding_cache.setdefault("always", spectral_embedding_from_affinity(a_always, spectral_dim, args.seed))
        elif v == "safe_rdg_heuristic":
            if graph_enabled:
                local_gate.update({"graph_enabled": True, "uses_graph_clustering": True, "uses_rdg_graph": True})
                emb = embedding_cache.setdefault("safe", spectral_embedding_from_affinity(a_safe, spectral_dim, args.seed))
            else:
                local_gate["fallback_to_pca"] = True
                emb = z_raw
        elif v == "neg_random_cell_graph":
            local_gate.update({"graph_enabled": True, "uses_graph_clustering": True, "uses_rdg_graph": False, "negative_control": "random_cell_graph"})
            emb = embedding_cache.setdefault(
                "neg_random_cell_graph",
                spectral_embedding_from_affinity((a_pca_bar + negative_graph_bar("random_cell")).tocsr(), spectral_dim, args.seed),
            )
        elif v == "neg_degree_shuffle_graph":
            local_gate.update({"graph_enabled": True, "uses_graph_clustering": True, "uses_rdg_graph": False, "negative_control": "degree_preserving_shuffled_cell_graph"})
            emb = embedding_cache.setdefault(
                "neg_degree_shuffle_graph",
                spectral_embedding_from_affinity((a_pca_bar + negative_graph_bar("degree_shuffle")).tocsr(), spectral_dim, args.seed),
            )
        elif v == "neg_shuffled_gene_cell_graph":
            local_gate.update({"graph_enabled": True, "uses_graph_clustering": True, "uses_rdg_graph": False, "negative_control": "shuffled_gene_module_cell_graph"})
            emb = embedding_cache.setdefault(
                "neg_shuffled_gene_cell_graph",
                spectral_embedding_from_affinity((a_pca_bar + negative_graph_bar("shuffled_gene_cell")).tocsr(), spectral_dim, args.seed),
            )
        else:
            raise ValueError(f"Unsupported variant: {v}")
        return emb, kmeans_labels(emb, n_clusters, args.seed, cfg.kmeans_n_init), local_gate

    diagnostics = {
        "dataset": dataset_name,
        "variant": variant,
        "runtime_seconds": float(time.time() - start),
        "config": asdict(cfg),
        "preprocess_profile": bundle.profile,
        "raw_pca_explained_variance_sum": float(np.sum(raw_var)),
        "gene_pca_explained_variance_sum": float(np.sum(gene_var)) if gene_var else 0.0,
        "effective_rank": {
            "Z_raw": effective_rank(z_raw),
            "Z_cell": effective_rank(z_cell),
            "Z_gene": effective_rank(z_gene),
            "Z_dual": effective_rank(z_dual),
        },
        "graphs": {
            "A_pca": final_pca_diag,
            "A_cell_reliable": {**cell_diag, **cell_info, "spectral_gap_proxy": cell_gap},
            "A_gene_graph": {**gene_info, **gene_module_info},
            "A_gene_cell": gene_cell_diag,
            "A_final_dual": graph_diagnostics(a_final_from_dual),
        },
        "gate": gate,
    }
    if variant == "stage_a_all":
        variants_to_write = STAGE_A_VARIANTS + (NEGATIVE_CONTROL_VARIANTS if args.include_negative_controls else [])
    elif variant == "negative_controls_all":
        variants_to_write = NEGATIVE_CONTROL_VARIANTS
    else:
        variants_to_write = [variant]
    all_metrics = {}
    for out_variant in variants_to_write:
        embedding, pred, local_gate = select_variant(out_variant)
        out_dir = save_dir / out_variant if variant == "stage_a_all" else save_dir
        out_dir.mkdir(parents=True, exist_ok=True)
        local_effective_rank = {**diagnostics["effective_rank"], "Z_final": effective_rank(embedding)}
        local_diag = {
            **diagnostics,
            "variant": out_variant,
            "effective_rank": local_effective_rank,
            "clusters": cluster_diagnostics(pred),
            "gate": local_gate,
        }
        if out_variant in NEGATIVE_CONTROL_VARIANTS:
            local_diag["graphs"] = {
                **diagnostics["graphs"],
                "negative_controls": {
                    "A_random_cell": graph_diagnostics(negative_graph("random_cell")),
                    "A_degree_shuffle_cell": graph_diagnostics(negative_graph("degree_shuffle")),
                    "A_shuffled_gene_cell": graph_diagnostics(negative_graph("shuffled_gene_cell")),
                },
            }
        metrics = write_outputs(
            out_dir,
            dataset_name,
            args.method_name,
            out_variant,
            args.seed,
            embedding,
            y,
            pred,
            n_clusters,
            local_diag,
            local_gate,
            args,
            {**bundle.preprocess_config, "safe_rdg_pca": asdict(cfg)},
        )
        all_metrics[out_variant] = metrics
        print(
            f"[RESULT] {dataset_name} {out_variant} seed={args.seed} "
            f"ARI={metrics.get('ari')} NMI={metrics.get('nmi')} ACC={metrics.get('acc')} "
            f"q={qs['q_total']:.3f} enabled={graph_enabled}"
        )
    if variant == "stage_a_all":
        save_json(all_metrics, str(save_dir / "stage_a_metrics.json"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
