#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import random
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import torch
import torch.nn.functional as F
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.metrics import (
    adjusted_rand_score,
    completeness_score,
    f1_score,
    fowlkes_mallows_score,
    homogeneity_score,
    normalized_mutual_info_score,
    v_measure_score,
)
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import LabelEncoder, StandardScaler, normalize
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

try:
    from .model import AutoEncoder
except ImportError:
    from model import AutoEncoder

from methods.DeepLearning.PlantSPADE_LGCL.utils import ensure_dir, save_json


LABEL_CANDIDATES = [
    "maintype",
    "cell_type",
    "Celltype",
    "celltype",
    "label",
    "labels",
    "cell_label",
    "Cluster",
    "cluster",
    "clusters",
    "Seurat_clusters",
]


@dataclass
class DataBundle:
    adata: sc.AnnData
    data: np.ndarray
    labels: np.ndarray
    label_names: np.ndarray
    label_key: str
    gene_names: np.ndarray
    profile: dict
    preprocess_config: dict


class IndexedExpressionDataset(Dataset):
    def __init__(self, data: np.ndarray, labels: np.ndarray):
        self.data = torch.as_tensor(data, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.data[idx], self.labels[idx]


def str2bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "t", "yes", "y"}:
        return True
    if value in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got {value!r}")


def parse_float_list(value: str):
    if value is None or str(value).strip() == "":
        return []
    return [float(item) for item in str(value).split(",") if item.strip()]


def parse_args():
    parser = argparse.ArgumentParser(description="NeighborMix-scMAE fixed-protocol runner")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default="NeighborMix_scMAE")
    parser.add_argument("--variant_name", default="nm_scmae_mid")
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.0)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_loss_weight", type=float, default=0.7)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("--neighbor_k", type=int, default=5)
    parser.add_argument("--mix_neighbors", type=int, default=4)
    parser.add_argument("--mix_weight", type=float, default=0.5)
    parser.add_argument("--consistency_weight", type=float, default=0.02)
    parser.add_argument("--target_mode", default="original", choices=["original", "mixed"])
    parser.add_argument("--tau", type=float, default=0.2)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--eval_neighbors", type=int, default=15)
    parser.add_argument("--leiden_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--louvain_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--leiden_resolutions", default="0.2,0.4,0.6,0.8,1.0,1.2")
    parser.add_argument("--include_louvain", type=str2bool, default=False)
    parser.add_argument("--run_oracle_sweep", type=str2bool, default=False)
    parser.add_argument("--sweep_max_cells", type=int, default=10000)
    parser.add_argument("--silhouette_sample_size", type=int, default=3000)
    parser.add_argument("--skip_eval", type=str2bool, default=False)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--no_save_h5ad", action="store_true")
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible:
        visible_ids = [item.strip() for item in visible.split(",") if item.strip()]
        if set(visible_ids).intersection({"0", "7"}):
            raise ValueError("CUDA_VISIBLE_DEVICES includes forbidden physical GPU 0 or 7.")
        if len(visible_ids) == 1:
            return torch.device("cuda:0")
        if str(gpu) in visible_ids:
            return torch.device(f"cuda:{visible_ids.index(str(gpu))}")
        if 0 <= gpu < len(visible_ids):
            return torch.device(f"cuda:{gpu}")
        raise ValueError(f"--gpu {gpu} is outside isolated CUDA_VISIBLE_DEVICES={visible!r}.")
    if gpu in {0, 7}:
        raise ValueError("Physical GPU 0 and GPU 7 are forbidden. Use 1,2,3,4,5,6 or --no_cuda.")
    return torch.device(f"cuda:{gpu}")


def _ensure_csr(matrix) -> sp.csr_matrix:
    if sp.issparse(matrix):
        out = matrix.tocsr().astype(np.float32)
    else:
        out = sp.csr_matrix(np.asarray(matrix, dtype=np.float32))
    out.data = np.nan_to_num(out.data, nan=0.0, posinf=0.0, neginf=0.0)
    out.data[out.data < 0.0] = 0.0
    out.eliminate_zeros()
    out.sort_indices()
    return out


def _sample_values(matrix, max_rows: int = 256) -> np.ndarray:
    sample = matrix[: min(max_rows, matrix.shape[0])]
    if sp.issparse(sample):
        return sample.data.astype(np.float32, copy=False) if sample.nnz else np.array([], dtype=np.float32)
    return np.asarray(sample, dtype=np.float32).ravel()


def _looks_like_raw_counts(matrix) -> bool:
    values = _sample_values(matrix)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return False
    values = values[: min(values.size, 100000)]
    return bool(np.all(values >= 0.0) and np.allclose(values, np.round(values), atol=1e-4))


def _var_gene_names(var) -> tuple[np.ndarray, object]:
    var = var.copy()
    gene_names = np.asarray(var.index).astype(str)
    for name_col in ["gene_name", "features", "gene_symbols", "symbol", "_index"]:
        if name_col in var.columns:
            candidate = var[name_col].astype(str).to_numpy()
            if len(candidate) == len(gene_names) and len(np.unique(candidate)) == len(candidate):
                gene_names = candidate
                break
    var.index = gene_names
    return gene_names, var


def _select_count_source(adata: sc.AnnData, input_mode: str):
    if input_mode in {"auto", "raw"} and "counts" in adata.layers:
        gene_names, var = _var_gene_names(adata.var)
        return adata.layers["counts"], gene_names, var, "layers[counts]", "raw"
    if input_mode in {"auto", "raw"} and adata.raw is not None:
        gene_names, var = _var_gene_names(adata.raw.var)
        return adata.raw.X, gene_names, var, "adata.raw.X", "raw"
    inferred = "raw" if _looks_like_raw_counts(adata.X) else "log1p"
    if input_mode == "raw" and inferred != "raw":
        raise ValueError("--input_mode raw was requested, but no raw-looking X/raw/layers[counts] source is available.")
    if input_mode == "log1p":
        inferred = "log1p"
    gene_names, var = _var_gene_names(adata.var)
    source = "adata.X" if inferred == "raw" else "adata.X_log1p_fallback"
    return adata.X, gene_names, var, source, inferred


def _resolve_labels(adata: sc.AnnData, label_key: str):
    if label_key and label_key != "auto":
        if label_key not in adata.obs.columns:
            raise KeyError(f"Configured label_key={label_key!r} is absent. Available obs columns: {list(adata.obs.columns)}")
        raw = adata.obs[label_key].astype(str).to_numpy()
        encoder = LabelEncoder()
        return encoder.fit_transform(raw).astype(np.int64), encoder.classes_, label_key
    for candidate in LABEL_CANDIDATES:
        if candidate in adata.obs.columns:
            raw = adata.obs[candidate].astype(str).to_numpy()
            encoder = LabelEncoder()
            return encoder.fit_transform(raw).astype(np.int64), encoder.classes_, candidate
    raise KeyError(f"No label column found. Available obs columns: {list(adata.obs.columns)}")


def _dense_float32(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        arr = matrix.toarray()
    else:
        arr = np.asarray(matrix)
    arr = np.asarray(arr, dtype=np.float32)
    return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def _normalize_total_log1p(matrix, target_sum: float) -> sp.csr_matrix:
    x = _ensure_csr(matrix).astype(np.float32)
    row_sum = np.asarray(x.sum(axis=1)).ravel().astype(np.float32)
    scale = np.divide(
        float(target_sum),
        row_sum,
        out=np.zeros_like(row_sum, dtype=np.float32),
        where=row_sum > 0.0,
    )
    x = x.multiply(scale[:, None]).tocsr()
    x.data = np.log1p(x.data).astype(np.float32, copy=False)
    x.eliminate_zeros()
    return x


def load_scmae_dataset(
    file_path: str,
    input_mode: str,
    n_top_genes: int,
    target_sum: float,
    scale_input: bool,
    label_key: str,
    seed: int,
) -> DataBundle:
    adata = sc.read_h5ad(file_path)
    source_x, gene_names, var, counts_source, inferred_mode = _select_count_source(adata, input_mode)
    counts = _ensure_csr(source_x)
    work = sc.AnnData(X=counts.copy(), obs=adata.obs.copy(), var=var.copy())
    work.obs_names = adata.obs_names.copy()
    work.var_names = gene_names.copy()

    if inferred_mode == "raw":
        work.X = _normalize_total_log1p(work.X, target_sum=target_sum)
    elif _sample_values(work.X).size and float(np.nanmax(_sample_values(work.X))) > 30.0:
        work.X = _normalize_total_log1p(work.X, target_sum=target_sum)

    if n_top_genes and n_top_genes > 0 and work.n_vars > n_top_genes:
        sc.pp.highly_variable_genes(work, flavor="seurat", n_top_genes=n_top_genes, subset=True)

    if scale_input:
        sc.pp.scale(work)

    data = _dense_float32(work.X)
    labels, label_names, resolved_label_key = _resolve_labels(work, label_key)
    gene_names = np.asarray(work.var_names).astype(str)

    label_counts = {
        str(key): int(value)
        for key, value in work.obs[resolved_label_key].astype(str).value_counts(dropna=False).sort_index().items()
    }
    profile = {
        "dataset_name": Path(file_path).stem,
        "n_cells": int(work.n_obs),
        "n_genes_original": int(adata.n_vars),
        "n_genes": int(work.n_vars),
        "label_key": resolved_label_key,
        "n_cell_types": int(len(label_names)),
        "cell_type_counts": label_counts,
        "counts_source": counts_source,
        "input_mode": inferred_mode,
        "has_raw": bool(adata.raw is not None),
        "has_layers_counts": bool("counts" in adata.layers),
        "scale_input": bool(scale_input),
    }
    preprocess_config = {
        "file_path": str(file_path),
        "counts_source": counts_source,
        "input_mode": inferred_mode,
        "normalization": f"normalize_total(target_sum={target_sum}) + log1p when raw",
        "hvg": {"n_top_genes": int(n_top_genes), "flavor": "seurat"},
        "scale_input": bool(scale_input),
        "selected_n_genes": int(work.n_vars),
        "seed": int(seed),
        "label_key": resolved_label_key,
    }
    return DataBundle(
        adata=work,
        data=data,
        labels=labels,
        label_names=label_names,
        label_key=resolved_label_key,
        gene_names=gene_names,
        profile=profile,
        preprocess_config=preprocess_config,
    )


def build_knn_distribution(
    data: np.ndarray,
    k: int,
    pca_dim: int,
    tau: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    n_cells, n_genes = data.shape
    k = max(1, min(int(k), n_cells - 1))
    n_components = max(2, min(int(pca_dim), n_cells - 1, n_genes - 1))
    if n_cells < 3:
        indices = np.zeros((n_cells, 1), dtype=np.int64)
        probs = np.ones((n_cells, 1), dtype=np.float32)
        return indices, probs, {"neighbor_k": 1, "pca_dim": 0, "note": "too_few_cells"}

    if n_genes > n_components:
        reducer = PCA(n_components=n_components, random_state=seed, svd_solver="randomized")
        emb = reducer.fit_transform(data).astype(np.float32)
    else:
        emb = data.astype(np.float32, copy=False)
    if emb.shape[1] > 1:
        emb = StandardScaler().fit_transform(emb).astype(np.float32)
    emb = normalize(emb, norm="l2", axis=1, copy=False).astype(np.float32)

    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine")
    nn.fit(emb)
    distances, neighbors = nn.kneighbors(emb, return_distance=True)
    neighbors = neighbors[:, 1:].astype(np.int64, copy=False)
    distances = distances[:, 1:].astype(np.float32, copy=False)
    sim = np.clip(1.0 - distances, a_min=0.0, a_max=None)
    logits = sim / max(float(tau), 1e-6)
    logits = logits - logits.max(axis=1, keepdims=True)
    probs = np.exp(logits).astype(np.float64)
    probs = probs / probs.sum(axis=1, keepdims=True).clip(min=1e-12)
    return neighbors, probs.astype(np.float32), {
        "enabled": True,
        "neighbor_k": int(k),
        "pca_dim": int(n_components),
        "tau": float(tau),
        "mean_max_neighbor_prob": float(probs.max(axis=1).mean()),
        "mean_neighbor_similarity": float(sim.mean()),
    }


def sample_mix(
    data_np: np.ndarray,
    batch_indices: np.ndarray,
    batch_x: torch.Tensor,
    alpha: float,
    mix_neighbors: int,
    rng: np.random.Generator,
    neighbor_indices: np.ndarray,
    neighbor_probs: np.ndarray,
) -> torch.Tensor:
    bsz = int(batch_indices.shape[0])
    mix_neighbors = max(1, int(mix_neighbors))
    sampled = np.empty((bsz, mix_neighbors), dtype=np.int64)
    weights = np.empty((bsz, mix_neighbors), dtype=np.float32)
    for pos, cell in enumerate(batch_indices):
        probs = neighbor_probs[cell]
        choices = rng.choice(neighbor_indices.shape[1], size=mix_neighbors, replace=True, p=probs)
        sampled[pos] = neighbor_indices[cell, choices]
        picked = probs[choices].astype(np.float32, copy=False)
        weights[pos] = picked / max(float(picked.sum()), 1e-12)

    neighbor_expr = data_np[sampled]
    neighbor_mean = np.sum(neighbor_expr * weights[:, :, None], axis=1).astype(np.float32)
    neighbor_t = torch.as_tensor(neighbor_mean, dtype=batch_x.dtype, device=batch_x.device)
    alpha = float(alpha)
    return alpha * batch_x + (1.0 - alpha) * neighbor_t


def apply_scmae_noise(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    should_swap = torch.bernoulli(float(mask_ratio) * torch.ones_like(x))
    if x.shape[0] <= 1:
        replacement = x
    else:
        replacement = x[torch.randperm(x.shape[0], device=x.device)]
    corrupted = torch.where(should_swap.bool(), replacement, x)
    mask = (corrupted != x).float()
    return corrupted, mask


@torch.no_grad()
def extract_embedding(model: AutoEncoder, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    embeddings = []
    labels = []
    for _, x, y in loader:
        z = model.feature(x.to(device))
        embeddings.append(z.detach().cpu().numpy())
        labels.append(y.numpy())
    emb = np.concatenate(embeddings, axis=0).astype(np.float32)
    labels_np = np.concatenate(labels, axis=0).astype(np.int64)
    emb = np.nan_to_num(emb, nan=0.0, posinf=0.0, neginf=0.0)
    return emb, labels_np


def save_embedding_h5(path: Path, embedding: np.ndarray, labels: np.ndarray) -> None:
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=embedding.astype(np.float32))
        handle.create_dataset("labels", data=labels.astype(np.int64))


def best_map(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    true_values = np.unique(y_true)
    pred_values = np.unique(y_pred)
    n = max(len(true_values), len(pred_values))
    counts = np.zeros((n, n), dtype=np.int64)
    for i, true_label in enumerate(true_values):
        for j, pred_label in enumerate(pred_values):
            counts[i, j] = int(np.sum((y_true == true_label) & (y_pred == pred_label)))
    rows, cols = linear_sum_assignment(-counts)
    mapped = np.zeros_like(y_pred, dtype=np.int64)
    for row, col in zip(rows, cols):
        if row < len(true_values) and col < len(pred_values):
            mapped[y_pred == pred_values[col]] = true_values[row]
    return mapped


def compute_kmeans_metrics(labels: np.ndarray, pred: np.ndarray) -> tuple[dict, np.ndarray]:
    mapped = best_map(labels, pred)
    metrics = {
        "acc": float(np.mean(mapped == labels)),
        "nmi": float(normalized_mutual_info_score(labels, pred)),
        "ari": float(adjusted_rand_score(labels, pred)),
        "f1_macro": float(f1_score(labels, mapped, average="macro", zero_division=0)),
        "fmi": float(fowlkes_mallows_score(labels, pred)),
        "v_measure": float(v_measure_score(labels, pred)),
        "homogeneity": float(homogeneity_score(labels, pred)),
        "completeness": float(completeness_score(labels, pred)),
        "n_pred_clusters": int(len(np.unique(pred))),
        "silhouette": float("nan"),
        "protocol": "fixed",
        "cluster_method": "kmeans_known_k",
        "uses_known_k": True,
    }
    return metrics, mapped.astype(np.int64)


def write_kmeans_known_k_outputs(
    output_dir: Path,
    dataset: str,
    method: str,
    seed: int,
    embedding: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    extra: dict,
) -> dict:
    pred = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    metrics, mapped = compute_kmeans_metrics(labels, pred.astype(np.int64))
    fixed = {"kmeans_known_k": metrics}
    row = {
        "dataset": dataset,
        "method": method,
        "seed": int(seed),
        **extra,
        **metrics,
    }
    pd.DataFrame([row]).to_csv(output_dir / "eval_fixed.csv", index=False)
    payload = {
        "dataset": dataset,
        "method": method,
        "seed": int(seed),
        "n_clusters": int(n_clusters),
        "fixed": fixed,
        "oracle": {},
        "sweep": [],
    }
    save_json(payload, str(output_dir / "eval_metrics.json"))
    np.save(output_dir / "eval_kmeans_known_k.npy", pred.astype(np.int64))
    np.save(output_dir / "eval_kmeans_known_k_mapped.npy", mapped)
    return {"fixed": fixed, "preds": {"kmeans_known_k": pred.astype(np.int64)}}


def main():
    args = parse_args()
    set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = get_device(args.gpu, args.no_cuda)
    print(f"Using device: {device}")

    bundle = load_scmae_dataset(
        file_path=args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        scale_input=args.scale_input,
        label_key=args.label_key,
        seed=args.seed,
    )
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    with open(save_dir / "selected_genes.txt", "w", encoding="utf-8") as handle:
        for gene in bundle.gene_names:
            handle.write(f"{gene}\n")

    data_np = bundle.data
    labels = bundle.labels
    n_clusters = int(args.n_clusters) if args.n_clusters and args.n_clusters > 0 else int(len(np.unique(labels)))
    dataset_name = args.dataset_name or Path(args.data_path).stem
    print(f"Cells={data_np.shape[0]} genes={data_np.shape[1]} clusters={n_clusters} variant={args.variant_name}")

    neighbor_indices, neighbor_probs, neighbor_profile = build_knn_distribution(
        data_np,
        k=args.neighbor_k,
        pca_dim=args.knn_pca_dim,
        tau=args.tau,
        seed=args.seed,
    )
    save_json(neighbor_profile, str(save_dir / "neighbor_graph_profile.json"))

    dataset = IndexedExpressionDataset(data_np, labels)
    generator = torch.Generator()
    generator.manual_seed(args.seed)
    train_loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True,
        generator=generator,
    )
    eval_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = AutoEncoder(
        num_genes=data_np.shape[1],
        hidden_size=args.hidden_size,
        dropout=args.dropout,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight,
    ).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)
    rng = np.random.default_rng(args.seed + 2027)
    history = {
        "loss": [],
        "self_loss": [],
        "mix_loss": [],
        "consistency_loss": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        totals = {key: 0.0 for key in history}
        n_batches = 0
        for idx_t, x_cpu, _ in train_loader:
            idx_np = idx_t.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)

            x_corrupt, self_mask = apply_scmae_noise(x, args.mask_ratio)
            z_self, loss_self = model.loss_mask(x_corrupt, x, self_mask)

            x_mix = sample_mix(
                data_np=data_np,
                batch_indices=idx_np,
                batch_x=x,
                alpha=args.alpha,
                mix_neighbors=args.mix_neighbors,
                rng=rng,
                neighbor_indices=neighbor_indices,
                neighbor_probs=neighbor_probs,
            )
            x_mix_corrupt, mix_mask = apply_scmae_noise(x_mix, args.mask_ratio)
            target = x if args.target_mode == "original" else x_mix
            z_mix, loss_mix = model.loss_mask(x_mix_corrupt, target, mix_mask)
            consistency_loss = (
                1.0 - F.cosine_similarity(F.normalize(z_mix, dim=1), F.normalize(z_self.detach(), dim=1), dim=1)
            ).mean()

            loss = loss_self + float(args.mix_weight) * loss_mix + float(args.consistency_weight) * consistency_loss
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            optimizer.step()

            totals["loss"] += float(loss.detach().cpu())
            totals["self_loss"] += float(loss_self.detach().cpu())
            totals["mix_loss"] += float(loss_mix.detach().cpu())
            totals["consistency_loss"] += float(consistency_loss.detach().cpu())
            n_batches += 1

        for key in totals:
            history[key].append(totals[key] / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"self={history['self_loss'][-1]:.4f} mix={history['mix_loss'][-1]:.4f} "
                f"cons={history['consistency_loss'][-1]:.4f}"
            )

    embedding, labels_out = extract_embedding(model, eval_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save(
        {
            "model_state": model.state_dict(),
            "args": vars(args),
            "gene_names": bundle.gene_names.astype(str),
            "neighbor_profile": neighbor_profile,
        },
        save_dir / "model.pt",
    )

    result = None
    eval_extra = {
        "variant": args.variant_name,
        "alpha": float(args.alpha),
        "neighbor_k": int(args.neighbor_k),
        "mix_neighbors": int(args.mix_neighbors),
        "mix_weight": float(args.mix_weight),
        "consistency_weight": float(args.consistency_weight),
        "target_mode": args.target_mode,
        "mask_ratio": float(args.mask_ratio),
    }
    if not args.skip_eval:
        result = write_kmeans_known_k_outputs(
            output_dir=save_dir,
            dataset=dataset_name,
            method=args.method_name,
            seed=args.seed,
            embedding=embedding,
            labels=labels_out,
            n_clusters=n_clusters,
            extra=eval_extra,
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))

    if not args.no_save_h5ad:
        bundle.adata.obsm["X_neighbormix_scmae"] = embedding
        bundle.adata.uns["neighbormix_scmae"] = {
            "method": args.method_name,
            "variant": args.variant_name,
            "alpha": float(args.alpha),
            "neighbor_k": int(args.neighbor_k),
            "mix_neighbors": int(args.mix_neighbors),
            "mix_weight": float(args.mix_weight),
            "consistency_weight": float(args.consistency_weight),
            "target_mode": args.target_mode,
            "mask_ratio": float(args.mask_ratio),
        }
        bundle.adata.write_h5ad(save_dir / "adata_neighbormix_scmae.h5ad", compression="gzip")

    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "variant": args.variant_name,
        "seed": int(args.seed),
        "n_cells": int(data_np.shape[0]),
        "n_genes": int(data_np.shape[1]),
        "n_clusters": int(n_clusters),
        "embedding_path": str((save_dir / "embedding_final.npy").resolve()),
        "fixed_metrics": result["fixed"] if result is not None else {},
        "note": "scMAE AutoEncoder is the base model; NeighborMix is used only as an auxiliary training branch.",
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(f"Results saved to: {save_dir}")


if __name__ == "__main__":
    main()
