"""Stage 1 scMAE-compatible runner for SRP182008.

This script has two modes:

1. diagnostics-only mode, which requires h5py/scipy/sklearn/numpy and validates
   HVG selection plus random gene-wise shuffle mask behavior.
2. training mode, reserved for environments with torch installed.

The current workspace Python has no torch, so diagnostics-only is the immediate
reproducible entry point.
"""

from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import h5py
import numpy as np
from scipy import sparse


def find_first_h5ad(root: Path) -> Path:
    for dirpath, _, filenames in os.walk(root):
        for filename in filenames:
            if filename.lower().endswith(".h5ad"):
                return Path(dirpath) / filename
    raise FileNotFoundError(f"No .h5ad file found under {root}")


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, obj: Dict) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def decode_array(values: np.ndarray) -> np.ndarray:
    if values.dtype.kind in ("S", "O"):
        decoded = []
        for value in values:
            if isinstance(value, bytes):
                decoded.append(value.decode("utf-8", "replace"))
            else:
                decoded.append(str(value))
        return np.asarray(decoded, dtype=object)
    return values


def read_obs_column(obs_group: h5py.Group, column: str) -> Optional[np.ndarray]:
    if column not in obs_group:
        return None
    node = obs_group[column]
    if isinstance(node, h5py.Dataset):
        return decode_array(np.asarray(node))
    if isinstance(node, h5py.Group) and "categories" in node and "codes" in node:
        categories = decode_array(np.asarray(node["categories"]))
        codes = np.asarray(node["codes"])
        values = np.empty(codes.shape[0], dtype=object)
        for i, code in enumerate(codes):
            values[i] = None if code < 0 else categories[int(code)]
        return values
    return None


def read_h5ad_csr(path: Path) -> Tuple[sparse.csr_matrix, Dict[str, np.ndarray], np.ndarray, Dict]:
    with h5py.File(path, "r") as f:
        x_group = f["X"]
        shape = tuple(int(x) for x in x_group.attrs["shape"])
        x = sparse.csr_matrix(
            (
                np.asarray(x_group["data"], dtype=np.float32),
                np.asarray(x_group["indices"], dtype=np.int32),
                np.asarray(x_group["indptr"], dtype=np.int64),
            ),
            shape=shape,
        )

        obs = {}
        for column in f["obs"].keys():
            values = read_obs_column(f["obs"], column)
            if values is not None:
                obs[column] = values

        var_index = decode_array(np.asarray(f["var"]["_index"]))
        metadata = {
            "root_keys": list(f.keys()),
            "x_encoding": {k: str(v) for k, v in x_group.attrs.items()},
            "obs_columns": sorted([k for k in f["obs"].keys()]),
            "var_columns": sorted([k for k in f["var"].keys()]),
            "obsm_keys": sorted(list(f["obsm"].keys())) if "obsm" in f else [],
        }
    return x, obs, var_index, metadata


def log1p_sparse(x: sparse.csr_matrix) -> sparse.csr_matrix:
    x_log = x.copy()
    x_log.data = np.log1p(x_log.data)
    return x_log


def select_hvg_by_log_variance(x: sparse.csr_matrix, n_top_genes: int) -> np.ndarray:
    if n_top_genes <= 0 or n_top_genes >= x.shape[1]:
        return np.arange(x.shape[1])
    x_log = log1p_sparse(x)
    mean = np.asarray(x_log.mean(axis=0)).ravel()
    mean_sq = np.asarray(x_log.multiply(x_log).mean(axis=0)).ravel()
    variance = np.maximum(mean_sq - mean * mean, 0.0)
    top = np.argpartition(-variance, n_top_genes - 1)[:n_top_genes]
    top = top[np.argsort(-variance[top])]
    return top.astype(np.int64)


def summarize_numeric(values: np.ndarray) -> Dict[str, float]:
    values = np.asarray(values, dtype=np.float64)
    return {
        "min": float(np.min(values)),
        "q25": float(np.quantile(values, 0.25)),
        "median": float(np.median(values)),
        "mean": float(np.mean(values)),
        "q75": float(np.quantile(values, 0.75)),
        "max": float(np.max(values)),
    }


def summarize_labels(values: np.ndarray, max_items: int = 30) -> Dict:
    counter = Counter(values.tolist())
    total = len(values)
    return {
        "n_unique": len(counter),
        "top": [
            {"label": str(label), "count": int(count), "fraction": float(count / total)}
            for label, count in counter.most_common(max_items)
        ],
    }


def dense_hvg_matrix(
    x: sparse.csr_matrix,
    hvg_idx: np.ndarray,
    max_cells: Optional[int],
    seed: int,
) -> Tuple[np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n_cells = x.shape[0]
    if max_cells is not None and max_cells > 0 and max_cells < n_cells:
        cell_idx = np.sort(rng.choice(n_cells, size=max_cells, replace=False))
    else:
        cell_idx = np.arange(n_cells)
    x_sub = x[cell_idx][:, hvg_idx].astype(np.float32)
    x_sub = log1p_sparse(x_sub)
    return np.asarray(x_sub.toarray(), dtype=np.float32), cell_idx


def random_mask(shape: Tuple[int, int], mask_ratio: float, rng: np.random.Generator) -> np.ndarray:
    return rng.random(shape, dtype=np.float32) < mask_ratio


def gene_wise_shuffle(x: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    x_prime = np.empty_like(x)
    for j in range(x.shape[1]):
        perm = rng.permutation(x.shape[0])
        x_prime[:, j] = x[perm, j]
    return x_prime


def mask_entropy(gene_frequency: np.ndarray, eps: float = 1e-12) -> float:
    probs = gene_frequency / max(float(np.sum(gene_frequency)), eps)
    probs = probs[probs > 0]
    if probs.size == 0:
        return 0.0
    return float(-(probs * np.log(probs + eps)).sum())


def compute_mask_diagnostics(
    x: np.ndarray,
    mask_ratio: float,
    seed: int,
) -> Tuple[Dict, np.ndarray]:
    rng = np.random.default_rng(seed)
    mask = random_mask(x.shape, mask_ratio, rng)
    x_prime = gene_wise_shuffle(x, rng)

    masked_count = int(mask.sum())
    total = int(mask.size)
    observed = x > 0
    masked_observed_count = int((mask & observed).sum())
    total_observed = int(observed.sum())

    changed = mask & (x_prime != x)
    zero_to_zero = mask & (x == 0) & (x_prime == 0)
    masked_zero = mask & (x == 0)
    masked_nonzero_to_zero = mask & (x > 0) & (x_prime == 0)
    masked_zero_to_nonzero = mask & (x == 0) & (x_prime > 0)

    gene_frequency = mask.mean(axis=0)
    top_freq = np.sort(gene_frequency)[::-1][:20]

    diagnostics = {
        "mask_ratio_target": float(mask_ratio),
        "n_cells": int(x.shape[0]),
        "n_genes": int(x.shape[1]),
        "total_entries": total,
        "masked_count": masked_count,
        "actual_mask_ratio_global": float(masked_count / total),
        "total_observed_entries": total_observed,
        "observed_fraction": float(total_observed / total),
        "masked_observed_count": masked_observed_count,
        "actual_mask_ratio_observed": float(masked_observed_count / max(total_observed, 1)),
        "masked_observed_fraction_among_masked": float(masked_observed_count / max(masked_count, 1)),
        "zero_to_zero_count": int(zero_to_zero.sum()),
        "zero_to_zero_fraction_among_masked": float(zero_to_zero.sum() / max(masked_count, 1)),
        "effective_changed_count": int(changed.sum()),
        "effective_changed_fraction_among_masked": float(changed.sum() / max(masked_count, 1)),
        "masked_zero_fraction_among_masked": float(masked_zero.sum() / max(masked_count, 1)),
        "masked_nonzero_to_zero_fraction_among_masked": float(
            masked_nonzero_to_zero.sum() / max(masked_count, 1)
        ),
        "masked_zero_to_nonzero_fraction_among_masked": float(
            masked_zero_to_nonzero.sum() / max(masked_count, 1)
        ),
        "gene_mask_frequency_mean": float(gene_frequency.mean()),
        "gene_mask_frequency_std": float(gene_frequency.std()),
        "gene_mask_frequency_max": float(gene_frequency.max()),
        "gene_mask_frequency_top20": [float(x) for x in top_freq],
        "mask_entropy": mask_entropy(gene_frequency),
        "mask_entropy_normalized": float(mask_entropy(gene_frequency) / np.log(x.shape[1])),
    }
    return diagnostics, gene_frequency


def build_data_summary(
    path: Path,
    x: sparse.csr_matrix,
    obs: Dict[str, np.ndarray],
    var_index: np.ndarray,
    metadata: Dict,
    hvg_idx: np.ndarray,
    selected_cells: np.ndarray,
) -> Dict:
    total = int(x.shape[0] * x.shape[1])
    summary = {
        "data_path": str(path),
        "shape": {"n_cells": int(x.shape[0]), "n_genes": int(x.shape[1])},
        "nnz": int(x.nnz),
        "density": float(x.nnz / total),
        "sparsity": float(1.0 - x.nnz / total),
        "metadata": metadata,
        "hvg": {
            "n_selected_genes": int(len(hvg_idx)),
            "selected_gene_indices_sample": [int(i) for i in hvg_idx[:20]],
            "selected_gene_names_sample": [str(var_index[i]) for i in hvg_idx[:20]],
        },
        "selected_cells": {
            "n_selected_cells": int(len(selected_cells)),
            "selected_cell_indices_sample": [int(i) for i in selected_cells[:20]],
        },
        "obs": {},
    }

    for column in [
        "Celltype",
        "Seurat_clusters",
        "Dataset",
        "Orig.ident",
        "ACE",
        "Condition",
        "Genotype",
        "Libraries",
        "Organ",
        "Tissue",
    ]:
        if column in obs:
            summary["obs"][column] = summarize_labels(obs[column])

    for column in ["nCount_RNA", "nFeature_RNA", "Percent.mt"]:
        if column in obs:
            summary["obs"][column] = summarize_numeric(obs[column])

    return summary


def require_torch():
    try:
        import torch  # noqa: F401
    except Exception as exc:  # pragma: no cover - environment-specific
        raise RuntimeError(
            "Torch is not installed in this environment. Run with --diagnostics-only, "
            "or install torch before training."
        ) from exc


def has_torch() -> bool:
    try:
        import torch  # noqa: F401
    except Exception:
        return False
    return True


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-path", type=Path, default=None)
    parser.add_argument("--save-dir", type=Path, required=True)
    parser.add_argument("--n-top-genes", type=int, default=2000)
    parser.add_argument("--max-cells", type=int, default=2048)
    parser.add_argument("--mask-ratio", type=float, default=0.3)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--diagnostics-only", action="store_true")
    parser.add_argument("--save-cache", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.time()
    save_dir = args.save_dir
    if save_dir.exists() and any(save_dir.iterdir()) and not args.force:
        raise FileExistsError(f"{save_dir} is not empty. Use --force to overwrite/add files.")
    ensure_dir(save_dir)

    data_path = args.data_path or find_first_h5ad(Path("."))
    x, obs, var_index, metadata = read_h5ad_csr(data_path)
    hvg_idx = select_hvg_by_log_variance(x, args.n_top_genes)
    x_dense, selected_cells = dense_hvg_matrix(x, hvg_idx, args.max_cells, args.seed)

    diagnostics, gene_frequency = compute_mask_diagnostics(
        x=x_dense,
        mask_ratio=args.mask_ratio,
        seed=args.seed,
    )

    config = {
        "mode": "diagnostics_only" if args.diagnostics_only else "training_requested",
        "data_path": str(data_path),
        "save_dir": str(save_dir),
        "n_top_genes": int(args.n_top_genes),
        "max_cells": int(args.max_cells) if args.max_cells is not None else None,
        "mask_ratio": float(args.mask_ratio),
        "seed": int(args.seed),
    }

    torch_available = has_torch()
    if not args.diagnostics_only:
        require_torch()
        raise NotImplementedError(
            "Training mode is intentionally not implemented in this lightweight "
            "diagnostics runner yet. Use this script to validate data/HVG/mask "
            "diagnostics, then implement the torch trainer as the next step."
        )

    data_summary = build_data_summary(
        path=data_path,
        x=x,
        obs=obs,
        var_index=var_index,
        metadata=metadata,
        hvg_idx=hvg_idx,
        selected_cells=selected_cells,
    )

    np.save(save_dir / "hvg_gene_indices.npy", hvg_idx)
    np.save(save_dir / "selected_cell_indices.npy", selected_cells)
    np.save(save_dir / "gene_mask_frequency.npy", gene_frequency)
    if args.save_cache:
        np.save(save_dir / "x_hvg_log1p.npy", x_dense.astype(np.float32))
        np.save(save_dir / "hvg_gene_names.npy", np.asarray([str(var_index[i]) for i in hvg_idx]))
        for label_name in ["Celltype", "Seurat_clusters"]:
            if label_name in obs:
                np.save(
                    save_dir / f"labels_{label_name}.npy",
                    np.asarray(obs[label_name][selected_cells], dtype=object),
                    allow_pickle=True,
                )

    write_json(save_dir / "config.json", config)
    write_json(save_dir / "data_summary.json", data_summary)
    write_json(save_dir / "mask_diagnostics.json", diagnostics)
    write_json(
        save_dir / "runtime.json",
        {
            "elapsed_seconds": float(time.time() - start),
            "torch_available": bool(torch_available),
            "training_executed": False,
        },
    )

    print(json.dumps({"config": config, "mask_diagnostics": diagnostics}, indent=2))


if __name__ == "__main__":
    main()
