#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import numpy as np
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader() -> None:
    try:
        import h5py
        from anndata._io.specs.registry import IOSpec, _REGISTRY

        def _read_null(*args, **kwargs):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


_register_null_h5ad_reader()


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--pca_dim", type=int, default=128)
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def effective_dimensionality(std: np.ndarray) -> dict:
    var = np.square(std.astype(np.float64))
    pr = float((var.sum() ** 2) / max(float(np.square(var).sum()), 1e-12))
    return {
        "std_min": float(std.min()),
        "std_median": float(np.median(std)),
        "std_max": float(std.max()),
        "effective_dim_pr": pr,
        "dims_std_gt_0p1": int((std > 0.1).sum()),
        "dims_std_gt_1p0": int((std > 1.0).sum()),
    }


def main() -> int:
    args = parse_args()
    start = time.time()
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    dataset_name = args.dataset_name or Path(args.data_path).stem
    bundle = family.load_scmae_dataset(
        args.data_path,
        args.input_mode,
        args.n_top_genes,
        args.target_sum,
        args.scale_input,
        args.label_key,
        args.seed,
    )
    data = np.asarray(bundle.data, dtype=np.float32)
    labels = np.asarray(bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    dim = max(2, min(int(args.pca_dim), data.shape[0] - 1, data.shape[1] - 1))
    embedding = PCA(n_components=dim, random_state=args.seed).fit_transform(data.astype(np.float64)).astype(np.float32)
    pred = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit_predict(embedding)
    metrics, mapped = family.compute_kmeans_metrics(labels, pred.astype(np.int64))
    metrics["pca_dim"] = int(dim)
    fixed = {
        "kmeans_known_k": metrics,
        "pca_kmeans_known_k": metrics,
        "direct_prototype_argmax": {},
    }
    np.save(save_dir / "embedding_final.npy", embedding)
    np.save(save_dir / "labels.npy", labels)
    np.save(save_dir / "eval_kmeans_known_k.npy", pred.astype(np.int64))
    np.save(save_dir / "eval_kmeans_known_k_mapped.npy", mapped.astype(np.int64))
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    save_json(fixed, str(save_dir / "metrics.json"))
    counts = np.bincount(pred.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    summary = {
        "dataset": dataset_name,
        "seed": int(args.seed),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "method": "pca_kmeans_baseline",
        "assignment_mode": "pca",
        "variance_weight": 0.0,
        "entropy_weight": 0.0,
        "confidence_weight": 0.0,
        "latent_dim": int(dim),
        "fixed_metrics": fixed,
        "std_profile": effective_dimensionality(embedding.std(axis=0)),
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(
        f"[RESULT] {dataset_name} seed={args.seed} PCA{dim}+KMeans "
        f"ARI={metrics['ari']:.4f} NMI={metrics['nmi']:.4f} ACC={metrics['acc']:.4f}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
