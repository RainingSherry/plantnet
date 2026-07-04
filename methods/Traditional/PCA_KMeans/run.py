#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
PCA_KMeans
==========

Traditional fixed-K baseline: HVG/log-normalized expression -> PCA -> KMeans.

This runner intentionally reports the known-K protocol: labels are not used for
training or PCA, but the number of ground-truth classes is supplied through
--n_clusters, as in many deep clustering benchmark tables.
"""

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

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
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
    parser.add_argument("--method_name", default="PCA+KMeans known-K")
    parser.add_argument("--variant_name", default="pca_kmeans_known_k")
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--pca_dim", type=int, default=128)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip_eval", type=family.str2bool, default=False)
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

    np.save(save_dir / "embedding_final.npy", embedding)
    np.save(save_dir / "labels.npy", labels)
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels)
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))

    eval_result = None
    preds = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            args.method_name,
            args.seed,
            embedding,
            labels,
            n_clusters,
            {
                "variant": args.variant_name,
                "pca_dim": int(dim),
                "preprocessing": "scMAE_family",
            },
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64) if preds is not None else np.zeros(n_clusters)
    frac = counts / max(1.0, counts.sum())
    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "method_raw": args.variant_name,
        "seed": int(args.seed),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "pca_dim": int(dim),
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "std_profile": effective_dimensionality(embedding.std(axis=0)),
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
    }
    save_json(summary, str(save_dir / "summary.json"))
    ari = summary["fixed_metrics"].get("kmeans_known_k", {}).get("ari")
    nmi = summary["fixed_metrics"].get("kmeans_known_k", {}).get("nmi")
    print(f"[RESULT] {dataset_name} PCA{dim}+KMeans known-K ARI={ari} NMI={nmi}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
