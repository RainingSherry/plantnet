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
import pandas as pd
import scanpy as sc
from anndata import AnnData
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score, silhouette_score

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
    parser.add_argument("--pca_dim", type=int, default=128)
    parser.add_argument("--n_neighbors", type=int, default=15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--resolutions", default="0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.8,1.0,1.2,1.5,2.0,2.5,3.0")
    parser.add_argument("--silhouette_sample", type=int, default=5000)
    return parser.parse_args()


def safe_unsup_scores(embedding: np.ndarray, pred: np.ndarray, sample_size: int, seed: int) -> dict:
    n_clusters = int(len(np.unique(pred)))
    if n_clusters <= 1 or n_clusters >= embedding.shape[0]:
        return {"silhouette": float("nan"), "davies_bouldin": float("nan"), "calinski_harabasz": float("nan")}
    rng = np.random.default_rng(seed)
    if embedding.shape[0] > sample_size:
        idx = rng.choice(embedding.shape[0], size=sample_size, replace=False)
        emb = embedding[idx]
        labels = pred[idx]
    else:
        emb = embedding
        labels = pred
    out = {}
    try:
        out["silhouette"] = float(silhouette_score(emb, labels))
    except Exception:
        out["silhouette"] = float("nan")
    try:
        out["davies_bouldin"] = float(davies_bouldin_score(emb, labels))
    except Exception:
        out["davies_bouldin"] = float("nan")
    try:
        out["calinski_harabasz"] = float(calinski_harabasz_score(emb, labels))
    except Exception:
        out["calinski_harabasz"] = float("nan")
    return out


def choose_label_blind(rows: list[dict], selector: str) -> dict:
    candidates = [row for row in rows if row["n_pred_clusters"] > 1]
    if not candidates:
        return rows[0]
    if selector == "silhouette":
        finite = [row for row in candidates if np.isfinite(row.get("silhouette", float("nan")))]
        return max(finite, key=lambda row: row["silhouette"]) if finite else candidates[0]
    if selector == "davies_bouldin":
        finite = [row for row in candidates if np.isfinite(row.get("davies_bouldin", float("nan")))]
        return min(finite, key=lambda row: row["davies_bouldin"]) if finite else candidates[0]
    finite = [row for row in candidates if np.isfinite(row.get("calinski_harabasz", float("nan")))]
    return max(finite, key=lambda row: row["calinski_harabasz"]) if finite else candidates[0]


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
    true_k = int(len(np.unique(labels)))
    dim = max(2, min(int(args.pca_dim), data.shape[0] - 1, data.shape[1] - 1))
    embedding = PCA(n_components=dim, random_state=args.seed).fit_transform(data.astype(np.float64)).astype(np.float32)
    np.save(save_dir / "embedding_final.npy", embedding)
    np.save(save_dir / "labels.npy", labels)
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))

    fixed_pred = KMeans(n_clusters=true_k, n_init=20, random_state=args.seed).fit_predict(embedding)
    fixed_metrics, _ = family.compute_kmeans_metrics(labels, fixed_pred.astype(np.int64))

    adata = AnnData(X=embedding)
    sc.pp.neighbors(adata, n_neighbors=args.n_neighbors, use_rep="X", random_state=args.seed)
    rows = []
    for resolution in [float(x) for x in args.resolutions.split(",") if x.strip()]:
        key = f"leiden_{resolution:g}"
        sc.tl.leiden(adata, resolution=resolution, key_added=key, random_state=args.seed)
        pred = adata.obs[key].astype(int).to_numpy().astype(np.int64)
        metrics, _ = family.compute_kmeans_metrics(labels, pred)
        unsup = safe_unsup_scores(embedding, pred, args.silhouette_sample, args.seed)
        rows.append(
            {
                "dataset": dataset_name,
                "seed": int(args.seed),
                "resolution": float(resolution),
                "n_pred_clusters": int(len(np.unique(pred))),
                **metrics,
                **unsup,
            }
        )
    best_oracle = max(rows, key=lambda row: row["ari"])
    selected = {
        "silhouette": choose_label_blind(rows, "silhouette"),
        "davies_bouldin": choose_label_blind(rows, "davies_bouldin"),
        "calinski_harabasz": choose_label_blind(rows, "calinski_harabasz"),
    }
    pd.DataFrame(rows).to_csv(save_dir / "leiden_resolution_sweep.csv", index=False)
    summary = {
        "dataset": dataset_name,
        "seed": int(args.seed),
        "true_k": true_k,
        "method": "pca_leiden_kfree",
        "pca_dim": int(dim),
        "n_neighbors": int(args.n_neighbors),
        "runtime_seconds": float(time.time() - start),
        "known_k_kmeans": fixed_metrics,
        "oracle_best_leiden": best_oracle,
        "label_blind_selected": selected,
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(
        f"[RESULT] {dataset_name} seed={args.seed} knownK_ARI={fixed_metrics['ari']:.4f} "
        f"oracleLeiden_ARI={best_oracle['ari']:.4f}@res={best_oracle['resolution']} k={best_oracle['n_pred_clusters']} "
        f"silhouette_ARI={selected['silhouette']['ari']:.4f}@res={selected['silhouette']['resolution']} "
        f"k={selected['silhouette']['n_pred_clusters']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
