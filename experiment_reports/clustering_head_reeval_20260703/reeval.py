#!/usr/bin/env python3
"""Re-cluster EXISTING saved embeddings with different clustering heads.

Diagnostic question (see README): across the whole scMAE search the clustering
head was ALWAYS KMeans (euclidean, known k). Is the ceiling the embedding, or the
head? We re-score frozen `embedding_final.npy` with KMeans / GMM / Leiden without
retraining anything.

Leakage guardrail: Leiden has a free resolution. We select it by SILHOUETTE
(label-free), then report ARI at that selection = the honest number. We also
report the ARI-oracle over the sweep, explicitly labeled as a leakage upper bound
(NOT a usable result), only to show headroom.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys
import warnings
from pathlib import Path

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "4")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import numpy as np
from sklearn.cluster import KMeans
from sklearn.mixture import GaussianMixture
from sklearn.metrics import silhouette_score

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
from methods.DeepLearning import scMAE_family as family  # compute_kmeans_metrics

warnings.filterwarnings("ignore")


def label_free_silhouette(emb: np.ndarray, pred: np.ndarray, seed: int, sample: int = 5000) -> float:
    """Silhouette on a subsample (full pairwise on 44k is too costly). Label-free."""
    if len(np.unique(pred)) < 2:
        return float("nan")
    n = emb.shape[0]
    if n > sample:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n, size=sample, replace=False)
        emb, pred = emb[idx], pred[idx]
        if len(np.unique(pred)) < 2:
            return float("nan")
    try:
        return float(silhouette_score(emb, pred, metric="euclidean"))
    except Exception:
        return float("nan")


def run_kmeans(emb, k, seed):
    return KMeans(n_clusters=k, n_init=20, random_state=seed).fit_predict(emb).astype(np.int64)


def run_gmm(emb, k, seed):
    gm = GaussianMixture(n_components=k, covariance_type="diag", n_init=3,
                         max_iter=200, random_state=seed)
    return gm.fit_predict(emb).astype(np.int64)


def build_leiden_adata(emb, seed, n_neighbors=15):
    """Build the kNN graph ONCE; reuse across resolutions."""
    import scanpy as sc
    import anndata as ad
    a = ad.AnnData(np.ascontiguousarray(emb.astype(np.float32)))
    sc.pp.neighbors(a, n_neighbors=n_neighbors, use_rep="X", random_state=seed)
    return a


def leiden_labels(adata, resolution, seed):
    import scanpy as sc
    sc.tl.leiden(adata, resolution=float(resolution), random_state=seed, key_added="leiden")
    return adata.obs["leiden"].to_numpy().astype(np.int64)


def eval_pred(labels, pred, emb, seed):
    metrics, _ = family.compute_kmeans_metrics(labels, pred)
    metrics["silhouette"] = label_free_silhouette(emb, pred, seed)
    return metrics


def process_run(run_dir: Path, args) -> list[dict]:
    emb_path = run_dir / "embedding_final.npy"
    lab_path = run_dir / "labels.npy"
    if not emb_path.exists() or not lab_path.exists():
        return []
    emb = np.load(emb_path).astype(np.float32)
    emb = np.nan_to_num(emb)
    labels = np.load(lab_path).astype(np.int64)
    k = int(len(np.unique(labels)))
    rows = []
    base = {"run": run_dir.name, "n_cells": int(emb.shape[0]), "dim": int(emb.shape[1]), "known_k": k}

    # --- fixed-k heads ---
    for head, fn in [("kmeans_known_k", run_kmeans), ("gmm_known_k", run_gmm)]:
        try:
            pred = fn(emb, k, args.seed)
            m = eval_pred(labels, pred, emb, args.seed)
            rows.append({**base, "head": head, "selection": "known_k",
                         "resolution": "", "n_pred": m["n_pred_clusters"],
                         "ari": m["ari"], "nmi": m["nmi"], "acc": m["acc"], "silhouette": m["silhouette"]})
        except Exception as e:
            rows.append({**base, "head": head, "selection": "known_k", "resolution": "",
                         "n_pred": 0, "ari": None, "nmi": None, "acc": None, "silhouette": None,
                         "error": str(e)[:120]})

    # --- Leiden: sweep resolution, select by silhouette (label-free), also report ARI-oracle ---
    if not args.no_leiden:
        sweep = []
        adata = build_leiden_adata(emb, args.seed, args.n_neighbors)
        for res in args.resolutions:
            try:
                pred = leiden_labels(adata, res, args.seed)
                m = eval_pred(labels, pred, emb, args.seed)
                sweep.append({"resolution": res, "n_pred": m["n_pred_clusters"],
                              "ari": m["ari"], "nmi": m["nmi"], "acc": m["acc"], "silhouette": m["silhouette"]})
            except Exception as e:
                sweep.append({"resolution": res, "error": str(e)[:120], "silhouette": None, "ari": None})
        valid = [s for s in sweep if s.get("silhouette") is not None and not np.isnan(s.get("silhouette", np.nan))]
        if valid:
            sel = max(valid, key=lambda s: s["silhouette"])  # label-free selection
            rows.append({**base, "head": "leiden", "selection": "silhouette",
                         "resolution": sel["resolution"], "n_pred": sel["n_pred"],
                         "ari": sel["ari"], "nmi": sel["nmi"], "acc": sel["acc"], "silhouette": sel["silhouette"]})
        oracle_pool = [s for s in sweep if s.get("ari") is not None]
        if oracle_pool:
            orc = max(oracle_pool, key=lambda s: s["ari"])  # LEAKAGE upper bound
            rows.append({**base, "head": "leiden", "selection": "ARI_ORACLE_leak",
                         "resolution": orc["resolution"], "n_pred": orc["n_pred"],
                         "ari": orc["ari"], "nmi": orc["nmi"], "acc": orc["acc"], "silhouette": orc.get("silhouette")})
        (run_dir_out := args.out / run_dir.name).mkdir(parents=True, exist_ok=True)
        with open(run_dir_out / "leiden_sweep.json", "w") as h:
            json.dump(sweep, h, indent=2)
    print(f"[done] {run_dir.name}: " + ", ".join(
        f"{r['head']}/{r['selection']}={r['ari']:.4f}" for r in rows if r.get("ari") is not None))
    return rows


def main() -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--run_dirs", nargs="+", required=True,
                   help="glob(s) of run directories containing embedding_final.npy + labels.npy")
    p.add_argument("--out", type=Path,
                   default=ROOT / "experiment_reports" / "clustering_head_reeval_20260703")
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--n_neighbors", type=int, default=15)
    p.add_argument("--resolutions", type=float, nargs="+",
                   default=[0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.8, 1.2])
    p.add_argument("--no_leiden", action="store_true")
    args = p.parse_args()

    import glob
    dirs = []
    for pattern in args.run_dirs:
        dirs.extend(sorted(Path(x) for x in glob.glob(pattern)))
    dirs = [d for d in dirs if d.is_dir()]
    if not dirs:
        print("No run dirs matched.")
        return 1
    print(f"Re-evaluating {len(dirs)} run(s).")

    all_rows = []
    for d in dirs:
        all_rows.extend(process_run(d, args))

    args.out.mkdir(parents=True, exist_ok=True)
    csv_path = args.out / "reeval.csv"
    fields = ["run", "head", "selection", "resolution", "n_cells", "dim", "known_k",
              "n_pred", "ari", "nmi", "acc", "silhouette"]
    with open(csv_path, "w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(all_rows)
    print(f"Wrote {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
