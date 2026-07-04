#!/usr/bin/env python
"""cuML (GPU) vs scikit-learn (CPU) benchmark on classic single-cell steps.

Steps benchmarked on the HVG-1000 expression matrix X (dense float32):
  1. PCA        (n_components=50)
  2. KMeans     (n_clusters=15, n_init=10)
  3. KNN graph  (n_neighbors=15, euclidean; fit + kneighbors query)

Design
------
Each (dataset, step, backend) unit runs in its OWN subprocess (worker mode),
launched by the orchestrator (--run-all). This gives:
  * clean GPU memory per unit (no cuML allocation carry-over),
  * a hard per-unit wall-clock timeout for the slow CPU side (>600s => infeasible),
  * identical, isolated timing conditions.

CPU threads are bound to 48 via OMP/OPENBLAS/MKL/NUMEXPR env vars, set in the
worker environment BEFORE numpy/sklearn import.

Usage
-----
  # single unit (prints one JSON line to stdout):
  python bench_cuml.py --worker --file F --step pca --backend sklearn
  CUDA_VISIBLE_DEVICES=1 python bench_cuml.py --worker --file F --step pca --backend cuml

  # full sweep (orchestrates all subprocesses, writes csv + markdown):
  python bench_cuml.py --run-all
"""
import argparse
import json
import os
import subprocess
import sys
import time

# ----------------------------------------------------------------------------
# Config
# ----------------------------------------------------------------------------
CPU_THREADS = 48
GPU_ID = "1"                 # only 1/2/3 usable; pick 1
CPU_TIMEOUT_S = 600          # sklearn per-step budget; over this => "CPU>600s"
GPU_TIMEOUT_S = 1200         # generous ceiling for GPU worker
PYTHON = sys.executable      # rapids_bench python (both cuml + sklearn live here)

PCA_COMPONENTS = 50
KMEANS_CLUSTERS = 15
KMEANS_NINIT = 10
KNN_NEIGHBORS = 15
KNN_METRIC = "euclidean"

DATA_DIR = "/home/luolie/.claude/jobs/23ea9258/tmp/scale_bench/hvg1000"
DATASETS = [
    ("real_10k",          f"{DATA_DIR}/real_10k_hvg1000.h5ad"),
    ("real_50k",          f"{DATA_DIR}/real_50k_hvg1000.h5ad"),
    ("real_86k",          f"{DATA_DIR}/real_86k_hvg1000.h5ad"),
    ("real_200k_xspecies", f"{DATA_DIR}/real_200k_hvg1000_xspecies.h5ad"),
]
STEPS = ["pca", "kmeans", "knn"]
BACKENDS = ["sklearn", "cuml"]


# ----------------------------------------------------------------------------
# Data loading (worker side)
# ----------------------------------------------------------------------------
def _register_null_h5ad_reader():
    """Tolerate uns entries stored with encoding-type 'null' (anndata can't)."""
    try:
        from anndata._io.specs.registry import _REGISTRY, IOSpec
        import h5py

        def _read_null(*a, **k):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


def load_X(path):
    """Return X as a C-contiguous dense float32 ndarray, plus (n, d)."""
    import numpy as np
    import scipy.sparse as sp

    _register_null_h5ad_reader()
    import anndata as ad

    adata = ad.read_h5ad(path)
    X = adata.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.ascontiguousarray(np.asarray(X, dtype=np.float32))
    return X, X.shape


# ----------------------------------------------------------------------------
# Timing primitives (worker side)
# ----------------------------------------------------------------------------
def _time_once(fn):
    t0 = time.perf_counter()
    fn()
    return time.perf_counter() - t0


def bench_sklearn(step, X):
    """Return wall-time (s) for one sklearn step. No warm-up needed (CPU)."""
    if step == "pca":
        from sklearn.decomposition import PCA

        def run():
            PCA(n_components=PCA_COMPONENTS, svd_solver="auto",
                random_state=0).fit_transform(X)
        return _time_once(run)

    if step == "kmeans":
        from sklearn.cluster import KMeans

        def run():
            KMeans(n_clusters=KMEANS_CLUSTERS, n_init=KMEANS_NINIT,
                   random_state=0).fit(X)
        return _time_once(run)

    if step == "knn":
        from sklearn.neighbors import NearestNeighbors

        def run():
            nn = NearestNeighbors(n_neighbors=KNN_NEIGHBORS, metric=KNN_METRIC,
                                  algorithm="auto", n_jobs=-1)
            nn.fit(X)
            nn.kneighbors(X)  # build the graph (query all points)
        return _time_once(run)

    raise ValueError(step)


def bench_cuml(step, X):
    """Return wall-time (s) for one cuML step, timed AFTER one warm-up call.

    Warm-up absorbs CUDA context / JIT / allocator init so the timed run
    reflects steady-state compute.
    """
    import cupy as cp

    Xg = cp.asarray(X)  # move to device once; transfer excluded from timing

    if step == "pca":
        from cuml import PCA

        def run():
            PCA(n_components=PCA_COMPONENTS).fit_transform(Xg)
            cp.cuda.Stream.null.synchronize()

    elif step == "kmeans":
        from cuml.cluster import KMeans

        def run():
            KMeans(n_clusters=KMEANS_CLUSTERS, n_init=KMEANS_NINIT,
                   random_state=0).fit(Xg)
            cp.cuda.Stream.null.synchronize()

    elif step == "knn":
        from cuml.neighbors import NearestNeighbors

        def run():
            nn = NearestNeighbors(n_neighbors=KNN_NEIGHBORS, metric=KNN_METRIC)
            nn.fit(Xg)
            nn.kneighbors(Xg)
            cp.cuda.Stream.null.synchronize()
    else:
        raise ValueError(step)

    run()          # warm-up (untimed)
    return _time_once(run)


def worker_main(args):
    X, (n, d) = load_X(args.file)
    if args.backend == "sklearn":
        secs = bench_sklearn(args.step, X)
    else:
        secs = bench_cuml(args.step, X)
    print(json.dumps({"step": args.step, "backend": args.backend,
                      "n": int(n), "d": int(d), "seconds": secs}))


# ----------------------------------------------------------------------------
# Orchestrator
# ----------------------------------------------------------------------------
def _worker_env(backend):
    env = dict(os.environ)
    env["OMP_NUM_THREADS"] = str(CPU_THREADS)
    env["OPENBLAS_NUM_THREADS"] = str(CPU_THREADS)
    env["MKL_NUM_THREADS"] = str(CPU_THREADS)
    env["NUMEXPR_NUM_THREADS"] = str(CPU_THREADS)
    if backend == "cuml":
        env["CUDA_VISIBLE_DEVICES"] = GPU_ID
    else:
        # keep sklearn honest & off the GPU
        env["CUDA_VISIBLE_DEVICES"] = ""
    return env


def run_unit(path, step, backend):
    """Launch one worker subprocess; return (seconds|None, status)."""
    timeout = CPU_TIMEOUT_S if backend == "sklearn" else GPU_TIMEOUT_S
    cmd = [PYTHON, os.path.abspath(__file__), "--worker",
           "--file", path, "--step", step, "--backend", backend]
    try:
        p = subprocess.run(cmd, env=_worker_env(backend), capture_output=True,
                           text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, f"{backend.upper()}>{timeout}s"
    if p.returncode != 0:
        tail = (p.stderr or "").strip().splitlines()
        tail = tail[-1] if tail else "unknown error"
        return None, f"ERROR: {tail}"
    try:
        rec = json.loads(p.stdout.strip().splitlines()[-1])
        return rec["seconds"], "ok"
    except Exception as e:
        return None, f"PARSE_ERR: {e}"


def run_all():
    out_dir = os.path.dirname(os.path.abspath(__file__))
    results = []
    for ds_name, path in DATASETS:
        if not os.path.exists(path):
            print(f"[skip] missing {path}", file=sys.stderr)
            continue
        for step in STEPS:
            row = {"dataset": ds_name, "step": step}
            for backend in BACKENDS:
                secs, status = run_unit(path, step, backend)
                row[f"{backend}_s"] = secs
                row[f"{backend}_status"] = status
                tag = f"{secs:.3f}s" if secs is not None else status
                print(f"[{ds_name:20s}] {step:6s} {backend:7s} -> {tag}",
                      flush=True)
            cpu, gpu = row.get("sklearn_s"), row.get("cuml_s")
            row["speedup"] = (cpu / gpu) if (cpu and gpu) else None
            results.append(row)

    csv_path = os.path.join(out_dir, "bench_cuml_results.csv")
    with open(csv_path, "w") as f:
        f.write("dataset,step,cpu_sklearn_s,gpu_cuml_s,speedup,cpu_status,gpu_status\n")
        for r in results:
            cpu = f"{r['sklearn_s']:.3f}" if r.get("sklearn_s") else ""
            gpu = f"{r['cuml_s']:.3f}" if r.get("cuml_s") else ""
            sp = f"{r['speedup']:.1f}" if r.get("speedup") else ""
            f.write(f"{r['dataset']},{r['step']},{cpu},{gpu},{sp},"
                    f"{r['sklearn_status']},{r['cuml_status']}\n")
    print(f"\n[done] wrote {csv_path}")
    print(json.dumps(results, indent=2))
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--worker", action="store_true")
    ap.add_argument("--run-all", action="store_true")
    ap.add_argument("--file")
    ap.add_argument("--step", choices=STEPS)
    ap.add_argument("--backend", choices=BACKENDS)
    args = ap.parse_args()
    if args.worker:
        worker_main(args)
    elif args.run_all:
        run_all()
    else:
        ap.error("pass --worker or --run-all")


if __name__ == "__main__":
    main()
