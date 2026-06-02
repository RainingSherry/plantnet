#!/usr/bin/env python
import argparse
import glob
import importlib.util
import os
import subprocess
import sys
import time
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mpl_config")
os.environ.setdefault("SCBENCH_N_TOP_GENES", "2000")
os.environ.setdefault("PYTHONUNBUFFERED", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "32")
os.environ.setdefault("OMP_NUM_THREADS", "32")

import h5py
import networkx as nx
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import coo_matrix
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import kneighbors_graph

from common import compute_metrics, ensure_dir, evaluate_embedding, labels_from_adata, leiden_labels, save_json


DATASETS = {
    "Mouse_Pancreas_1": ("benchmarks/unified_protocol/preprocessed/Mouse_Pancreas_1_hvg2000.h5ad", 13),
    "SRP182008": ("benchmarks/unified_protocol/preprocessed/SRP182008_hvg2000.h5ad", 15),
    "SRP235541": ("data/SRP235541.h5ad", 18),
    "SRP171040": ("benchmarks/unified_protocol/preprocessed/SRP171040_hvg2000.h5ad", 12),
}

_PCA_CACHE = {}
METHOD_DEPS = {
    "desc": ("tensorflow",),
    "scDeepCluster": ("tensorflow", "keras"),
    "scNAME": ("tensorflow", "keras"),
    "scziDesk": ("tensorflow", "keras"),
}


def parse_args():
    parser = argparse.ArgumentParser(description="Run unified PlantNet benchmark jobs.")
    parser.add_argument("--datasets", default="Mouse_Pancreas_1,SRP182008,SRP171040")
    parser.add_argument(
        "--methods",
        default=(
            "traditional_pca,traditional_leiden,traditional_louvain,traditional_sc3,"
            "scMAE,scVI,PhytoCluster,scCDCG,plantspade_lgcl,"
            "dec,scDCC,desc,scDeepCluster,scNAME,scziDesk,"
            "scDSC,AttentionAE-sc,scGNN"
        ),
    )
    parser.add_argument("--root", default="benchmarks/unified_protocol")
    parser.add_argument("--gpus", default="1,2,3,4,5,6")
    parser.add_argument("--max_parallel", type=int, default=6)
    parser.add_argument("--deep_epochs", type=int, default=100)
    parser.add_argument("--phyto_pretrain_iter", type=int, default=8000)
    parser.add_argument("--phyto_cluster_iter", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--skip_existing", action="store_true")
    parser.add_argument("--force_eval_existing", action="store_true")
    parser.add_argument("--fail_on_error", action="store_true")
    return parser.parse_args()


def py():
    return sys.executable


def missing_dependencies(method):
    return [dep for dep in METHOD_DEPS.get(method, ()) if importlib.util.find_spec(dep) is None]


def build_command(method, data_path, save_dir, n_clusters, gpu, args):
    standard_h5_methods = {
        "dec": ("methods/DeepLearning/dec/run.py", ["--pretrain_epochs", str(args.deep_epochs)]),
        "scDCC": ("methods/DeepLearning/scDCC/run.py", ["--pretrain_epochs", str(args.deep_epochs)]),
        "desc": ("methods/DeepLearning/desc/run.py", ["--pretrain_epochs", str(args.deep_epochs), "--max_iter", str(args.deep_epochs)]),
        "scDeepCluster": ("methods/DeepLearning/scDeepCluster/run.py", ["--pretrain_epochs", str(args.deep_epochs), "--maxiter", str(args.deep_epochs)]),
        "scNAME": ("methods/DeepLearning/scNAME/run.py", ["--pretrain_epochs", str(args.deep_epochs)]),
        "scziDesk": ("methods/DeepLearning/scziDesk/run.py", ["--pretrain_epochs", str(args.deep_epochs)]),
    }
    if method in standard_h5_methods:
        script, extra = standard_h5_methods[method]
        return [
            py(), script,
            "--data_path", data_path, "--save_dir", save_dir, "--n_clusters", str(n_clusters),
            "--epochs", str(args.deep_epochs), "--batch_size", "256", "--gpu", str(gpu),
            "--seed", str(args.seed), *extra,
        ]
    if method in {"scDSC", "AttentionAE-sc"}:
        script = {
            "scDSC": "methods/GNN/scDSC/run.py",
            "AttentionAE-sc": "methods/GNN/AttentionAE-sc/run.py",
        }[method]
        return [
            py(), script,
            "--data_path", data_path, "--save_dir", save_dir, "--n_clusters", str(n_clusters),
            "--epochs", str(args.deep_epochs), "--pretrain_epochs", str(args.deep_epochs),
            "--gpu", str(gpu), "--seed", str(args.seed),
        ]
    if method == "scGNN":
        return [
            py(), "methods/GNN/scGNN/run.py",
            "--data_path", data_path, "--save_dir", save_dir, "--n_clusters", str(n_clusters),
            "--Regu_epochs", str(args.deep_epochs), "--EM_epochs", str(args.deep_epochs),
            "--cluster_epochs", str(args.deep_epochs), "--EM_iteration", "1",
            "--quickmode", "--batch_size", "4096", "--cores_usage", "1",
            "--gpu", str(gpu), "--seed", str(args.seed),
        ]
    if method == "scMAE":
        return [
            py(), "methods/DeepLearning/scMAE/run.py",
            "--data_path", data_path, "--save_dir", save_dir, "--n_clusters", str(n_clusters),
            "--epochs", str(args.deep_epochs), "--batch_size", "256", "--gpu", str(gpu),
            "--seed", str(args.seed), "--eval_interval", str(args.deep_epochs),
        ]
    if method == "scCDCG":
        return [
            py(), "methods/GNN/scCDCG/run.py",
            "--data_path", data_path, "--save_dir", save_dir, "--n_clusters", str(n_clusters),
            "--epochs", str(args.deep_epochs), "--gpu", str(gpu), "--seed", str(args.seed),
        ]
    if method == "plantspade_lgcl":
        return [
            py(), "methods/DeepLearning/PlantSPADE_LGCL/run_plantspade.py",
            "--data_path", data_path, "--save_dir", save_dir, "--input_mode", "auto",
            "--n_top_genes", "2000", "--n_clusters", str(n_clusters),
            "--latent_dim", "32", "--layers", "2", "--epochs", str(args.deep_epochs),
            "--pairs_per_epoch", "262144", "--contrastive_batch_size", "2048",
            "--lr", "1e-3", "--weight_decay", "1e-5", "--edge_dropout", "0.1",
            "--contrastive_weight", "0.05", "--module_weight", "0.001",
            "--num_modules", "16", "--module_top_k", "30",
            "--gpu", str(gpu), "--seed", str(args.seed),
        ]
    if method == "PhytoCluster":
        return [
            py(), "methods/DeepLearning/PhytoCluster/run.py",
            "--data_path", data_path, "--save_dir", save_dir, "--n_clusters", str(n_clusters),
            "--n_top_genes", "2000", "--pretrain_max_iter", str(args.phyto_pretrain_iter),
            "--cluster_max_iter", str(args.phyto_cluster_iter), "--batch_size", "256", "--gpu", str(gpu),
            "--seed", str(args.seed),
        ]
    if method == "scVI":
        return [
            py(), "methods/DeepLearning/scVI/run.py",
            "--data_path", data_path, "--save_dir", save_dir, "--n_clusters", str(n_clusters),
            "--n_top_genes", "2000", "--epochs", str(args.deep_epochs),
            "--batch_size", "256", "--gpu", str(gpu), "--seed", str(args.seed),
        ]
    raise ValueError(f"Unknown subprocess method: {method}")


def load_embedding(path):
    if path.endswith(".h5") or path.endswith(".hdf5"):
        with h5py.File(path, "r") as handle:
            return np.asarray(handle["X"])
    return np.load(path)


def find_embeddings(method, save_dir):
    candidates = []
    if method in {
        "scMAE",
        "scCDCG",
        "PhytoCluster",
        "dec",
        "scDCC",
        "desc",
        "scDeepCluster",
        "scNAME",
        "scziDesk",
        "scDSC",
        "AttentionAE-sc",
        "scGNN",
    }:
        candidates.append((method, os.path.join(save_dir, "embedding.h5")))
    elif method == "plantspade_lgcl":
        candidates.append(("plantspade_lgcl", os.path.join(save_dir, "embedding_final.npy")))
    elif method == "scVI":
        hits = sorted(glob.glob(os.path.join(save_dir, "**", "embedding.h5"), recursive=True))
        if hits:
            candidates.append((method, hits[-1]))
    return [(name, path) for name, path in candidates if os.path.exists(path)]


def evaluate_and_write(dataset, method_name, embedding_path, data_path, out_dir, seed):
    adata = sc.read_h5ad(data_path)
    labels, label_key = labels_from_adata(adata)
    embedding = load_embedding(embedding_path)
    if embedding.shape[0] != labels.shape[0]:
        raise ValueError(f"{method_name}: embedding cells {embedding.shape[0]} != labels {labels.shape[0]}")
    metrics, preds = evaluate_embedding(embedding, labels, n_clusters=len(np.unique(labels)), seed=seed)
    ensure_dir(out_dir)
    payload = {
        "dataset": dataset,
        "method": method_name,
        "embedding_path": os.path.abspath(embedding_path),
        "label_key": label_key,
        "metrics": metrics,
    }
    save_json(payload, os.path.join(out_dir, f"{method_name}.json"))
    rows = []
    for cluster_method, vals in metrics.items():
        row = {"dataset": dataset, "method": method_name, "cluster_method": cluster_method}
        row.update(vals)
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, f"{method_name}.csv"), index=False)
    for pred_name, pred in preds.items():
        np.save(os.path.join(out_dir, f"{method_name}_{pred_name}.npy"), pred)
    return rows


def get_pca_embedding(data_path, seed, n_components=50):
    cache_key = (data_path, seed, n_components)
    if cache_key in _PCA_CACHE:
        return _PCA_CACHE[cache_key]
    adata = sc.read_h5ad(data_path)
    x = adata.X.toarray() if hasattr(adata.X, "toarray") else np.asarray(adata.X)
    n_components = min(n_components, x.shape[0] - 1, x.shape[1])
    emb = PCA(n_components=n_components, random_state=seed).fit_transform(x)
    _PCA_CACHE[cache_key] = emb
    return emb


def run_pca(dataset, data_path, out_dir, seed, method_name="traditional_pca"):
    emb = get_pca_embedding(data_path, seed)
    emb_path = os.path.join(out_dir, "pca_embedding.npy")
    ensure_dir(out_dir)
    np.save(emb_path, emb.astype(np.float32))
    return evaluate_and_write(dataset, method_name, emb_path, data_path, out_dir, seed)


def evaluate_labels_and_write(dataset, method_name, cluster_method, y_pred, data_path, out_dir, seed, embedding=None, extra=None):
    adata = sc.read_h5ad(data_path)
    labels, label_key = labels_from_adata(adata)
    if y_pred.shape[0] != labels.shape[0]:
        raise ValueError(f"{method_name}: pred cells {y_pred.shape[0]} != labels {labels.shape[0]}")
    metrics, mapped = compute_metrics(labels, y_pred, embedding=embedding)
    if extra:
        metrics.update(extra)
    ensure_dir(out_dir)
    payload = {
        "dataset": dataset,
        "method": method_name,
        "cluster_method": cluster_method,
        "label_key": label_key,
        "metrics": metrics,
    }
    save_json(payload, os.path.join(out_dir, f"{method_name}.json"))
    row = {"dataset": dataset, "method": method_name, "cluster_method": cluster_method}
    row.update(metrics)
    pd.DataFrame([row]).to_csv(os.path.join(out_dir, f"{method_name}.csv"), index=False)
    np.save(os.path.join(out_dir, f"{method_name}_{cluster_method}.npy"), y_pred.astype(np.int64))
    np.save(os.path.join(out_dir, f"{method_name}_{cluster_method}_mapped.npy"), mapped.astype(np.int64))
    return [row]


def run_traditional_leiden(dataset, data_path, out_dir, n_clusters, seed):
    emb = get_pca_embedding(data_path, seed)
    pca_leiden_path = os.path.join(out_dir, "traditional_pca_leiden.npy")
    if os.path.exists(pca_leiden_path):
        pred = np.load(pca_leiden_path)
        res = float("nan")
    else:
        pred, res = leiden_labels(emb, n_clusters, seed=seed)
    return evaluate_labels_and_write(
        dataset,
        "traditional_leiden",
        "leiden",
        pred.astype(np.int64),
        data_path,
        out_dir,
        seed,
        embedding=emb,
        extra={"resolution": float(res)},
    )


def louvain_labels(embedding, n_clusters, seed=42, n_neighbors=15):
    graph = kneighbors_graph(
        np.asarray(embedding, dtype=np.float32),
        n_neighbors=n_neighbors,
        mode="distance",
        include_self=False,
    )
    graph = coo_matrix(graph)
    nx_graph = nx.Graph()
    nx_graph.add_nodes_from(range(embedding.shape[0]))
    for row, col, dist in zip(graph.row, graph.col, graph.data):
        if row != col:
            nx_graph.add_edge(int(row), int(col), weight=float(1.0 / (1.0 + dist)))

    best_pred = None
    best_res = 1.0
    best_diff = float("inf")
    for res in [1.0]:
        communities = nx.algorithms.community.louvain_communities(
            nx_graph,
            resolution=float(res),
            seed=seed,
            weight="weight",
        )
        pred = np.zeros(embedding.shape[0], dtype=np.int64)
        for cluster_id, nodes in enumerate(communities):
            pred[list(nodes)] = cluster_id
        diff = abs(len(communities) - n_clusters)
        if diff < best_diff:
            best_diff = diff
            best_res = float(res)
            best_pred = pred
        if diff == 0:
            break
    return best_pred, best_res


def run_traditional_louvain(dataset, data_path, out_dir, n_clusters, seed, max_cells=8000):
    emb = get_pca_embedding(data_path, seed)
    if emb.shape[0] > max_cells:
        row = {
            "dataset": dataset,
            "method": "traditional_louvain",
            "cluster_method": "louvain",
            "status": f"skipped_n_cells_gt_{max_cells}",
            "n_cells": int(emb.shape[0]),
        }
        ensure_dir(out_dir)
        save_json(row, os.path.join(out_dir, "traditional_louvain.SKIPPED.json"))
        pd.DataFrame([row]).to_csv(os.path.join(out_dir, "traditional_louvain.csv"), index=False)
        return [row]
    pred, res = louvain_labels(emb, n_clusters, seed=seed)
    return evaluate_labels_and_write(
        dataset,
        "traditional_louvain",
        "louvain",
        pred.astype(np.int64),
        data_path,
        out_dir,
        seed,
        embedding=emb,
        extra={"resolution": float(res)},
    )


def run_traditional_sc3(dataset, data_path, out_dir, n_clusters, seed, max_cells=8000):
    emb = get_pca_embedding(data_path, seed, n_components=20)
    if emb.shape[0] > max_cells:
        row = {
            "dataset": dataset,
            "method": "traditional_sc3",
            "cluster_method": "sc3_consensus",
            "status": f"skipped_n_cells_gt_{max_cells}",
            "n_cells": int(emb.shape[0]),
        }
        ensure_dir(out_dir)
        save_json(row, os.path.join(out_dir, "traditional_sc3.SKIPPED.json"))
        pd.DataFrame([row]).to_csv(os.path.join(out_dir, "traditional_sc3.csv"), index=False)
        return [row]

    n_cells = emb.shape[0]
    consensus = np.zeros((n_cells, n_cells), dtype=np.float32)
    configs = [
        (n_clusters, 20, seed),
        (n_clusters, 10, seed + 1),
        (n_clusters, 10, seed + 2),
        (n_clusters + 1, 20, seed + 3),
        (max(2, n_clusters - 1), 20, seed + 4),
    ]
    used = 0
    for k, n_init, random_state in configs:
        if k < 2 or k > n_cells:
            continue
        labels = KMeans(n_clusters=k, n_init=n_init, random_state=random_state).fit_predict(emb)
        for label in np.unique(labels):
            idx = np.flatnonzero(labels == label)
            consensus[np.ix_(idx, idx)] += 1.0
        used += 1
    consensus /= max(1, used)
    np.fill_diagonal(consensus, 1.0)
    dist_matrix = np.nan_to_num(1.0 - consensus, nan=1.0)
    try:
        pred = AgglomerativeClustering(
            n_clusters=n_clusters,
            metric="precomputed",
            linkage="average",
        ).fit_predict(dist_matrix)
    except TypeError:
        pred = AgglomerativeClustering(
            n_clusters=n_clusters,
            affinity="precomputed",
            linkage="average",
        ).fit_predict(dist_matrix)
    return evaluate_labels_and_write(
        dataset,
        "traditional_sc3",
        "sc3_consensus",
        pred.astype(np.int64),
        data_path,
        out_dir,
        seed,
        embedding=emb,
    )


def launch_job(job, env):
    ensure_dir(job["save_dir"])
    log_path = os.path.join(job["save_dir"], "train.log")
    log = open(log_path, "w", encoding="utf-8")
    print("LAUNCH", job["dataset"], job["method"], "gpu", job["gpu"], "log", log_path)
    proc = subprocess.Popen(job["cmd"], cwd=job["cwd"], env=env, stdout=log, stderr=subprocess.STDOUT)
    job["proc"] = proc
    job["log_handle"] = log
    job["log_path"] = log_path
    job["start"] = time.time()
    return job


def main():
    args = parse_args()
    root = Path(args.root)
    runs_root = root / "runs"
    eval_root = root / "eval"
    ensure_dir(runs_root)
    ensure_dir(eval_root)

    datasets = [item.strip() for item in args.datasets.split(",") if item.strip()]
    methods = [item.strip() for item in args.methods.split(",") if item.strip()]
    gpus = [int(item) for item in args.gpus.split(",") if item.strip()]
    if any(gpu in {0, 7} for gpu in gpus):
        raise ValueError("GPU 0 and GPU 7 are forbidden by request.")

    rows_all = []
    jobs = []
    direct_methods = {"pca", "traditional_pca", "traditional_leiden", "traditional_louvain", "traditional_sc3"}
    for dataset in datasets:
        data_path, n_clusters = DATASETS[dataset]
        for method in methods:
            if method in {"pca", "traditional_pca"}:
                rows_all.extend(run_pca(dataset, data_path, str(eval_root / dataset), args.seed, method_name="traditional_pca"))
                continue
            if method == "traditional_leiden":
                rows_all.extend(run_traditional_leiden(dataset, data_path, str(eval_root / dataset), n_clusters, args.seed))
                continue
            if method == "traditional_louvain":
                rows_all.extend(run_traditional_louvain(dataset, data_path, str(eval_root / dataset), n_clusters, args.seed))
                continue
            if method == "traditional_sc3":
                rows_all.extend(run_traditional_sc3(dataset, data_path, str(eval_root / dataset), n_clusters, args.seed))
                continue
            if method in direct_methods:
                continue
            save_dir = str(runs_root / dataset / method)
            missing = missing_dependencies(method)
            if missing:
                ensure_dir(save_dir)
                payload = {
                    "dataset": dataset,
                    "method": method,
                    "status": "skipped_missing_dependencies",
                    "missing_dependencies": missing,
                }
                save_json(payload, os.path.join(save_dir, "SKIPPED.json"))
                rows_all.append({
                    "dataset": dataset,
                    "method": method,
                    "cluster_method": "skipped",
                    "status": "missing_dependencies:" + ",".join(missing),
                })
                continue
            done_marker = os.path.join(save_dir, "UNIFIED_DONE")
            if args.skip_existing and os.path.exists(done_marker):
                dataset_eval_dir = str(eval_root / dataset)
                for name, emb_path in find_embeddings(method, save_dir):
                    csv_path = os.path.join(dataset_eval_dir, f"{name}.csv")
                    if os.path.exists(csv_path) and not args.force_eval_existing:
                        rows_all.extend(pd.read_csv(csv_path).to_dict("records"))
                    else:
                        if args.skip_existing and not args.force_eval_existing:
                            rows_all.append({
                                "dataset": dataset,
                                "method": name,
                                "cluster_method": "eval_missing",
                                "status": "trained_embedding_exists_but_eval_csv_missing",
                                "embedding_path": os.path.abspath(emb_path),
                            })
                        else:
                            rows_all.extend(evaluate_and_write(dataset, name, emb_path, data_path, dataset_eval_dir, args.seed))
                continue
            jobs.append({
                "dataset": dataset,
                "method": method,
                "data_path": data_path,
                "save_dir": save_dir,
                "n_clusters": n_clusters,
                "cwd": os.getcwd(),
            })

    env = os.environ.copy()
    env["NUMBA_DISABLE_JIT"] = "1"
    env["MPLCONFIGDIR"] = "/tmp/mpl_config"
    env["SCBENCH_N_TOP_GENES"] = "2000"
    env["PYTHONUNBUFFERED"] = "1"
    env["OPENBLAS_NUM_THREADS"] = "32"
    env["OMP_NUM_THREADS"] = "32"

    pending = jobs[:]
    running = []
    failed = []
    while pending or running:
        while pending and len(running) < args.max_parallel:
            used_gpus = {job["gpu"] for job in running if "gpu" in job}
            free_gpus = [gpu for gpu in gpus if gpu not in used_gpus]
            if not free_gpus:
                break
            job = pending.pop(0)
            job["gpu"] = free_gpus[0]
            job["cmd"] = build_command(
                job["method"],
                job["data_path"],
                job["save_dir"],
                job["n_clusters"],
                job["gpu"],
                args,
            )
            running.append(launch_job(job, env))
        time.sleep(5)
        still = []
        for job in running:
            ret = job["proc"].poll()
            if ret is None:
                still.append(job)
                continue
            job["log_handle"].close()
            elapsed = time.time() - job["start"]
            print("FINISH", job["dataset"], job["method"], "ret", ret, f"{elapsed/60:.1f} min")
            if ret != 0:
                failed.append(job)
                save_json({"returncode": ret, "log": job["log_path"], "cmd": job["cmd"]}, os.path.join(job["save_dir"], "FAILED.json"))
                rows_all.append({
                    "dataset": job["dataset"],
                    "method": job["method"],
                    "cluster_method": "failed",
                    "status": "failed",
                    "log": job["log_path"],
                })
                continue
            Path(os.path.join(job["save_dir"], "UNIFIED_DONE")).write_text("ok\n", encoding="utf-8")
            try:
                for name, emb_path in find_embeddings(job["method"], job["save_dir"]):
                    rows_all.extend(
                        evaluate_and_write(job["dataset"], name, emb_path, job["data_path"], str(eval_root / job["dataset"]), args.seed)
                    )
            except Exception as exc:
                failed.append(job)
                save_json({"error": str(exc), "log": job["log_path"], "cmd": job["cmd"]}, os.path.join(job["save_dir"], "EVAL_FAILED.json"))
                rows_all.append({
                    "dataset": job["dataset"],
                    "method": job["method"],
                    "cluster_method": "eval_failed",
                    "status": str(exc),
                    "log": job["log_path"],
                })
        running = still

    if rows_all:
        out = pd.DataFrame(rows_all)
        out = out.drop_duplicates(subset=["dataset", "method", "cluster_method"], keep="last")
        out = out.sort_values(["dataset", "method", "cluster_method"])
        out.to_csv(root / "metrics_long.csv", index=False)
        print(out.to_string(index=False))
    if failed:
        print("FAILED JOBS:")
        for job in failed:
            print(job["dataset"], job["method"], job.get("log_path"))
        if args.fail_on_error:
            raise SystemExit(1)


if __name__ == "__main__":
    main()
