#!/usr/bin/env python
import argparse
import os

import h5py
import numpy as np
import scanpy as sc

from common import ensure_dir, evaluate_embedding, labels_from_adata, save_json


def parse_args():
    parser = argparse.ArgumentParser(description="Unified KMeans/Leiden evaluation for an embedding.")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--embedding", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--out_dir", required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--n_neighbors", type=int, default=15)
    return parser.parse_args()


def load_embedding(path):
    if path.endswith(".h5") or path.endswith(".hdf5"):
        with h5py.File(path, "r") as handle:
            return np.asarray(handle["X"])
    return np.load(path)


def main():
    args = parse_args()
    ensure_dir(args.out_dir)
    adata = sc.read_h5ad(args.data_path)
    labels, label_key = labels_from_adata(adata)
    embedding = load_embedding(args.embedding)
    if embedding.shape[0] != labels.shape[0]:
        raise ValueError(f"Embedding cells {embedding.shape[0]} != labels {labels.shape[0]}")

    metrics, preds = evaluate_embedding(
        embedding,
        labels,
        n_clusters=len(np.unique(labels)),
        seed=args.seed,
        n_neighbors=args.n_neighbors,
    )
    payload = {
        "dataset": args.dataset,
        "method": args.method,
        "embedding_path": os.path.abspath(args.embedding),
        "label_key": label_key,
        "n_cells": int(embedding.shape[0]),
        "embedding_dim": int(embedding.shape[1]),
        "metrics": metrics,
    }
    save_json(payload, os.path.join(args.out_dir, f"{args.method}.json"))
    rows = []
    for cluster_method, vals in metrics.items():
        row = {
            "dataset": args.dataset,
            "method": args.method,
            "cluster_method": cluster_method,
        }
        row.update(vals)
        rows.append(row)
    import pandas as pd
    pd.DataFrame(rows).to_csv(os.path.join(args.out_dir, f"{args.method}.csv"), index=False)
    for name, pred in preds.items():
        np.save(os.path.join(args.out_dir, f"{args.method}_{name}.npy"), pred)
    print(payload)


if __name__ == "__main__":
    main()

