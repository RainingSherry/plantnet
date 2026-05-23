#!/usr/bin/env python
import argparse
import os

import numpy as np
import scanpy as sc
import scipy.sparse as sp

from common import ensure_dir, find_label_key, to_dense


def parse_args():
    parser = argparse.ArgumentParser(description="Create unified HVG2000 h5ad files.")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--out_path", required=True)
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--target_sum", type=float, default=1e4)
    return parser.parse_args()


def main():
    args = parse_args()
    ensure_dir(os.path.dirname(args.out_path))

    adata = sc.read_h5ad(args.data_path)
    label_key = find_label_key(adata)

    counts = adata.raw.X.copy() if adata.raw is not None else adata.X.copy()
    counts = counts.tocsr().astype(np.float32) if sp.issparse(counts) else np.asarray(counts, dtype=np.float32)
    counts = counts.copy()
    counts_dense_probe = to_dense(counts[: min(256, adata.n_obs)])
    if not np.allclose(counts_dense_probe, np.round(counts_dense_probe), atol=1e-4):
        print("Warning: input does not look like integer raw counts; preserving nonzero support anyway.")

    raw_var_names = (adata.raw.var_names if adata.raw is not None else adata.var_names).to_numpy()
    raw_var = adata.raw.var.copy() if adata.raw is not None else adata.var.copy()

    if sp.issparse(counts):
        counts_csr = counts.tocsr().astype(np.float32)
        lib = np.asarray(counts_csr.sum(axis=1)).ravel()
        scale = np.divide(args.target_sum, lib, out=np.zeros_like(lib, dtype=np.float32), where=lib > 0)
        values = counts_csr.multiply(scale[:, None]).tocsr()
        values.data = np.log1p(values.data)
        mean = np.asarray(values.mean(axis=0)).ravel()
        sq_mean = np.asarray(values.power(2).mean(axis=0)).ravel()
    else:
        lib = counts.sum(axis=1)
        scale = np.divide(args.target_sum, lib, out=np.zeros_like(lib, dtype=np.float32), where=lib > 0)
        values = np.log1p(counts * scale[:, None]).astype(np.float32)
        mean = values.mean(axis=0)
        sq_mean = np.square(values).mean(axis=0)
    var = np.nan_to_num(sq_mean - np.square(mean), nan=0.0, posinf=0.0, neginf=0.0)
    n_select = min(args.n_top_genes, var.shape[0])
    selected_idx = np.argsort(var)[::-1][:n_select]
    selected_idx = np.sort(selected_idx)
    selected = raw_var_names[selected_idx]

    values_selected = values[:, selected_idx] if sp.issparse(values) else values[:, selected_idx]
    counts_selected = counts[:, selected_idx] if sp.issparse(counts) else counts[:, selected_idx]
    support = (counts_selected > 0).astype(np.float32)
    if sp.issparse(support):
        support = support.tocsr()

    obs = adata.obs.copy()
    var = raw_var.iloc[selected_idx].copy()
    for frame in (obs, var):
        if "_index" in frame.columns:
            frame.drop(columns=["_index"], inplace=True)

    out = sc.AnnData(
        X=values_selected.copy(),
        obs=obs,
        var=var,
    )
    out.obs_names = adata.obs_names.copy()
    out.var_names = selected.copy()
    out.layers["counts"] = counts_selected.copy()
    out.layers["support"] = support.copy()
    out.uns["unified_protocol"] = {
        "source_path": os.path.abspath(args.data_path),
        "n_top_genes": int(args.n_top_genes),
        "target_sum": float(args.target_sum),
        "label_key": label_key,
        "x_is": "normalize_total_log1p_hvg",
        "layers.counts": "raw_counts_hvg",
        "layers.support": "observed_support_from_counts_hvg",
    }
    out.write_h5ad(args.out_path)
    print(f"Wrote {args.out_path}: {out.n_obs} cells x {out.n_vars} genes, support density={float(to_dense(support).mean()):.4f}")


if __name__ == "__main__":
    main()
