#!/usr/bin/env python
import argparse
import os
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
BENCHMARK_DIR = ROOT / "benchmarks" / "unified_protocol"
for path in [str(ROOT), str(BENCHMARK_DIR)]:
    if path not in sys.path:
        sys.path.insert(0, path)

from common import evaluate_embedding, save_json as save_benchmark_json  # noqa: E402
from methods.DeepLearning.PlantSPADE_LGCL.data import load_lgcl_dataset  # noqa: E402
from methods.DeepLearning.PlantSPADE_LGCL.train import (  # noqa: E402
    LGCLTrainConfig,
    PlantSPADELGCL,
    normalized_bipartite_support,
    scipy_to_torch_sparse,
    train_lgcl,
)
from methods.DeepLearning.PlantSPADE_LGCL.utils import ensure_dir, save_json, set_seed  # noqa: E402


def parse_args():
    parser = argparse.ArgumentParser(description="PlantSPADE-LGCL: LightGCN and SVD contrastive cell-gene clustering")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--latent_dim", type=int, default=32)
    parser.add_argument("--svd_dim", type=int, default=0)
    parser.add_argument("--svd_iter", type=int, default=7)
    parser.add_argument("--layers", type=int, default=2)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--pairs_per_epoch", type=int, default=262144)
    parser.add_argument("--contrastive_batch_size", type=int, default=2048)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--edge_dropout", type=float, default=0.1)
    parser.add_argument("--temperature", type=float, default=0.2)
    parser.add_argument("--contrastive_weight", type=float, default=0.05)
    parser.add_argument("--module_weight", type=float, default=0.001)
    parser.add_argument("--num_modules", type=int, default=16)
    parser.add_argument("--module_top_k", type=int, default=30)
    parser.add_argument("--target_sum", type=float, default=1e4)
    parser.add_argument("--global_blend", type=float, default=0.0)
    parser.add_argument("--eval_neighbors", type=int, default=15)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--no_normalize_embedding", action="store_true")
    parser.add_argument("--no_save_h5ad", action="store_true")
    parser.add_argument("--method_name", default="plantspade_lgcl")
    return parser.parse_args()


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    if gpu in {0, 7}:
        raise ValueError("GPU 0 and GPU 7 are forbidden for this run. Use --gpu 1,2,3,4,5,6 or --no_cuda.")
    return torch.device(f"cuda:{gpu}")


def save_embedding_h5(path: str, embedding: np.ndarray, labels: np.ndarray = None, pred_labels: np.ndarray = None) -> None:
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=embedding.astype(np.float32))
        if labels is not None:
            handle.create_dataset("labels", data=labels.astype(np.int64))
        if pred_labels is not None:
            handle.create_dataset("Y", data=pred_labels.astype(np.int64))


def sanitize_anndata_for_write(adata) -> None:
    for frame in (adata.obs, adata.var):
        if "_index" in frame.columns:
            replacement = "reserved_index"
            suffix = 1
            while replacement in frame.columns:
                replacement = f"reserved_index_{suffix}"
                suffix += 1
            frame.rename(columns={"_index": replacement}, inplace=True)


def write_unified_eval(
    save_dir: str,
    dataset_name: str,
    method_name: str,
    data_path: str,
    embedding_path: str,
    embedding: np.ndarray,
    labels: np.ndarray,
    label_key: str,
    n_clusters: int,
    seed: int,
    n_neighbors: int,
) -> dict:
    metrics, preds = evaluate_embedding(
        embedding,
        labels,
        n_clusters=n_clusters,
        seed=seed,
        n_neighbors=n_neighbors,
    )
    payload = {
        "dataset": dataset_name,
        "method": method_name,
        "data_path": os.path.abspath(data_path),
        "embedding_path": os.path.abspath(embedding_path),
        "label_key": label_key,
        "n_cells": int(embedding.shape[0]),
        "embedding_dim": int(embedding.shape[1]),
        "metrics": metrics,
    }
    save_benchmark_json(payload, os.path.join(save_dir, f"{method_name}.json"))
    save_benchmark_json(metrics, os.path.join(save_dir, "metrics.json"))
    rows = []
    for cluster_method, vals in metrics.items():
        row = {"dataset": dataset_name, "method": method_name, "cluster_method": cluster_method}
        row.update(vals)
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(save_dir, f"{method_name}.csv"), index=False)
    for name, pred in preds.items():
        np.save(os.path.join(save_dir, f"{method_name}_{name}.npy"), pred.astype(np.int64))
        np.save(os.path.join(save_dir, f"pred_labels_{name}.npy"), pred.astype(np.int64))
    return {"metrics": metrics, "preds": preds, "payload": payload}


def main():
    args = parse_args()
    set_seed(args.seed)
    save_dir = ensure_dir(args.save_dir)
    save_json(vars(args), os.path.join(save_dir, "args.json"))

    device = get_device(args.gpu, args.no_cuda)
    print(f"Using device: {device}")

    svd_dim = args.svd_dim if args.svd_dim > 0 else args.latent_dim
    print("=" * 72)
    print("Step 1: load h5ad, build raw-count support M and log1p amplitude A")
    print("=" * 72)
    bundle = load_lgcl_dataset(
        args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        svd_dim=svd_dim,
        svd_iter=args.svd_iter,
        seed=args.seed,
    )
    if bundle.labels is None:
        raise ValueError("No labels found in h5ad obs; requested benchmark evaluation requires labels.")
    n_clusters = args.n_clusters if args.n_clusters > 0 else int(len(np.unique(bundle.labels)))
    dataset_name = Path(args.data_path).stem
    print(
        f"Cells={bundle.support.shape[0]} genes={bundle.support.shape[1]} "
        f"edges={bundle.support.nnz} density={bundle.support_density:.4f} clusters={n_clusters}"
    )

    print("=" * 72)
    print("Step 2: build normalized bipartite graph and PlantSPADE-LGCL model")
    print("=" * 72)
    adj_norm = normalized_bipartite_support(bundle.support)
    adj_torch = scipy_to_torch_sparse(adj_norm, device)
    model = PlantSPADELGCL(
        n_cells=bundle.support.shape[0],
        n_genes=bundle.support.shape[1],
        latent_dim=args.latent_dim,
        adj_norm=adj_torch,
        global_cell_embedding=bundle.global_embedding,
        num_layers=args.layers,
        edge_dropout=args.edge_dropout,
        temperature=args.temperature,
        num_modules=args.num_modules,
        module_top_k=args.module_top_k,
    ).to(device)
    print(f"Parameters: {sum(p.numel() for p in model.parameters()):,}")

    print("=" * 72)
    print("Step 3: train with BPR ranking loss plus SVD InfoNCE alignment")
    print("=" * 72)
    train_config = LGCLTrainConfig(
        epochs=args.epochs,
        pairs_per_epoch=args.pairs_per_epoch,
        contrastive_batch_size=args.contrastive_batch_size,
        lr=args.lr,
        weight_decay=args.weight_decay,
        contrastive_weight=args.contrastive_weight,
        module_weight=args.module_weight,
        seed=args.seed,
    )
    history = train_lgcl(model, bundle.support, train_config, device)
    save_json(history, os.path.join(save_dir, "training_history.json"))

    print("=" * 72)
    print("Step 4: extract direct embedding and run unified KMeans/Leiden evaluation")
    print("=" * 72)
    normalize_embedding = not args.no_normalize_embedding
    local_embedding, gene_embedding, projected_global = model.get_embeddings(normalize_output=normalize_embedding)
    if args.global_blend != 0.0:
        local_t = torch.as_tensor(local_embedding)
        global_t = torch.as_tensor(projected_global)
        blended = local_t + float(args.global_blend) * global_t
        if normalize_embedding:
            blended = F.normalize(blended, dim=1)
        embedding = blended.numpy().astype(np.float32)
    else:
        embedding = local_embedding.astype(np.float32)

    bundle.adata.obsm["X_plantspade_lgcl"] = embedding
    bundle.adata.uns["plantspade_lgcl"] = {
        "method": args.method_name,
        "input_mode": bundle.input_mode,
        "support_density": bundle.support_density,
        "n_edges": int(bundle.support.nnz),
        "n_top_genes": int(bundle.support.shape[1]),
    }

    embedding_path = os.path.join(save_dir, "embedding_final.npy")
    np.save(embedding_path, embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "embeddings_direct.npy"), embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "global_embedding_svd_projected.npy"), projected_global.astype(np.float32))
    np.save(os.path.join(save_dir, "gene_embedding.npy"), gene_embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "labels.npy"), bundle.labels.astype(np.int64))

    eval_result = write_unified_eval(
        save_dir=save_dir,
        dataset_name=dataset_name,
        method_name=args.method_name,
        data_path=args.data_path,
        embedding_path=embedding_path,
        embedding=embedding,
        labels=bundle.labels,
        label_key=bundle.label_key,
        n_clusters=n_clusters,
        seed=args.seed,
        n_neighbors=args.eval_neighbors,
    )
    pred_for_h5 = eval_result["preds"].get("kmeans")
    save_embedding_h5(os.path.join(save_dir, "embedding.h5"), embedding, labels=bundle.labels, pred_labels=pred_for_h5)

    module_rows = model.module_top_genes(bundle.gene_names, gene_embedding=gene_embedding, top_k=args.module_top_k)
    save_json({"top_genes": module_rows}, os.path.join(save_dir, "module_top_genes.json"))
    pd.DataFrame(module_rows).to_csv(os.path.join(save_dir, "module_top_genes.csv"), index=False)

    torch.save(
        {
            "model_state": model.state_dict(),
            "args": vars(args),
            "gene_names": bundle.gene_names.astype(str),
        },
        os.path.join(save_dir, "model.pt"),
    )
    if not args.no_save_h5ad:
        sanitize_anndata_for_write(bundle.adata)
        bundle.adata.write_h5ad(os.path.join(save_dir, "adata_plantspade_lgcl.h5ad"), compression="gzip")

    summary = {
        "method": args.method_name,
        "dataset": dataset_name,
        "n_cells": int(bundle.support.shape[0]),
        "n_genes": int(bundle.support.shape[1]),
        "n_edges": int(bundle.support.nnz),
        "support_density": float(bundle.support_density),
        "n_clusters": int(n_clusters),
        "input_mode": bundle.input_mode,
        "metrics": eval_result["metrics"],
    }
    save_json(summary, os.path.join(save_dir, "summary.json"))
    print("Final metrics:", eval_result["metrics"])
    print(f"Results saved to: {save_dir}")


if __name__ == "__main__":
    main()
