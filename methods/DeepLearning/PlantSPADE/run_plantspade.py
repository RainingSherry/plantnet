import argparse
import os
import sys

import numpy as np

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(CURRENT_DIR)
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(PARENT)))
METHODS_DIR = os.path.join(ROOT, "methods")
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
if METHODS_DIR not in sys.path:
    sys.path.insert(0, METHODS_DIR)
if PARENT not in sys.path:
    sys.path.insert(0, PARENT)
if CURRENT_DIR not in sys.path:
    sys.path.insert(0, CURRENT_DIR)

from methods.utils import save as save_benchmark
from maskdiffusion.data import load_sc_dataset
from maskdiffusion.eval.cluster_eval import cluster_and_evaluate
from maskdiffusion.eval.marker_eval import marker_gene_enrichment
from maskdiffusion.eval.sparsity_eval import evaluate_support_predictions
from maskdiffusion.train.train_joint import run_two_stage_training
from maskdiffusion.utils import ensure_dir, save_json, set_seed


def parse_args():
    parser = argparse.ArgumentParser(description="MaskDiffusion for scRNA-seq clustering")
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, required=True)
    parser.add_argument("--input_mode", type=str, default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--latent_dim", type=int, default=16)
    parser.add_argument("--hidden_dim", type=int, default=256)
    parser.add_argument("--diffusion_hidden_dim", type=int, default=128)
    parser.add_argument("--diffusion_steps", type=int, default=100)
    parser.add_argument("--mask_epochs", type=int, default=20)
    parser.add_argument("--embedding_epochs", type=int, default=50)
    parser.add_argument("--joint_epochs", type=int, default=10)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--mask_threshold", type=float, default=0.5)
    parser.add_argument("--cluster_method", type=str, default="leiden", choices=["kmeans", "leiden"])
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=0)
    parser.add_argument("--no_cuda", action="store_true")
    return parser.parse_args()


def main():
    args = parse_args()
    set_seed(args.seed)
    save_dir = ensure_dir(args.save_dir)

    bundle = load_sc_dataset(
        file_path=args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
    )

    result = run_two_stage_training(args=args, bundle=bundle)
    embedding = result["embedding"]
    mask_probs = result["mask_probs"]
    denoised_values = result["denoised_values"]

    labels = bundle.labels
    if labels is None:
        raise ValueError("Dataset lacks label column, so clustering evaluation cannot run.")
    n_clusters = len(np.unique(labels))

    cluster_metrics = cluster_and_evaluate(
        embedding=embedding,
        labels=labels,
        n_clusters=n_clusters,
        method=args.cluster_method,
    )

    bundle.adata.obsm["X_scspade"] = embedding.astype(np.float32)
    bundle.adata.layers["scspade_mask"] = mask_probs.astype(np.float32)
    bundle.adata.layers["scspade_denoised"] = denoised_values.astype(np.float32)
    bundle.adata.obs["scspade_cluster"] = cluster_metrics["pred_labels"].astype(str)

    marker_metrics = marker_gene_enrichment(
        adata=bundle.adata,
        pred_key="scspade_cluster",
        label_key="Celltype" if "Celltype" in bundle.adata.obs.columns else None,
        top_n=20,
    )
    sparsity_metrics = evaluate_support_predictions(bundle.support, mask_probs, threshold=args.mask_threshold)

    eval_summary = {
        **cluster_metrics["metrics"],
        **marker_metrics,
        **sparsity_metrics,
        "n_cells": int(bundle.values.shape[0]),
        "n_genes": int(bundle.values.shape[1]),
        "input_mode": bundle.input_mode,
    }
    save_json(eval_summary, os.path.join(save_dir, "summary.json"))
    save_benchmark(save_dir, labels, cluster_metrics["pred_labels"], epoch="final", embedding=embedding)

    adata_to_save = bundle.adata.copy()
    if "_index" in adata_to_save.var.columns:
        adata_to_save.var = adata_to_save.var.drop(columns=["_index"])
    if "_index" in adata_to_save.obs.columns:
        adata_to_save.obs = adata_to_save.obs.drop(columns=["_index"])
    dataset_name = os.path.splitext(os.path.basename(args.data_path))[0]
    adata_to_save.write_h5ad(os.path.join(save_dir, f"{dataset_name}_scspade.h5ad"))
    np.save(os.path.join(save_dir, "mask_probs.npy"), mask_probs)
    np.save(os.path.join(save_dir, "denoised_values.npy"), denoised_values)

    print("Finished MaskDiffusion clustering.")
    print(eval_summary)


if __name__ == "__main__":
    main()
