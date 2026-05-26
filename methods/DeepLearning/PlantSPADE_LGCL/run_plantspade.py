#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning.PlantSPADE_LGCL.data import load_lgcl_dataset, write_dataset_artifacts
from methods.DeepLearning.PlantSPADE_LGCL.eval import write_evaluation_outputs
from methods.DeepLearning.PlantSPADE_LGCL.eval.marker_analysis import write_marker_outputs
from methods.DeepLearning.PlantSPADE_LGCL.support_gene_attention import (
    SparseAttentionWeights,
    SupportGeneAttention,
    TrainableSupportGeneAttentionRefiner,
)
from methods.DeepLearning.PlantSPADE_LGCL.train import (
    LGCLTrainConfig,
    PlantSPADELGCL,
    normalized_bipartite_support,
    scipy_to_torch_sparse,
    train_lgcl,
)
from methods.DeepLearning.PlantSPADE_LGCL.utils import ensure_dir, save_json, set_seed


def str2bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "t", "yes", "y"}:
        return True
    if value in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got {value!r}")


def parse_float_list(value: str):
    if value is None or str(value).strip() == "":
        return []
    return [float(item) for item in str(value).split(",") if item.strip()]


def parse_int_list(value: str):
    if value is None or str(value).strip() == "":
        return []
    return [int(item) for item in str(value).split(",") if item.strip()]


def parse_args():
    parser = argparse.ArgumentParser(description="PlantSPADE-LGCL fixed-protocol runner")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--label_key", default="auto")
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
    parser.add_argument("--negative_sampler", default="random_zero", choices=["random_zero", "idf_weighted_zero", "neighbor_conflict_zero"])
    parser.add_argument("--negative_neighbor_k", type=int, default=15)
    parser.add_argument("--eval_neighbors", type=int, default=15)
    parser.add_argument("--leiden_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--louvain_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--leiden_resolutions", default="0.2,0.4,0.6,0.8,1.0,1.2")
    parser.add_argument("--sweep_max_cells", type=int, default=10000)
    parser.add_argument("--include_louvain", type=str2bool, default=True)
    parser.add_argument("--use_support_attention", type=str2bool, default=False)
    parser.add_argument("--attention_topk_genes", type=int, default=128)
    parser.add_argument("--attention_beta", type=float, default=0.1)
    parser.add_argument("--attention_gamma", type=float, default=0.1)
    parser.add_argument("--attention_eta", type=float, default=0.5)
    parser.add_argument("--attention_dropout", type=float, default=0.1)
    parser.add_argument("--run_attention_ablations", type=str2bool, default=False)
    parser.add_argument("--attention_topk_sweep", default="64,128,256")
    parser.add_argument("--use_trainable_attention_refiner", type=str2bool, default=False)
    parser.add_argument("--attention_refiner_epochs", type=int, default=20)
    parser.add_argument("--attention_refiner_lr", type=float, default=1e-2)
    parser.add_argument("--attention_refiner_batch_size", type=int, default=2048)
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
        for reserved in ("_index", "reserved_index"):
            if reserved in frame.columns:
                replacement = reserved + "_renamed"
                suffix = 1
                while replacement in frame.columns:
                    replacement = f"{reserved}_renamed_{suffix}"
                    suffix += 1
                frame.rename(columns={reserved: replacement}, inplace=True)
        if "_index" in frame.index:
            frame.index.name = "cell_name" if frame is adata.obs else "gene_name"


def compute_gene_idf(support) -> np.ndarray:
    df = np.diff(support.tocsc().indptr).astype(np.float32)
    n_cells = float(support.shape[0])
    return np.log1p(n_cells / (1.0 + df)).astype(np.float32)


@torch.no_grad()
def make_support_attention_embedding(
    base_embedding: np.ndarray,
    gene_embedding: np.ndarray,
    support,
    amplitude,
    gene_idf: np.ndarray,
    device: torch.device,
    top_k: int,
    beta: float,
    gamma: float,
    eta: float,
    dropout: float,
    normalize_output: bool,
    return_attention: bool = False,
):
    attention = SupportGeneAttention(
        support=support,
        amplitude=amplitude,
        gene_idf=torch.as_tensor(gene_idf, dtype=torch.float32, device=device),
        top_k_genes=top_k,
        beta=beta,
        gamma=gamma,
        eta=eta,
        dropout=dropout,
    ).to(device)
    attention.eval()
    cell_t = torch.as_tensor(base_embedding, dtype=torch.float32, device=device)
    gene_t = torch.as_tensor(gene_embedding, dtype=torch.float32, device=device)
    if return_attention:
        refined_t, weights = attention(cell_t, gene_t, return_attention=True)
    else:
        refined_t = attention(cell_t, gene_t, return_attention=False)
        weights = None
    if normalize_output:
        refined_t = F.normalize(refined_t, dim=1)
    return refined_t.detach().cpu().numpy().astype(np.float32), weights


def train_attention_refiner(
    base_embedding: np.ndarray,
    gene_embedding: np.ndarray,
    projected_global: np.ndarray,
    support,
    amplitude,
    gene_idf: np.ndarray,
    device: torch.device,
    top_k: int,
    beta: float,
    gamma: float,
    eta: float,
    dropout: float,
    epochs: int,
    lr: float,
    batch_size: int,
    temperature: float,
    seed: int,
    normalize_output: bool,
) -> tuple[np.ndarray, Optional[SparseAttentionWeights], dict]:
    refiner = TrainableSupportGeneAttentionRefiner(
        support=support,
        amplitude=amplitude,
        gene_idf=torch.as_tensor(gene_idf, dtype=torch.float32, device=device),
        top_k_genes=top_k,
        beta=beta,
        gamma=gamma,
        eta=eta,
        dropout=dropout,
    ).to(device)
    cell_t = torch.as_tensor(base_embedding, dtype=torch.float32, device=device)
    gene_t = torch.as_tensor(gene_embedding, dtype=torch.float32, device=device)
    global_t = F.normalize(torch.as_tensor(projected_global, dtype=torch.float32, device=device), dim=1)
    optimizer = torch.optim.AdamW(refiner.parameters(), lr=lr, weight_decay=0.0)
    rng = np.random.default_rng(seed + 777)
    n_cells = cell_t.shape[0]
    history = {"loss": [], "contrastive": [], "consistency": [], "beta": [], "gamma": [], "eta": []}

    for epoch in range(1, max(1, epochs) + 1):
        refiner.train()
        csz = min(batch_size, n_cells)
        batch = torch.as_tensor(rng.choice(n_cells, size=csz, replace=False), dtype=torch.long, device=device)
        refined = refiner(cell_t, gene_t, batch_cells=batch, return_attention=False)
        refined_norm = F.normalize(refined, dim=1)
        target = global_t[batch]
        logits = refined_norm @ target.T / temperature
        labels = torch.arange(logits.shape[0], device=device)
        contrastive = 0.5 * (F.cross_entropy(logits, labels) + F.cross_entropy(logits.T, labels))
        consistency = (1.0 - F.cosine_similarity(refined_norm, F.normalize(cell_t[batch], dim=1), dim=1)).mean()
        loss = contrastive + 0.05 * consistency
        optimizer.zero_grad(set_to_none=True)
        loss.backward()
        optimizer.step()
        attn = refiner.attention
        history["loss"].append(float(loss.detach().cpu()))
        history["contrastive"].append(float(contrastive.detach().cpu()))
        history["consistency"].append(float(consistency.detach().cpu()))
        history["beta"].append(float(attn.beta_param.detach().cpu()))
        history["gamma"].append(float(attn.gamma_param.detach().cpu()))
        history["eta"].append(float(attn.eta_param.detach().cpu()))

    refiner.eval()
    with torch.no_grad():
        refined_t, weights = refiner(cell_t, gene_t, return_attention=True)
        if normalize_output:
            refined_t = F.normalize(refined_t, dim=1)
    return refined_t.detach().cpu().numpy().astype(np.float32), weights, history


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
        label_key=args.label_key,
    )
    write_dataset_artifacts(bundle, save_dir)
    if bundle.labels is None:
        raise ValueError("No labels found in h5ad obs; fixed benchmark evaluation requires labels.")
    n_clusters = args.n_clusters if args.n_clusters > 0 else int(len(np.unique(bundle.labels)))
    dataset_name = args.dataset_name or Path(args.data_path).stem
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
    print(f"Step 3: train with BPR + SVD InfoNCE, negative_sampler={args.negative_sampler}")
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
        negative_sampler=args.negative_sampler,
        negative_neighbor_k=args.negative_neighbor_k,
    )
    history = train_lgcl(model, bundle.support, train_config, device)
    save_json(history, os.path.join(save_dir, "training_history.json"))

    print("=" * 72)
    print("Step 4: extract embeddings and evaluate fixed/oracle/sweep protocols")
    print("=" * 72)
    normalize_embedding = not args.no_normalize_embedding
    local_embedding, gene_embedding, projected_global = model.get_embeddings(normalize_output=normalize_embedding)
    if args.global_blend != 0.0:
        blended = torch.as_tensor(local_embedding) + float(args.global_blend) * torch.as_tensor(projected_global)
        if normalize_embedding:
            blended = F.normalize(blended, dim=1)
        base_embedding = blended.numpy().astype(np.float32)
    else:
        base_embedding = local_embedding.astype(np.float32)

    resolutions = parse_float_list(args.leiden_resolutions)
    topk_sweep = parse_int_list(args.attention_topk_sweep)
    gene_idf = compute_gene_idf(bundle.support)
    variant_embeddings = {"baseline": base_embedding}
    attention_weights_for_markers = None

    should_compute_attention = args.use_support_attention or args.run_attention_ablations or args.use_trainable_attention_refiner
    if should_compute_attention:
        print("Computing sparse SupportGeneAttention variants")
        attn_embedding, attention_weights_for_markers = make_support_attention_embedding(
            base_embedding,
            gene_embedding,
            bundle.support,
            bundle.amplitude,
            gene_idf,
            device,
            top_k=args.attention_topk_genes,
            beta=args.attention_beta,
            gamma=args.attention_gamma,
            eta=args.attention_eta,
            dropout=args.attention_dropout,
            normalize_output=normalize_embedding,
            return_attention=True,
        )
        variant_embeddings["support_attention"] = attn_embedding

    if args.use_trainable_attention_refiner:
        trainable_embedding, trainable_weights, refiner_history = train_attention_refiner(
            base_embedding,
            gene_embedding,
            projected_global,
            bundle.support,
            bundle.amplitude,
            gene_idf,
            device,
            top_k=args.attention_topk_genes,
            beta=args.attention_beta,
            gamma=args.attention_gamma,
            eta=args.attention_eta,
            dropout=args.attention_dropout,
            epochs=args.attention_refiner_epochs,
            lr=args.attention_refiner_lr,
            batch_size=args.attention_refiner_batch_size,
            temperature=args.temperature,
            seed=args.seed,
            normalize_output=normalize_embedding,
        )
        variant_embeddings["support_attention_trainable"] = trainable_embedding
        attention_weights_for_markers = trainable_weights
        save_json(refiner_history, os.path.join(save_dir, "attention_refiner_history.json"))

    if args.run_attention_ablations:
        ablation_configs = [
            ("attention_no_amplitude", args.attention_topk_genes, 0.0, args.attention_gamma),
            ("attention_no_idf", args.attention_topk_genes, args.attention_beta, 0.0),
        ]
        for top_k in topk_sweep:
            ablation_configs.append((f"attention_topk_{top_k}", top_k, args.attention_beta, args.attention_gamma))
        for name, top_k, beta, gamma in ablation_configs:
            emb, _ = make_support_attention_embedding(
                base_embedding,
                gene_embedding,
                bundle.support,
                bundle.amplitude,
                gene_idf,
                device,
                top_k=top_k,
                beta=beta,
                gamma=gamma,
                eta=args.attention_eta,
                dropout=args.attention_dropout,
                normalize_output=normalize_embedding,
                return_attention=False,
            )
            variant_embeddings[name] = emb

    if args.use_trainable_attention_refiner:
        primary_variant = "support_attention_trainable"
    else:
        primary_variant = "support_attention" if args.use_support_attention else "baseline"
    embedding = variant_embeddings[primary_variant]

    bundle.adata.obsm["X_plantspade_lgcl_base"] = base_embedding
    if "support_attention" in variant_embeddings:
        bundle.adata.obsm["X_plantspade_lgcl_attn"] = variant_embeddings["support_attention"]
    bundle.adata.obsm["X_plantspade_lgcl"] = embedding
    bundle.adata.uns["plantspade_lgcl"] = {
        "method": args.method_name,
        "input_mode": bundle.input_mode,
        "counts_source": bundle.counts_source,
        "support_density": bundle.support_density,
        "n_edges": int(bundle.support.nnz),
        "n_top_genes": int(bundle.support.shape[1]),
        "negative_sampler": args.negative_sampler,
        "primary_variant": primary_variant,
        "use_support_attention": bool(args.use_support_attention),
        "use_trainable_attention_refiner": bool(args.use_trainable_attention_refiner),
        "attention_topk_genes": int(args.attention_topk_genes),
        "attention_beta": float(args.attention_beta),
        "attention_gamma": float(args.attention_gamma),
        "attention_eta": float(args.attention_eta),
    }

    embedding_path = os.path.join(save_dir, "embedding_final.npy")
    np.save(embedding_path, embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "embeddings_base.npy"), base_embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "global_embedding_svd_projected.npy"), projected_global.astype(np.float32))
    np.save(os.path.join(save_dir, "gene_embedding.npy"), gene_embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "labels.npy"), bundle.labels.astype(np.int64))

    variant_eval_results = {}
    for variant_name, variant_embedding in variant_embeddings.items():
        variant_path = os.path.join(save_dir, f"embedding_{variant_name}.npy")
        np.save(variant_path, variant_embedding.astype(np.float32))
        variant_method = f"{args.method_name}_{variant_name}"
        variant_eval_results[variant_name] = write_evaluation_outputs(
            output_dir=save_dir,
            dataset=dataset_name,
            method=variant_method,
            seed=args.seed,
            embedding=variant_embedding,
            labels=bundle.labels,
            n_clusters=n_clusters,
            n_neighbors=args.eval_neighbors,
            leiden_fixed_resolution=args.leiden_fixed_resolution,
            louvain_fixed_resolution=args.louvain_fixed_resolution,
            leiden_sweep_resolutions=resolutions,
            sweep_max_cells=args.sweep_max_cells,
            include_louvain=args.include_louvain,
            prefix=f"eval_{variant_name}",
            extra={"variant": variant_name, "negative_sampler": args.negative_sampler},
        )

    primary_result = variant_eval_results[primary_variant]
    primary_payload = {
        "dataset": dataset_name,
        "method": args.method_name,
        "seed": int(args.seed),
        "primary_variant": primary_variant,
        "embedding_path": os.path.abspath(embedding_path),
        "fixed": primary_result["fixed"],
        "oracle": primary_result["oracle"],
        "note": "Main results use only fixed protocol. leiden_oracle_best is supplementary upper bound.",
    }
    save_json(primary_payload, os.path.join(save_dir, f"{args.method_name}.json"))
    save_json(primary_result["fixed"], os.path.join(save_dir, "metrics.json"))
    pd.DataFrame(
        [
            {"dataset": dataset_name, "method": args.method_name, "seed": args.seed, "cluster_method": key, **vals}
            for key, vals in primary_result["fixed"].items()
        ]
    ).to_csv(os.path.join(save_dir, f"{args.method_name}_fixed.csv"), index=False)

    pred_for_h5 = primary_result["preds"].get("kmeans_known_k")
    save_embedding_h5(os.path.join(save_dir, "embedding.h5"), embedding, labels=bundle.labels, pred_labels=pred_for_h5)

    if attention_weights_for_markers is not None:
        cluster_for_markers = primary_result["preds"].get("leiden_fixed")
        if cluster_for_markers is None:
            cluster_for_markers = primary_result["preds"].get("kmeans_known_k")
        write_marker_outputs(
            save_dir,
            attention_weights_for_markers,
            cluster_for_markers,
            bundle.support,
            bundle.amplitude,
            bundle.gene_names,
            top_n_attention=20,
            top_n_deg=50,
        )

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
        "counts_source": bundle.counts_source,
        "negative_sampler": args.negative_sampler,
        "primary_variant": primary_variant,
        "fixed_metrics": primary_result["fixed"],
        "oracle_metrics": primary_result["oracle"],
        "variant_names": sorted(variant_embeddings.keys()),
    }
    save_json(summary, os.path.join(save_dir, "summary.json"))
    print("Fixed-protocol metrics:", primary_result["fixed"])
    print(f"Results saved to: {save_dir}")


if __name__ == "__main__":
    main()
