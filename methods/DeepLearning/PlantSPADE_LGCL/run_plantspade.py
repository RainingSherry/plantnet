#!/usr/bin/env python
import argparse
import os
import sys
import warnings
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import torch.nn.functional as F
from sklearn.cluster import KMeans

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
BENCHMARK_DIR = ROOT / "benchmarks" / "unified_protocol"
for path in [str(ROOT), str(BENCHMARK_DIR)]:
    if path not in sys.path:
        sys.path.insert(0, path)

from common import compute_metrics, save_json as save_benchmark_json  # noqa: E402
from methods.DeepLearning.PlantSPADE_LGCL.data import load_lgcl_dataset  # noqa: E402
from methods.DeepLearning.PlantSPADE_LGCL.support_gene_attention import (  # noqa: E402
    SparseAttentionWeights,
    SupportGeneAttention,
)
from methods.DeepLearning.PlantSPADE_LGCL.train import (  # noqa: E402
    LGCLTrainConfig,
    PlantSPADELGCL,
    normalized_bipartite_support,
    scipy_to_torch_sparse,
    train_lgcl,
)
from methods.DeepLearning.PlantSPADE_LGCL.utils import ensure_dir, save_json, set_seed  # noqa: E402


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
    return [float(item) for item in str(value).split(",") if item.strip()]


def parse_int_list(value: str):
    return [int(item) for item in str(value).split(",") if item.strip()]


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
    parser.add_argument("--leiden_resolutions", default="0.2,0.4,0.6,0.8,1.0,1.2")
    parser.add_argument("--use_support_attention", type=str2bool, default=False)
    parser.add_argument("--attention_topk_genes", type=int, default=128)
    parser.add_argument("--attention_beta", type=float, default=0.1)
    parser.add_argument("--attention_gamma", type=float, default=0.1)
    parser.add_argument("--attention_eta", type=float, default=0.5)
    parser.add_argument("--attention_dropout", type=float, default=0.1)
    parser.add_argument("--run_attention_ablations", type=str2bool, default=False)
    parser.add_argument("--attention_topk_sweep", default="64,128,256")
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


def format_resolution(resolution: float) -> str:
    return str(float(resolution)).replace(".", "p")


def evaluate_embedding_with_leiden_sweep(
    embedding: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    seed: int,
    n_neighbors: int,
    resolutions,
) -> dict:
    embedding = np.asarray(embedding, dtype=np.float32)
    metrics = {}
    preds = {}
    mapped_preds = {}

    pred_kmeans = KMeans(n_clusters=n_clusters, n_init=1, random_state=seed).fit_predict(embedding)
    metrics["kmeans"], mapped = compute_metrics(labels, pred_kmeans, embedding=embedding)
    preds["kmeans"] = pred_kmeans.astype(np.int64)
    mapped_preds["kmeans"] = mapped.astype(np.int64)

    adata = sc.AnnData(embedding)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X")
    best_key = None
    best_nmi = -np.inf
    for resolution in resolutions:
        key = f"leiden_res_{format_resolution(resolution)}"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", FutureWarning)
            try:
                sc.tl.leiden(
                    adata,
                    random_state=seed,
                    resolution=float(resolution),
                    key_added=key,
                    flavor="igraph",
                    n_iterations=2,
                    directed=False,
                )
            except TypeError:
                sc.tl.leiden(adata, random_state=seed, resolution=float(resolution), key_added=key)
        pred = adata.obs[key].astype(int).to_numpy()
        vals, mapped = compute_metrics(labels, pred, embedding=embedding)
        vals["resolution"] = float(resolution)
        metrics[key] = vals
        preds[key] = pred.astype(np.int64)
        mapped_preds[key] = mapped.astype(np.int64)
        if vals["nmi"] > best_nmi:
            best_nmi = vals["nmi"]
            best_key = key

    if best_key is not None:
        best_vals = dict(metrics[best_key])
        best_vals["selected_from"] = best_key
        metrics["leiden_best"] = best_vals
        metrics["leiden"] = dict(best_vals)
        preds["leiden_best"] = preds[best_key]
        preds["leiden"] = preds[best_key]
        mapped_preds["leiden_best"] = mapped_preds[best_key]
        mapped_preds["leiden"] = mapped_preds[best_key]
    return {"metrics": metrics, "preds": preds, "mapped_preds": mapped_preds, "best_leiden_key": best_key}


def write_eval_outputs(
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
    resolutions,
    file_prefix: Optional[str] = None,
) -> dict:
    result = evaluate_embedding_with_leiden_sweep(
        embedding,
        labels,
        n_clusters=n_clusters,
        seed=seed,
        n_neighbors=n_neighbors,
        resolutions=resolutions,
    )
    metrics = result["metrics"]
    preds = result["preds"]
    mapped_preds = result["mapped_preds"]
    prefix = file_prefix or method_name
    payload = {
        "dataset": dataset_name,
        "method": method_name,
        "data_path": os.path.abspath(data_path),
        "embedding_path": os.path.abspath(embedding_path),
        "label_key": label_key,
        "n_cells": int(embedding.shape[0]),
        "embedding_dim": int(embedding.shape[1]),
        "best_leiden_key": result["best_leiden_key"],
        "metrics": metrics,
    }
    save_benchmark_json(payload, os.path.join(save_dir, f"{prefix}.json"))
    rows = []
    for cluster_method, vals in metrics.items():
        row = {"dataset": dataset_name, "method": method_name, "cluster_method": cluster_method}
        row.update(vals)
        rows.append(row)
    pd.DataFrame(rows).to_csv(os.path.join(save_dir, f"{prefix}.csv"), index=False)
    for name, pred in preds.items():
        np.save(os.path.join(save_dir, f"{prefix}_{name}.npy"), pred.astype(np.int64))
        np.save(os.path.join(save_dir, f"{prefix}_{name}_mapped.npy"), mapped_preds[name].astype(np.int64))
    result["payload"] = payload
    return result


def compute_gene_idf(support) -> np.ndarray:
    df = np.diff(support.tocsc().indptr).astype(np.float32)
    n_cells = float(support.shape[0])
    return np.log(n_cells / (1.0 + df)).astype(np.float32)


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


def attention_top_genes_per_cluster(
    attention_weights: SparseAttentionWeights,
    cluster_labels: np.ndarray,
    support,
    amplitude,
    gene_names: np.ndarray,
    top_n: int = 20,
) -> pd.DataFrame:
    rows = []
    n_genes = support.shape[1]
    labels = np.asarray(cluster_labels, dtype=np.int64)
    for cluster_id in sorted(np.unique(labels)):
        cluster_cells = np.flatnonzero(labels == cluster_id)
        if cluster_cells.size == 0:
            continue
        edge_mask = labels[attention_weights.cells] == cluster_id
        if np.any(edge_mask):
            mean_attention = np.bincount(
                attention_weights.genes[edge_mask],
                weights=attention_weights.weights[edge_mask],
                minlength=n_genes,
            ).astype(np.float64)
            mean_attention /= float(cluster_cells.size)
        else:
            mean_attention = np.zeros(n_genes, dtype=np.float64)
        top_genes = np.argsort(-mean_attention, kind="stable")[:top_n]
        support_block = support[cluster_cells][:, top_genes]
        amplitude_block = amplitude[cluster_cells][:, top_genes]
        support_rate = np.asarray(support_block.mean(axis=0)).ravel()
        mean_expression = np.asarray(amplitude_block.mean(axis=0)).ravel()
        for local_idx, gene_idx in enumerate(top_genes):
            rows.append(
                {
                    "cluster_id": int(cluster_id),
                    "gene_name": str(gene_names[gene_idx]),
                    "mean_attention": float(mean_attention[gene_idx]),
                    "support_rate": float(support_rate[local_idx]),
                    "mean_expression": float(mean_expression[local_idx]),
                }
            )
    return pd.DataFrame(rows)


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
    print("Step 4: extract base/attention embeddings and run unified KMeans/Leiden evaluation")
    print("=" * 72)
    normalize_embedding = not args.no_normalize_embedding
    local_embedding, gene_embedding, projected_global = model.get_embeddings(normalize_output=normalize_embedding)
    if args.global_blend != 0.0:
        local_t = torch.as_tensor(local_embedding)
        global_t = torch.as_tensor(projected_global)
        blended = local_t + float(args.global_blend) * global_t
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

    should_compute_attention = args.use_support_attention or args.run_attention_ablations
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

    primary_variant = "support_attention" if args.use_support_attention else "baseline"
    embedding = variant_embeddings[primary_variant]

    bundle.adata.obsm["X_plantspade_lgcl_base"] = base_embedding
    if "support_attention" in variant_embeddings:
        bundle.adata.obsm["X_plantspade_lgcl_attn"] = variant_embeddings["support_attention"]
    bundle.adata.obsm["X_plantspade_lgcl"] = embedding
    bundle.adata.uns["plantspade_lgcl"] = {
        "method": args.method_name,
        "input_mode": bundle.input_mode,
        "support_density": bundle.support_density,
        "n_edges": int(bundle.support.nnz),
        "n_top_genes": int(bundle.support.shape[1]),
        "primary_variant": primary_variant,
        "use_support_attention": bool(args.use_support_attention),
        "attention_topk_genes": int(args.attention_topk_genes),
        "attention_beta": float(args.attention_beta),
        "attention_gamma": float(args.attention_gamma),
        "attention_eta": float(args.attention_eta),
    }

    embedding_path = os.path.join(save_dir, "embedding_final.npy")
    np.save(embedding_path, embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "embeddings_base.npy"), base_embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "embeddings_direct.npy"), base_embedding.astype(np.float32))
    if "support_attention" in variant_embeddings:
        np.save(os.path.join(save_dir, "embeddings_attn.npy"), variant_embeddings["support_attention"].astype(np.float32))
    np.save(os.path.join(save_dir, "global_embedding_svd_projected.npy"), projected_global.astype(np.float32))
    np.save(os.path.join(save_dir, "gene_embedding.npy"), gene_embedding.astype(np.float32))
    np.save(os.path.join(save_dir, "labels.npy"), bundle.labels.astype(np.int64))

    variant_eval_results = {}
    for variant_name, variant_embedding in variant_embeddings.items():
        variant_path = os.path.join(save_dir, f"embedding_{variant_name}.npy")
        np.save(variant_path, variant_embedding.astype(np.float32))
        variant_eval_results[variant_name] = write_eval_outputs(
            save_dir=save_dir,
            dataset_name=dataset_name,
            method_name=f"{args.method_name}_{variant_name}",
            data_path=args.data_path,
            embedding_path=variant_path,
            embedding=variant_embedding,
            labels=bundle.labels,
            label_key=bundle.label_key,
            n_clusters=n_clusters,
            seed=args.seed,
            n_neighbors=args.eval_neighbors,
            resolutions=resolutions,
            file_prefix=f"eval_{variant_name}",
        )

    eval_result = variant_eval_results[primary_variant]
    primary_payload = dict(eval_result["payload"])
    primary_payload["method"] = args.method_name
    primary_payload["primary_variant"] = primary_variant
    primary_payload["embedding_path"] = os.path.abspath(embedding_path)
    save_benchmark_json(primary_payload, os.path.join(save_dir, f"{args.method_name}.json"))
    save_benchmark_json(eval_result["metrics"], os.path.join(save_dir, "metrics.json"))
    primary_rows = []
    for cluster_method, vals in eval_result["metrics"].items():
        row = {"dataset": dataset_name, "method": args.method_name, "cluster_method": cluster_method}
        row.update(vals)
        primary_rows.append(row)
    pd.DataFrame(primary_rows).to_csv(os.path.join(save_dir, f"{args.method_name}.csv"), index=False)
    for name, pred in eval_result["preds"].items():
        np.save(os.path.join(save_dir, f"{args.method_name}_{name}.npy"), pred.astype(np.int64))
        np.save(os.path.join(save_dir, f"pred_labels_{name}.npy"), pred.astype(np.int64))
        np.save(os.path.join(save_dir, f"{args.method_name}_{name}_mapped.npy"), eval_result["mapped_preds"][name].astype(np.int64))
        np.save(os.path.join(save_dir, f"pred_labels_{name}_mapped.npy"), eval_result["mapped_preds"][name].astype(np.int64))

    ablation_rows = []
    for variant_name, result in variant_eval_results.items():
        for cluster_method, vals in result["metrics"].items():
            row = {"dataset": dataset_name, "variant": variant_name, "cluster_method": cluster_method}
            row.update(vals)
            ablation_rows.append(row)
    pd.DataFrame(ablation_rows).to_csv(os.path.join(save_dir, "attention_ablation_metrics.csv"), index=False)
    save_json(
        {name: result["metrics"] for name, result in variant_eval_results.items()},
        os.path.join(save_dir, "attention_ablation_metrics.json"),
    )

    pred_for_h5 = eval_result["preds"].get("kmeans")
    save_embedding_h5(os.path.join(save_dir, "embedding.h5"), embedding, labels=bundle.labels, pred_labels=pred_for_h5)

    if attention_weights_for_markers is not None and "leiden_best" in eval_result["preds"]:
        top_gene_df = attention_top_genes_per_cluster(
            attention_weights_for_markers,
            eval_result["preds"]["leiden_best"],
            bundle.support,
            bundle.amplitude,
            bundle.gene_names,
            top_n=20,
        )
        top_gene_df.to_csv(os.path.join(save_dir, "attention_top_genes_per_cluster.csv"), index=False)

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
        "primary_variant": primary_variant,
        "metrics": eval_result["metrics"],
        "ablation_metrics": {name: result["metrics"] for name, result in variant_eval_results.items()},
    }
    save_json(summary, os.path.join(save_dir, "summary.json"))
    print("Final metrics:", eval_result["metrics"])
    print(f"Results saved to: {save_dir}")


if __name__ == "__main__":
    main()
