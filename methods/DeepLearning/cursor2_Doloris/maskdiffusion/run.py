#!/usr/bin/env python
"""ScSpade: Support-Masked Latent Diffusion Autoencoder for scRNA-seq Clustering.

This is the unified cursor2_Doloris maskdiffusion implementation that integrates:
  1. SupportMaskNet: predicts per-gene activation probability from expression.
  2. LatentDiffusionAE: DDPM denoising in latent space.
  3. ClusterLoss: DEC-style soft clustering (KMeans-initialized).

Training strategy (matching design doc):
  Phase 1 (epoch < warmup): mask + recon only (stable structure learning)
  Phase 2 (warmup to warmup+50): mask + recon + gentle diffusion (0.05)
  Phase 3 (epoch >= warmup+50): all losses with configured weights

Evaluation: Both direct embedding (raw encoder) and diffusion embedding are
extracted and evaluated separately. This ablation is critical for understanding
whether the diffusion component contributes to clustering quality.
"""

import os
import sys
import json
import time
import random
import argparse
from pathlib import Path

import numpy as np
import torch

# Ensure the package is importable
sys.path.insert(0, str(Path(__file__).parent))

from data import load_sc_dataset
from train import run_scspade_training
from eval import cluster_and_evaluate, evaluate_support_predictions, marker_gene_enrichment

# PlantNet evaluation integration
try:
    from methods.evaluation import evaluation as plantnet_evaluation
    from methods.utils import save as plantnet_save_benchmark
    HAS_PLANTNET = True
except ImportError:
    plantnet_evaluation = None
    plantnet_save_benchmark = None
    HAS_PLANTNET = False


# ── Utilities ───────────────────────────────────────────────────────────────────


def set_seed(seed: int):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True


def get_device(gpu: int = 0, no_cuda: bool = False) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    return torch.device(f"cuda:{gpu}")


def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def save_json(obj, path: str):
    """Save a JSON file, handling numpy types and circular references."""
    # Strip numpy arrays (pred_labels) which can cause issues
    def strip_arrays(o):
        if isinstance(o, dict):
            return {k: strip_arrays(v) for k, v in o.items()}
        elif isinstance(o, (list, tuple)):
            return [strip_arrays(v) for v in o]
        elif isinstance(o, np.ndarray):
            return None  # strip arrays
        elif isinstance(o, (np.floating,)):
            return float(o)
        elif isinstance(o, (np.integer,)):
            return int(o)
        else:
            return o

    clean = strip_arrays(obj)
    with open(path, "w") as f:
        json.dump(clean, f, indent=2)


# ── Evaluation callback ─────────────────────────────────────────────────────────


def make_eval_fn(n_clusters: int, cluster_methods: tuple = ("kmeans", "leiden")):
    """Create the evaluation function passed to training for per-epoch eval."""
    def eval_fn(embeddings: np.ndarray, labels: np.ndarray, **kwargs) -> dict:
        result = cluster_and_evaluate(
            embeddings=embeddings,
            labels=labels,
            n_clusters=n_clusters,
            methods=cluster_methods,
            metric="nmi",
            return_all=True,
        )
        best = result["best_metrics"]
        return {k: v for k, v in best.items() if k != "pred_labels"}
    return eval_fn


# ── Main evaluation and saving ─────────────────────────────────────────────────


def evaluate_and_save(
    embeddings_direct: np.ndarray,
    embeddings_diffusion: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    save_dir: str,
    cluster_methods: tuple = ("kmeans", "leiden"),
) -> dict:
    """Run full evaluation on both embedding types using PlantNet interface.

    Returns a dict with all metrics for both direct and diffusion embeddings.
    """
    all_metrics = {}

    for emb_name, emb in [("direct", embeddings_direct), ("diffusion", embeddings_diffusion)]:
        print(f"\n[Evaluate] {emb_name} embedding: shape={emb.shape}")
        result = cluster_and_evaluate(
            embeddings=emb,
            labels=labels,
            n_clusters=n_clusters,
            methods=cluster_methods,
            metric="nmi",
            save_dir=save_dir,
            save_key=f"scspade_{emb_name}",
        )

        best = result["best_metrics"]
        best_method = result["best_method"]

        print(f"  Best method: {best_method}")
        print(f"  ACC={best.get('acc', 0):.4f}  NMI={best.get('nmi', 0):.4f}  "
              f"ARI={best.get('ari', 0):.4f}  F1={best.get('f1_macro', 0):.4f}")

        all_metrics[emb_name] = {
            "best_method": best_method,
            **best,
        }

        # Save using PlantNet benchmark interface
        if HAS_PLANTNET and plantnet_save_benchmark:
            best_pred = result[best_method]["pred_labels"]
            bench_dir = os.path.join(save_dir, f"benchmark_{emb_name}")
            ensure_dir(bench_dir)
            plantnet_save_benchmark(
                bench_dir,
                labels,
                best_pred,
                epoch="final",
                embedding=emb,
            )

    return all_metrics


# ── CLI ─────────────────────────────────────────────────────────────────────────


def parse_args():
    parser = argparse.ArgumentParser(
        description="ScSpade: Support-Masked Latent Diffusion for scRNA-seq Clustering"
    )

    # Data
    parser.add_argument("--data_path", type=str, required=True,
                        help="Path to .h5ad file")
    parser.add_argument("--save_dir", type=str, required=True,
                        help="Directory to save results")
    parser.add_argument("--label_key", type=str, default=None,
                        help="obs column name for cell-type labels (auto-detected if None)")
    parser.add_argument("--n_clusters", type=int, default=0,
                        help="Expected number of clusters (auto-detected if 0)")

    # Preprocessing
    parser.add_argument("--n_top_genes", type=int, default=2000,
                        help="Number of highly variable genes")
    parser.add_argument("--input_mode", type=str, default="auto",
                        choices=["auto", "raw", "log1p"])
    parser.add_argument("--min_genes", type=int, default=0,
                        help="Minimum genes per cell before HVG; 0 disables filtering for unified benchmark")
    parser.add_argument("--min_cells", type=int, default=0,
                        help="Minimum cells per gene before HVG; 0 disables filtering for unified benchmark")

    # Model
    parser.add_argument("--latent_dim", type=int, default=32,
                        help="Latent embedding dimension")
    parser.add_argument("--hidden_dim", type=int, default=256,
                        help="Encoder/decoder hidden dimension")
    parser.add_argument("--diffusion_hidden_dim", type=int, default=256,
                        help="Denoiser hidden dimension")
    parser.add_argument("--diffusion_steps", type=int, default=100,
                        help="Number of diffusion steps")
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--mask_dropout", type=float, default=0.1)

    # Training
    parser.add_argument("--epochs", type=int, default=150,
                        help="Total training epochs")
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--warmup_epochs", type=int, default=30,
                        help="Phase 1 warmup epochs (mask + recon only)")
    parser.add_argument("--eval_interval", type=int, default=10,
                        help="Evaluate every N epochs")

    # Loss weights
    parser.add_argument("--mask_loss_weight", type=float, default=0.2)
    parser.add_argument("--recon_loss_weight", type=float, default=0.8,
                        help="Reconstruction loss weight (should be dominant)")
    parser.add_argument("--diffusion_loss_weight", type=float, default=0.1)
    parser.add_argument("--cluster_loss_weight", type=float, default=0.0,
                        help="DEC clustering loss weight (start at 0, increase after warmup)")

    # Evaluation
    parser.add_argument("--cluster_methods", type=str, default="kmeans,leiden",
                        help="Comma-separated clustering methods")
    parser.add_argument("--eval_both_embeddings", action="store_true", default=True,
                        help="Evaluate both direct and diffusion embeddings")

    # System
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=2,
                        help="GPU device ID (avoid 0 per user request)")
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--progress_bar", action="store_true", default=True,
                        help="Show tqdm progress bar")
    parser.add_argument("--resume", type=str, default=None,
                        help="Path to checkpoint .pt to resume from")
    parser.add_argument("--skip_train", action="store_true",
                        help="Skip training and load embeddings from save_dir")

    return parser.parse_args()


def main():
    args = parse_args()

    # ── Setup ─────────────────────────────────────────────────────────────────
    set_seed(args.seed)
    ensure_dir(args.save_dir)
    device = get_device(gpu=args.gpu, no_cuda=args.no_cuda)

    print("=" * 70)
    print("ScSpade — Support-Masked Latent Diffusion Autoencoder")
    print("=" * 70)
    print(f"Data:       {args.data_path}")
    print(f"Save dir:   {args.save_dir}")
    print(f"Device:     {device}")
    print(f"Latent dim: {args.latent_dim}")
    print(f"Epochs:     {args.epochs}")
    print(f"Weights:    mask={args.mask_loss_weight}, recon={args.recon_loss_weight}, "
          f"diff={args.diffusion_loss_weight}, cluster={args.cluster_loss_weight}")

    cluster_methods = tuple(m.strip() for m in args.cluster_methods.split(","))

    # ── Load data ─────────────────────────────────────────────────────────────
    t0 = time.time()
    bundle = load_sc_dataset(
        file_path=args.data_path,
        n_top_genes=args.n_top_genes,
        input_mode=args.input_mode,
        label_key=args.label_key,
        n_clusters=args.n_clusters,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
    )
    n_clusters = bundle.values.shape[0]  # Will be overridden by detected value
    # Re-detect n_clusters
    from collections import Counter
    n_clusters = len(Counter(bundle.labels))
    print(f"\n[Data] Loaded in {time.time()-t0:.1f}s: "
          f"{bundle.values.shape[0]} cells × {bundle.values.shape[1]} genes, "
          f"{n_clusters} clusters")

    # Store in adata
    bundle.adata.obs["scspade_true_label"] = bundle.labels.astype(str)

    # ── Training ─────────────────────────────────────────────────────────────
    if args.skip_train and os.path.exists(os.path.join(args.save_dir, "best_embeddings_direct.npy")):
        print("\n[Main] Loading saved embeddings (--skip_train)...")
        embeddings_direct = np.load(os.path.join(args.save_dir, "best_embeddings_direct.npy"))
        emb_diff_path = os.path.join(args.save_dir, "best_embeddings_diffusion.npy")
        embeddings_diffusion = (np.load(emb_diff_path) if os.path.exists(emb_diff_path)
                               else embeddings_direct)
        model = None
    else:
        t1 = time.time()

        train_result = run_scspade_training(
            bundle=bundle,
            n_clusters=n_clusters,
            latent_dim=args.latent_dim,
            hidden_dim=args.hidden_dim,
            diffusion_hidden_dim=args.diffusion_hidden_dim,
            diffusion_steps=args.diffusion_steps,
            dropout=args.dropout,
            mask_dropout=args.mask_dropout,
            epochs=args.epochs,
            batch_size=args.batch_size,
            lr=args.lr,
            weight_decay=args.weight_decay,
            warmup_epochs=args.warmup_epochs,
            mask_loss_weight=args.mask_loss_weight,
            recon_loss_weight=args.recon_loss_weight,
            diffusion_loss_weight=args.diffusion_loss_weight,
            cluster_loss_weight=args.cluster_loss_weight,
            eval_interval=args.eval_interval,
            save_dir=args.save_dir,
            device=device,
            seed=args.seed,
            progress_bar=args.progress_bar,
            eval_fn=make_eval_fn(n_clusters, cluster_methods),
        )

        model = train_result["model"]
        embeddings_direct = train_result["embeddings_direct"]
        embeddings_diffusion = train_result["embeddings_diffusion"]
        np.save(os.path.join(args.save_dir, "embeddings_direct.npy"), embeddings_direct.astype(np.float32))
        np.save(os.path.join(args.save_dir, "embeddings_diffusion.npy"), embeddings_diffusion.astype(np.float32))

        print(f"\n[Main] Training completed in {time.time()-t1:.1f}s")

    # ── Final evaluation ─────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("FINAL EVALUATION — Both Embedding Types")
    print("=" * 70)

    all_metrics = evaluate_and_save(
        embeddings_direct=embeddings_direct,
        embeddings_diffusion=embeddings_diffusion,
        labels=bundle.labels,
        n_clusters=n_clusters,
        save_dir=args.save_dir,
        cluster_methods=cluster_methods,
    )

    # ── Save adata ───────────────────────────────────────────────────────────
    bundle.adata.obsm["X_scspade_direct"] = embeddings_direct.astype(np.float32)
    bundle.adata.obsm["X_scspade_diffusion"] = embeddings_diffusion.astype(np.float32)

    best_emb_name = max(all_metrics.keys(), key=lambda k: all_metrics[k].get("nmi", 0))
    bundle.adata.obsm["X_scspade"] = (
        embeddings_direct if best_emb_name == "direct" else embeddings_diffusion
    ).astype(np.float32)

    adata_path = os.path.join(args.save_dir, "scspade_result.h5ad")
    bundle.adata.write_h5ad(adata_path)
    print(f"\n[Saved] AnnData to {adata_path}")

    # ── Save summary ────────────────────────────────────────────────────────
    summary = {
        "method": "ScSpade (cursor2_Doloris maskdiffusion)",
        "data_path": args.data_path,
        "n_cells": int(bundle.values.shape[0]),
        "n_genes": int(bundle.values.shape[1]),
        "n_clusters": int(n_clusters),
        "latent_dim": args.latent_dim,
        "epochs": args.epochs,
        "device": str(device),
        "weights": {
            "mask": args.mask_loss_weight,
            "recon": args.recon_loss_weight,
            "diffusion": args.diffusion_loss_weight,
            "cluster": args.cluster_loss_weight,
        },
        "metrics": all_metrics,
    }
    save_json(summary, os.path.join(args.save_dir, "summary.json"))

    # Also save a PlantNet-compatible metrics table
    if HAS_PLANTNET and plantnet_evaluation:
        for emb_name, emb in [("direct", embeddings_direct), ("diffusion", embeddings_diffusion)]:
            best_method = all_metrics[emb_name]["best_method"]
            pred = np.load(
                os.path.join(args.save_dir, f"scspade_{emb_name}_{best_method}_pred_labels.npy")
            )
            try:
                acc, nmi, ari, f1, fmi, vm, hom, com, _ = plantnet_evaluation(
                    bundle.labels, pred
                )
                row = {
                    "method": f"ScSpade_{emb_name}",
                    "ACC": float(acc), "NMI": float(nmi), "ARI": float(ari),
                    "F1-macro": float(f1), "FMI": float(fmi),
                    "V-measure": float(vm), "Homogeneity": float(hom),
                    "Completeness": float(com),
                }
                save_json(row, os.path.join(args.save_dir, f"metrics_{emb_name}.json"))
            except Exception as e:
                print(f"  Warning: PlantNet evaluation failed for {emb_name}: {e}")

    print(f"\n[Done] Results saved to {args.save_dir}")
    return all_metrics


if __name__ == "__main__":
    main()
