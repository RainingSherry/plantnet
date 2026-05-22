"""
Main pipeline runner for Dual-Domain Diffusion Bridge.

Runs the full three-stage pipeline on a given dataset:
    Stage 1: Train source diffusion (raw -> shared latent)
    Stage 2: Train target diffusion (shared latent -> cluster domain)
    Stage 3: Train bridge + cluster head jointly

Usage:
    python run_pipeline.py \
        --data-path /path/to/data.h5ad \
        --output-dir ./results/SRP182008 \
        --n-hvg 2000 \
        --latent-dim 64 \
        --source-epochs 20 \
        --target-epochs 20 \
        --bridge-epochs 30 \
        --batch-size 256
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str]):
    """Run a command and stream output."""
    print(f"\n{'='*70}")
    print(f"RUNNING: {' '.join(cmd)}")
    print(f"{'='*70}")
    result = subprocess.run(cmd, check=False, cwd=str(Path(__file__).resolve().parent))
    if result.returncode != 0:
        print(f"FAILED with exit code {result.returncode}")
        sys.exit(result.returncode)
    else:
        print(f"COMPLETED successfully")


def run_pipeline(args):
    """Run the full three-stage pipeline."""
    root = Path(__file__).resolve().parent
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Stage 1: Source diffusion
    print("\n" + "="*70)
    print("STAGE 1: Training Source Diffusion")
    print("="*70)
    run([
        sys.executable,
        str(root / "train" / "train_source.py"),
        "--data-path", args.data_path,
        "--output-dir", str(output_dir),
        "--n-hvg", str(args.n_hvg),
        "--latent-dim", str(args.latent_dim),
        "--hidden-dim", str(args.hidden_dim),
        "--time-embed-dim", str(args.time_embed_dim),
        "--diffusion-steps", str(args.diffusion_steps),
        "--dropout", str(args.dropout),
        "--batch-size", str(args.batch_size),
        "--epochs", str(args.source_epochs),
        "--lr", str(args.lr),
        "--weight-decay", str(args.weight_decay),
        "--recon-weight", str(args.recon_weight),
        "--prior-weight", str(args.prior_weight),
        "--zero-weight", str(args.zero_weight),
        "--teacher-mode", str(args.teacher_mode),
    ])

    # Stage 2: Target diffusion
    print("\n" + "="*70)
    print("STAGE 2: Training Target Diffusion")
    print("="*70)
    run([
        sys.executable,
        str(root / "train" / "train_target.py"),
        "--data-path", args.data_path,
        "--output-dir", str(output_dir),
        "--n-hvg", str(args.n_hvg),
        "--latent-dim", str(args.latent_dim),
        "--hidden-dim", str(args.hidden_dim),
        "--time-embed-dim", str(args.time_embed_dim),
        "--diffusion-steps", str(args.diffusion_steps),
        "--dropout", str(args.dropout),
        "--batch-size", str(args.batch_size),
        "--epochs", str(args.target_epochs),
        "--lr", str(args.lr),
        "--weight-decay", str(args.weight_decay),
        "--recon-weight", str(args.recon_weight),
        "--prior-weight", str(args.prior_weight),
        "--teacher-mode", str(args.teacher_mode),
    ])

    # Stage 3: Bridge + cluster head
    print("\n" + "="*70)
    print("STAGE 3: Training Bridge + Cluster Head")
    print("="*70)
    run([
        sys.executable,
        str(root / "train" / "train_bridge_cluster.py"),
        "--data-path", args.data_path,
        "--output-dir", str(output_dir),
        "--n-hvg", str(args.n_hvg),
        "--latent-dim", str(args.latent_dim),
        "--hidden-dim", str(args.hidden_dim),
        "--time-embed-dim", str(args.time_embed_dim),
        "--diffusion-steps", str(args.diffusion_steps),
        "--dropout", str(args.dropout),
        "--cluster-hidden-dim", str(args.cluster_hidden_dim),
        "--batch-size", str(args.batch_size),
        "--epochs", str(args.bridge_epochs),
        "--lr", str(args.lr),
        "--weight-decay", str(args.weight_decay),
        "--warmup-epochs", str(args.warmup_epochs),
        "--teacher-weight", str(args.teacher_weight),
        "--cluster-weight", str(args.cluster_weight),
        "--gaussian-weight", str(args.gaussian_weight),
        "--support-weight", str(args.support_weight),
        "--entropy-weight", str(args.entropy_weight),
        "--support-topk", str(args.support_topk),
        "--support-blend", str(args.support_blend),
        "--teacher-mode", str(args.teacher_mode),
    ])

    print("\n" + "="*70)
    print("PIPELINE COMPLETE")
    print(f"Results saved to: {output_dir}")
    print("="*70)
    for f in sorted(output_dir.glob("*")):
        size_mb = f.stat().st_size / 1024 / 1024
        print(f"  {f.name} ({size_mb:.1f} MB)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Run Dual-Domain Diffusion Bridge pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Data
    parser.add_argument("--data-path", type=str, required=True, help="Path to h5ad file")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--n-hvg", type=int, default=2000, help="Number of highly variable genes")

    # Architecture
    parser.add_argument("--latent-dim", type=int, default=64, help="Shared latent dimension")
    parser.add_argument("--hidden-dim", type=int, default=256, help="Hidden dimension")
    parser.add_argument("--time-embed-dim", type=int, default=128, help="Time embedding dimension")
    parser.add_argument("--diffusion-steps", type=int, default=50, help="Number of diffusion steps")
    parser.add_argument("--dropout", type=float, default=0.0, help="Dropout rate")
    parser.add_argument("--cluster-hidden-dim", type=int, default=256, help="Cluster head hidden dimension")

    # Training
    parser.add_argument("--source-epochs", type=int, default=20, help="Source model epochs")
    parser.add_argument("--target-epochs", type=int, default=20, help="Target model epochs")
    parser.add_argument("--bridge-epochs", type=int, default=30, help="Bridge+cluster epochs")
    parser.add_argument("--batch-size", type=int, default=256, help="Batch size")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="Weight decay")
    parser.add_argument("--warmup-epochs", type=int, default=5, help="Warmup epochs for cluster loss")

    # Loss weights
    parser.add_argument("--recon-weight", type=float, default=1.0, help="Reconstruction loss weight (source/target)")
    parser.add_argument("--prior-weight", type=float, default=1e-3, help="Prior loss weight")
    parser.add_argument("--zero-weight", type=float, default=0.25, help="Zero-value weight in reconstruction")
    parser.add_argument("--teacher-weight", type=float, default=1.0, help="Teacher loss weight (bridge)")
    parser.add_argument("--cluster-weight", type=float, default=1.0, help="Cluster loss weight")
    parser.add_argument("--gaussian-weight", type=float, default=1e-3, help="Gaussian prior loss weight")
    parser.add_argument("--support-weight", type=float, default=0.5, help="Support loss weight")
    parser.add_argument("--entropy-weight", type=float, default=1e-3, help="Entropy loss weight")

    # Support mask
    parser.add_argument("--support-topk", type=int, default=256, help="Top-k genes in support mask")
    parser.add_argument("--support-blend", type=float, default=0.2, help="Support blend weight")

    # Teacher
    parser.add_argument("--teacher-mode", type=str, default="pca_graph", choices=["pca", "pca_graph"])

    args = parser.parse_args()
    run_pipeline(args)
