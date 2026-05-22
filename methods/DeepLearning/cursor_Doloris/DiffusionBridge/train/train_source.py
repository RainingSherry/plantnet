"""
Train source diffusion model.

The source diffusion learns the raw sparse expression distribution.
It learns to encode raw gene expression vectors into the shared Gaussian
latent space via a denoising objective.

Architecture:
    Input (n_genes) -> domain_encoder -> shared latent (shared_dim)
                        -> denoiser (time-conditioned) -> output (shared_dim)

Training objective:
    L = L_diffusion + recon_weight * L_reconstruction + prior_weight * L_prior

Usage:
    python train_source.py --data-path DATA.h5ad --output-dir ./output --n-hvg 2000 --latent-dim 64
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch
from torch.utils.data import DataLoader, TensorDataset

# Add parent to path for imports
_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from source_diffusion import SourceDiffusion


def read_h5ad_compat(path: str | Path):
    """Read h5ad with numpy compatibility fixes."""
    if not hasattr(np, "string_"):
        np.string_ = np.bytes_
    if not hasattr(np, "unicode_"):
        np.unicode_ = np.str_
    return sc.read_h5ad(str(path))


def preprocess_adata(adata, n_hvg: int = 2000):
    """
    Preprocess scRNA-seq data for the bridge model.

    Steps:
    1. Normalize total counts per cell (CPM-style, target_sum=1e4)
    2. log1p transformation
    3. Highly variable gene selection
    4. Z-score scaling per gene

    Returns:
        Preprocessed AnnData object (in-place modification)
    """
    adata = adata.copy()
    if adata.raw is None:
        adata.raw = adata.copy()

    # Per-cell normalization
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Log1p transformation
    sc.pp.log1p(adata)

    # HVG selection
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat",
        n_top_genes=n_hvg,
        subset=False,
    )
    if "highly_variable_rank" in adata.var.columns:
        adata = adata[:, adata.var["highly_variable_rank"] < n_hvg].copy()
    else:
        adata = adata[:, adata.var["highly_variable"]].copy()

    # Z-score scaling
    sc.pp.scale(adata)

    return adata


def build_teacher_target(adata, latent_dim: int, mode: str = "pca_graph"):
    """
    Build a teacher target embedding for the target diffusion.

    Three modes:
    - 'pca': raw PCA on normalized data
    - 'pca_graph': PCA followed by graph smoothing (recommended)
    - 'raw': use log1p normalized data directly

    The teacher target provides a "clean" reference that the bridge tries to match.
    """
    x = adata.X.toarray().astype(np.float32) if hasattr(adata.X, "toarray") else np.asarray(adata.X, dtype=np.float32)

    from sklearn.decomposition import PCA

    pca = PCA(n_components=latent_dim, random_state=42)
    z = pca.fit_transform(x)

    if mode == "pca_graph":
        # Build neighbor graph on PCA embedding and smooth
        adata_tmp = adata.copy()
        adata_tmp.obsm["X_pca"] = z
        adata_tmp.X = x
        sc.pp.neighbors(adata_tmp, n_neighbors=15, use_rep="X_pca")
        graph = adata_tmp.obsp["connectivities"]
        if sp.issparse(graph):
            graph = graph.toarray()
        z = graph @ z

    return np.asarray(z, dtype=np.float32)


def train_source(args):
    """Train the source diffusion model."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[Source] Device: {device}")

    # Load and preprocess data
    print(f"[Source] Loading data from {args.data_path}")
    adata = preprocess_adata(read_h5ad_compat(args.data_path), n_hvg=args.n_hvg)

    # Extract expression matrix
    x = adata.X.toarray().astype(np.float32) if hasattr(adata.X, "toarray") else np.asarray(
        adata.X, dtype=np.float32
    )
    print(f"[Source] Data shape: {x.shape}, sparsity: {(x == 0).mean():.4f}")

    # Handle NaN/Inf from scaling
    x = np.nan_to_num(x, nan=0.0, posinf=10.0, neginf=-10.0)

    dataset = TensorDataset(torch.from_numpy(x))
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True)

    # Build source diffusion model
    model = SourceDiffusion(
        domain_dim=x.shape[1],
        shared_dim=args.latent_dim,
        hidden_dim=args.hidden_dim,
        time_embed_dim=args.time_embed_dim,
        num_steps=args.diffusion_steps,
        dropout=args.dropout,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    total_params = sum(p.numel() for p in model.parameters())
    print(f"[Source] Model params: {total_params:,}")

    # Training loop
    for epoch in range(args.epochs):
        model.train()
        epoch_loss = 0.0
        epoch_diff = 0.0
        epoch_recon = 0.0
        epoch_prior = 0.0
        n_batches = 0

        for (batch_x,) in loader:
            batch_x = batch_x.to(device)
            loss_dict = model.training_loss(
                batch_x,
                recon_weight=args.recon_weight,
                prior_weight=args.prior_weight,
                zero_weight=args.zero_weight,
            )
            loss = loss_dict["loss"]

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            epoch_loss += loss.item()
            epoch_diff += loss_dict["diffusion_loss"].item()
            epoch_recon += loss_dict["recon_loss"].item()
            epoch_prior += loss_dict["prior_loss"].item()
            n_batches += 1

        avg_loss = epoch_loss / n_batches
        avg_diff = epoch_diff / n_batches
        avg_recon = epoch_recon / n_batches
        avg_prior = epoch_prior / n_batches

        print(
            f"[Source] Epoch {epoch+1}/{args.epochs} | "
            f"loss={avg_loss:.6f} | diff={avg_diff:.6f} | "
            f"recon={avg_recon:.6f} | prior={avg_prior:.6f}"
        )

    # Save model
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ckpt = {
        "state_dict": model.state_dict(),
        "input_dim": x.shape[1],
        "n_genes": x.shape[1],
        "latent_dim": args.latent_dim,
        "hidden_dim": args.hidden_dim,
        "time_embed_dim": args.time_embed_dim,
        "diffusion_steps": args.diffusion_steps,
        "args": vars(args),
    }
    torch.save(ckpt, output_dir / "source.pt")
    print(f"[Source] Saved to {output_dir / 'source.pt'}")

    # Also save teacher target for downstream use
    print(f"[Source] Building teacher target (mode={args.teacher_mode})...")
    teacher_z = build_teacher_target(adata, latent_dim=args.latent_dim, mode=args.teacher_mode)
    np.save(output_dir / "teacher_target.npy", teacher_z)
    print(f"[Source] Teacher target saved: shape={teacher_z.shape}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train source diffusion model")
    # Data
    parser.add_argument("--data-path", type=str, required=True, help="Path to h5ad file")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--n-hvg", type=int, default=2000, help="Number of highly variable genes")
    # Model
    parser.add_argument("--latent-dim", type=int, default=64, help="Shared latent dimension")
    parser.add_argument("--hidden-dim", type=int, default=256, help="Hidden dimension")
    parser.add_argument("--time-embed-dim", type=int, default=128, help="Time embedding dimension")
    parser.add_argument("--diffusion-steps", type=int, default=50, help="Number of diffusion steps")
    parser.add_argument("--dropout", type=float, default=0.0, help="Dropout rate")
    # Training
    parser.add_argument("--batch-size", type=int, default=256, help="Batch size")
    parser.add_argument("--epochs", type=int, default=20, help="Number of epochs")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="Weight decay")
    # Loss weights
    parser.add_argument("--recon-weight", type=float, default=1.0, help="Reconstruction loss weight")
    parser.add_argument("--prior-weight", type=float, default=1e-3, help="Prior loss weight")
    parser.add_argument("--zero-weight", type=float, default=0.25, help="Zero-value weight in reconstruction")
    # Teacher
    parser.add_argument("--teacher-mode", type=str, default="pca_graph", choices=["pca", "pca_graph"])

    args = parser.parse_args()
    train_source(args)
