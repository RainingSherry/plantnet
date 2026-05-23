"""
Train target diffusion model.

The target diffusion learns the cluster-friendly denoised representation distribution.
It trains on teacher embeddings (from PCA or PCA+graph smoothing) and learns
to sample from the cluster-separable manifold in the shared latent space.

Training objective:
    L = L_diffusion + recon_weight * L_reconstruction + prior_weight * L_prior

The target diffusion only needs to be trained once per dataset.
It provides the "bridge head" for the cluster-friendly domain.
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

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from target_diffusion import TargetDiffusion


def read_h5ad_compat(path: str | Path):
    if not hasattr(np, "string_"):
        np.string_ = np.bytes_
    if not hasattr(np, "unicode_"):
        np.unicode_ = np.str_
    return sc.read_h5ad(str(path))


def preprocess_adata(adata, n_hvg: int = 2000):
    """Preprocess scRNA-seq data. Must be consistent with train_source.py."""
    adata = adata.copy()
    if adata.raw is None:
        adata.raw = adata.copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

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

    sc.pp.scale(adata)
    return adata


def build_teacher_target(adata, latent_dim: int, mode: str = "pca_graph"):
    """Build teacher target embedding. Must be consistent with train_source.py."""
    x = adata.X.toarray().astype(np.float32) if hasattr(adata.X, "toarray") else np.asarray(
        adata.X, dtype=np.float32
    )

    from sklearn.decomposition import PCA

    pca = PCA(n_components=latent_dim, random_state=42)
    z = pca.fit_transform(x)

    if mode == "pca_graph":
        adata_tmp = adata.copy()
        adata_tmp.X = x  # restore original expression for graph building
        adata_tmp.obsm["X_pca"] = z
        sc.pp.neighbors(adata_tmp, n_neighbors=15, use_rep="X_pca")
        graph = adata_tmp.obsp["connectivities"]
        if sp.issparse(graph):
            graph = graph.toarray()
        z = graph @ z

    return np.asarray(z, dtype=np.float32)


def train_target(args):
    """Train the target diffusion model."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[Target] Device: {device}")

    # Load and preprocess
    print(f"[Target] Loading data from {args.data_path}")
    adata = preprocess_adata(read_h5ad_compat(args.data_path), n_hvg=args.n_hvg)

    # Build teacher target
    teacher_z = build_teacher_target(adata, latent_dim=args.latent_dim, mode=args.teacher_mode)
    print(f"[Target] Teacher target shape: {teacher_z.shape}")

    # Handle NaN/Inf
    teacher_z = np.nan_to_num(teacher_z, nan=0.0, posinf=10.0, neginf=-10.0)

    dataset = TensorDataset(torch.from_numpy(teacher_z))
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True)

    # Build target diffusion model
    model = TargetDiffusion(
        domain_dim=args.latent_dim,  # domain_dim = latent_dim (teacher target dim)
        shared_dim=args.latent_dim,  # shared_dim = latent_dim
        hidden_dim=args.hidden_dim,
        time_embed_dim=args.time_embed_dim,
        cond_dim=args.latent_dim,  # condition on support anchor of same dim
        num_steps=args.diffusion_steps,
        dropout=args.dropout,
    ).to(device)

    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    total_params = sum(p.numel() for p in model.parameters())
    print(f"[Target] Model params: {total_params:,}")

    # Training loop
    for epoch in range(args.epochs):
        model.train()
        epoch_loss = 0.0
        epoch_diff = 0.0
        epoch_recon = 0.0
        n_batches = 0

        for (batch_z,) in loader:
            batch_z = batch_z.to(device)
            loss_dict = model.training_loss(
                batch_z,
                recon_weight=args.recon_weight,
                prior_weight=args.prior_weight,
            )
            loss = loss_dict["loss"]

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
            optimizer.step()

            epoch_loss += loss.item()
            epoch_diff += loss_dict["diffusion_loss"].item()
            epoch_recon += loss_dict["recon_loss"].item()
            n_batches += 1

        avg_loss = epoch_loss / n_batches
        avg_diff = epoch_diff / n_batches
        avg_recon = epoch_recon / n_batches

        print(
            f"[Target] Epoch {epoch+1}/{args.epochs} | "
            f"loss={avg_loss:.6f} | diff={avg_diff:.6f} | recon={avg_recon:.6f}"
        )

    # Save model
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    ckpt = {
        "state_dict": model.state_dict(),
        "latent_dim": args.latent_dim,
        "hidden_dim": args.hidden_dim,
        "time_embed_dim": args.time_embed_dim,
        "diffusion_steps": args.diffusion_steps,
        "args": vars(args),
    }
    torch.save(ckpt, output_dir / "target.pt")
    print(f"[Target] Saved to {output_dir / 'target.pt'}")

    # Ensure teacher target is also saved (in case source didn't run)
    np.save(output_dir / "teacher_target.npy", teacher_z)
    print(f"[Target] Teacher target saved: shape={teacher_z.shape}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train target diffusion model")
    parser.add_argument("--data-path", type=str, required=True, help="Path to h5ad file")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--n-hvg", type=int, default=2000, help="Number of highly variable genes")
    parser.add_argument("--latent-dim", type=int, default=64, help="Latent dimension (also teacher dim)")
    parser.add_argument("--hidden-dim", type=int, default=256, help="Hidden dimension")
    parser.add_argument("--time-embed-dim", type=int, default=128, help="Time embedding dimension")
    parser.add_argument("--diffusion-steps", type=int, default=50, help="Number of diffusion steps")
    parser.add_argument("--dropout", type=float, default=0.0, help="Dropout rate")
    parser.add_argument("--batch-size", type=int, default=256, help="Batch size")
    parser.add_argument("--epochs", type=int, default=20, help="Number of epochs")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="Weight decay")
    parser.add_argument("--recon-weight", type=float, default=1.0, help="Reconstruction loss weight")
    parser.add_argument("--prior-weight", type=float, default=1e-3, help="Prior loss weight")
    parser.add_argument("--teacher-mode", type=str, default="pca_graph", choices=["pca", "pca_graph"])

    args = parser.parse_args()
    train_target(args)
