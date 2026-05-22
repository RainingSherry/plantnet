"""
Train bridge + cluster head jointly.

This is the core paper-like stage. The bridge connects:
    Raw Sparse Count Domain  --[DDIM reverse]-->  Shared Gaussian Latent  --[DDIM sample]-->  Cluster-Separable Denoised Domain

The joint training optimizes:
    L_total = L_teacher + L_cluster + L_latent_prior + L_support + L_entropy

Where:
    - L_teacher: MSE between predicted embedding and teacher target (from PCA/graph smoothing)
    - L_cluster: KL-divergence between DEC soft assignments and target distribution
    - L_latent_prior: encourages shared latent to be near isotropic Gaussian
    - L_support: MSE between support anchor and teacher target
    - L_entropy: entropy regularization to prevent premature cluster collapse

Mathematical narrative:
    Traditional dimensionality reduction is a lossy compression -- it projects all cells
    into a shared manifold without distinguishing cluster boundaries. Our bridge-based
    representation learning instead constructs a path through the high-dimensional
    space: raw sparse observations are first encoded into a Gaussian latent space
    (which preserves distributional information), then decoded into a cluster-separable
    denoised space (which separates cell types). The shared Gaussian intermediate
    serves as an information bottleneck that prevents the model from memorizing
    expression patterns, while the target diffusion learns to reconstruct the
    cluster manifold.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.metrics import (
    adjusted_rand_score,
    normalized_mutual_info_score,
    silhouette_score,
)
from torch.utils.data import DataLoader, TensorDataset

_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_ROOT))

from bridge import BridgeSampler, DiffusionBridge
from cluster_head import ClusterHead
from source_diffusion import SourceDiffusion, TargetDiffusion
from support_mask import GeneSupportMask, build_support_mask


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


def extract_labels(adata):
    """Extract cell type labels from AnnData.obs."""
    candidates = [
        "cell_type", "Celltype", "celltype", "cell_label",
        "label", "Seurat_clusters", "celltype_after", "Cluster", "cluster", "cell_type_label",
    ]
    for key in candidates:
        if key in adata.obs.columns:
            return adata.obs[key].astype(str).to_numpy(), key
    raise KeyError(f"No label column found. Available columns: {list(adata.obs.columns)}")


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


def load_bridge(args, input_dim: int, device: torch.device):
    """Load pre-trained source and target diffusion models, build bridge."""
    source = SourceDiffusion(
        domain_dim=input_dim,
        shared_dim=args.latent_dim,
        hidden_dim=args.hidden_dim,
        time_embed_dim=args.time_embed_dim,
        num_steps=args.diffusion_steps,
        dropout=args.dropout,
    )
    target = TargetDiffusion(
        domain_dim=args.latent_dim,
        shared_dim=args.latent_dim,
        hidden_dim=args.hidden_dim,
        time_embed_dim=args.time_embed_dim,
        cond_dim=args.latent_dim,
        num_steps=args.diffusion_steps,
        dropout=args.dropout,
    )

    # Load pre-trained weights
    source_ckpt = torch.load(
        Path(args.output_dir) / "source.pt", map_location="cpu", weights_only=True
    )
    target_ckpt = torch.load(
        Path(args.output_dir) / "target.pt", map_location="cpu", weights_only=True
    )
    source.load_state_dict(source_ckpt["state_dict"])
    target.load_state_dict(target_ckpt["state_dict"])

    bridge = DiffusionBridge(
        source=source,
        target=target,
        support_mask=GeneSupportMask(blend=args.support_blend),
        support_hidden_dim=args.hidden_dim,
    ).to(device)

    return bridge


def dec_target_distribution(q: torch.Tensor) -> torch.Tensor:
    """DEC target distribution: sharpen confident assignments."""
    weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
    return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)


def train_bridge_cluster(args):
    """Joint training of bridge and cluster head."""
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"[Bridge] Device: {device}")

    # Load and preprocess
    print(f"[Bridge] Loading data from {args.data_path}")
    adata = preprocess_adata(read_h5ad_compat(args.data_path), n_hvg=args.n_hvg)
    labels, label_key = extract_labels(adata)

    # Encode labels
    unique_labels = sorted(np.unique(labels).tolist())
    label_to_id = {label: idx for idx, label in enumerate(unique_labels)}
    y = np.array([label_to_id[label] for label in labels], dtype=np.int64)
    n_clusters = len(unique_labels)
    print(f"[Bridge] Cells: {len(y)}, Genes: {adata.n_vars}, Clusters: {n_clusters}")
    print(f"[Bridge] Label column: {label_key}, values: {unique_labels}")

    # Raw expression matrix
    x = adata.X.toarray().astype(np.float32) if hasattr(adata.X, "toarray") else np.asarray(
        adata.X, dtype=np.float32
    )
    x = np.nan_to_num(x, nan=0.0, posinf=10.0, neginf=-10.0)
    print(f"[Bridge] X shape: {x.shape}, sparsity: {(x == 0).mean():.4f}")

    # Load teacher target
    teacher_path = Path(args.output_dir) / "teacher_target.npy"
    if teacher_path.exists():
        teacher_target = np.load(teacher_path).astype(np.float32)
        teacher_target = teacher_target[: x.shape[0]]
    else:
        print("[Bridge] Teacher target not found, building on-the-fly...")
        teacher_target = build_teacher_target(adata, latent_dim=args.latent_dim, mode=args.teacher_mode)
    print(f"[Bridge] Teacher target shape: {teacher_target.shape}")

    # Load bridge
    bridge = load_bridge(args, input_dim=x.shape[1], device=device)
    bridge.train()
    sampler = BridgeSampler(bridge)

    # Cluster head
    cluster_head = ClusterHead(
        input_dim=args.latent_dim,
        n_clusters=n_clusters,
        hidden_dim=args.cluster_hidden_dim,
        alpha=1.0,
    ).to(device)

    # Initialize prototypes from teacher target via KMeans
    print(f"[Bridge] Initializing cluster prototypes from teacher target...")
    teacher_tensor = torch.from_numpy(teacher_target).to(device).float()
    teacher_kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=20).fit(
        teacher_target
    )

    with torch.no_grad():
        encoded_teacher = cluster_head.encoder(teacher_tensor)
        encoded_centers = []
        for cluster_id in range(n_clusters):
            mask = teacher_kmeans.labels_ == cluster_id
            if mask.sum() > 0:
                encoded_centers.append(encoded_teacher[mask].mean(dim=0))
            else:
                encoded_centers.append(encoded_teacher.mean(dim=0))
        cluster_head.prototypes.data.copy_(torch.stack(encoded_centers, dim=0))

    print(f"[Bridge] Cluster prototypes initialized. Sample assignment distribution:")
    with torch.no_grad():
        q_init = cluster_head(teacher_tensor[:1000].to(device))
        print(f"  Mean max prob: {q_init.max(dim=1)[0].mean():.4f}")
        print(f"  Mean entropy: {(-q_init * q_init.clamp_min(1e-8).log()).sum(dim=1).mean():.4f}")

    # Dataset
    dataset = TensorDataset(
        torch.from_numpy(x),
        torch.from_numpy(y),
        torch.from_numpy(teacher_target),
    )
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=False)

    # Optimizer
    optimizer = torch.optim.AdamW(
        list(bridge.parameters()) + list(cluster_head.parameters()),
        lr=args.lr,
        weight_decay=args.weight_decay,
    )

    # Learning rate scheduler
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer, T_max=args.epochs, eta_min=args.lr * 0.01
    )

    # Training loop
    print(f"[Bridge] Starting joint training for {args.epochs} epochs...")
    for epoch in range(args.epochs):
        bridge.train()
        cluster_head.train()

        epoch_loss = 0.0
        epoch_teacher = 0.0
        epoch_cluster = 0.0
        epoch_prior = 0.0
        epoch_support = 0.0
        epoch_entropy = 0.0
        n_batches = 0

        for batch_x, batch_y, batch_teacher in loader:
            batch_x = batch_x.to(device)
            batch_teacher = batch_teacher.to(device)

            # Build support mask from raw expression
            raw_mask = build_support_mask(batch_x, topk=args.support_topk)

            # --- Bridge forward ---
            # Step 1: DDIM reverse (encode raw -> shared latent)
            latent = sampler.encode(batch_x, raw_mask=raw_mask)

            # Step 2: DDIM sample (decode latent -> target embedding)
            output = bridge(batch_x, raw_mask=raw_mask)
            pred_z = output.target_embedding
            support_anchor = output.support_anchor

            # --- Cluster head ---
            q = cluster_head(pred_z)
            p = dec_target_distribution(q.detach())

            # --- Losses ---
            # Teacher loss: match bridge output to clean teacher target
            teacher_loss = F.mse_loss(pred_z, batch_teacher)

            # Cluster loss: DEC-style self-training
            cluster_loss = F.kl_div(
                q.clamp_min(1e-8).log(), p, reduction="batchmean"
            )

            # Latent prior: keep shared latent near Gaussian
            latent_prior_loss = latent.pow(2).mean()

            # Support loss: align support anchor with teacher
            support_loss = F.mse_loss(support_anchor, batch_teacher)

            # Entropy regularization: prevent premature cluster collapse
            entropy = -(q * q.clamp_min(1e-8).log()).sum(dim=1).mean()

            # Warmup: gradually increase cluster weight
            cluster_warmup = min(1.0, (epoch + 1) / max(1, args.warmup_epochs))
            entropy_warmup = min(1.0, (epoch + 1) / max(1, args.warmup_epochs // 2))

            loss = (
                args.teacher_weight * teacher_loss
                + cluster_warmup * args.cluster_weight * cluster_loss
                + args.gaussian_weight * latent_prior_loss
                + args.support_weight * support_loss
                + entropy_warmup * args.entropy_weight * entropy
            )

            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(
                list(bridge.parameters()) + list(cluster_head.parameters()),
                max_norm=1.0,
            )
            optimizer.step()

            epoch_loss += loss.item()
            epoch_teacher += teacher_loss.item()
            epoch_cluster += cluster_loss.item()
            epoch_prior += latent_prior_loss.item()
            epoch_support += support_loss.item()
            epoch_entropy += entropy.item()
            n_batches += 1

        scheduler.step()

        avg_loss = epoch_loss / n_batches
        avg_teacher = epoch_teacher / n_batches
        avg_cluster = epoch_cluster / n_batches
        avg_prior = epoch_prior / n_batches
        avg_support = epoch_support / n_batches
        avg_entropy = epoch_entropy / n_batches

        print(
            f"[Bridge] Epoch {epoch+1}/{args.epochs} | "
            f"loss={avg_loss:.4f} | teacher={avg_teacher:.4f} | "
            f"cluster={avg_cluster:.4f} | prior={avg_prior:.4f} | "
            f"support={avg_support:.4f} | entropy={avg_entropy:.4f} | "
            f"lr={scheduler.get_last_lr()[0]:.6f}"
        )

    # --- Evaluation ---
    print("\n[Bridge] Evaluating...")
    bridge.eval()
    cluster_head.eval()

    with torch.no_grad():
        all_x = torch.from_numpy(x).to(device)
        eval_mask = build_support_mask(all_x, topk=args.support_topk)
        output = bridge(all_x, raw_mask=eval_mask)
        pred_embedding = output.target_embedding.cpu().numpy()

        # DEC predictions
        probs = cluster_head(output.target_embedding).cpu().numpy()
        pred_dec = probs.argmax(axis=1)

        # KMeans on bridge embeddings
        kmeans_pred = KMeans(n_clusters=n_clusters, random_state=42, n_init=20).fit_predict(
            pred_embedding
        )

    # Compute metrics
    ari_dec = adjusted_rand_score(y, pred_dec)
    nmi_dec = normalized_mutual_info_score(y, pred_dec)
    ari_kmeans = adjusted_rand_score(y, kmeans_pred)
    nmi_kmeans = normalized_mutual_info_score(y, kmeans_pred)
    sil = silhouette_score(pred_embedding, pred_dec) if len(np.unique(pred_dec)) > 1 else float("nan")

    print(f"\n{'='*60}")
    print(f"[Bridge] Final Results ({label_key})")
    print(f"{'='*60}")
    print(f"  DEC-based:   ARI={ari_dec:.4f}, NMI={nmi_dec:.4f}")
    print(f"  KMeans:      ARI={ari_kmeans:.4f}, NMI={nmi_kmeans:.4f}")
    print(f"  Silhouette:  {sil:.4f}")
    print(f"{'='*60}\n")

    # Save outputs
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    torch.save(
        {
            "bridge": bridge.state_dict(),
            "cluster_head": cluster_head.state_dict(),
            "args": vars(args),
        },
        output_dir / "bridge_cluster.pt",
    )

    np.save(output_dir / "bridge_embedding.npy", pred_embedding)
    np.save(output_dir / "bridge_dec_labels.npy", pred_dec)
    np.save(output_dir / "bridge_kmeans_labels.npy", kmeans_pred)
    np.save(output_dir / "bridge_probs.npy", probs)

    print(f"[Bridge] Results saved to {output_dir}/")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Train bridge + cluster head")
    parser.add_argument("--data-path", type=str, required=True, help="Path to h5ad file")
    parser.add_argument("--output-dir", type=str, required=True, help="Output directory")
    parser.add_argument("--n-hvg", type=int, default=2000, help="Number of highly variable genes")
    # Model architecture
    parser.add_argument("--latent-dim", type=int, default=64, help="Shared latent dimension")
    parser.add_argument("--hidden-dim", type=int, default=256, help="Hidden dimension")
    parser.add_argument("--time-embed-dim", type=int, default=128, help="Time embedding dimension")
    parser.add_argument("--diffusion-steps", type=int, default=50, help="Number of diffusion steps")
    parser.add_argument("--dropout", type=float, default=0.0, help="Dropout rate")
    parser.add_argument("--cluster-hidden-dim", type=int, default=256, help="Cluster head hidden dim")
    # Support mask
    parser.add_argument("--support-topk", type=int, default=256, help="Top-k genes in support mask")
    parser.add_argument("--support-blend", type=float, default=0.2, help="Support blend weight")
    # Training
    parser.add_argument("--batch-size", type=int, default=256, help="Batch size")
    parser.add_argument("--epochs", type=int, default=30, help="Number of epochs")
    parser.add_argument("--lr", type=float, default=1e-3, help="Learning rate")
    parser.add_argument("--weight-decay", type=float, default=1e-4, help="Weight decay")
    parser.add_argument("--warmup-epochs", type=int, default=5, help="Warmup epochs for cluster loss")
    # Loss weights
    parser.add_argument("--teacher-weight", type=float, default=1.0, help="Teacher loss weight")
    parser.add_argument("--cluster-weight", type=float, default=1.0, help="Cluster loss weight")
    parser.add_argument("--gaussian-weight", type=float, default=1e-3, help="Gaussian prior loss weight")
    parser.add_argument("--support-weight", type=float, default=0.5, help="Support loss weight")
    parser.add_argument("--entropy-weight", type=float, default=1e-3, help="Entropy loss weight")
    # Teacher
    parser.add_argument("--teacher-mode", type=str, default="pca_graph", choices=["pca", "pca_graph"])

    args = parser.parse_args()
    train_bridge_cluster(args)
