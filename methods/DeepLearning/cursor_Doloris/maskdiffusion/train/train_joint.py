# -*- coding: utf-8 -*-
"""
Stage C: Joint Training of SupportMaskNet + LatentDiffusionAE

This script trains both models jointly, where:
1. SupportMaskNet learns gene activation patterns
2. LatentDiffusionAE learns denoised representations
3. The two are coupled: mask guides diffusion, and diffusion improves mask prediction

Usage:
    python train_joint.py --data_path /path/to/data.h5ad --save_dir ./results_joint
"""

import os
import sys
import argparse
import random
import math
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from sklearn.mixture import GaussianMixture

# Add project root to path
PROJECT_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
sys.path.insert(0, PROJECT_ROOT)
from preprocess import prepare_data_for_model
from models.support_mask import SupportMaskNet
from models.latent_diffusion import LatentDiffusionAE


def set_seed(seed: int = 42):
    """Set random seeds for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


class AverageMeter:
    """Computes and stores the average and current value."""

    def __init__(self):
        self.reset()

    def reset(self):
        self.val = 0
        self.avg = 0
        self.sum = 0
        self.count = 0

    def update(self, val, n=1):
        self.val = val
        self.sum += val * n
        self.count += n
        self.avg = self.sum / self.count


def parse_args():
    parser = argparse.ArgumentParser(
        description='Joint training of SupportMaskNet + LatentDiffusionAE',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Data
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results_joint',
                        help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, default=0,
                        help='Number of clusters (0 = auto-detect from data)')

    # Mask Model
    parser.add_argument('--mask_hidden_dims', type=str, default='512,256,128',
                        help='Comma-separated hidden layer dimensions for mask model')
    parser.add_argument('--mask_dropout', type=float, default=0.1,
                        help='Dropout rate for mask model')

    # Diffusion Model
    parser.add_argument('--latent_dim', type=int, default=32,
                        help='Latent space dimension')
    parser.add_argument('--diffusion_hidden_dims', type=str, default='512,256',
                        help='Comma-separated encoder/decoder hidden dimensions')
    parser.add_argument('--diffusion_steps', type=int, default=100,
                        help='Number of diffusion timesteps')
    parser.add_argument('--diffusion_type', type=str, default='ddpm',
                        choices=['ddpm', 'ddim'],
                        help='Diffusion model type')
    parser.add_argument('--diffusion_dropout', type=float, default=0.1,
                        help='Dropout rate for diffusion model')

    # Training
    parser.add_argument('--epochs', type=int, default=200,
                        help='Number of training epochs')
    parser.add_argument('--batch_size', type=int, default=128,
                        help='Batch size')
    parser.add_argument('--lr', type=float, default=1e-3,
                        help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-4,
                        help='Weight decay')

    # Loss weights
    parser.add_argument('--mask_loss_weight', type=float, default=1.0,
                        help='Weight for mask loss')
    parser.add_argument('--diffusion_loss_weight', type=float, default=1.0,
                        help='Weight for diffusion loss')
    parser.add_argument('--recon_loss_weight', type=float, default=0.1,
                        help='Weight for reconstruction loss')
    parser.add_argument('--cluster_loss_weight', type=float, default=0.01,
                        help='Weight for clustering loss')

    # Warmup
    parser.add_argument('--warmup_epochs', type=int, default=20,
                        help='Number of warmup epochs (mask-only training)')

    # Evaluation
    parser.add_argument('--eval_interval', type=int, default=10,
                        help='Evaluation interval (epochs)')

    # Other
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--gpu', type=int, default=0,
                        help='GPU device number')
    parser.add_argument('--no_cuda', action='store_true',
                        help='Disable CUDA')

    return parser.parse_args()


class ClusterLoss(nn.Module):
    """
    Enhanced clustering loss with DEC-style + contrastive learning.

    Combines:
    1. DEC-style soft assignment loss
    2. Contrastive loss to pull similar cells together
    3. Separability loss to push different clusters apart
    """

    def __init__(self, n_clusters: int, latent_dim: int, temperature: float = 0.1):
        super().__init__()
        self.n_clusters = n_clusters
        self.temperature = temperature

        # Learnable cluster centers
        self.register_parameter(
            'cluster_centers',
            nn.Parameter(torch.randn(n_clusters, latent_dim) * 0.1)
        )

        # Learnable temperature for contrastive loss
        self.log_temp = nn.Parameter(torch.tensor(0.0))

    def forward(self, z: torch.Tensor, return_assignments: bool = False) -> tuple:
        """
        Args:
            z: Latent embeddings (batch, latent_dim)
            return_assignments: Whether to return soft assignments

        Returns:
            (soft_assignments, total_cluster_loss)
        """
        # 1. DEC-style soft assignment loss
        dist = torch.cdist(z, self.cluster_centers)  # (batch, n_clusters)
        q = torch.softmax(-dist / (self.temperature + 1e-8), dim=-1)

        # Target distribution P
        p = q ** 2 / (q.sum(dim=0, keepdim=True) + 1e-8)
        p = p / (p.sum(dim=1, keepdim=True) + 1e-8)

        # KL divergence loss
        kl_loss = torch.nn.functional.kl_div(
            torch.log(q + 1e-10),
            p,
            reduction='batchmean'
        )

        # 2. Center contrastive loss - pull embeddings toward their assigned center
        assigned_centers = torch.index_select(
            self.cluster_centers, 0, q.argmax(dim=1)
        )
        center_loss = ((z - assigned_centers) ** 2).mean()

        # 3. Separability loss - push cluster centers apart
        center_dist = torch.cdist(self.cluster_centers, self.cluster_centers)
        # Only consider upper triangle (i < j)
        mask = torch.triu(torch.ones_like(center_dist), diagonal=1) > 0
        if mask.sum() > 0:
            sep_loss = -center_dist[mask].mean()  # Negative because we want to maximize
        else:
            sep_loss = torch.tensor(0.0, device=z.device)

        # Combined loss with higher weight on center alignment
        total_loss = kl_loss + 0.5 * center_loss + 0.1 * sep_loss

        if return_assignments:
            return q, total_loss
        return None, total_loss

    def get_assignments(self, z: torch.Tensor) -> torch.Tensor:
        """Get soft cluster assignments."""
        dist = torch.cdist(z, self.cluster_centers)
        return torch.softmax(-dist / (self.temperature + 1e-8), dim=-1)


class ScSpade(nn.Module):
    """
    Complete ScSpade model: Support-Masked Diffusion Autoencoder.

    Combines:
    1. SupportMaskNet: Predicts gene activation
    2. LatentDiffusionAE: Denoises latent representations
    3. ClusterLoss: Guides clustering in latent space
    """

    def __init__(
        self,
        num_genes: int,
        n_clusters: int,
        latent_dim: int = 32,
        mask_hidden_dims: list = None,
        diffusion_hidden_dims: list = None,
        diffusion_steps: int = 100,
        dropout: float = 0.1,
    ):
        super().__init__()

        # Support Mask Network
        if mask_hidden_dims is None:
            mask_hidden_dims = [512, 256, 128]
        self.mask_net = SupportMaskNet(
            num_genes=num_genes,
            hidden_dims=mask_hidden_dims,
            dropout=dropout,
        )

        # Latent Diffusion Autoencoder
        if diffusion_hidden_dims is None:
            diffusion_hidden_dims = [512, 256]
        self.diffusion_ae = LatentDiffusionAE(
            num_genes=num_genes,
            latent_dim=latent_dim,
            hidden_dims=diffusion_hidden_dims,
            diffusion_steps=diffusion_steps,
            dropout=dropout,
        )

        # Clustering head
        self.cluster_loss = ClusterLoss(n_clusters, latent_dim)

    def forward(
        self,
        x: torch.Tensor,
        return_recon: bool = True,
        sample_diffusion: bool = False,
    ) -> dict:
        """
        Full forward pass.

        Args:
            x: Input expression (batch, n_genes)
            return_recon: Whether to return reconstruction
            sample_diffusion: Whether to run diffusion sampling

        Returns:
            dict with keys: mask_prob, z, z_denoised, x_recon, losses
        """
        losses = {}

        # 1. Mask prediction
        mask_output = self.mask_net(x)
        mask_prob = torch.nan_to_num(mask_output['gene_activation_prob'], nan=0.5, posinf=1.0, neginf=0.0)
        mask_prob = mask_prob.clamp(1e-6, 1 - 1e-6)

        # Binary mask: 1 = expressed, 0 = zero
        mask = (x > 0).float()

        # Mask loss
        losses['mask'] = nn.functional.binary_cross_entropy(
            mask_prob, mask, reduction='mean'
        )

        # 2. Latent diffusion
        diffusion_result = self.diffusion_ae(
            x,
            mask=mask,
            return_recon=return_recon,
            sample_diffusion=sample_diffusion,
        )

        z = diffusion_result['z']
        z_denoised = diffusion_result['z_denoised']

        # Diffusion and reconstruction losses
        if 'losses' in diffusion_result:
            losses['diffusion'] = diffusion_result['losses'].get(
                'diffusion', torch.tensor(0.0).to(x.device)
            )
            losses['recon'] = diffusion_result['losses'].get(
                'recon', torch.tensor(0.0).to(x.device)
            )

        # 3. Clustering loss - use z for clustering (more stable than z_denoised during training)
        if self.training:
            # Use z_denoised when available for better clustering, otherwise use z
            z_for_cluster = z_denoised if z_denoised is not None else z
            _, cluster_loss = self.cluster_loss(z_for_cluster)
            losses['cluster'] = cluster_loss
        else:
            losses['cluster'] = torch.tensor(0.0, device=x.device)

        return {
            'mask_prob': mask_prob,
            'z': z,
            'z_denoised': z_denoised,
            'x_recon': diffusion_result.get('x_recon'),
            'losses': losses,
        }

    def get_cluster_assignments(self, z: torch.Tensor) -> torch.Tensor:
        """Get soft cluster assignments for given embeddings."""
        return self.cluster_loss.get_assignments(z)

    def get_embedding(self, x: torch.Tensor, use_diffusion: bool = False) -> torch.Tensor:
        """
        Get embedding for clustering.

        Args:
            x: Input expression
            use_diffusion: If True, use diffusion denoised embedding;
                          If False (recommended), use direct encoder output
        """
        self.eval()
        with torch.no_grad():
            if use_diffusion:
                result = self.forward(x, return_recon=False, sample_diffusion=True)
                return result['z_denoised']
            else:
                # Use direct encoder output - faster and often better for clustering
                z = self.diffusion_ae.get_direct_embedding(x)
                return z


def train_epoch(
    model: nn.Module,
    dataloader: DataLoader,
    optimizer: optim.Optimizer,
    device: torch.device,
    epoch: int,
    warmup_epochs: int = 20,
    mask_weight: float = 1.0,
    diffusion_weight: float = 1.0,
    recon_weight: float = 0.1,
    cluster_weight: float = 0.01,
) -> dict:
    """Train for one epoch."""
    model.train()

    loss_meter = AverageMeter()
    mask_loss_meter = AverageMeter()
    diffusion_loss_meter = AverageMeter()
    recon_loss_meter = AverageMeter()
    cluster_loss_meter = AverageMeter()

    is_warmup = epoch < warmup_epochs

    for batch_idx, (x, _) in enumerate(dataloader):
        x = x.to(device)

        optimizer.zero_grad()

        # Forward pass
        result = model(x, return_recon=True, sample_diffusion=False)
        losses = result['losses']

        # Compute total loss
        if is_warmup:
            # Warmup: only train mask model
            total_loss = mask_weight * losses['mask']
        else:
            # Full joint training
            total_loss = (
                mask_weight * losses['mask'] +
                diffusion_weight * losses['diffusion'] +
                recon_weight * losses['recon'] +
                cluster_weight * losses['cluster']
            )

        # Backward pass
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        # Metrics
        loss_meter.update(total_loss.item(), x.size(0))
        mask_loss_meter.update(losses['mask'].item(), x.size(0))
        diffusion_loss_meter.update(losses['diffusion'].item(), x.size(0))
        recon_loss_meter.update(losses['recon'].item(), x.size(0))
        cluster_loss_meter.update(losses['cluster'].item(), x.size(0))

    return {
        'loss': loss_meter.avg,
        'mask_loss': mask_loss_meter.avg,
        'diffusion_loss': diffusion_loss_meter.avg,
        'recon_loss': recon_loss_meter.avg,
        'cluster_loss': cluster_loss_meter.avg,
        'is_warmup': is_warmup,
    }


@torch.no_grad()
def extract_embeddings(
    model: nn.Module,
    dataloader: DataLoader,
    device: torch.device,
    has_labels: bool = True,
) -> tuple:
    """
    Extract embeddings from the model.

    Args:
        model: ScSpade model
        dataloader: DataLoader providing data
        device: Computation device
        has_labels: Whether the dataloader provides labels

    Returns:
        (embeddings, labels) tuple
    """
    model.eval()

    embeddings = []
    labels = []

    for batch in dataloader:
        if has_labels:
            x, y = batch
        else:
            x = batch

        x = x.to(device)
        z = model.get_embedding(x)
        # Handle NaN
        z = torch.nan_to_num(z, nan=0.0, posinf=0.0, neginf=0.0)
        embeddings.append(z.cpu())

        if has_labels:
            labels.append(y)

    embeddings = torch.cat(embeddings, dim=0).numpy()
    if has_labels:
        labels = torch.cat(labels, dim=0).numpy()
    else:
        labels = None

    return embeddings, labels


def initialize_clusters(embeddings: np.ndarray, n_clusters: int) -> np.ndarray:
    """Initialize cluster assignments using K-Means."""
    from sklearn.cluster import KMeans

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    return kmeans.fit_predict(embeddings)


def main():
    args = parse_args()

    # Setup
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    set_seed(args.seed)
    os.makedirs(args.save_dir, exist_ok=True)

    # Load and preprocess data
    print('=' * 60)
    print('Loading and preprocessing data...')
    print('=' * 60)

    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=True,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True,
    )

    X_np = np.array(X).astype(np.float32)
    Y_np = np.array(Y)

    # Determine number of clusters
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y_np))

    print(f'Number of cells: {X_np.shape[0]}')
    print(f'Number of genes: {X_np.shape[1]}')
    print(f'Number of clusters: {n_clusters}')
    print(f'Number of cell types in data: {len(np.unique(Y_np))}')

    # Create dataloader
    dataset = TensorDataset(
        torch.from_numpy(X_np),
        torch.from_numpy(Y_np),
    )
    dataloader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True,
    )

    # Build model
    print('=' * 60)
    print('Building ScSpade model...')
    print('=' * 60)

    mask_hidden_dims = [int(d) for d in args.mask_hidden_dims.split(',')]
    diffusion_hidden_dims = [int(d) for d in args.diffusion_hidden_dims.split(',')]

    model = ScSpade(
        num_genes=X_np.shape[1],
        n_clusters=n_clusters,
        latent_dim=args.latent_dim,
        mask_hidden_dims=mask_hidden_dims,
        diffusion_hidden_dims=diffusion_hidden_dims,
        diffusion_steps=args.diffusion_steps,
        dropout=max(args.mask_dropout, args.diffusion_dropout),
    ).to(device)

    optimizer = optim.Adam(
        model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
    )

    total_params = sum(p.numel() for p in model.parameters())
    trainable_params = sum(p.numel() for p in model.parameters() if p.requires_grad)
    print(f'Total parameters: {total_params:,}')
    print(f'Trainable parameters: {trainable_params:,}')

    # Training loop
    print('=' * 60)
    print('Starting training...')
    print('=' * 60)

    best_loss = float('inf')

    for epoch in range(args.epochs):
        train_metrics = train_epoch(
            model,
            dataloader,
            optimizer,
            device,
            epoch,
            warmup_epochs=args.warmup_epochs,
            mask_weight=args.mask_loss_weight,
            diffusion_weight=args.diffusion_loss_weight,
            recon_weight=args.recon_loss_weight,
            cluster_weight=args.cluster_loss_weight,
        )

        # Evaluate
        if (epoch + 1) % args.eval_interval == 0 or epoch == args.epochs - 1:
            phase = 'WARMUP' if train_metrics['is_warmup'] else 'JOINT'
            print(
                f'Epoch {epoch + 1}/{args.epochs} [{phase}] | '
                f'Loss: {train_metrics["loss"]:.4f} | '
                f'Mask: {train_metrics["mask_loss"]:.4f} | '
                f'Diff: {train_metrics["diffusion_loss"]:.4f} | '
                f'Recon: {train_metrics["recon_loss"]:.4f} | '
                f'Cluster: {train_metrics["cluster_loss"]:.4f}'
            )

            # Extract and save embeddings periodically
            embeddings, labels = extract_embeddings(model, dataloader, device)

            # Quick clustering evaluation
            from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score
            from sklearn.cluster import KMeans

            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
            pred_labels = kmeans.fit_predict(embeddings)

            nmi = normalized_mutual_info_score(labels, pred_labels)
            ari = adjusted_rand_score(labels, pred_labels)

            print(f'  -> NMI: {nmi:.4f}, ARI: {ari:.4f}')

            # Save best model
            if train_metrics['loss'] < best_loss:
                best_loss = train_metrics['loss']
                torch.save(
                    model.state_dict(),
                    os.path.join(args.save_dir, 'best_scspade_model.pt')
                )
                print(f'  -> Saved best model (loss={best_loss:.4f})')
        else:
            print(f'Epoch {epoch + 1}/{args.epochs} | Loss: {train_metrics["loss"]:.4f}')

    # Save final model
    torch.save(
        {
            'model_state_dict': model.state_dict(),
            'args': vars(args),
            'n_clusters': n_clusters,
        },
        os.path.join(args.save_dir, 'scspade_model_final.pt')
    )

    # Final embedding extraction
    print('=' * 60)
    print('Extracting final embeddings...')
    print('=' * 60)

    embeddings, labels = extract_embeddings(model, dataloader, device)

    # Save embeddings
    np.save(os.path.join(args.save_dir, 'embeddings.npy'), embeddings)
    np.save(os.path.join(args.save_dir, 'labels.npy'), labels)

    # Final clustering
    print('=' * 60)
    print('Final clustering evaluation...')
    print('=' * 60)

    from sklearn.metrics import normalized_mutual_info_score, adjusted_rand_score, accuracy_score
    from sklearn.cluster import KMeans

    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    pred_labels = kmeans.fit_predict(embeddings)

    nmi = normalized_mutual_info_score(labels, pred_labels)
    ari = adjusted_rand_score(labels, pred_labels)

    print(f'Final NMI: {nmi:.4f}')
    print(f'Final ARI: {ari:.4f}')

    # Save results
    results = {
        'nmi': nmi,
        'ari': ari,
        'n_clusters': n_clusters,
        'embedding_dim': args.latent_dim,
        'epochs': args.epochs,
        'batch_size': args.batch_size,
    }

    import json
    with open(os.path.join(args.save_dir, 'results.json'), 'w') as f:
        json.dump(results, f, indent=2)

    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
