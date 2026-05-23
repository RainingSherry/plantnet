# -*- coding: utf-8 -*-
"""
Stage B: Train LatentDiffusionAE

This script trains the LatentDiffusionAE to learn denoised latent representations.
The model encodes cells to a low-dimensional space, applies diffusion-based denoising,
and decodes back to expression space.

Usage:
    python train_embedding.py --data_path /path/to/data.h5ad --save_dir ./results_embedding
"""

import os
import sys
import argparse
import random
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset

# Add project root to path
PROJECT_ROOT = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', '..', '..'))
sys.path.insert(0, PROJECT_ROOT)
from preprocess import prepare_data_for_model
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
        description='Train LatentDiffusionAE for denoised embedding',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Data
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results_embedding',
                        help='Directory to save results')
    parser.add_argument('--mask_model_path', type=str, default=None,
                        help='Path to pretrained mask model')

    # Model
    parser.add_argument('--latent_dim', type=int, default=32,
                        help='Latent space dimension')
    parser.add_argument('--hidden_dims', type=str, default='512,256',
                        help='Comma-separated encoder/decoder hidden dimensions')
    parser.add_argument('--diffusion_steps', type=int, default=100,
                        help='Number of diffusion timesteps')
    parser.add_argument('--diffusion_type', type=str, default='ddpm',
                        choices=['ddpm', 'ddim'],
                        help='Diffusion model type')
    parser.add_argument('--dropout', type=float, default=0.1,
                        help='Dropout rate')

    # Training
    parser.add_argument('--epochs', type=int, default=100,
                        help='Number of training epochs')
    parser.add_argument('--batch_size', type=int, default=128,
                        help='Batch size')
    parser.add_argument('--lr', type=float, default=1e-3,
                        help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-4,
                        help='Weight decay')
    parser.add_argument('--diffusion_loss_weight', type=float, default=1.0,
                        help='Weight for diffusion loss')
    parser.add_argument('--recon_loss_weight', type=float, default=0.1,
                        help='Weight for reconstruction loss')

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


def train_epoch(
    model: nn.Module,
    dataloader: DataLoader,
    optimizer: optim.Optimizer,
    device: torch.device,
    epoch: int,
    diffusion_weight: float = 1.0,
    recon_weight: float = 0.1,
) -> dict:
    """Train for one epoch."""
    model.train()

    loss_meter = AverageMeter()
    diffusion_loss_meter = AverageMeter()
    recon_loss_meter = AverageMeter()

    for batch_idx, (x, _) in enumerate(dataloader):
        x = x.to(device)

        # Compute mask: 1 = expressed (X > 0), 0 = zero
        mask = (x > 0).float()

        # Forward pass
        optimizer.zero_grad()
        result = model(x, mask=mask, return_recon=True, sample_diffusion=False)

        # Compute total loss
        losses = result['losses']
        total_loss = (
            diffusion_weight * losses.get('diffusion', torch.tensor(0.0).to(device)) +
            recon_weight * losses.get('recon', torch.tensor(0.0).to(device))
        )

        # Backward pass
        total_loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
        optimizer.step()

        # Metrics
        loss_meter.update(total_loss.item(), x.size(0))
        diffusion_loss_meter.update(
            losses.get('diffusion', torch.tensor(0.0)).item(), x.size(0)
        )
        recon_loss_meter.update(
            losses.get('recon', torch.tensor(0.0)).item(), x.size(0)
        )

    return {
        'loss': loss_meter.avg,
        'diffusion_loss': diffusion_loss_meter.avg,
        'recon_loss': recon_loss_meter.avg,
    }


@torch.no_grad()
def extract_embeddings(
    model: nn.Module,
    dataloader: DataLoader,
    device: torch.device,
) -> np.ndarray:
    """Extract denoised embeddings."""
    model.eval()

    embeddings = []
    labels = []

    for x, y in dataloader:
        x = x.to(device)

        # Get denoised embedding
        result = model(x, sample_diffusion=True)
        z_denoised = result['z_denoised']

        embeddings.append(z_denoised.cpu())
        labels.append(y)

    embeddings = torch.cat(embeddings, dim=0).numpy()
    labels = torch.cat(labels, dim=0).numpy()

    return embeddings, labels


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

    print(f'Number of cells: {X_np.shape[0]}')
    print(f'Number of genes: {X_np.shape[1]}')
    print(f'Number of cell types: {len(np.unique(Y_np))}')

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
    print('Building model...')
    print('=' * 60)

    hidden_dims = [int(d) for d in args.hidden_dims.split(',')]
    model = LatentDiffusionAE(
        num_genes=X_np.shape[1],
        latent_dim=args.latent_dim,
        hidden_dims=hidden_dims,
        diffusion_steps=args.diffusion_steps,
        diffusion_type=args.diffusion_type,
        dropout=args.dropout,
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
            diffusion_weight=args.diffusion_loss_weight,
            recon_weight=args.recon_loss_weight,
        )

        # Evaluate
        if (epoch + 1) % args.eval_interval == 0 or epoch == args.epochs - 1:
            print(
                f'Epoch {epoch + 1}/{args.epochs} | '
                f'Loss: {train_metrics["loss"]:.4f} | '
                f'Diff: {train_metrics["diffusion_loss"]:.4f} | '
                f'Recon: {train_metrics["recon_loss"]:.4f}'
            )

            # Save best model
            if train_metrics['loss'] < best_loss:
                best_loss = train_metrics['loss']
                torch.save(
                    model.state_dict(),
                    os.path.join(args.save_dir, 'best_diffusion_model.pt')
                )
                print(f'  -> Saved best model (loss={best_loss:.4f})')
        else:
            print(f'Epoch {epoch + 1}/{args.epochs} | Loss: {train_metrics["loss"]:.4f}')

    # Save final model
    torch.save(
        {
            'model_state_dict': model.state_dict(),
            'args': args,
        },
        os.path.join(args.save_dir, 'diffusion_model_final.pt')
    )

    # Extract embeddings
    print('=' * 60)
    print('Extracting denoised embeddings...')
    print('=' * 60)

    embeddings, labels = extract_embeddings(model, dataloader, device)

    np.save(os.path.join(args.save_dir, 'diffusion_embeddings.npy'), embeddings)
    np.save(os.path.join(args.save_dir, 'labels.npy'), labels)

    print(f'Embeddings shape: {embeddings.shape}')
    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
