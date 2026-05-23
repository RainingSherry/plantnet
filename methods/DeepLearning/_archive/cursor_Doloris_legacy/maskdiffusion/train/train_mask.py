# -*- coding: utf-8 -*-
"""
Stage A: Train SupportMaskNet

This script trains the SupportMaskNet to predict gene activation probabilities.
The mask model learns the observed support structure (which genes are expressed)
rather than distinguishing true zeros from dropouts.

Usage:
    python train_mask.py --data_path /path/to/data.h5ad --save_dir ./results_mask
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
from models.support_mask import SupportMaskNet


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
        description='Train SupportMaskNet for gene activation prediction',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Data
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results_mask',
                        help='Directory to save results')

    # Model
    parser.add_argument('--hidden_dims', type=str, default='512,256,128',
                        help='Comma-separated hidden layer dimensions')
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
) -> dict:
    """Train for one epoch."""
    model.train()

    loss_meter = AverageMeter()
    acc_meter = AverageMeter()
    bce_meter = AverageMeter()

    for batch_idx, (x, _) in enumerate(dataloader):
        x = x.to(device)

        # Compute mask: 1 = expressed (X > 0), 0 = zero
        mask = (x > 0).float()

        # Forward pass
        optimizer.zero_grad()
        loss, gene_activation_prob = model.get_gene_activation_loss(
            x, mask=mask, pos_weight=1.0
        )

        # Backward pass
        loss.backward()
        optimizer.step()

        # Compute accuracy
        pred = (gene_activation_prob > 0.5).float()
        acc = (pred == mask).float().mean()

        # Metrics
        loss_meter.update(loss.item(), x.size(0))
        acc_meter.update(acc.item(), x.size(0))

    return {
        'loss': loss_meter.avg,
        'acc': acc_meter.avg,
    }


@torch.no_grad()
def evaluate(
    model: nn.Module,
    dataloader: DataLoader,
    device: torch.device,
) -> dict:
    """Evaluate model."""
    model.eval()

    total_loss = 0
    total_correct = 0
    total_count = 0

    all_preds = []
    all_masks = []

    for x, _ in dataloader:
        x = x.to(device)
        mask = (x > 0).float()

        loss, gene_activation_prob = model.get_gene_activation_loss(
            x, mask=mask
        )

        pred = (gene_activation_prob > 0.5).float()
        correct = (pred == mask).float().sum().item()

        total_loss += loss.item() * x.size(0)
        total_correct += correct
        total_count += mask.numel()

        all_preds.append(gene_activation_prob.cpu())
        all_masks.append(mask.cpu())

    all_preds = torch.cat(all_preds, dim=0)
    all_masks = torch.cat(all_masks, dim=0)

    # Compute per-class metrics
    expressed_correct = ((all_preds > 0.5) == (all_masks > 0.5))[all_masks > 0.5].float().mean().item()
    zero_correct = ((all_preds > 0.5) == (all_masks > 0.5))[all_masks == 0].float().mean().item()

    return {
        'loss': total_loss / len(dataloader.dataset),
        'accuracy': total_correct / total_count,
        'expressed_accuracy': expressed_correct,
        'zero_accuracy': zero_correct,
    }


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
    model = SupportMaskNet(
        num_genes=X_np.shape[1],
        hidden_dims=hidden_dims,
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
            model, dataloader, optimizer, device, epoch
        )

        # Evaluate
        if (epoch + 1) % args.eval_interval == 0 or epoch == args.epochs - 1:
            eval_metrics = evaluate(model, dataloader, device)

            print(
                f'Epoch {epoch + 1}/{args.epochs} | '
                f'Loss: {train_metrics["loss"]:.4f} | '
                f'Acc: {eval_metrics["accuracy"]:.4f} | '
                f'Expressed Acc: {eval_metrics["expressed_accuracy"]:.4f} | '
                f'Zero Acc: {eval_metrics["zero_accuracy"]:.4f}'
            )

            # Save best model
            if eval_metrics['loss'] < best_loss:
                best_loss = eval_metrics['loss']
                torch.save(
                    model.state_dict(),
                    os.path.join(args.save_dir, 'best_mask_model.pt')
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
        os.path.join(args.save_dir, 'mask_model_final.pt')
    )

    # Extract embeddings
    print('=' * 60)
    print('Extracting embeddings...')
    print('=' * 60)

    model.eval()
    embeddings = []
    with torch.no_grad():
        for x, _ in dataloader:
            x = x.to(device)
            output = model(x)
            embeddings.append(output['cell_embedding'].cpu())

    embeddings = torch.cat(embeddings, dim=0).numpy()
    np.save(os.path.join(args.save_dir, 'mask_embeddings.npy'), embeddings)

    print(f'Embeddings shape: {embeddings.shape}')
    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
