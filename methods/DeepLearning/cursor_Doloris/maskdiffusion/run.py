# -*- coding: utf-8 -*-
"""
Main run script for ScSpade (Support-Masked Diffusion Autoencoder for scRNA-seq Clustering).

This script provides an end-to-end pipeline:
1. Data loading and preprocessing
2. Model training (joint SupportMaskNet + LatentDiffusionAE)
3. Embedding extraction
4. Clustering evaluation
5. Results saving

Usage:
    python run.py --data_path /path/to/data.h5ad --save_dir ./results_scspade
"""

import os
import sys
import argparse
import random
import json
import numpy as np
import torch
from sklearn.cluster import KMeans
import scanpy as sc

# Add project root to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from data import load_and_preprocess, create_dataloader
from models.support_mask import SupportMaskNet
from models.latent_diffusion import LatentDiffusionAE
from eval.cluster_eval import evaluate_clustering, run_all_evaluations
from eval.sparsity_eval import compute_sparsity_stats, compare_sparsity_patterns
from train.train_joint import ScSpade, set_seed, train_epoch, extract_embeddings


def parse_args():
    parser = argparse.ArgumentParser(
        description='ScSpade: Support-Masked Diffusion Autoencoder for scRNA-seq Clustering',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Data
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results_scspade',
                        help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, default=0,
                        help='Number of clusters (0 = auto-detect from data)')

    # Preprocessing
    parser.add_argument('--n_top_genes', type=int, default=1000,
                        help='Number of highly variable genes')

    # Mask Model
    parser.add_argument('--mask_hidden_dims', type=str, default='512,256,128',
                        help='Comma-separated hidden layer dimensions for mask model')
    parser.add_argument('--mask_dropout', type=float, default=0.1,
                        help='Dropout rate for mask model')

    # Diffusion Model
    parser.add_argument('--latent_dim', type=int, default=64,
                        help='Latent space dimension (64 recommended for better clustering)')
    parser.add_argument('--diffusion_hidden_dims', type=str, default='1024,512',
                        help='Comma-separated encoder/decoder hidden dimensions')
    parser.add_argument('--diffusion_steps', type=int, default=100,
                        help='Number of diffusion timesteps')
    parser.add_argument('--dropout', type=float, default=0.2,
                        help='Dropout rate')

    # Training
    parser.add_argument('--epochs', type=int, default=100,
                        help='Number of training epochs')
    parser.add_argument('--batch_size', type=int, default=256,
                        help='Batch size')
    parser.add_argument('--lr', type=float, default=1e-3,
                        help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-3,
                        help='Weight decay')
    parser.add_argument('--warmup_epochs', type=int, default=10,
                        help='Number of warmup epochs (mask-only training)')

    # Loss weights - OPTIMIZED for better clustering
    parser.add_argument('--mask_loss_weight', type=float, default=0.2,
                        help='Weight for mask loss')
    parser.add_argument('--diffusion_loss_weight', type=float, default=0.2,
                        help='Weight for diffusion loss')
    parser.add_argument('--recon_loss_weight', type=float, default=0.8,
                        help='Weight for reconstruction loss')
    parser.add_argument('--cluster_loss_weight', type=float, default=1.0,
                        help='Weight for clustering loss')

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


def main():
    args = parse_args()

    # Setup
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    set_seed(args.seed)
    os.makedirs(args.save_dir, exist_ok=True)

    # Save arguments
    with open(os.path.join(args.save_dir, 'args.json'), 'w') as f:
        json.dump(vars(args), f, indent=2)

    # =========================================================================
    # Step 1: Load and preprocess data
    # =========================================================================
    print('=' * 60)
    print('Step 1: Loading and preprocessing data...')
    print('=' * 60)

    X, Y, adata = load_and_preprocess(
        args.data_path,
        n_top_genes=args.n_top_genes,
        log_transform=True,
        normalize=True,
        scale=False,  # Don't scale to preserve sparsity
        preserve_sparsity=True,  # Use [0,1] normalization instead
    )

    # Determine number of clusters
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))

    print(f'Number of cells: {X.shape[0]}')
    print(f'Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')
    print(f'Ground truth cell types: {len(np.unique(Y))}')

    # Sparsity analysis
    sparsity_stats = compute_sparsity_stats(X)
    print(f'Data sparsity: {sparsity_stats["total_sparsity"]:.4f}')

    # =========================================================================
    # Step 2: Create dataloader
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 2: Creating dataloader...')
    print('=' * 60)

    # Training dataloader (with drop_last=True for stable training)
    dataloader = create_dataloader(
        X, Y,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True,
    )

    # Evaluation dataloader (drop_last=False to process all samples)
    eval_dataloader = create_dataloader(
        X, Y,
        batch_size=args.batch_size,
        shuffle=False,
        drop_last=False,
    )

    print(f'Number of training batches: {len(dataloader)}')
    print(f'Number of evaluation batches: {len(eval_dataloader)}')

    # =========================================================================
    # Step 3: Build model
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 3: Building ScSpade model...')
    print('=' * 60)

    mask_hidden_dims = [int(d) for d in args.mask_hidden_dims.split(',')]
    diffusion_hidden_dims = [int(d) for d in args.diffusion_hidden_dims.split(',')]

    model = ScSpade(
        num_genes=X.shape[1],
        n_clusters=n_clusters,
        latent_dim=args.latent_dim,
        mask_hidden_dims=mask_hidden_dims,
        diffusion_hidden_dims=diffusion_hidden_dims,
        diffusion_steps=args.diffusion_steps,
        dropout=args.dropout,
    ).to(device)

    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=args.lr,
        weight_decay=args.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=args.epochs)

    total_params = sum(p.numel() for p in model.parameters())
    print(f'Total parameters: {total_params:,}')

    # =========================================================================
    # Step 4: Training loop with early stopping
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 4: Training...')
    print('=' * 60)

    best_nmi = 0
    best_metrics = None
    patience = 20
    patience_counter = 0

    from train.train_joint import train_epoch, extract_embeddings

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

        scheduler.step()

        # Evaluate periodically
        if (epoch + 1) % args.eval_interval == 0 or epoch == args.epochs - 1:
            embeddings, labels = extract_embeddings(model, eval_dataloader, device, has_labels=True)

            # KMeans clustering with more inits
            kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=20)
            pred_labels = kmeans.fit_predict(embeddings)

            # Evaluate
            metrics = evaluate_clustering(labels, pred_labels, embeddings, verbose=False)

            phase = 'WARMUP' if train_metrics['is_warmup'] else 'JOINT'
            print(
                f'Epoch {epoch + 1}/{args.epochs} [{phase}] | '
                f'Loss: {train_metrics["loss"]:.4f} | '
                f'NMI: {metrics["nmi"]:.4f} | '
                f'ARI: {metrics["ari"]:.4f}'
            )

            # Save best model
            if metrics['nmi'] > best_nmi:
                best_nmi = metrics['nmi']
                best_metrics = metrics.copy()
                best_metrics['epoch'] = epoch + 1
                torch.save(
                    model.state_dict(),
                    os.path.join(args.save_dir, 'best_model.pt')
                )
                np.save(
                    os.path.join(args.save_dir, 'best_embeddings.npy'),
                    embeddings
                )
                np.save(
                    os.path.join(args.save_dir, 'best_labels.npy'),
                    labels
                )
                patience_counter = 0
            else:
                patience_counter += 1

            # Early stopping
            if patience_counter >= patience:
                print(f'Early stopping at epoch {epoch + 1}')
                break

    print(f'\nBest NMI: {best_nmi:.4f} (epoch {best_metrics["epoch"]})')

    # =========================================================================
    # Step 5: Final evaluation
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 5: Final evaluation...')
    print('=' * 60)

    # Load best model
    model.load_state_dict(
        torch.load(os.path.join(args.save_dir, 'best_model.pt'), weights_only=False)
    )

    # Extract final embeddings
    embeddings, labels = extract_embeddings(model, eval_dataloader, device, has_labels=True)

    # KMeans clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    pred_labels = kmeans.fit_predict(embeddings)

    # Comprehensive evaluation
    results = run_all_evaluations(
        labels,
        pred_labels,
        X=embeddings,
        min_rare_size=10,
        verbose=True,
    )

    # =========================================================================
    # Step 6: Save results
    # =========================================================================
    print('\n' + '=' * 60)
    print('Step 6: Saving results...')
    print('=' * 60)

    # Save embeddings
    np.save(os.path.join(args.save_dir, 'embeddings.npy'), embeddings)
    np.save(os.path.join(args.save_dir, 'labels.npy'), labels)
    np.save(os.path.join(args.save_dir, 'pred_labels.npy'), pred_labels)

    # Save results summary
    results_summary = {
        'clustering_metrics': results['clustering'],
        'rare_cell_metrics': results['rare_cells'],
        'sparsity_stats': sparsity_stats,
        'best_epoch': best_metrics.get('epoch', 'N/A'),
        'best_nmi': best_nmi,
        'n_clusters': n_clusters,
        'latent_dim': args.latent_dim,
    }

    with open(os.path.join(args.save_dir, 'results.json'), 'w') as f:
        json.dump(results_summary, f, indent=2, default=str)

    print(f'Results saved to: {args.save_dir}')
    print('\nDone!')


if __name__ == '__main__':
    main()
