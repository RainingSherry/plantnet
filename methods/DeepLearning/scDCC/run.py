# -*- coding: utf-8 -*-
"""
Unified scDCC Model Interface for scBench
=========================================

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
import torch
import random

# Add methods/ root and local model dir to path
# Local dir must come before methods/ root so local utils.py takes precedence
_LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _LOCAL_DIR)
sys.path.insert(0, os.path.dirname(os.path.dirname(_LOCAL_DIR)))

from preprocess import prepare_data_for_model
from utils import save, save_json

from scDCC import scDCC


def set_seed(seed):
    """Set random seed for reproducibility."""
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='scDCC: Deep Constrained Clustering for scRNA-seq Data',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results',
                        help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, required=True,
                        help='Number of clusters')
    parser.add_argument('--n_z', type=int, default=32,
                        help='Dimension of latent space')
    parser.add_argument('--encode_dims', type=int, nargs='+', default=[256, 64],
                        help='Encoder hidden layer dimensions')
    parser.add_argument('--decode_dims', type=int, nargs='+', default=[64, 256],
                        help='Decoder hidden layer dimensions')
    parser.add_argument('--sigma', type=float, default=2.5,
                        help='Standard deviation of Gaussian noise')
    parser.add_argument('--alpha', type=float, default=1.0,
                        help='Alpha parameter for clustering layer')
    parser.add_argument('--gamma', type=float, default=1.0,
                        help='Gamma parameter for clustering loss')
    parser.add_argument('--pretrain_epochs', type=int, default=400,
                        help='Number of pretraining epochs')
    parser.add_argument('--epochs', type=int, default=200,
                        help='Number of training epochs')
    parser.add_argument('--batch_size', type=int, default=256,
                        help='Batch size for training')
    parser.add_argument('--lr', type=float, default=1.0,
                        help='Learning rate for clustering')
    parser.add_argument('--pretrain_lr', type=float, default=0.001,
                        help='Learning rate for pretraining')
    parser.add_argument('--tol', type=float, default=0.001,
                        help='Tolerance for stopping criterion')
    parser.add_argument('--update_interval', type=int, default=1,
                        help='Update interval for target distribution')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--gpu', type=int, default=0,
                        help='GPU device ID')
    parser.add_argument('--no_cuda', action='store_true',
                        help='Disable CUDA')
    parser.add_argument('--ae_weights', type=str, default=None,
                        help='Path to pretrained autoencoder weights')

    return parser.parse_args()


def main():
    args = parse_args()

    if torch.cuda.is_available():
        torch.cuda.set_device(args.gpu)
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    set_seed(args.seed)

    os.makedirs(args.save_dir, exist_ok=True)
    save_json(vars(args), os.path.join(args.save_dir, 'run_config.json'))

    print('Loading data...')
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=True,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True
    )

    X = np.array(X).astype(np.float32)
    Y = np.array(Y)
    sf = np.array(sf).astype(np.float32).reshape(-1, 1)

    from sklearn.preprocessing import LabelEncoder
    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    # Use norm_log layer (normalized+log1p HVG data) for ZINB loss; shape must match X
    if 'norm_log' in adata.layers:
        raw_counts = adata.layers['norm_log']
        if hasattr(raw_counts, 'toarray'):
            raw_counts = raw_counts.toarray()
        raw_counts = np.array(raw_counts).astype(np.float32)
    elif adata.raw is not None:
        # Fallback: use raw data but HVG subset only (match X columns)
        raw_counts = adata.raw.X
        if hasattr(raw_counts, 'toarray'):
            raw_counts = raw_counts.toarray()
        raw_counts = np.array(raw_counts).astype(np.float32)
    else:
        raw_counts = X.copy()

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    print('Initializing scDCC model...')
    model = scDCC(
        input_dim=X.shape[1],
        z_dim=args.n_z,
        n_clusters=n_clusters,
        encodeLayer=args.encode_dims,
        decodeLayer=args.decode_dims,
        sigma=args.sigma,
        alpha=args.alpha,
        gamma=args.gamma
    )

    ae_weights_path = os.path.join(args.save_dir, 'ae_weights.pth.tar')
    if args.ae_weights is None:
        print('Pretraining autoencoder...')
        model.pretrain_autoencoder(
            x=X,
            X_raw=raw_counts,
            size_factor=sf,
            batch_size=args.batch_size,
            lr=args.pretrain_lr,
            epochs=args.pretrain_epochs,
            ae_save=True,
            ae_weights=ae_weights_path,
            device=device
        )
    else:
        print(f'Loading pretrained weights from {args.ae_weights}')
        checkpoint = torch.load(args.ae_weights, map_location='cpu')
        model.load_state_dict(checkpoint['ae_state_dict'])

    print('Training clustering model...')
    y_pred, final_acc, final_nmi, final_ari, final_epoch = model.fit(
        X=X,
        X_raw=raw_counts,
        sf=sf,
        y=Y,
        lr=args.lr,
        batch_size=args.batch_size,
        num_epochs=args.epochs,
        update_interval=args.update_interval,
        tol=args.tol,
        save_dir=args.save_dir,
        device=device
    )

    X_tensor = torch.tensor(X).to(device)
    embedding = model.encodeBatch(X_tensor, device=device).cpu().numpy()

    save(args.save_dir, Y, y_pred, args.epochs, embedding, args=vars(args))

    print(f'Training completed.')
    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
