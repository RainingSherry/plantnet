# -*- coding: utf-8 -*-
"""
Unified scNAME Model Interface for scBench
==========================================

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd

# Add parent directory to path for imports
_repo_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, _repo_root)

from preprocess import prepare_data_for_model
from utils import save
from scNAME_preprocess import normalize as scNAME_normalize

_TF_READY = False

def _ensure_tf():
    """Lazy import TensorFlow only when actually needed (not for --help)."""
    global _TF_READY, tf, tf, autoencoder
    if _TF_READY:
        return
    import tensorflow.compat.v1 as tf
    tf.disable_v2_behavior()
    from scNAME_network import autoencoder
    _TF_READY = True

def set_seed(seed):
    """Set random seed for reproducibility."""
    np.random.seed(seed)
    try:
        tf.set_random_seed(seed)
    except NameError:
        pass  # TF not imported yet


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='scNAME: Self-supervised Contrastive Learning for scRNA-seq Clustering',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Data arguments
    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results',
                        help='Directory to save results')

    # Model arguments
    parser.add_argument('--n_clusters', type=int, required=True,
                        help='Number of clusters')
    parser.add_argument('--dims', type=int, nargs='+', default=[256, 64, 32],
                        help='Hidden layer dimensions (excluding input dim)')
    parser.add_argument('--k', type=int, default=10,
                        help='Top k most similar features to be neighborhoods')
    parser.add_argument('--temperature', type=float, default=0.7,
                        help='Temperature parameter in contrastive loss')
    parser.add_argument('--p_m', type=float, default=0.3,
                        help='Corruption probability for self-supervised learning')
    parser.add_argument('--noise_sd', type=float, default=1.5,
                        help='Standard deviation of Gaussian noise')

    # Loss weights
    parser.add_argument('--alpha', type=float, default=1.0,
                        help='Weight for mask loss')
    parser.add_argument('--beta', type=float, default=0.1,
                        help='Weight for neighbor loss')
    parser.add_argument('--gamma', type=float, default=0.1,
                        help='Weight for kmeans loss')

    # Training arguments
    parser.add_argument('--pretrain_epochs', type=int, default=500,
                        help='Number of pretraining epochs')
    parser.add_argument('--epochs', type=int, default=1000,
                        help='Number of fine-tuning epochs')
    parser.add_argument('--batch_size', type=int, default=128,
                        help='Batch size for training')
    parser.add_argument('--lr', type=float, default=0.0001,
                        help='Learning rate')
    parser.add_argument('--update_epoch', type=int, default=50,
                        help='Update interval for checking convergence')
    parser.add_argument('--tol', type=float, default=0.001,
                        help='Tolerance for stopping criterion')

    # Other arguments
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--gpu', type=int, default=1,
                        help='GPU device ID')

    parser.add_argument('--no_cuda', action='store_true',
                        help='Disable CUDA')
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    # Handle no_cuda
    if args.no_cuda:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    # Lazy import TF only when actually running
    _ensure_tf()
    set_seed(args.seed)

    # Create save directory
    os.makedirs(args.save_dir, exist_ok=True)

    # Load and preprocess data using standard interface
    # NOTE: scNAME_normalize() requires RAW COUNTS for its assert check.
    # Do NOT normalize/log-transform here — let scNAME_normalize handle it.
    print('Loading data...')
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=False,
        filter_min_counts=True,
        logtrans_input=False,
        normalize_input=False
    )

    # Apply scNAME-specific normalization (ZINB-aware, without additional HVG filtering).
    # prepare_data_for_model() already did HVG filtering → adata.X has 2000 genes.
    # Do NOT call highly_variable_genes here again (would fail on already-HVG'd data).
    print('Applying scNAME normalization...')
    adata = scNAME_normalize(adata, copy=True, highly_genes=None,
                             size_factors=False, normalize_input=False, logtrans_input=False)

    # Convert to numpy arrays
    X = np.array(X).astype(np.float32)
    Y = np.array(Y)

    # Compute size factors from raw counts if not provided by prepare_data_for_model
    if sf is None:
        total_counts = np.array(adata.X.sum(axis=1)).flatten()
        median_count = np.median(total_counts)
        sf = (total_counts / median_count).astype(np.float32).reshape(-1, 1)

    # Encode labels to integers if needed
    from sklearn.preprocessing import LabelEncoder
    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    # Get raw counts for ZINB loss — must match adata.X shape (same genes as HVG subset)
    # adata.raw contains ALL genes; adata.X only has HVG genes after prepare_data_for_model
    # Use adata.to_df() which matches adata.X
    raw_counts = np.array(adata.to_df()).astype(np.float32)

    # Shuffle data
    n = X.shape[0]
    shuffle_ix = np.random.permutation(np.arange(n))
    X = X[shuffle_ix]
    Y = Y[shuffle_ix]
    raw_counts = raw_counts[shuffle_ix]
    sf = sf[shuffle_ix]

    # Get number of clusters from data if not specified
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    # Build model dimensions
    dims = [X.shape[1]] + args.dims

    # Reset TensorFlow graph
    tf.reset_default_graph()

    # Initialize model
    print('Initializing scNAME model...')
    model = autoencoder(
        dataname='scNAME',
        n=n,
        batch_size=args.batch_size,
        k=args.k,
        temperature=args.temperature,
        dims=dims,
        cluster_num=n_clusters,
        alpha=args.alpha,
        beta=args.beta,
        gamma=args.gamma,
        learning_rate=args.lr,
        noise_sd=args.noise_sd
    )

    # Pretrain
    print('Pretraining...')
    gpu_option = "" if args.no_cuda else str(args.gpu)
    model.pretrain(
        X=X,
        count_X=raw_counts,
        p_m=args.p_m,
        size_factor=sf,
        batch_size=args.batch_size,
        pretrain_epoch=args.pretrain_epochs,
        gpu_option=gpu_option
    )

    # Fine-tune
    print('Fine-tuning...')
    y_pred = model.funetrain(
        dataname='scNAME',
        X=X,
        Y=Y,
        count_X=raw_counts,
        p_m=args.p_m,
        size_factor=sf,
        batch_size=args.batch_size,
        funetrain_epoch=args.epochs,
        update_epoch=args.update_epoch,
        error=args.tol
    )

    # Get embeddings (use bank_current as embedding)
    embedding = model.bank_current

    # Restore original order for saving
    restore_ix = np.argsort(shuffle_ix)
    y_pred_restored = y_pred[restore_ix]
    Y_restored = Y[restore_ix]
    embedding_restored = embedding[restore_ix]

    # Save results using standard interface
    save(args.save_dir, Y_restored, y_pred_restored, args.epochs, embedding_restored)

    print(f'Training completed.')
    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
