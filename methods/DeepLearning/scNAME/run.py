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


def _to_dense_float32(matrix):
    """Convert sparse/dense AnnData matrices to dense float32 arrays."""
    if hasattr(matrix, 'toarray'):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def _validate_scname_count_matrix(matrix, expected_shape, source_name):
    if matrix.shape != expected_shape:
        raise ValueError(
            f"scNAME count_X from {source_name} has shape {matrix.shape}, "
            f"expected {expected_shape}."
        )
    if not np.isfinite(matrix).all():
        raise ValueError(f"scNAME count_X from {source_name} contains NaN or inf values.")
    if matrix.size and float(np.min(matrix)) < 0:
        raise ValueError(
            f"scNAME count_X from {source_name} contains negative values; "
            "scaled adata.X/adata.to_df() cannot be used as count_X."
        )
    return matrix


def _get_scname_count_matrix(adata):
    """Return nonnegative count/count-like input for scNAME NB/ZINB loss."""
    expected_shape = adata.X.shape

    if adata.raw is not None:
        try:
            raw_counts = _to_dense_float32(adata.raw[:, adata.var_names].X)
            return _validate_scname_count_matrix(
                raw_counts, expected_shape, "adata.raw[:, adata.var_names].X"
            )
        except (KeyError, ValueError, IndexError) as exc:
            raw_error = exc
        else:
            raw_error = None
    else:
        raw_error = None

    if 'norm_log' in adata.layers:
        norm_log = _to_dense_float32(adata.layers['norm_log'])
        return _validate_scname_count_matrix(norm_log, expected_shape, "adata.layers['norm_log']")

    detail = f" Raw count lookup failed: {raw_error}" if raw_error is not None else ""
    raise ValueError(
        "scNAME requires nonnegative count_X from adata.raw[:, adata.var_names].X "
        "or adata.layers['norm_log']; scaled adata.X/adata.to_df() cannot be used "
        f"as count_X.{detail}"
    )


def _size_factors_from_counts(raw_counts):
    total_counts = raw_counts.sum(axis=1).astype(np.float32)
    if not np.isfinite(total_counts).all():
        raise ValueError("Cannot compute scNAME size factors: raw_counts row sums are not finite.")
    median_count = np.median(total_counts)
    if not np.isfinite(median_count) or median_count <= 0:
        raise ValueError("Cannot compute scNAME size factors: raw_counts median row sum is not positive.")
    return (total_counts / median_count).astype(np.float32)


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
    # NOTE: prepare_data_for_model does ALL preprocessing (normalize, log1p, scale, HVG).
    # Do NOT call scNAME_normalize here — it would redo normalization on already-normalized
    # data, causing NaN/inf corruption (filter_genes on float sparse matrix is destructive).
    print('Loading data...')
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=True,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True
    )
    print('Data ready (normalized, log-transformed, scaled, HVG-selected).')

    # Convert to numpy arrays
    X = np.array(X).astype(np.float32)
    Y = np.array(Y)

    raw_counts = _get_scname_count_matrix(adata)

    # Compute size factors — always reshape to (n_cells, 1) for TF placeholder
    if sf is None:
        sf = _size_factors_from_counts(raw_counts)
    else:
        sf = np.array(sf).astype(np.float32)
        if not np.isfinite(sf).all():
            sf = _size_factors_from_counts(raw_counts)
    sf = sf.reshape(-1, 1)

    # Encode labels to integers if needed
    from sklearn.preprocessing import LabelEncoder
    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

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
    save(args.save_dir, Y_restored, y_pred_restored, args.epochs, embedding_restored, args=vars(args))

    print(f'Training completed.')
    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
