# -*- coding: utf-8 -*-
"""
Unified scDeepCluster Model Interface for scBench
=================================================

[TensorFlow-gated] Requires TensorFlow/Keras. If you see this error,
please install: pip install tensorflow

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
from time import time

# Add methods/ root and local model dir to path
_LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _LOCAL_DIR)
# Insert methods/ root (two levels up from model dir)
sys.path.insert(0, os.path.dirname(os.path.dirname(_LOCAL_DIR)))

from preprocess import prepare_data_for_model
from utils import save, save_json


def _check_tf():
    """Lazily import TensorFlow; raise clear error if not available."""
    try:
        import tensorflow as tf
        from tensorflow import keras
        from keras.models import Model
        from keras.layers import Dense, Input, GaussianNoise, Layer, Activation
        from keras.optimizers import SGD, Adam
        from keras.callbacks import EarlyStopping
        from tensorflow.keras.optimizers import Adam as keras_Adam
        return tf, keras, keras_Adam
    except ImportError as e:
        raise ImportError(
            "scDeepCluster requires TensorFlow/Keras but it is not installed.\n"
            "Please install: pip install tensorflow\n"
            f"Original error: {e}"
        )


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='scDeepCluster: Deep Clustering for scRNA-seq Data [TensorFlow]',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--data_path', type=str, required=True,
                        help='Path to input h5ad file')
    parser.add_argument('--save_dir', type=str, default='./results',
                        help='Directory to save results')
    parser.add_argument('--n_clusters', type=int, required=True,
                        help='Number of clusters')
    parser.add_argument('--dims', type=int, nargs='+', default=[256, 64, 32],
                        help='Hidden layer dimensions (excluding input dim)')
    parser.add_argument('--noise_sd', type=float, default=2.5,
                        help='Standard deviation of Gaussian noise')
    parser.add_argument('--alpha', type=float, default=1.0,
                        help='Alpha parameter for clustering layer')
    parser.add_argument('--pretrain_epochs', type=int, default=400,
                        help='Number of pretraining epochs')
    parser.add_argument('--epochs', type=int, default=None,
                       help='Unified CLI alias for --pretrain_epochs (smoke test shortcut)')
    parser.add_argument('--maxiter', type=int, default=20000,
                        help='Maximum number of iterations')
    parser.add_argument('--batch_size', type=int, default=256,
                        help='Batch size for training')
    parser.add_argument('--gamma', type=float, default=1.0,
                        help='Coefficient of clustering loss')
    parser.add_argument('--tol', type=float, default=0.001,
                        help='Tolerance for stopping criterion')
    parser.add_argument('--update_interval', type=int, default=0,
                        help='Update interval (0 for auto)')
    parser.add_argument('--ae_weights', type=str, default=None,
                        help='Path to pretrained autoencoder weights')

    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--gpu', type=int, default=1,
                        help='GPU device ID')
    parser.add_argument('--no_cuda', action='store_true',
                        help='Disable CUDA')
    return parser.parse_args()


def main():
    args = parse_args()

    # Forward --epochs -> --pretrain_epochs (unified CLI convenience)
    if args.epochs is not None:
        args.pretrain_epochs = args.epochs

    # Handle no_cuda
    if args.no_cuda:
        os.environ["CUDA_VISIBLE_DEVICES"] = ""

    # Import TF lazily (only when actually running, not just for --help)
    tf, keras, keras_Adam = _check_tf()

    from numpy.random import seed
    seed(args.seed)
    tf.random.set_seed(args.seed)

    from sklearn.cluster import KMeans
    from sklearn.preprocessing import LabelEncoder

    from scdeepcluster_src.scDeepCluster import SCDeepCluster

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

    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    if adata.raw is not None:
        raw_counts = adata.raw.X
        if hasattr(raw_counts, 'toarray'):
            raw_counts = raw_counts.toarray()
        raw_counts = np.array(raw_counts).astype(np.float32)
    else:
        raw_counts = X.copy()

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    if args.update_interval == 0:
        args.update_interval = max(1, int(X.shape[0] / args.batch_size))

    dims = [X.shape[1]] + args.dims

    print('Initializing scDeepCluster model...')
    model = SCDeepCluster(
        dims=dims,
        n_clusters=n_clusters,
        noise_sd=args.noise_sd,
        alpha=args.alpha
    )

    print('Model summary:')
    model.autoencoder.summary()

    t0 = time()

    if args.ae_weights is None:
        print('Pretraining autoencoder...')
        optimizer1 = keras_Adam(amsgrad=True)
        ae_weight_file = os.path.join(args.save_dir, 'ae_weights.h5')
        model.pretrain(
            x=[X, sf],
            y=raw_counts,
            batch_size=args.batch_size,
            epochs=args.pretrain_epochs,
            optimizer=optimizer1,
            ae_file=ae_weight_file
        )

    print('Training clustering model...')
    y_pred = model.fit(
        x_counts=X,
        sf=sf,
        y=Y,
        raw_counts=raw_counts,
        batch_size=args.batch_size,
        tol=args.tol,
        maxiter=args.maxiter,
        update_interval=args.update_interval,
        ae_weights=args.ae_weights,
        save_dir=args.save_dir,
        loss_weights=[args.gamma, 1],
        optimizer='adadelta'
    )

    embedding = model.extract_feature([X, sf])

    save(args.save_dir, Y, y_pred, args.maxiter, embedding, args=vars(args))

    print(f'Training completed in {int(time() - t0)} seconds.')
    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
