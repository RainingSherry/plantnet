# -*- coding: utf-8 -*-
"""
Unified scMAE Model Interface for scBench
==========================================

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --epochs 100 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
import torch
import random
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans
from torch.utils.data import DataLoader, Dataset

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from preprocess import prepare_data_for_model
from utils import save

# Import model components
from model import AutoEncoder


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


def apply_noise(X, p):
    """Apply noise mask to input data."""
    p = torch.tensor(p)
    should_swap = torch.bernoulli(p.to(X.device) * torch.ones((X.shape)).to(X.device))
    corrupted_X = torch.where(should_swap == 1, X[torch.randperm(X.shape[0])], X)
    masked = (corrupted_X != X).float()
    return corrupted_X, masked


class scRNADataset(Dataset):
    """Simple dataset class for scRNA-seq data."""
    def __init__(self, data, labels):
        self.data = torch.FloatTensor(data)
        self.labels = torch.LongTensor(labels)

    def __len__(self):
        return len(self.data)

    def __getitem__(self, idx):
        return self.data[idx], self.labels[idx]


def res_search_fixed_clus(adata, fixed_clus_count, increment=0.02):
    """Search for resolution parameter to get fixed number of clusters."""
    dis = []
    resolutions = sorted(list(np.arange(0.01, 2.5, increment)), reverse=True)
    for res in resolutions:
        sc.tl.leiden(adata, random_state=0, resolution=res)
        count_unique_leiden = len(pd.DataFrame(adata.obs['leiden']).leiden.unique())
        dis.append(abs(count_unique_leiden - fixed_clus_count))
        if count_unique_leiden == fixed_clus_count:
            break
    reso = resolutions[np.argmin(dis)]
    return reso


def inference(net, data_loader, device):
    """Extract features from the model."""
    net.eval()
    feature_vector = []
    labels_vector = []
    with torch.no_grad():
        for x, y in data_loader:
            x = x.to(device)
            feature_vector.extend(net.feature(x).detach().cpu().numpy())
            labels_vector.extend(y.numpy())
    return np.array(feature_vector), np.array(labels_vector)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='scMAE: Masked Autoencoder for scRNA-seq Clustering',
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
    parser.add_argument('--hidden_size', type=int, default=128,
                        help='Hidden layer size')
    parser.add_argument('--mask_prob', type=float, default=0.4,
                        help='Mask probability')
    parser.add_argument('--masked_data_weight', type=float, default=0.75,
                        help='Weight for masked data reconstruction')
    parser.add_argument('--mask_loss_weight', type=float, default=0.7,
                        help='Weight for mask prediction loss')

    # Training arguments
    parser.add_argument('--epochs', type=int, default=100,
                        help='Number of training epochs')
    parser.add_argument('--batch_size', type=int, default=256,
                        help='Batch size for training')
    parser.add_argument('--lr', type=float, default=1e-3,
                        help='Learning rate')

    # Other arguments
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    parser.add_argument('--gpu', type=int, default=0,
                        help='GPU device ID')
    parser.add_argument('--no_cuda', action='store_true',
                        help='Disable CUDA')
    parser.add_argument('--eval_interval', type=int, default=10,
                        help='Evaluation interval (epochs)')

    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    # Set device
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    # Set random seed
    set_seed(args.seed)

    # Create save directory
    os.makedirs(args.save_dir, exist_ok=True)

    # Load and preprocess data using standard interface
    print('Loading data...')
    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=True,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True
    )

    # Convert to numpy arrays
    X = np.array(X).astype(np.float32)
    Y = np.array(Y)

    # Encode labels to integers if needed
    from sklearn.preprocessing import LabelEncoder
    if Y.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y = le.fit_transform(Y)

    # Get number of clusters from data if not specified
    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    # Create dataset and dataloader
    dataset = scRNADataset(X, Y)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True)
    test_loader = DataLoader(dataset, batch_size=args.batch_size * 5, shuffle=False, drop_last=False)

    # Initialize model
    model = AutoEncoder(
        num_genes=X.shape[1],
        hidden_size=args.hidden_size,
        masked_data_weight=args.masked_data_weight,
        mask_loss_weight=args.mask_loss_weight
    ).to(device)

    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    # Mask probabilities
    mask_probas = [args.mask_prob] * X.shape[1]

    # Training loop
    print('Starting training...')
    best_acc = 0
    for epoch in range(args.epochs):
        model.train()
        total_loss = 0
        for x, y in train_loader:
            x = x.to(device)
            x_corrupted, mask = apply_noise(x, mask_probas)

            optimizer.zero_grad()
            _, loss = model.loss_mask(x_corrupted, x, mask)
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        avg_loss = total_loss / len(train_loader)

        # Evaluation
        if (epoch + 1) % args.eval_interval == 0 or epoch == args.epochs - 1:
            # Extract embeddings
            embedding, true_labels = inference(model, test_loader, device)

            # Clustering
            if embedding.shape[0] < 10000:
                kmeans = KMeans(n_clusters=n_clusters, random_state=args.seed, n_init=20)
                pred_labels = kmeans.fit_predict(embedding)
            else:
                # Use Leiden for large datasets
                adata_emb = sc.AnnData(embedding)
                sc.pp.neighbors(adata_emb, n_neighbors=10, use_rep="X")
                reso = res_search_fixed_clus(adata_emb, n_clusters)
                sc.tl.leiden(adata_emb, resolution=reso)
                pred_labels = np.array([int(x) for x in adata_emb.obs['leiden'].to_list()])

            # Save results using standard interface
            save(args.save_dir, true_labels, pred_labels, epoch + 1, embedding)

            print(f'Epoch {epoch + 1}/{args.epochs}, Loss: {avg_loss:.4f}')

    # Save final model
    torch.save({
        'model': model.state_dict(),
        'optimizer': optimizer.state_dict(),
        'args': vars(args)
    }, os.path.join(args.save_dir, 'model_checkpoint.pth'))

    print(f'Training completed. Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
