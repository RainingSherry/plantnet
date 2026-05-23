# -*- coding: utf-8 -*-
"""
Unified scCDCG Model Interface for scBench
==========================================

Usage:
    python run.py --data_path /path/to/data.h5ad --n_clusters 10 --save_dir ./results
"""

import os
import sys
import argparse
import numpy as np
import torch
import random
import pickle

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from preprocess import prepare_data_for_model
from utils import save

from sklearn.cluster import KMeans
from torchmetrics.functional import pairwise_cosine_similarity

# Import model components
from model import AE_NN, FULL_NN, ClusterAssignment
from utils import get_laplace_matrix
import torch.nn as nn


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


def sinkhorn(pred, lambdas, row, col):
    """Sinkhorn algorithm for optimal transport."""
    num_node = pred.shape[0]
    num_class = pred.shape[1]
    p = np.power(pred, lambdas)

    u = np.ones(num_node)
    v = np.ones(num_class)

    for index in range(1000):
        u = row * np.power(np.dot(p, v), -1)
        u[np.isinf(u)] = -9e-15
        v = col * np.power(np.dot(u, p), -1)
        v[np.isinf(v)] = -9e-15
    u = row * np.power(np.dot(p, v), -1)
    target = np.dot(np.dot(np.diag(u), p), np.diag(v))
    return target


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='scCDCG: Graph Neural Network for scRNA-seq Clustering',
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
    parser.add_argument('--embedding_dim', type=int, default=16,
                        help='Dimension of embedding space')
    parser.add_argument('--hidden_dim', type=int, default=256,
                        help='Hidden layer dimension')

    # Loss weights
    parser.add_argument('--factor_construct', type=float, default=0.23,
                        help='Weight for reconstruction loss')
    parser.add_argument('--factor_ort', type=float, default=0.65,
                        help='Weight for orthogonality loss')
    parser.add_argument('--factor_corvar', type=float, default=0.17,
                        help='Weight for covariance loss')
    parser.add_argument('--factor_KL', type=float, default=0.12,
                        help='Weight for KL divergence loss')
    parser.add_argument('--balancer', type=float, default=0.55,
                        help='Balancer for Laplacian matrices')
    parser.add_argument('--lambdas', type=float, default=5,
                        help='Lambda for sinkhorn algorithm')

    # Training arguments
    parser.add_argument('--epochs', type=int, default=200,
                        help='Number of training epochs')
    parser.add_argument('--lr', type=float, default=1e-3,
                        help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=5e-3,
                        help='Weight decay')

    # Other arguments
    parser.add_argument('--seed', type=int, default=3047,
                        help='Random seed')
    parser.add_argument('--gpu', type=int, default=0,
                        help='GPU device ID')
    parser.add_argument('--no_cuda', action='store_true',
                        help='Disable CUDA')

    return parser.parse_args()


def main():
    """Main function."""
    args = parse_args()

    # Set device
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    if args.cuda:
        torch.cuda.set_device(args.gpu)

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

    # Convert to tensors
    x = torch.tensor(X, dtype=torch.float)
    y = torch.tensor(Y, dtype=torch.float)

    # Normalize and compute adjacency matrices
    x_ = torch.nn.functional.normalize(x, p=2, dim=1)
    adj_self_loop = torch.mm(x_, x_.T)
    adj_f = torch.abs(pairwise_cosine_similarity(x_, x_))
    adj_f = torch.mm(adj_f, adj_f.T)
    L_1 = get_laplace_matrix(adj_self_loop)
    L_2 = get_laplace_matrix(adj_f)

    # Model dimensions
    dims_encoder = [args.hidden_dim, args.embedding_dim]
    dims_decoder = [args.embedding_dim, args.hidden_dim]

    # Paths for saving pretrained model
    pretrain_model_path = os.path.join(args.save_dir, 'pretrain_model.pkl')
    pretrain_centers_path = os.path.join(args.save_dir, 'pretrain_centers.pkl')
    pretrain_labels_path = os.path.join(args.save_dir, 'pretrain_labels.pkl')

    # ==================== PRE-TRAINING ====================
    print('Pre-training autoencoder...')
    Model = AE_NN(dim_input=x.shape[1], dims_encoder=dims_encoder, dims_decoder=dims_decoder)
    if args.cuda:
        Model = Model.cuda()

    optimizer = torch.optim.Adam(Model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    acc_max = 0
    for epoch in range(1, args.epochs + 1):
        Model.train()

        if args.cuda:
            h, x_hat = Model.forward(x.cuda(), adj_self_loop.cuda())
        else:
            h, x_hat = Model.forward(x, adj_self_loop)

        z = torch.nn.functional.normalize(h, p=2, dim=0)

        if args.cuda:
            loss_x = torch.nn.functional.mse_loss(x_hat, x.cuda())
            loss_corvariates = -torch.mm(torch.mm(z.T, (args.balancer * L_1.cuda() + (1-args.balancer) * L_2.cuda())), z).trace() / len(z.T)
            loss_ort = torch.nn.functional.mse_loss(torch.mm(z.T, z).view(-1).cuda(), torch.eye(len(z.T)).view(-1).cuda())
        else:
            loss_x = torch.nn.functional.mse_loss(x_hat, x)
            loss_corvariates = -torch.mm(torch.mm(z.T, (args.balancer * L_1 + (1-args.balancer) * L_2)), z).trace() / len(z.T)
            loss_ort = torch.nn.functional.mse_loss(torch.mm(z.T, z).view(-1), torch.eye(len(z.T)).view(-1))

        loss = args.factor_construct * loss_x + args.factor_ort * loss_ort + args.factor_corvar * loss_corvariates

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        with torch.no_grad():
            kmeans = KMeans(n_clusters=n_clusters, random_state=args.seed, n_init=20).fit(z.cpu().numpy())
            centers = torch.tensor(kmeans.cluster_centers_)

            if acc_max == 0 or epoch == args.epochs:
                acc_max = 1
                torch.save(Model.state_dict(), pretrain_model_path)
                with open(pretrain_centers_path, 'wb') as f:
                    pickle.dump(centers, f, protocol=pickle.HIGHEST_PROTOCOL)
                pseudo_labels = torch.LongTensor(kmeans.labels_)
                with open(pretrain_labels_path, 'wb') as f:
                    pickle.dump(pseudo_labels, f, protocol=pickle.HIGHEST_PROTOCOL)

        if epoch % 50 == 0:
            print(f'Pre-train Epoch {epoch}/{args.epochs}, Loss: {loss.item():.6f}')

    # ==================== FINE-TUNING ====================
    print('Fine-tuning with clustering...')
    Model = FULL_NN(
        dim_input=x.shape[1],
        dims_encoder=dims_encoder,
        dims_decoder=dims_decoder,
        num_class=n_clusters,
        pretrain_model_load_path=pretrain_model_path
    )
    if args.cuda:
        Model = Model.cuda()

    optimizer = torch.optim.Adam(Model.parameters(), lr=args.lr)

    with open(pretrain_centers_path, 'rb') as f:
        centers = pickle.load(f)
        if args.cuda:
            centers = centers.cuda()
    with open(pretrain_labels_path, 'rb') as f:
        pseudo_labels = pickle.load(f)
        if args.cuda:
            pseudo_labels = pseudo_labels.cuda()

    best_embedding = None
    best_y_pred = None
    acc_max = 0

    for epoch in range(1, args.epochs + 1):
        Model.train()

        if args.cuda:
            z, x_hat = Model.forward(x.cuda(), adj_self_loop.cuda())
        else:
            z, x_hat = Model.forward(x, adj_self_loop)

        z = torch.nn.functional.normalize(z, p=2, dim=0)
        centers = centers.detach()

        if args.cuda:
            loss_x = torch.nn.functional.mse_loss(x_hat, x.cuda())
            loss_corvariates = -torch.mm(torch.mm(z.T, (args.balancer * L_1.cuda() + (1-args.balancer) * L_2.cuda())), z).trace() / len(z.T)
            loss_ort = torch.nn.functional.mse_loss(torch.mm(z.T, z).view(-1).cuda(), torch.eye(len(z.T)).view(-1).cuda())
        else:
            loss_x = torch.nn.functional.mse_loss(x_hat, x)
            loss_corvariates = -torch.mm(torch.mm(z.T, (args.balancer * L_1 + (1-args.balancer) * L_2)), z).trace() / len(z.T)
            loss_ort = torch.nn.functional.mse_loss(torch.mm(z.T, z).view(-1), torch.eye(len(z.T)).view(-1))

        # DEC clustering
        if args.cuda:
            class_assign_model = ClusterAssignment(n_clusters, len(z.T), 1, centers).cuda()
            temp_class = class_assign_model(z.cuda())
        else:
            class_assign_model = ClusterAssignment(n_clusters, len(z.T), 1, centers)
            temp_class = class_assign_model(z)

        # Sinkhorn for target distribution
        if epoch == 1 or epoch % 10 == 0:
            p_distribution = torch.tensor(
                sinkhorn(
                    temp_class.cpu().detach().numpy(),
                    args.lambdas,
                    torch.ones(x.shape[0]).numpy(),
                    torch.tensor([torch.sum(pseudo_labels.cpu() == i) for i in range(n_clusters)]).numpy()
                )
            ).float()
            if args.cuda:
                p_distribution = p_distribution.cuda()
            p_distribution = p_distribution.detach()

        KL_loss_function = nn.KLDivLoss(reduction='sum')
        if args.cuda:
            loss_KL = KL_loss_function(temp_class.cuda(), p_distribution.cuda())
        else:
            loss_KL = KL_loss_function(temp_class, p_distribution)

        loss = args.factor_construct * loss_x + args.factor_ort * loss_ort + args.factor_corvar * loss_corvariates + args.factor_KL * loss_KL

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        with torch.no_grad():
            kmeans = KMeans(n_clusters=n_clusters, random_state=args.seed, n_init=20).fit(z.cpu().numpy())
            y_pred = kmeans.labels_
            pseudo_labels = torch.LongTensor(kmeans.labels_)
            if args.cuda:
                pseudo_labels = pseudo_labels.cuda()
            centers = torch.tensor(kmeans.cluster_centers_)
            if args.cuda:
                centers = centers.cuda()

            # Save best results
            best_embedding = z.cpu().numpy()
            best_y_pred = y_pred

        if epoch % 50 == 0:
            print(f'Fine-tune Epoch {epoch}/{args.epochs}, Loss: {loss.item():.6f}')

    # Save results using standard interface
    save(args.save_dir, Y, best_y_pred, args.epochs, best_embedding)

    print(f'Training completed.')
    print(f'Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
