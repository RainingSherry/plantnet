# -*- coding: utf-8 -*-
"""
PhytoCluster: VAE + GMM for Plant Single-Cell RNA-seq Clustering

Based on:
  Wang et al. (2025). PhytoCluster: a generative deep learning model
  for clustering plant single-cell RNA-seq data. aBIOTECH.
  GitHub: https://github.com/Llana-168/PhytoCluster
"""

import os
import sys
import argparse
import math
import random
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from torch.nn import init
from torch.utils.data import DataLoader, TensorDataset
from sklearn.mixture import GaussianMixture

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from preprocess import prepare_data_for_model
from utils import save


def set_seed(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


class EarlyStopping:
    def __init__(self, patience=10, verbose=False, outdir=None):
        self.patience = patience
        self.verbose = verbose
        self.counter = 0
        self.best_score = None
        self.early_stop = False
        self.loss_min = np.inf
        self.model_file = os.path.join(outdir, 'model.pt') if outdir else None

    def __call__(self, loss, model):
        loss_val = loss.item() if torch.is_tensor(loss) else loss
        if np.isnan(loss_val):
            self.early_stop = True
        score = -loss_val
        if self.best_score is None:
            self.best_score = score
            self.save_checkpoint(loss, model)
        elif score < self.best_score:
            self.counter += 1
            if self.verbose:
                print(f'EarlyStopping counter: {self.counter} out of {self.patience}')
            if self.counter >= self.patience:
                self.early_stop = True
                if self.model_file and os.path.exists(self.model_file):
                    sd = torch.load(self.model_file, map_location='cpu', weights_only=False)
                    model.load_state_dict(sd, strict=False)
        else:
            self.best_score = score
            self.save_checkpoint(loss, model)
            self.counter = 0

    def save_checkpoint(self, loss, model):
        if self.verbose:
            print(f'Loss decreased ({self.loss_min:.6f} --> {loss:.6f}). Saving model ...')
        if self.model_file:
            torch.save(model.state_dict(), self.model_file)
        self.loss_min = loss.item() if torch.is_tensor(loss) else loss


def create_mlp(layers, activation=nn.ReLU(), bn=False, dropout=0):
    net = []
    for i in range(1, len(layers)):
        net.append(nn.Linear(layers[i - 1], layers[i]))
        if bn:
            net.append(nn.BatchNorm1d(layers[i - 1]))
        net.append(activation)
        if dropout > 0:
            net.append(nn.Dropout(dropout))
    return nn.Sequential(*net)


class DeterministicWarmup:
    def __init__(self, n=100, t_max=1):
        self.t = 0
        self.t_max = t_max
        self.inc = t_max / n

    def __iter__(self):
        return self

    def __next__(self):
        self.t = min(self.t + self.inc, self.t_max)
        return self.t


class Stochastic(nn.Module):
    def reparametrize(self, mu, logvar):
        epsilon = torch.randn(mu.size(), requires_grad=False, device=mu.device)
        std = logvar.mul(0.5).exp_()
        return mu + std * epsilon


class GaussianSample(Stochastic):
    def __init__(self, in_features, out_features):
        super(GaussianSample, self).__init__()
        self.mu = nn.Linear(in_features, out_features)
        self.log_var = nn.Linear(in_features, out_features)

    def forward(self, x):
        mu = self.mu(x)
        log_var = self.log_var(x)
        return self.reparametrize(mu, log_var), mu / 10, log_var


class Encoder(nn.Module):
    def __init__(self, dims, bn=False, dropout=0):
        super(Encoder, self).__init__()
        x_dim, h_dim, z_dim = dims[0], dims[1], dims[2]
        hidden_layers = [x_dim] + h_dim
        self.hidden = create_mlp(hidden_layers, bn=bn, dropout=dropout)
        self.sample = GaussianSample(hidden_layers[-1], z_dim)

    def forward(self, x):
        return self.sample(self.hidden(x))


class Decoder(nn.Module):
    def __init__(self, dims, bn=False, dropout=0):
        super(Decoder, self).__init__()
        z_dim, h_dim, x_dim = dims[0], dims[1], dims[2]
        hidden_layers = [z_dim] + h_dim
        self.hidden = create_mlp(hidden_layers, bn=bn, dropout=dropout)
        self.reconstruction = nn.Linear(hidden_layers[-1], x_dim)

    def forward(self, x):
        return self.reconstruction(self.hidden(x))


class VAE(nn.Module):
    """Variational Autoencoder for latent feature extraction."""

    def __init__(self, dims, bn=False, dropout=0):
        super(VAE, self).__init__()
        x_dim, z_dim, encode_dim, decode_dim = dims[0], dims[1], dims[2], dims[3]
        self.encoder = Encoder([x_dim, encode_dim, z_dim], bn=bn, dropout=dropout)
        self.decoder = Decoder([z_dim, decode_dim, x_dim], bn=bn, dropout=dropout)
        for m in self.modules():
            if isinstance(m, nn.Linear):
                init.xavier_normal_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.zero_()

    def forward(self, x):
        z, mu, logvar = self.encoder(x)
        return self.decoder(z), mu, logvar

    def compute_loss(self, x):
        recon_x, mu, logvar = self.forward(x)
        # MSE reconstruction loss
        recon_loss = F.mse_loss(recon_x, x, reduction='mean') * x.size(1)
        # KL divergence with clamping for stability
        logvar = torch.clamp(logvar, min=-10, max=10)
        kld = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1).mean()
        return recon_loss, kld

    def encode_batch(self, dataloader, device='cpu', output_type='z'):
        self.eval()
        output = []
        with torch.no_grad():
            for item in dataloader:
                if isinstance(item, (tuple, list)):
                    x = item[0]
                else:
                    x = item
                x = x.view(x.size(0), -1).float().to(device)
                recon_x, mu, logvar = self.forward(x)
                if output_type == 'z':
                    z, _, _ = self.encoder(x)
                    output.append(z.detach().cpu())
                elif output_type == 'mu':
                    output.append(mu.cpu().detach())
                elif output_type == 'x':
                    output.append(recon_x.cpu().detach())
        return torch.cat(output, dim=0).numpy()


class PhytoCluster(VAE):
    """
    PhytoCluster: VAE + GMM for joint clustering.

    Adds GMM clustering parameters (pi, mu_c, var_c) to the VAE,
    and jointly optimizes the clustering objective function.
    """

    def __init__(self, dims, n_centroids):
        super(PhytoCluster, self).__init__(dims)
        self.n_centroids = n_centroids
        z_dim = dims[1]
        self.pi = nn.Parameter(torch.ones(n_centroids) / n_centroids)
        self.mu_c = nn.Parameter(torch.zeros(z_dim, n_centroids))
        self.var_c = nn.Parameter(torch.ones(z_dim, n_centroids))

    def infer_clusters(self, z):
        """Compute posterior cluster assignments q(c|z) using Bayes rule."""
        N = z.size(0)
        K = self.n_centroids
        z_exp = z.unsqueeze(2).expand(N, z.size(1), K)
        pi = self.pi.unsqueeze(0).expand(N, K)
        mu_c = self.mu_c.unsqueeze(0).expand(N, z.size(1), K)
        var_c = self.var_c.unsqueeze(0).expand(N, z.size(1), K) + 1e-6

        log_p_z_given_c = -0.5 * (
            torch.log(2 * math.pi * var_c) + (z_exp - mu_c) ** 2 / var_c
        ).sum(dim=1)

        log_joint = torch.log(pi + 1e-10) + log_p_z_given_c
        log_joint_max = log_joint.max(dim=1, keepdim=True)[0]
        log_evidence = torch.log(torch.exp(log_joint - log_joint_max).sum(dim=1, keepdim=True)) + log_joint_max
        gamma = torch.exp(log_joint - log_evidence)

        return gamma, mu_c, var_c, pi

    def compute_loss(self, x):
        """
        PhytoCluster loss = reconstruction + clustering KL divergence.

        Clustering KL includes:
          -E_q[log p(z|c)]: expected log-likelihood under cluster posterior
          -E_q[log p(c)]: cluster prior
          -H(q(c|z)): entropy of cluster posterior
          -H(q(z|x)): VAE KL entropy term
        """
        recon_x, mu, logvar = self.forward(x)
        gamma, _, _, _ = self.infer_clusters(mu)

        # Reconstruction loss
        recon_loss = F.mse_loss(recon_x, x, reduction='mean') * x.size(1)

        # Clustering KL divergence
        N = mu.size(0)
        z_dim = mu.size(1)
        K = self.n_centroids

        logvar_clamped = torch.clamp(logvar, min=-10, max=10)
        mu_c_expanded = self.mu_c.unsqueeze(0).expand(N, z_dim, K)
        var_c_expanded = (self.var_c + 1e-6).unsqueeze(0).expand(N, z_dim, K)
        logvar_expanded = logvar_clamped.unsqueeze(2).expand(N, z_dim, K)
        mu_expanded = mu.unsqueeze(2).expand(N, z_dim, K)

        log_p_z_given_c = -0.5 * (
            torch.log(2 * math.pi * var_c_expanded) +
            torch.exp(logvar_expanded) / var_c_expanded +
            (mu_expanded - mu_c_expanded) ** 2 / var_c_expanded
        ).sum(dim=1)

        logpzc = (gamma * log_p_z_given_c).sum(dim=1)
        logpc = (gamma * torch.log(self.pi.unsqueeze(0) + 1e-10)).sum(dim=1)
        logqcx = -(gamma * torch.log(gamma + 1e-10)).sum(dim=1)
        kl_entropy = -0.5 * torch.sum(1 + logvar_clamped + math.log(2 * math.pi), dim=1)

        kl_loss = (-logpzc - logpc + kl_entropy + logqcx).mean()

        return recon_loss, kl_loss

    def initialize_gmm_params(self, dataloader, device='cpu'):
        """Initialize GMM parameters using pretrained VAE latent features."""
        gmm = GaussianMixture(
            n_components=self.n_centroids,
            covariance_type='diag',
            random_state=42,
            n_init=10,
            reg_covar=1e-6
        )
        z = self.encode_batch(dataloader, device)
        mask = ~np.isnan(z).any(axis=1)
        z_clean = z[mask] if mask.sum() >= self.n_centroids else z
        gmm.fit(z_clean)
        with torch.no_grad():
            self.mu_c.copy_(torch.from_numpy(gmm.means_.T.astype(np.float32)))
            self.var_c.copy_(torch.from_numpy(gmm.covariances_.T.astype(np.float32)))


def fit(model, dataloader,
        lr=0.001,
        weight_decay=1e-4,
        device='cuda',
        beta=1.0,
        n=2000,
        max_iter=30000,
        verbose=True,
        patience=30,
        outdir=None):
    """Train VAE or PhytoCluster model."""
    model.to(device)

    encoder_params = (list(model.encoder.hidden.parameters()) +
                      list(model.encoder.sample.parameters()))
    decoder_params = list(model.decoder.parameters())
    all_params = [{'params': encoder_params, 'lr': lr},
                  {'params': decoder_params, 'lr': lr}]

    if isinstance(model, PhytoCluster):
        all_params += [{'params': [model.pi, model.mu_c, model.var_c], 'lr': lr * 0.1}]

    optimizer = torch.optim.Adam(all_params, weight_decay=weight_decay)

    Beta = DeterministicWarmup(n=n, t_max=beta)
    n_epoch = int(np.ceil(max_iter / len(dataloader)))
    early_stopping = EarlyStopping(patience=patience, outdir=outdir, verbose=verbose)

    for epoch in range(n_epoch):
        model.train()
        epoch_loss = 0.0

        for item in dataloader:
            if isinstance(item, (tuple, list)):
                x = item[0]
            else:
                x = item
            x = x.float().to(device)
            optimizer.zero_grad()

            recon_loss, kl_loss = model.compute_loss(x)
            b = next(Beta)
            loss = recon_loss + b * kl_loss

            if torch.isnan(loss):
                continue

            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            epoch_loss += loss.item()

        avg_loss = epoch_loss / len(dataloader)

        if verbose and (epoch + 1) % 5 == 0:
            print(f"Epoch {epoch + 1}/{n_epoch}: loss={avg_loss:.4f} "
                  f"(rec={recon_loss.item():.4f}, kl={kl_loss.item():.4f}, beta={b:.4f})")

        early_stopping(avg_loss, model)
        if early_stopping.early_stop:
            if verbose:
                print(f"Early stopping triggered at epoch {epoch + 1}")
            break

    return model


def parse_args():
    parser = argparse.ArgumentParser(
        description='PhytoCluster: VAE + GMM for plant scRNA-seq clustering',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('--data_path', type=str, required=True,
                       help='Input h5ad file path')
    parser.add_argument('--save_dir', type=str, default='./results',
                       help='Output directory')
    parser.add_argument('--n_clusters', type=int, required=True,
                       help='Number of clusters (K)')
    parser.add_argument('--latent_dim', type=int, default=10,
                       help='Latent space dimension')
    parser.add_argument('--encode_dim', type=str, default='1024,128',
                       help='Encoder hidden dims (comma-separated)')
    parser.add_argument('--decode_dim', type=str, default='128,1024',
                       help='Decoder hidden dims (comma-separated)')
    parser.add_argument('--pretrain_iter', type=int, default=30000,
                       help='VAE pretrain iterations')
    parser.add_argument('--finetune_iter', type=int, default=300,
                       help='PhytoCluster finetune iterations')
    parser.add_argument('--batch_size', type=int, default=64,
                       help='Batch size')
    parser.add_argument('--lr', type=float, default=0.001,
                       help='Learning rate')
    parser.add_argument('--weight_decay', type=float, default=1e-4,
                       help='Weight decay')
    parser.add_argument('--patience', type=int, default=30,
                       help='EarlyStopping patience')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed')
    parser.add_argument('--gpu', type=int, default=0,
                       help='GPU device')
    parser.add_argument('--no_cuda', action='store_true',
                       help='Disable CUDA')
    return parser.parse_args()


def main():
    args = parse_args()

    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    set_seed(args.seed)
    os.makedirs(args.save_dir, exist_ok=True)

    # Step 1: Load and preprocess data
    print('=' * 60)
    print('Step 1: Loading and preprocessing data...')
    print('=' * 60)

    X, Y, sf, adata = prepare_data_for_model(
        args.data_path,
        size_factors=True,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True
    )

    X_np = np.array(X).astype(np.float32)
    Y_np = np.array(Y)

    from sklearn.preprocessing import LabelEncoder
    if Y_np.dtype.kind not in ['i', 'u']:
        le = LabelEncoder()
        Y_np = le.fit_transform(Y_np)

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y_np))
    print(f'Number of cells: {X_np.shape[0]}, Number of genes: {X_np.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    # Step 2: Prepare data loaders
    print('\n' + '=' * 60)
    print('Step 2: Preparing data loaders...')
    print('=' * 60)

    all_data_tensor = torch.from_numpy(X_np).float()
    all_labels_tensor = torch.from_numpy(Y_np).long()
    full_dataset = TensorDataset(all_data_tensor, all_labels_tensor)

    train_loader = DataLoader(
        full_dataset, batch_size=args.batch_size, shuffle=True, drop_last=False
    )
    eval_loader = DataLoader(
        full_dataset, batch_size=args.batch_size * 4, shuffle=False, drop_last=False
    )

    # Step 3: Build model
    print('\n' + '=' * 60)
    print('Step 3: Building model...')
    print('=' * 60)

    input_dim = X_np.shape[1]
    encode_dim = [int(d.strip()) for d in args.encode_dim.split(',')]
    decode_dim = [int(d.strip()) for d in args.decode_dim.split(',')]
    dims = [input_dim, args.latent_dim, encode_dim, decode_dim]

    print(f'Input dim: {input_dim}, Latent dim: {args.latent_dim}')
    print(f'Encode dims: {encode_dim}')
    print(f'Decode dims: {decode_dim}')
    print(f'Number of centroids: {n_clusters}')

    # Step 4: Stage 1 - VAE Pretraining
    print('\n' + '=' * 60)
    print('Step 4: Stage 1 - VAE Pretraining...')
    print('=' * 60)

    pretrain_model = VAE(dims, bn=False, dropout=0).to(device)
    pretrain_model = fit(
        pretrain_model, train_loader,
        lr=args.lr,
        weight_decay=args.weight_decay,
        device=device,
        beta=1.0,
        n=args.pretrain_iter,
        max_iter=args.pretrain_iter,
        verbose=True,
        patience=args.patience,
        outdir=args.save_dir
    )

    # Step 5: GMM initialization
    print('\n' + '=' * 60)
    print('Step 5: Initializing GMM parameters...')
    print('=' * 60)

    pretrain_feat = pretrain_model.encode_batch(eval_loader, device=device, output_type='mu')
    print(f'Pretrained latent features shape: {pretrain_feat.shape}')
    print(f'NaN in features: {np.isnan(pretrain_feat).sum()}')

    # Step 6: Stage 2 - PhytoCluster Joint Training
    print('\n' + '=' * 60)
    print('Step 6: Stage 2 - PhytoCluster Joint Training...')
    print('=' * 60)

    phyto_model = PhytoCluster(dims, n_clusters).to(device)
    # Load VAE weights carefully to avoid name collision with GMM params
    pretrain_sd = pretrain_model.state_dict()
    phyto_sd = phyto_model.state_dict()
    loaded_keys = []
    for key in pretrain_sd:
        if key in phyto_sd:
            phyto_sd[key] = pretrain_sd[key]
            loaded_keys.append(key)
    phyto_model.load_state_dict(phyto_sd, strict=False)
    print(f'Loaded {len(loaded_keys)} VAE weights into PhytoCluster')
    phyto_model.initialize_gmm_params(eval_loader, device=device)

    warmup_n = max(args.finetune_iter // 5, 100)
    phyto_model = fit(
        phyto_model, train_loader,
        lr=args.lr * 0.5,
        weight_decay=args.weight_decay,
        device=device,
        beta=1.0,
        n=warmup_n,
        max_iter=args.finetune_iter,
        verbose=True,
        patience=args.patience,
        outdir=args.save_dir
    )

    # Step 7: Clustering
    print('\n' + '=' * 60)
    print('Step 7: Extracting features and clustering...')
    print('=' * 60)

    phyto_model.eval()
    final_feat = phyto_model.encode_batch(eval_loader, device=device, output_type='mu')
    print(f'Final latent features shape: {final_feat.shape}')

    # GMM final clustering
    gmm = GaussianMixture(
        n_components=n_clusters, covariance_type='diag',
        random_state=args.seed, n_init=10, reg_covar=1e-6
    )
    pred_labels = gmm.fit_predict(final_feat)

    # Step 8: Save results
    print('\n' + '=' * 60)
    print('Step 8: Saving results...')
    print('=' * 60)

    save(args.save_dir, Y_np, pred_labels, 1, final_feat)
    print(f'\nPhytoCluster completed. Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
