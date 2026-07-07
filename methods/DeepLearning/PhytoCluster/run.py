# -*- coding: utf-8 -*-
"""
PhytoCluster — 植物单细胞RNA-seq变分自编码器聚类方法

两阶段流程：
1. 预训练 VAE 学习低维潜在表示
2. 用 GMM 初始化并联合优化 PhytoCluster

参考论文：Wang et al. (2025), aBIOTECH
"""

import os
import sys
import argparse
import math
import json
import random

import numpy as np

if not hasattr(np, 'string_'):
    np.string_ = np.bytes_
if not hasattr(np, 'unicode_'):
    np.unicode_ = np.str_

import pandas as pd
import scanpy as sc
import torch
import torch.nn as nn
import torch.nn.functional as F
import scipy.sparse as sp
from torch.nn import init
from torch.utils.data import DataLoader, TensorDataset
from sklearn.mixture import GaussianMixture
from sklearn.preprocessing import LabelEncoder

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


def parse_args():
    parser = argparse.ArgumentParser(
        description='PhytoCluster: VAE + GMM for Plant scRNA-seq Clustering',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('--data_path', type=str, required=True)
    parser.add_argument('--save_dir', type=str, default='./results')
    parser.add_argument('--n_clusters', type=int, required=True)
    parser.add_argument('--n_top_genes', type=int, default=2000)
    parser.add_argument('--latent_dim', type=int, default=10)
    parser.add_argument('--encode_dim', type=int, nargs=2, default=[1024, 128])
    parser.add_argument('--decode_dim', type=int, nargs=2, default=[128, 1024])
    parser.add_argument('--pretrain_max_iter', type=int, default=30000)
    parser.add_argument('--cluster_max_iter', type=int, default=300)
    parser.add_argument('--batch_size', type=int, default=64)
    parser.add_argument('--lr', type=float, default=1e-4)
    parser.add_argument('--var_lr', type=float, default=1e-4)
    parser.add_argument('--weight_decay', type=float, default=5e-4)
    parser.add_argument('--pretrain_warmup_mult', type=int, default=200)
    parser.add_argument('--cluster_warmup_mult', type=int, default=300)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--gpu', type=int, default=0)
    parser.add_argument('--no_cuda', action='store_true')
    return parser.parse_args()


def binary_cross_entropy_loss(recon_x, x):
    return -torch.sum(x * torch.log(recon_x + 1e-8) + (1 - x) * torch.log(1 - recon_x + 1e-8), dim=-1)


def compute_kl_divergence(mu, logvar):
    return -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1)


def compute_elbo_loss(recon_x, x, z_params, binary=True):
    mu, logvar = z_params
    kld = compute_kl_divergence(mu, logvar)
    if binary:
        likelihood = -binary_cross_entropy_loss(recon_x, x)
    else:
        likelihood = -F.mse_loss(recon_x, x)
    return torch.sum(likelihood), torch.sum(kld)


def compute_elbo(recon_x, x, gamma, c_params, z_params, binary=True):
    mu_c, var_c, pi = c_params
    var_c = var_c + 1e-8
    n_centroids = pi.size(1)
    mu, logvar = z_params
    mu_expand = mu.unsqueeze(2).expand(mu.size(0), mu.size(1), n_centroids)
    logvar_expand = logvar.unsqueeze(2).expand(logvar.size(0), logvar.size(1), n_centroids)

    if binary:
        likelihood = -binary_cross_entropy_loss(recon_x, x)
    else:
        likelihood = -F.mse_loss(recon_x, x)

    logpzc = -0.5 * torch.sum(
        gamma * torch.sum(
            math.log(2 * math.pi) + torch.log(var_c) + torch.exp(logvar_expand) / var_c + (mu_expand - mu_c) ** 2 / var_c,
            dim=1,
        ),
        dim=1,
    )
    logpc = torch.sum(gamma * torch.log(pi), 1)
    qentropy = -0.5 * torch.sum(1 + logvar + math.log(2 * math.pi), 1)
    logqcx = torch.sum(gamma * torch.log(gamma), 1)
    kld = -logpzc - logpc + qentropy + logqcx
    return torch.sum(likelihood), torch.sum(kld)


def create_mlp(layers, activation=nn.ReLU(), bn=False, dropout=0):
    net = []
    for i in range(1, len(layers)):
        net.append(nn.Linear(layers[i - 1], layers[i]))
        if bn:
            net.append(nn.BatchNorm1d(layers[i]))
        net.append(activation)
        if dropout > 0:
            net.append(nn.Dropout(dropout))
    return nn.Sequential(*net)


class DeterministicWarmup:
    def __init__(self, n=100, t_max=1.0):
        self.t = 0
        self.t_max = t_max
        self.inc = t_max / n

    def __iter__(self):
        return self

    def __next__(self):
        t = self.t + self.inc
        self.t = self.t_max if t > self.t_max else t
        return self.t


class Stochastic(nn.Module):
    def reparametrize(self, mu, logvar):
        epsilon = torch.randn(mu.size(), requires_grad=False, device=mu.device)
        std = logvar.mul(0.5).exp_()
        return mu.addcmul(std, epsilon)


class GaussianSample(Stochastic):
    def __init__(self, in_features, out_features):
        super().__init__()
        self.mu = nn.Linear(in_features, out_features)
        self.log_var = nn.Linear(in_features, out_features)

    def forward(self, x):
        mu = self.mu(x)
        log_var = self.log_var(x)
        return self.reparametrize(mu, log_var), mu / 10, log_var


class Encoder(nn.Module):
    def __init__(self, dims, bn=False, dropout=0):
        super().__init__()
        x_dim, h_dim, z_dim = dims
        self.hidden = create_mlp([x_dim] + h_dim, bn=bn, dropout=dropout)
        self.sample = GaussianSample(([x_dim] + h_dim)[-1], z_dim)

    def forward(self, x):
        x = self.hidden(x)
        return self.sample(x)


class Decoder(nn.Module):
    def __init__(self, dims, bn=False, dropout=0, output_activation=None):
        super().__init__()
        z_dim, h_dim, x_dim = dims
        self.hidden = create_mlp([z_dim, *h_dim], bn=bn, dropout=dropout)
        self.reconstruction = nn.Linear([z_dim, *h_dim][-1], x_dim)
        self.output_activation = output_activation

    def forward(self, x):
        x = self.hidden(x)
        x = self.reconstruction(x)
        return self.output_activation(x) if self.output_activation is not None else x


class VAE(nn.Module):
    def __init__(self, dims, bn=False, dropout=0, binary=False):
        super().__init__()
        x_dim, z_dim, encode_dim, decode_dim = dims
        self.binary = binary
        decode_activation = nn.Sigmoid() if binary else None
        self.encoder = Encoder([x_dim, encode_dim, z_dim], bn=bn, dropout=dropout)
        self.decoder = Decoder([z_dim, decode_dim, x_dim], bn=bn, dropout=dropout, output_activation=decode_activation)
        self.initialize_weights()

    def initialize_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                init.xavier_normal_(m.weight.data)
                if m.bias is not None:
                    m.bias.data.zero_()

    def forward(self, x):
        z, mu, logvar = self.encoder(x)
        return self.decoder(z)

    def compute_loss(self, x):
        z, mu, logvar = self.encoder(x)
        recon_x = self.decoder(z)
        likelihood, kl_loss = compute_elbo_loss(recon_x, x, (mu, logvar), binary=self.binary)
        return -likelihood, kl_loss

    def encode_batch(self, dataloader, device='cpu', output_type='z'):
        output = []
        self.eval()
        with torch.no_grad():
            for batch in dataloader:
                x = batch[0] if isinstance(batch, (list, tuple)) else batch
                x = x.view(x.size(0), -1).float().to(device)
                z, mu, logvar = self.encoder(x)
                if output_type == 'z':
                    output.append(z.detach().cpu())
                elif output_type == 'x':
                    output.append(self.decoder(z).detach().cpu())
                elif output_type == 'mu':
                    output.append(mu.detach().cpu())
                elif output_type == 'log_var':
                    output.append(logvar.detach().cpu())
        return torch.cat(output).numpy()


class PhytoCluster(VAE):
    def __init__(self, dims, n_centroids, bn=False, dropout=0, binary=False):
        super().__init__(dims, bn=bn, dropout=dropout, binary=binary)
        self.n_centroids = n_centroids
        z_dim = dims[1]
        self.pi = nn.Parameter(torch.ones(n_centroids) / n_centroids, requires_grad=False)
        self.mu_c = nn.Parameter(torch.zeros(z_dim, n_centroids), requires_grad=False)
        self.var_c = nn.Parameter(torch.ones(z_dim, n_centroids), requires_grad=False)

    def infer_clusters(self, z):
        n_centroids = self.n_centroids
        N = z.size(0)
        z_expanded = z.unsqueeze(2).expand(z.size(0), z.size(1), n_centroids)
        pi = self.pi.repeat(N, 1)
        mu_c = self.mu_c.repeat(N, 1, 1)
        var_c = self.var_c.repeat(N, 1, 1) + 1e-8
        p_c_z = torch.exp(
            torch.log(pi)
            - torch.sum(
                0.5 * torch.log(2 * math.pi * var_c) + (z_expanded - mu_c) ** 2 / (2 * var_c),
                dim=1,
            )
        ) + 1e-10
        gamma = p_c_z / torch.sum(p_c_z, dim=1, keepdim=True)
        return gamma, mu_c, var_c, pi

    def compute_loss(self, x):
        z, mu, logvar = self.encoder(x)
        recon_x = self.decoder(z)
        gamma, mu_c, var_c, pi = self.infer_clusters(z)
        likelihood, kl_loss = compute_elbo(recon_x, x, gamma, (mu_c, var_c, pi), (mu, logvar), binary=self.binary)
        return -likelihood, kl_loss

    def initialize_gmm_params(self, dataloader, device='cpu'):
        gmm = GaussianMixture(n_components=self.n_centroids, covariance_type='diag', random_state=42)
        z = self.encode_batch(dataloader, device, output_type='z')
        gmm.fit(z)
        self.mu_c.data.copy_(torch.from_numpy(gmm.means_.T.astype(np.float32)))
        self.var_c.data.copy_(torch.from_numpy(gmm.covariances_.T.astype(np.float32)))

    def predict(self, dataloader, device):
        self.eval()
        embedding_mu = self.encode_batch(dataloader, device, output_type='mu')
        gamma_input = torch.from_numpy(embedding_mu).float().to(device)
        with torch.no_grad():
            gamma, _, _, _ = self.infer_clusters(gamma_input)
        return embedding_mu, gamma.cpu().numpy()


def train_model(model, dataloader, device, args, phase='pretrain'):
    model.to(device)
    if phase == 'pretrain':
        optimizer = torch.optim.Adam([
            {'params': model.encoder.hidden.parameters(), 'lr': args.lr},
            {'params': model.encoder.sample.mu.parameters(), 'lr': args.var_lr},
            {'params': model.encoder.sample.log_var.parameters(), 'lr': args.var_lr},
            {'params': model.decoder.parameters(), 'lr': args.lr},
        ], weight_decay=args.weight_decay)
        max_iter = args.pretrain_max_iter
        warmup_steps = max(1, args.pretrain_max_iter * args.pretrain_warmup_mult)
    else:
        optimizer = torch.optim.Adam([
            {'params': model.encoder.hidden.parameters(), 'lr': args.lr},
            {'params': model.encoder.sample.mu.parameters(), 'lr': args.var_lr},
            {'params': model.encoder.sample.log_var.parameters(), 'lr': args.var_lr},
            {'params': model.decoder.parameters(), 'lr': args.lr},
        ], weight_decay=args.weight_decay)
        max_iter = args.cluster_max_iter
        warmup_steps = max(1, args.cluster_max_iter * args.cluster_warmup_mult)

    beta_scheduler = DeterministicWarmup(n=warmup_steps, t_max=1.0)
    n_epochs = int(np.ceil(max_iter / len(dataloader)))

    for epoch in range(n_epochs):
        model.train()
        epoch_loss = 0.0
        for batch in dataloader:
            x = batch[0] if isinstance(batch, (list, tuple)) else batch
            x = x.float().to(device)
            optimizer.zero_grad()
            recon_loss, kl_loss = model.compute_loss(x)
            beta = next(beta_scheduler)
            loss = recon_loss + beta * kl_loss
            if torch.isnan(loss):
                raise FloatingPointError('NaN encountered during training')
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 10)
            optimizer.step()
            epoch_loss += loss.item() / len(x)
        if (epoch + 1) % 20 == 0 or epoch == n_epochs - 1:
            print(f'[{phase}] Epoch {epoch + 1}/{n_epochs}: loss={epoch_loss:.4f}')
    return model


def main():
    args = parse_args()
    args.cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f'cuda:{args.gpu}' if args.cuda else 'cpu')
    print(f'Using device: {device}')

    set_seed(args.seed)
    os.makedirs(args.save_dir, exist_ok=True)

    print('Loading data...')
    X_df, Y, _, adata = prepare_data_for_model(
        args.data_path,
        size_factors=False,
        filter_min_counts=True,
        logtrans_input=True,
        normalize_input=True,
        n_top_genes=args.n_top_genes,
    )
    Y = np.array(Y)
    if Y.dtype.kind not in ['i', 'u']:
        Y = LabelEncoder().fit_transform(Y)

    X = X_df.to_numpy(dtype=np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    n_clusters = args.n_clusters if args.n_clusters > 0 else len(np.unique(Y))
    print(f'Number of cells: {X.shape[0]}, Number of genes: {X.shape[1]}')
    print(f'Number of clusters: {n_clusters}')

    input_dim = X.shape[1]
    dims = [input_dim, args.latent_dim, args.encode_dim, args.decode_dim]

    all_data = torch.from_numpy(X).float()
    dataset = TensorDataset(all_data, torch.zeros(len(X)).long())
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True)
    full_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, drop_last=False)

    print('\nPhase 1: Pretraining VAE')
    pretrain_model = VAE(dims, binary=False)
    pretrain_model = train_model(pretrain_model, train_loader, device, args, phase='pretrain')

    print('\nPhase 2: PhytoCluster joint training')
    model = PhytoCluster(dims, n_clusters, binary=False)
    model.load_state_dict(pretrain_model.state_dict(), strict=False)
    model.to(device)
    model.initialize_gmm_params(full_loader, device)
    model = train_model(model, train_loader, device, args, phase='cluster')

    print('\nPhase 3: Final prediction')
    embedding, gamma = model.predict(full_loader, device)

    final_gmm = GaussianMixture(
        n_components=n_clusters,
        covariance_type='diag',
        random_state=args.seed,
    )
    pred_labels = final_gmm.fit_predict(embedding)

    save(args.save_dir, Y, pred_labels, 0, embedding)
    np.save(os.path.join(args.save_dir, 'gamma_0.npy'), gamma)
    np.save(os.path.join(args.save_dir, 'gmm_labels_0.npy'), pred_labels)

    with open(os.path.join(args.save_dir, 'gmm_params.json'), 'w') as f:
        json.dump({
            'pi': model.pi.detach().cpu().numpy().tolist(),
            'mu_c': model.mu_c.detach().cpu().numpy().tolist(),
            'var_c': model.var_c.detach().cpu().numpy().tolist(),
            'final_gmm_weights': final_gmm.weights_.tolist(),
            'final_gmm_means': final_gmm.means_.tolist(),
            'final_gmm_covariances': final_gmm.covariances_.tolist(),
            'n_centroids': n_clusters,
        }, f)

    torch.save(model.state_dict(), os.path.join(args.save_dir, 'phytocluster_model.pt'))
    print(f'PhytoCluster completed. Results saved to: {args.save_dir}')


if __name__ == '__main__':
    main()
