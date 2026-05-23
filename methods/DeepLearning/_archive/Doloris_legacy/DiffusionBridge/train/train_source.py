from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import torch
from sklearn.decomposition import PCA
from torch.utils.data import DataLoader, TensorDataset

from models.source_diffusion import SourceDiffusion


def read_h5ad_compat(path: str | Path):
    import numpy as np

    if not hasattr(np, 'string_'):
        np.string_ = np.bytes_
    return sc.read_h5ad(str(path))


def preprocess_adata(adata, n_hvg: int = 2000):
    adata = adata.copy()
    if adata.raw is None:
        adata.raw = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=True)
    return adata


def build_teacher_target(adata, latent_dim: int, mode: str = 'pca_graph'):
    x = adata.X.toarray() if hasattr(adata.X, 'toarray') else np.asarray(adata.X)
    pca = PCA(n_components=latent_dim, random_state=42)
    z = pca.fit_transform(x)
    if mode == 'pca_graph':
        sc.pp.neighbors(adata, n_neighbors=15, use_rep='X')
        graph = adata.obsp['connectivities']
        z = graph @ z
    return np.asarray(z, dtype=np.float32)


def train_source(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    adata = preprocess_adata(read_h5ad_compat(args.data_path), n_hvg=args.n_hvg)
    x = adata.X.toarray().astype(np.float32) if hasattr(adata.X, 'toarray') else np.asarray(adata.X, dtype=np.float32)
    dataset = TensorDataset(torch.from_numpy(x))
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True)

    model = SourceDiffusion(domain_dim=x.shape[1], shared_dim=args.latent_dim, hidden_dim=args.hidden_dim, num_steps=args.diffusion_steps).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr)

    for epoch in range(args.epochs):
        total_loss = 0.0
        for (batch_x,) in loader:
            batch_x = batch_x.to(device)
            loss_dict = model.training_loss(batch_x, recon_weight=args.recon_weight, prior_weight=args.prior_weight, zero_weight=args.zero_weight)
            loss = loss_dict['loss']
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * batch_x.size(0)
        print(f"source epoch={epoch+1} loss={total_loss / len(dataset):.6f}")

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    torch.save({'state_dict': model.state_dict(), 'input_dim': x.shape[1]}, output_dir / 'source.pt')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-path', type=str, required=True)
    parser.add_argument('--output-dir', type=str, required=True)
    parser.add_argument('--n-hvg', type=int, default=2000)
    parser.add_argument('--latent-dim', type=int, default=64)
    parser.add_argument('--hidden-dim', type=int, default=256)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--epochs', type=int, default=20)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--diffusion-steps', type=int, default=50)
    parser.add_argument('--recon-weight', type=float, default=1.0)
    parser.add_argument('--prior-weight', type=float, default=1e-3)
    parser.add_argument('--zero-weight', type=float, default=0.25)
    args = parser.parse_args()
    train_source(args)
