from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import scanpy as sc
import torch
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score, silhouette_score
from torch.utils.data import DataLoader, TensorDataset

from models.bridge import BridgeSampler
from models.cluster_head import ClusterHead
from models.source_diffusion import DiffusionBridge, SourceDiffusion, TargetDiffusion
from models.support_mask import GeneSupportMask, build_support_mask


# The bridge stage is the paper-like part.
# Raw sparse counts are inverted into a shared Gaussian latent, then decoded into a cluster-separable target domain.
# We supervise the target domain with a teacher embedding and use DEC-style self-training to sharpen cluster boundaries.

def read_h5ad_compat(path: str | Path):
    import numpy as np

    if not hasattr(np, 'string_'):
        np.string_ = np.bytes_
    return sc.read_h5ad(str(path))


def preprocess_adata(adata, n_hvg: int = 2000):
    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_hvg, subset=True)
    return adata


def extract_labels(adata):
    for key in ['cell_type', 'Celltype', 'celltype', 'celltype_after', 'cell_label', 'label', 'Seurat_clusters']:
        if key in adata.obs.columns:
            return adata.obs[key].astype(str).to_numpy(), key
    raise KeyError('No label column found in adata.obs')


def load_teacher_target(output_dir: Path, fallback_shape: tuple[int, int]) -> np.ndarray:
    teacher_path = output_dir / 'teacher_target.npy'
    if teacher_path.exists():
        teacher_target = np.load(teacher_path)
        return teacher_target.astype(np.float32)
    return np.zeros(fallback_shape, dtype=np.float32)


def load_bridge(args, input_dim: int, device: torch.device):
    source = SourceDiffusion(domain_dim=input_dim, shared_dim=args.latent_dim, hidden_dim=args.hidden_dim, num_steps=args.diffusion_steps)
    target = TargetDiffusion(domain_dim=args.latent_dim, shared_dim=args.latent_dim, hidden_dim=args.hidden_dim, num_steps=args.diffusion_steps)
    source_ckpt = torch.load(Path(args.output_dir) / 'source.pt', map_location='cpu', weights_only=True)
    target_ckpt = torch.load(Path(args.output_dir) / 'target.pt', map_location='cpu', weights_only=True)
    source.load_state_dict(source_ckpt['state_dict'])
    target.load_state_dict(target_ckpt['state_dict'])
    bridge = DiffusionBridge(source=source, target=target, support_mask=GeneSupportMask(blend=args.support_blend))
    return bridge.to(device)


def cluster_distribution(q: torch.Tensor) -> torch.Tensor:
    # DEC target distribution: sharpen confident assignments and downweight diffuse ones.
    weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
    return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)


def train_bridge_cluster(args):
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    adata = preprocess_adata(read_h5ad_compat(args.data_path), n_hvg=args.n_hvg)
    x = adata.X.toarray().astype(np.float32) if hasattr(adata.X, 'toarray') else np.asarray(adata.X, dtype=np.float32)
    labels, label_key = extract_labels(adata)
    unique_labels = sorted(np.unique(labels).tolist())
    label_to_id = {label: idx for idx, label in enumerate(unique_labels)}
    y = np.array([label_to_id[label] for label in labels], dtype=np.int64)

    bridge = load_bridge(args, input_dim=x.shape[1], device=device)
    bridge.train()
    sampler = BridgeSampler(bridge)
    teacher_target = load_teacher_target(Path(args.output_dir), fallback_shape=(x.shape[0], args.latent_dim))
    teacher_target = teacher_target[: x.shape[0]]

    cluster_head = ClusterHead(input_dim=args.latent_dim, n_clusters=len(unique_labels), hidden_dim=args.hidden_dim).to(device)
    teacher_kmeans = KMeans(n_clusters=len(unique_labels), random_state=42, n_init=20).fit(teacher_target)
    with torch.no_grad():
        teacher_tensor = torch.from_numpy(teacher_target).to(device).float()
        encoded_teacher = cluster_head.encoder(teacher_tensor)
        encoded_centers = []
        for cluster_id in range(len(unique_labels)):
            cluster_points = encoded_teacher[teacher_kmeans.labels_ == cluster_id]
            if cluster_points.shape[0] == 0:
                encoded_centers.append(encoded_teacher.mean(dim=0))
            else:
                encoded_centers.append(cluster_points.mean(dim=0))
        cluster_head.prototypes.copy_(torch.stack(encoded_centers, dim=0))

    dataset = TensorDataset(torch.from_numpy(x), torch.from_numpy(y), torch.from_numpy(teacher_target))
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=False)
    optimizer = torch.optim.Adam(list(bridge.parameters()) + list(cluster_head.parameters()), lr=args.lr)

    for epoch in range(args.epochs):
        total_loss = 0.0
        for batch_x, batch_y, batch_teacher in loader:
            batch_x = batch_x.to(device)
            batch_teacher = batch_teacher.to(device)
            raw_mask = build_support_mask(batch_x, topk=args.support_topk)

            # 1) reverse: raw sparse counts -> shared Gaussian latent.
            latent = sampler.ddim_reverse_sample_loop(image=batch_x, raw_mask=raw_mask)

            # 2) forward: shared latent -> denoised cluster-friendly target embedding.
            output = bridge(batch_x, raw_mask=raw_mask)
            pred_z = output.target_embedding
            support_anchor = output.support_anchor

            q = cluster_head(pred_z)
            p = cluster_distribution(q).detach()
            cluster_loss = F.kl_div(q.clamp_min(1e-8).log(), p, reduction='batchmean')

            # The teacher target is the clean manifold anchor: scVI/DCA/PhytoCluster/PCA+graph smoothing.
            teacher_loss = F.mse_loss(pred_z, batch_teacher)

            # Keep the Gaussian bridge faithful to the source latent while preventing collapse.
            latent_prior_loss = latent.pow(2).mean()
            support_loss = F.mse_loss(support_anchor, batch_teacher)

            # Entropy regularization discourages trivial cluster collapse at the beginning of training.
            entropy = -(q * q.clamp_min(1e-8).log()).sum(dim=1).mean()

            loss = (
                args.teacher_weight * teacher_loss
                + args.cluster_weight * cluster_loss
                + args.gaussian_weight * latent_prior_loss
                + args.support_weight * support_loss
                + args.entropy_weight * entropy
            )

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            total_loss += loss.item() * batch_x.size(0)

        print(f'bridge epoch={epoch + 1} loss={total_loss / len(dataset):.6f}')

    bridge.eval()
    cluster_head.eval()
    with torch.no_grad():
        all_x = torch.from_numpy(x).to(device)
        eval_mask = build_support_mask(all_x, topk=args.support_topk)
        output = bridge(all_x, raw_mask=eval_mask)
        pred_embedding = output.target_embedding.cpu().numpy()
        probs = cluster_head(output.target_embedding).cpu().numpy()
        pred = probs.argmax(axis=1)
        kmeans_pred = KMeans(n_clusters=len(unique_labels), random_state=42, n_init=20).fit_predict(pred_embedding)

    ari_dec = adjusted_rand_score(y, pred)
    nmi_dec = normalized_mutual_info_score(y, pred)
    ari_kmeans = adjusted_rand_score(y, kmeans_pred)
    nmi_kmeans = normalized_mutual_info_score(y, kmeans_pred)
    sil = silhouette_score(pred_embedding, pred) if len(np.unique(pred)) > 1 else float('nan')
    print(f'label={label_key} dec_ari={ari_dec:.4f} dec_nmi={nmi_dec:.4f} kmeans_ari={ari_kmeans:.4f} kmeans_nmi={nmi_kmeans:.4f} sil={sil:.4f}')

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    torch.save({'bridge': bridge.state_dict(), 'cluster_head': cluster_head.state_dict()}, output_dir / 'bridge_cluster.pt')
    np.save(output_dir / 'bridge_embedding.npy', pred_embedding)
    np.save(output_dir / 'bridge_dec_labels.npy', pred)
    np.save(output_dir / 'bridge_kmeans_labels.npy', kmeans_pred)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--data-path', type=str, required=True)
    parser.add_argument('--output-dir', type=str, required=True)
    parser.add_argument('--n-hvg', type=int, default=2000)
    parser.add_argument('--latent-dim', type=int, default=64)
    parser.add_argument('--hidden-dim', type=int, default=256)
    parser.add_argument('--batch-size', type=int, default=256)
    parser.add_argument('--epochs', type=int, default=30)
    parser.add_argument('--lr', type=float, default=1e-3)
    parser.add_argument('--diffusion-steps', type=int, default=50)
    parser.add_argument('--support-topk', type=int, default=256)
    parser.add_argument('--support-blend', type=float, default=0.2)
    parser.add_argument('--teacher-weight', type=float, default=1.0)
    parser.add_argument('--cluster-weight', type=float, default=1.0)
    parser.add_argument('--gaussian-weight', type=float, default=1e-3)
    parser.add_argument('--support-weight', type=float, default=0.5)
    parser.add_argument('--entropy-weight', type=float, default=1e-3)
    parser.add_argument('--zero-weight', type=float, default=0.25)
    args = parser.parse_args()
    train_bridge_cluster(args)
