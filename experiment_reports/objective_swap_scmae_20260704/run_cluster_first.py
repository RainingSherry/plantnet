#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import random
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from torch.utils.data import DataLoader, Dataset

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader() -> None:
    try:
        import h5py
        from anndata._io.specs.registry import IOSpec, _REGISTRY

        def _read_null(*args, **kwargs):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


_register_null_h5ad_reader()


class ExprDataset(Dataset):
    def __init__(self, data: np.ndarray, labels: np.ndarray):
        self.data = torch.as_tensor(data, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.data[idx], self.labels[idx]


class PrototypeEncoder(nn.Module):
    def __init__(
        self,
        input_dim: int,
        latent_dim: int,
        hidden_dim: int,
        n_clusters: int,
        dropout: float,
        temperature: float,
    ):
        super().__init__()
        self.temperature = float(temperature)
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.GELU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dim, latent_dim),
        )
        self.norm = nn.LayerNorm(latent_dim)
        self.prototypes = nn.Parameter(torch.empty(n_clusters, latent_dim))
        nn.init.xavier_uniform_(self.prototypes)

    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.norm(self.encoder(x))

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        z_raw = self.feature(x)
        z = F.normalize(z_raw, dim=1)
        proto = F.normalize(self.prototypes, dim=1)
        scores = z @ proto.t()
        logits = scores / max(self.temperature, 1e-6)
        return {"z_raw": z_raw, "z": z, "scores": scores, "logits": logits}


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--epochs", type=int, default=100)
    parser.add_argument("--batch_size", type=int, default=512)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--hidden_dim", type=int, default=512)
    parser.add_argument("--latent_dim", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.05)
    parser.add_argument("--view_dropout", type=float, default=0.15)
    parser.add_argument("--gaussian_std", type=float, default=0.02)
    parser.add_argument("--temperature", type=float, default=0.2)
    parser.add_argument("--target_temperature", type=float, default=0.1)
    parser.add_argument("--target_power", type=float, default=2.0)
    parser.add_argument("--assignment_mode", default="sharpen", choices=["sharpen", "sinkhorn"])
    parser.add_argument("--sinkhorn_iters", type=int, default=3)
    parser.add_argument("--consistency_weight", type=float, default=1.0)
    parser.add_argument("--entropy_weight", type=float, default=0.05)
    parser.add_argument("--confidence_weight", type=float, default=0.02)
    parser.add_argument("--variance_weight", type=float, default=0.02)
    parser.add_argument("--cov_weight", type=float, default=0.0)
    parser.add_argument("--pca_dim", type=int, default=128)
    parser.add_argument("--skip_eval", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    if gpu in {0, 7}:
        raise ValueError("Physical GPU 0 and GPU 7 are intentionally avoided.")
    return torch.device(f"cuda:{gpu}")


def augment_view(x: torch.Tensor, dropout: float, gaussian_std: float) -> torch.Tensor:
    out = x
    if dropout > 0.0:
        keep = torch.rand_like(out).ge(float(dropout))
        out = out * keep.float()
    if gaussian_std > 0.0:
        out = out + torch.randn_like(out) * float(gaussian_std)
    return out


@torch.no_grad()
def sharpen_targets(logits: torch.Tensor, power: float) -> torch.Tensor:
    q = F.softmax(logits, dim=1)
    q = torch.pow(q, float(power))
    q = q / q.sum(dim=1, keepdim=True).clamp_min(1e-12)
    return q


@torch.no_grad()
def sinkhorn_targets(scores: torch.Tensor, epsilon: float, n_iters: int) -> torch.Tensor:
    q = torch.exp(scores / max(float(epsilon), 1e-6)).t()
    q = q / q.sum().clamp_min(1e-12)
    k, b = q.shape
    for _ in range(max(1, int(n_iters))):
        q = q / q.sum(dim=1, keepdim=True).clamp_min(1e-12)
        q = q / float(k)
        q = q / q.sum(dim=0, keepdim=True).clamp_min(1e-12)
        q = q / float(b)
    q = q * float(b)
    return q.t().contiguous()


def assignment_targets(out: dict[str, torch.Tensor], args) -> torch.Tensor:
    if args.assignment_mode == "sinkhorn":
        return sinkhorn_targets(out["scores"].detach(), args.target_temperature, args.sinkhorn_iters)
    return sharpen_targets(out["logits"].detach(), args.target_power)


def cross_view_loss(logits: torch.Tensor, target: torch.Tensor) -> torch.Tensor:
    return -(target * F.log_softmax(logits, dim=1)).sum(dim=1).mean()


def variance_floor_loss(z: torch.Tensor) -> torch.Tensor:
    std = torch.sqrt(z.var(dim=0, unbiased=False) + 1e-4)
    return F.relu(1.0 - std).mean()


def covariance_loss(z: torch.Tensor) -> torch.Tensor:
    zc = z - z.mean(dim=0, keepdim=True)
    zc = zc / torch.sqrt(zc.var(dim=0, unbiased=False, keepdim=True) + 1e-4)
    cov = (zc.t() @ zc) / max(1, zc.shape[0] - 1)
    off = cov - torch.diag(torch.diag(cov))
    return (off.pow(2).sum() / z.shape[1]).float()


def entropy_terms(logits: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
    probs = F.softmax(logits, dim=1)
    sample_entropy = -(probs * torch.log(probs.clamp_min(1e-12))).sum(dim=1).mean()
    mean_probs = probs.mean(dim=0)
    mean_entropy_loss = (mean_probs * torch.log(mean_probs.clamp_min(1e-12))).sum()
    return sample_entropy, mean_entropy_loss


def effective_dimensionality(std: np.ndarray) -> dict:
    var = np.square(std.astype(np.float64))
    pr = float((var.sum() ** 2) / max(float(np.square(var).sum()), 1e-12))
    return {
        "std_min": float(std.min()),
        "std_median": float(np.median(std)),
        "std_max": float(std.max()),
        "effective_dim_pr": pr,
        "dims_std_gt_0p1": int((std > 0.1).sum()),
        "dims_std_gt_1p0": int((std > 1.0).sum()),
    }


@torch.no_grad()
def extract_all(model: PrototypeEncoder, loader: DataLoader, device: torch.device):
    model.eval()
    emb, logits, labels = [], [], []
    for _, x, y in loader:
        out = model(x.to(device))
        emb.append(out["z_raw"].detach().cpu().numpy())
        logits.append(out["logits"].detach().cpu().numpy())
        labels.append(y.numpy())
    return (
        np.nan_to_num(np.concatenate(emb).astype(np.float32)),
        np.nan_to_num(np.concatenate(logits).astype(np.float32)),
        np.concatenate(labels).astype(np.int64),
    )


def pca_kmeans_baseline(data: np.ndarray, labels: np.ndarray, n_clusters: int, pca_dim: int, seed: int) -> dict:
    dim = max(2, min(int(pca_dim), data.shape[0] - 1, data.shape[1] - 1))
    emb = PCA(n_components=dim, random_state=seed).fit_transform(data.astype(np.float64)).astype(np.float32)
    pred = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(emb)
    metrics, _ = family.compute_kmeans_metrics(labels, pred.astype(np.int64))
    metrics["pca_dim"] = int(dim)
    return metrics


def direct_assignment_metrics(labels: np.ndarray, logits: np.ndarray) -> dict:
    pred = logits.argmax(axis=1).astype(np.int64)
    metrics, _ = family.compute_kmeans_metrics(labels, pred)
    metrics["cluster_method"] = "direct_prototype_argmax"
    return metrics


def main() -> int:
    args = parse_args()
    if args.smoke:
        args.epochs = min(args.epochs, 3)
        args.batch_size = min(args.batch_size, 512)
    set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem

    bundle = family.load_scmae_dataset(
        args.data_path,
        args.input_mode,
        args.n_top_genes,
        args.target_sum,
        args.scale_input,
        args.label_key,
        args.seed,
    )
    data = np.asarray(bundle.data, dtype=np.float32)
    labels = np.asarray(bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))

    dataset = ExprDataset(data, labels)
    generator = torch.Generator().manual_seed(args.seed)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 1024), shuffle=False, drop_last=False)

    pca_metrics = pca_kmeans_baseline(data, labels, n_clusters, args.pca_dim, args.seed)
    save_json({"kmeans_known_k_on_pca": pca_metrics}, str(save_dir / "pca_baseline_metrics.json"))

    model = PrototypeEncoder(
        input_dim=data.shape[1],
        latent_dim=args.latent_dim,
        hidden_dim=args.hidden_dim,
        n_clusters=n_clusters,
        dropout=args.dropout,
        temperature=args.temperature,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {k: [] for k in ["loss", "consistency", "variance", "cov", "sample_entropy", "mean_entropy", "usage_min", "usage_max"]}
    start = time.time()
    print(
        f"Device={device} dataset={dataset_name} cells={data.shape[0]} genes={data.shape[1]} "
        f"clusters={n_clusters} mode={args.assignment_mode} latent={args.latent_dim}"
    )
    print(f"PCA baseline ARI={pca_metrics['ari']:.4f} NMI={pca_metrics['nmi']:.4f} dim={pca_metrics['pca_dim']}")

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {k: 0.0 for k in history}
        batches = 0
        for _, x_cpu, _ in loader:
            x = x_cpu.to(device)
            x1 = augment_view(x, args.view_dropout, args.gaussian_std)
            x2 = augment_view(x, args.view_dropout, args.gaussian_std)
            out1 = model(x1)
            out2 = model(x2)
            t1 = assignment_targets(out1, args)
            t2 = assignment_targets(out2, args)
            consistency = 0.5 * (cross_view_loss(out2["logits"], t1) + cross_view_loss(out1["logits"], t2))
            var_loss = 0.5 * (variance_floor_loss(out1["z_raw"]) + variance_floor_loss(out2["z_raw"]))
            cov_loss = 0.5 * (covariance_loss(out1["z_raw"]) + covariance_loss(out2["z_raw"]))
            sample_entropy1, mean_entropy1 = entropy_terms(out1["logits"])
            sample_entropy2, mean_entropy2 = entropy_terms(out2["logits"])
            sample_entropy = 0.5 * (sample_entropy1 + sample_entropy2)
            mean_entropy = 0.5 * (mean_entropy1 + mean_entropy2)
            loss = (
                float(args.consistency_weight) * consistency
                + float(args.variance_weight) * var_loss
                + float(args.cov_weight) * cov_loss
                + float(args.confidence_weight) * sample_entropy
                + float(args.entropy_weight) * mean_entropy
            )
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            with torch.no_grad():
                usage = F.softmax(out1["logits"], dim=1).mean(dim=0)
            sums["loss"] += float(loss.detach().cpu())
            sums["consistency"] += float(consistency.detach().cpu())
            sums["variance"] += float(var_loss.detach().cpu())
            sums["cov"] += float(cov_loss.detach().cpu())
            sums["sample_entropy"] += float(sample_entropy.detach().cpu())
            sums["mean_entropy"] += float((-mean_entropy).detach().cpu())
            sums["usage_min"] += float(usage.min().detach().cpu())
            sums["usage_max"] += float(usage.max().detach().cpu())
            batches += 1
        for key, value in sums.items():
            history[key].append(value / max(1, batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"cons={history['consistency'][-1]:.4f} var={history['variance'][-1]:.4f} "
                f"Hs={history['sample_entropy'][-1]:.4f} Hmean={history['mean_entropy'][-1]:.4f} "
                f"usage=[{history['usage_min'][-1]:.4f},{history['usage_max'][-1]:.4f}]"
            )

    embedding, logits, labels_out = extract_all(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "logits_final.npy", logits.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    save_json(history, str(save_dir / "training_history.json"))
    std_profile = effective_dimensionality(embedding.std(axis=0))

    fixed = {}
    preds = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            "cluster_first_prototype",
            args.seed,
            embedding,
            labels_out,
            n_clusters,
            {
                "assignment_mode": args.assignment_mode,
                "variance_weight": float(args.variance_weight),
                "latent_dim": int(args.latent_dim),
            },
        )
        fixed = eval_result["fixed"]
        preds = eval_result["preds"]["kmeans_known_k"]
        fixed["direct_prototype_argmax"] = direct_assignment_metrics(labels_out, logits)
        fixed["pca_kmeans_known_k"] = pca_metrics
        save_json(fixed, str(save_dir / "metrics.json"))

    if preds is None:
        pred_for_mass = logits.argmax(axis=1).astype(np.int64)
    else:
        pred_for_mass = preds.astype(np.int64)
    counts = np.bincount(pred_for_mass, minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    summary = {
        "dataset": dataset_name,
        "seed": int(args.seed),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "method": "cluster_first_prototype",
        "assignment_mode": args.assignment_mode,
        "variance_weight": float(args.variance_weight),
        "entropy_weight": float(args.entropy_weight),
        "confidence_weight": float(args.confidence_weight),
        "latent_dim": int(args.latent_dim),
        "fixed_metrics": fixed,
        "std_profile": std_profile,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
    }
    save_json(summary, str(save_dir / "summary.json"))
    k_metrics = fixed.get("kmeans_known_k", {})
    d_metrics = fixed.get("direct_prototype_argmax", {})
    print(
        f"[RESULT] {dataset_name} mode={args.assignment_mode} "
        f"kmeans_ARI={k_metrics.get('ari')} direct_ARI={d_metrics.get('ari')} "
        f"pca_ARI={pca_metrics.get('ari')}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
