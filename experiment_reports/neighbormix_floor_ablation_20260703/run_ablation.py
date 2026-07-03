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
import torch.nn.functional as F
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import normalize
from torch.utils.data import DataLoader, Dataset

ROOT = Path(__file__).resolve().parents[2]
ADAPTIVE_SWITCH_DIR = ROOT / "experimental_retired_models" / "Granularity_scMAE_experiments" / "AdaptiveSwitch_scMAE"
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(ADAPTIVE_SWITCH_DIR) not in sys.path:
    sys.path.insert(0, str(ADAPTIVE_SWITCH_DIR))

from clusterability import compute_clusterability
from loss import AdaptiveSwitchLoss, compute_gate
from model import AdaptiveSwitchScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader() -> None:
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec

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
    def __init__(self, enc: np.ndarray, log: np.ndarray, labels: np.ndarray):
        self.enc = torch.as_tensor(enc, dtype=torch.float32)
        self.log = torch.as_tensor(log, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.enc.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.enc[idx], self.log[idx], self.labels[idx]


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
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--skip_eval", action="store_true")
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.05)
    parser.add_argument("--mask_prob", type=float, default=0.4)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_weight", type=float, default=0.65)
    parser.add_argument("--cluster_weight", type=float, default=0.35)
    parser.add_argument("--consistency_weight", type=float, default=0.05)
    parser.add_argument("--variance_weight", type=float, default=0.0)
    parser.add_argument("--var_mode", default="hinge", choices=["hinge", "cov", "koleo", "both"])
    parser.add_argument("--entropy_weight", type=float, default=0.10)
    parser.add_argument("--confidence_threshold", type=float, default=0.35)
    parser.add_argument("--warmup_epochs", type=int, default=20)
    parser.add_argument("--target_update_interval", type=int, default=5)
    parser.add_argument("--neighbor_k", type=int, default=15)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--gate_kappa", type=float, default=0.15)
    parser.add_argument("--force_gate", type=float, default=1.0)
    parser.add_argument("--gate_ema", type=float, default=0.5)
    parser.add_argument("--mix_mode", default="none", choices=["none", "neighbor"])
    parser.add_argument("--mix_alpha", type=float, default=0.8)
    parser.add_argument("--mix_neighbors", type=int, default=4)
    parser.add_argument("--mix_tau", type=float, default=0.2)
    parser.add_argument("--denoise_weight", type=float, default=1.0)
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


def build_neighbor_distribution(data: np.ndarray, k: int, pca_dim: int, tau: float, seed: int):
    n_cells, n_genes = data.shape
    k = max(1, min(int(k), n_cells - 1))
    dim = max(2, min(int(pca_dim), n_cells - 1, n_genes - 1))
    emb = PCA(n_components=dim, random_state=seed).fit_transform(data.astype(np.float64))
    emb = normalize(emb, norm="l2", axis=1)
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(emb)
    distances, neighbors = nn.kneighbors(emb, return_distance=True)
    neighbors = neighbors[:, 1:].astype(np.int64)
    distances = distances[:, 1:].astype(np.float32)
    sim = np.clip(1.0 - distances, a_min=0.0, a_max=None)
    logits = sim / max(float(tau), 1e-6)
    logits -= logits.max(axis=1, keepdims=True)
    probs = np.exp(logits).astype(np.float64)
    probs /= probs.sum(axis=1, keepdims=True).clip(min=1e-12)
    profile = {
        "neighbor_k": int(k),
        "pca_dim": int(dim),
        "tau": float(tau),
        "mean_neighbor_similarity": float(sim.mean()),
        "mean_max_neighbor_prob": float(probs.max(axis=1).mean()),
    }
    return neighbors, probs.astype(np.float32), profile


def sample_neighbor_mix(
    data_np: np.ndarray,
    batch_indices: np.ndarray,
    batch_x: torch.Tensor,
    neighbor_indices: np.ndarray,
    neighbor_probs: np.ndarray,
    alpha: float,
    mix_neighbors: int,
    rng: np.random.Generator,
) -> torch.Tensor:
    batch_size = int(batch_indices.shape[0])
    mix_neighbors = max(1, int(mix_neighbors))
    sampled = np.empty((batch_size, mix_neighbors), dtype=np.int64)
    weights = np.empty((batch_size, mix_neighbors), dtype=np.float32)
    for row, cell in enumerate(batch_indices):
        probs = neighbor_probs[cell]
        choices = rng.choice(neighbor_indices.shape[1], size=mix_neighbors, replace=True, p=probs)
        sampled[row] = neighbor_indices[cell, choices]
        picked = probs[choices].astype(np.float32, copy=False)
        weights[row] = picked / max(float(picked.sum()), 1e-12)
    neighbor_expr = data_np[sampled]
    neighbor_mean = np.sum(neighbor_expr * weights[:, :, None], axis=1).astype(np.float32)
    neighbor_t = torch.as_tensor(neighbor_mean, dtype=batch_x.dtype, device=batch_x.device)
    alpha = float(alpha)
    return alpha * batch_x + (1.0 - alpha) * neighbor_t


def scmae_reconstruction_part(out, target, mask, masked_data_weight: float, mask_weight: float):
    w = mask * float(masked_data_weight) + (1.0 - mask) * (1.0 - float(masked_data_weight))
    rec = (w * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
    mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
    return (1.0 - float(mask_weight)) * rec + float(mask_weight) * mask_loss


@torch.no_grad()
def extract_all(model, loader, device):
    model.eval()
    emb, q_all, labels = [], [], []
    for _, x, _, y in loader:
        out = model(x.to(device))
        emb.append(out["latent"].detach().cpu().numpy())
        q_all.append(out["cluster_q"].detach().cpu().numpy())
        labels.append(y.numpy())
    return (
        np.nan_to_num(np.concatenate(emb).astype(np.float32)),
        np.concatenate(q_all).astype(np.float32),
        np.concatenate(labels).astype(np.int64),
    )


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


def main() -> int:
    args = parse_args()
    if args.smoke:
        args.epochs = min(args.epochs, 3)
        args.warmup_epochs = min(args.warmup_epochs, 1)
    set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem

    target_bundle = family.load_scmae_dataset(
        args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed
    )
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(
            args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed
        )
        encoder_data = np.asarray(encoder_bundle.data, dtype=np.float32)
    else:
        encoder_data = np.asarray(target_bundle.data, dtype=np.float32)
    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))

    neighbor_indices = None
    neighbor_probs = None
    neighbor_profile = {"enabled": False}
    nb_indices_for_clusterability = None
    if args.mix_mode == "neighbor":
        neighbor_indices, neighbor_probs, neighbor_profile = build_neighbor_distribution(
            encoder_data, args.neighbor_k, args.knn_pca_dim, args.mix_tau, args.seed
        )
        neighbor_profile["enabled"] = True
        nb_indices_for_clusterability = neighbor_indices
    else:
        nb_indices_for_clusterability = build_neighbor_distribution(
            encoder_data, args.neighbor_k, args.knn_pca_dim, args.mix_tau, args.seed
        )[0]
    save_json(neighbor_profile, str(save_dir / "neighbor_profile.json"))

    dataset = ExprDataset(encoder_data, log_expr, labels)
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = AdaptiveSwitchScMAE(encoder_data.shape[1], n_clusters, args.hidden_size, args.dropout).to(device)
    criterion = AdaptiveSwitchLoss(
        args.masked_data_weight,
        args.mask_weight,
        args.cluster_weight,
        args.consistency_weight,
        args.variance_weight,
        args.entropy_weight,
        args.confidence_threshold,
        args.var_mode,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    rng = np.random.default_rng(args.seed + 20260703)

    p_targets = None
    clusterab = np.ones(encoder_data.shape[0], dtype=np.float32)
    gate = float(args.force_gate)
    kl_ref = 0.0
    centers_initialized = False
    history = {k: [] for k in ["loss", "base_loss", "denoise_loss", "sharp_loss", "variance_loss", "gate", "kl_ref"]}
    start = time.time()
    print(
        f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} "
        f"clusters={n_clusters} mix={args.mix_mode} varw={args.variance_weight}"
    )

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch > args.warmup_epochs and (
            (epoch - args.warmup_epochs - 1) % max(1, args.target_update_interval) == 0 or p_targets is None
        ):
            emb, q_full, _ = extract_all(model, full_loader, device)
            if not centers_initialized:
                km = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit(emb)
                model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=device))
                emb, q_full, _ = extract_all(model, full_loader, device)
                centers_initialized = True
            sharp_p = AdaptiveSwitchScMAE.sharpen(torch.as_tensor(q_full)).numpy().astype(np.float32)
            p_targets = sharp_p
            clusterab, _ = compute_clusterability(emb, q_full, nb_indices_for_clusterability, k=args.neighbor_k)
            g_new, kl_ref = compute_gate(sharp_p, q_full, args.gate_kappa)
            gate = args.gate_ema * gate + (1.0 - args.gate_ema) * g_new if args.force_gate < 0 else float(args.force_gate)

        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        model.train()
        sums = {k: 0.0 for k in ["loss", "base_loss", "denoise_loss", "sharp_loss", "variance_loss"]}
        batches = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            idx_np = idx.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            c_batch = torch.as_tensor(clusterab[idx_np], dtype=torch.float32, device=device)
            strong, mask = model.random_mask(x, args.mask_prob)
            weak, _ = model.random_mask(x, max(0.05, args.mask_prob * 0.5))
            out = model(strong)
            weak_out = model(weak)
            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)
            base_loss, parts = criterion(out, weak_out, target, mask, p_batch, c_batch, cluster_scale, gate)

            denoise_loss = torch.zeros((), dtype=torch.float32, device=device)
            if args.mix_mode == "neighbor" and args.denoise_weight > 0.0:
                mixed = sample_neighbor_mix(
                    encoder_data,
                    idx_np,
                    x,
                    neighbor_indices,
                    neighbor_probs,
                    args.mix_alpha,
                    args.mix_neighbors,
                    rng,
                )
                mixed_corrupt, mixed_mask = model.random_mask(mixed, args.mask_prob)
                mixed_out = model(mixed_corrupt)
                denoise_loss = scmae_reconstruction_part(
                    mixed_out, target, mixed_mask, args.masked_data_weight, args.mask_weight
                )
            loss = base_loss + float(args.denoise_weight) * denoise_loss
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            sums["loss"] += float(loss.detach().cpu())
            sums["base_loss"] += float(base_loss.detach().cpu())
            sums["denoise_loss"] += float(denoise_loss.detach().cpu())
            sums["sharp_loss"] += float(parts["sharp_loss"])
            sums["variance_loss"] += float(parts["variance_loss"])
            batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, batches))
        history["gate"].append(float(gate))
        history["kl_ref"].append(float(kl_ref))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"base={history['base_loss'][-1]:.4f} denoise={history['denoise_loss'][-1]:.4f} "
                f"sharp={history['sharp_loss'][-1]:.4f} var={history['variance_loss'][-1]:.4f} "
                f"gate={gate:.3f} kl_ref={kl_ref:.4f}"
            )

    embedding, q_out, labels_out = extract_all(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    save_json(history, str(save_dir / "training_history.json"))
    std_profile = effective_dimensionality(embedding.std(axis=0))

    eval_result = None
    preds = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            "NeighborMix x std-floor ablation",
            args.seed,
            embedding,
            labels_out,
            n_clusters,
            {
                "mix_mode": args.mix_mode,
                "variance_weight": float(args.variance_weight),
                "var_mode": args.var_mode,
                "denoise_weight": float(args.denoise_weight),
            },
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64) if preds is not None else np.zeros(n_clusters)
    frac = counts / max(1.0, counts.sum())
    summary = {
        "dataset": dataset_name,
        "seed": int(args.seed),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "mix_mode": args.mix_mode,
        "variance_weight": float(args.variance_weight),
        "denoise_weight": float(args.denoise_weight),
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "std_profile": std_profile,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
        "final_gate": float(gate),
        "final_kl_ref": float(kl_ref),
    }
    save_json(summary, str(save_dir / "summary.json"))
    ari = summary["fixed_metrics"].get("kmeans_known_k", {}).get("ari")
    nmi = summary["fixed_metrics"].get("kmeans_known_k", {}).get("nmi")
    print(f"[RESULT] {dataset_name} mix={args.mix_mode} varw={args.variance_weight} ARI={ari} NMI={nmi}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
