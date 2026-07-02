#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path

for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "8")

import numpy as np
import torch
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = next(p for p in [CURRENT_DIR, *CURRENT_DIR.parents] if (p / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

from loss import AdaptiveSwitchLoss, compute_gate
from model import AdaptiveSwitchScMAE
from clusterability import compute_clusterability
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader() -> None:
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec

        def _n(*a, **k):
            return None
        for t in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(t, IOSpec("null", "0.1.0"))(_n)
            except Exception:
                pass
    except Exception:
        pass


_register_null_h5ad_reader()

METHOD_NAME = "AdaptiveSwitch_scMAE"
DISPLAY_NAME = "scMAE + adaptive sharp/soft switch"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class ExprDataset(Dataset):
    def __init__(self, enc, log, labels):
        self.enc = torch.as_tensor(enc, dtype=torch.float32)
        self.log = torch.as_tensor(log, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self):
        return int(self.enc.shape[0])

    def __getitem__(self, i):
        return int(i), self.enc[i], self.log[i], self.labels[i]


def parse_args():
    p = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--data_path", required=True)
    p.add_argument("--save_dir", required=True)
    p.add_argument("--dataset_name", default=None)
    p.add_argument("--label_key", default="auto")
    p.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    p.add_argument("--n_top_genes", type=int, default=1000)
    p.add_argument("--target_sum", type=float, default=10000.0)
    p.add_argument("--scale_input", type=family.str2bool, default=True)
    p.add_argument("--n_clusters", type=int, required=True)
    p.add_argument("--epochs", type=int, default=80)
    p.add_argument("--batch_size", type=int, default=256)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight_decay", type=float, default=1e-5)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--gpu", type=int, default=1)
    p.add_argument("--no_cuda", action="store_true")
    p.add_argument("--skip_eval", action="store_true")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--hidden_size", type=int, default=128)
    p.add_argument("--dropout", type=float, default=0.05)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--masked_data_weight", type=float, default=0.75)
    p.add_argument("--mask_weight", type=float, default=0.65)
    p.add_argument("--cluster_weight", type=float, default=0.35)
    p.add_argument("--consistency_weight", type=float, default=0.05)
    p.add_argument("--variance_weight", type=float, default=0.02)
    p.add_argument("--var_mode", default="hinge", choices=["hinge", "cov", "koleo", "both"], help="anti-collapse mechanism")
    p.add_argument("--entropy_weight", type=float, default=0.10)
    p.add_argument("--confidence_threshold", type=float, default=0.35)
    p.add_argument("--warmup_epochs", type=int, default=20)
    p.add_argument("--target_update_interval", type=int, default=5)
    p.add_argument("--neighbor_k", type=int, default=15)
    p.add_argument("--knn_pca_dim", type=int, default=50)
    p.add_argument("--gate_kappa", type=float, default=0.15, help="gate=1/(1+(kl_ref/kappa)^2)")
    p.add_argument("--force_gate", type=float, default=-1.0, help="-1=auto; else fix gate in [0,1] for ablation")
    p.add_argument("--gate_ema", type=float, default=0.5, help="EMA smoothing of gate across refreshes")
    return p.parse_args()


# ===HELPERS===
def build_knn_graph(data_np, k, pca_dim, seed):
    n = int(data_np.shape[0]); max_k = min(int(k), max(1, n - 1))
    if max_k <= 0:
        return None
    dim = min(int(pca_dim), min(data_np.shape) - 1)
    emb = PCA(n_components=max(2, dim), random_state=seed).fit_transform(data_np.astype(np.float64)) if dim >= 2 else data_np.astype(np.float64)
    emb = emb / np.clip(np.linalg.norm(emb, axis=1, keepdims=True), 1e-12, None)
    nn = NearestNeighbors(n_neighbors=max_k + 1, metric="cosine").fit(emb)
    return nn.kneighbors(emb, return_distance=False)[:, 1:max_k + 1].astype(np.int64)


@torch.no_grad()
def extract_all(model, loader, device):
    model.eval()
    emb, q_all, labels = [], [], []
    for _, x, _, y in loader:
        out = model(x.to(device))
        emb.append(out["latent"].detach().cpu().numpy())
        q_all.append(out["cluster_q"].detach().cpu().numpy())
        labels.append(y.numpy())
    return (np.nan_to_num(np.concatenate(emb).astype(np.float32)),
            np.concatenate(q_all).astype(np.float32),
            np.concatenate(labels).astype(np.int64))


def metric(res, name):
    if not res:
        return None
    v = res.get("fixed", {}).get("kmeans_known_k", {}).get(name)
    return None if v is None else float(v)


def diagnostics(embedding, labels, n_clusters, seed, preds, cdiag, gate, kl_ref):
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    nn = NearestNeighbors(n_neighbors=min(11, embedding.shape[0])).fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    d = {
        "neighbor_purity_proxy": float(np.mean(labels[idx] == labels[:, None])),
        "mixed_cell_fraction": 0.0,
        "embedding_variance": var,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
        "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or (frac.min() if frac.size else 0) < 0.001 or (frac.max() if frac.size else 1) > 0.95),
        "final_gate": float(gate),
        "final_kl_ref": float(kl_ref),
    }
    d.update(cdiag)
    return d


# ===MAIN===
def main() -> int:
    args = parse_args()
    if args.gpu in {0, 7} and not args.no_cuda:
        raise ValueError("GPU 0 and GPU 7 are forbidden.")
    if args.smoke:
        args.epochs = min(args.epochs, 3); args.warmup_epochs = min(args.warmup_epochs, 1)
    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem
    stage = "smoke" if args.smoke else "screen"

    target_bundle = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed)
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed)
        encoder_data = encoder_bundle.data
    else:
        encoder_data = target_bundle.data
    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))

    print("Building PCA-KNN graph ...")
    nb_indices = build_knn_graph(encoder_data, args.neighbor_k, args.knn_pca_dim, args.seed)

    dataset = ExprDataset(encoder_data, log_expr, labels)
    gen = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=gen)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = AdaptiveSwitchScMAE(encoder_data.shape[1], n_clusters, args.hidden_size, args.dropout).to(device)
    criterion = AdaptiveSwitchLoss(args.masked_data_weight, args.mask_weight, args.cluster_weight,
                                   args.consistency_weight, args.variance_weight, args.entropy_weight, args.confidence_threshold, args.var_mode)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    p_targets = None
    clusterab = np.ones(encoder_data.shape[0], dtype=np.float32)
    cdiag = {"clusterability_mean": 1.0, "core_fraction": 1.0}
    gate = 1.0 if args.force_gate < 0 else float(args.force_gate)
    kl_ref = 0.0
    centers_initialized = False
    hist = {k: [] for k in ["loss", "scmae_loss", "sharp_loss", "soft_loss", "variance_loss", "gate", "kl_ref"]}
    hist["stage"] = stage
    start = time.time()
    print(f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} clusters={n_clusters} force_gate={args.force_gate} kappa={args.gate_kappa}")

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch > args.warmup_epochs and ((epoch - args.warmup_epochs - 1) % max(1, args.target_update_interval) == 0 or p_targets is None):
            emb, q_full, _ = extract_all(model, full_loader, device)
            if not centers_initialized:
                km = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit(emb)
                model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=device))
                emb, q_full, _ = extract_all(model, full_loader, device)
                centers_initialized = True
            sharp_p = AdaptiveSwitchScMAE.sharpen(torch.as_tensor(q_full)).numpy().astype(np.float32)
            p_targets = sharp_p
            clusterab, cdiag = compute_clusterability(emb, q_full, nb_indices, k=args.neighbor_k)
            g_new, kl_ref = compute_gate(sharp_p, q_full, args.gate_kappa)
            if args.force_gate < 0:
                gate = args.gate_ema * gate + (1.0 - args.gate_ema) * g_new if centers_initialized else g_new
            else:
                gate = float(args.force_gate)

        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        model.train()
        sums = {k: 0.0 for k in ["loss", "scmae_loss", "sharp_loss", "soft_loss", "variance_loss"]}
        nb = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            idx_np = idx.numpy()
            x = x_cpu.to(device); target = log_cpu.to(device)
            c_batch = torch.as_tensor(clusterab[idx_np], dtype=torch.float32, device=device)
            strong, mask = model.random_mask(x, args.mask_prob)
            weak, _ = model.random_mask(x, max(0.05, args.mask_prob * 0.5))
            out = model(strong); weak_out = model(weak)
            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)
            loss, parts = criterion(out, weak_out, target, mask, p_batch, c_batch, cluster_scale, gate)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss epoch {epoch}: {parts}")
            optimizer.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0); optimizer.step()
            for k in sums:
                sums[k] += parts[k]
            nb += 1
        for k in sums:
            hist[k].append(sums[k] / max(1, nb))
        hist["gate"].append(float(gate)); hist["kl_ref"].append(float(kl_ref))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={hist['loss'][-1]:.4f} scmae={hist['scmae_loss'][-1]:.4f} sharp={hist['sharp_loss'][-1]:.4f} soft={hist['soft_loss'][-1]:.4f} gate={gate:.3f} kl_ref={kl_ref:.4f} scale={cluster_scale:.2f}")

    embedding, q_out, labels_out = extract_all(model, full_loader, device)
    clusterab, cdiag = compute_clusterability(embedding, q_out, nb_indices, k=args.neighbor_k)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    save_json(hist, str(save_dir / "training_history.json"))

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    diag = diagnostics(embedding, labels_out, n_clusters, args.seed, preds, cdiag, gate, kl_ref)
    save_json(diag, str(save_dir / "diagnostics.json"))
    baseline = BASELINES.get(dataset_name, {})
    nmi, ari, acc = metric(eval_result, "nmi"), metric(eval_result, "ari"), metric(eval_result, "acc")
    meets = bool((nmi is not None and nmi >= baseline.get("nmi", np.inf)) or (ari is not None and ari >= baseline.get("ari", np.inf)))
    save_json({"dataset": dataset_name, "method": DISPLAY_NAME, "stage": stage, "seed": int(args.seed), "n_clusters": int(n_clusters), "runtime_seconds": float(time.time() - start), "fixed_metrics": eval_result["fixed"] if eval_result else {}, "diagnostics": diag, "baseline": baseline, "meets_screen_baseline_any": meets, "final_gate": float(gate), "final_kl_ref": float(kl_ref)}, str(save_dir / "summary.json"))
    print(f"[RESULT] {dataset_name} NMI={nmi} ARI={ari} meets={meets} collapse={diag['collapse_warning']} gate={gate:.3f} kl_ref={kl_ref:.4f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

