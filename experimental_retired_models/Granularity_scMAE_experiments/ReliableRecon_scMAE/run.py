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

from loss import ReliableReconLoss
from model import ReliableReconScMAE
from reliability_recon import compute_local_reliability
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


def _register_null_h5ad_reader():
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

METHOD_NAME = "ReliableRecon_scMAE"
DISPLAY_NAME = "scMAE + local-reliability precision-weighted reconstruction"
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
    p.add_argument("--reliability_lambda", type=float, default=1.0, help="0=vanilla scMAE, 1=full precision weighting")
    p.add_argument("--reliability_k", type=int, default=15)
    p.add_argument("--reliability_pca_dim", type=int, default=50)
    p.add_argument("--reliability_floor", type=float, default=0.2)
    return p.parse_args()


@torch.no_grad()
def extract_all(model, loader, device):
    model.eval()
    emb, labels = [], []
    for _, x, _, y in loader:
        emb.append(model.feature(x.to(device)).detach().cpu().numpy())
        labels.append(y.numpy())
    return np.nan_to_num(np.concatenate(emb).astype(np.float32)), np.concatenate(labels).astype(np.int64)


def metric(res, name):
    if not res:
        return None
    v = res.get("fixed", {}).get("kmeans_known_k", {}).get(name)
    return None if v is None else float(v)


def main() -> int:
    args = parse_args()
    if args.gpu in {0, 7} and not args.no_cuda:
        raise ValueError("GPU 0 and GPU 7 are forbidden.")
    if args.smoke:
        args.epochs = min(args.epochs, 3)
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

    reliability, rdiag = None, {"reliability_mean": 1.0}
    if args.reliability_lambda > 0:
        print("Estimating local per-gene reliability from decoupled raw-KNN graph ...")
        reliability, rdiag = compute_local_reliability(
            log_expr, encoder_data, k=args.reliability_k,
            pca_dim=args.reliability_pca_dim, floor=args.reliability_floor, seed=args.seed)
        np.save(save_dir / "reliability_field.npy", reliability.astype(np.float32))

    dataset = ExprDataset(encoder_data, log_expr, labels)
    gen = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=gen)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = ReliableReconScMAE(encoder_data.shape[1], args.hidden_size, args.dropout).to(device)
    criterion = ReliableReconLoss(args.masked_data_weight, args.mask_weight, args.reliability_lambda)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    rel_t = None if reliability is None else torch.as_tensor(reliability, dtype=torch.float32)

    hist = {k: [] for k in ["loss", "reconstruction_loss", "mask_loss"]}
    hist["stage"] = stage
    start = time.time()
    print(f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters} lambda={args.reliability_lambda}")

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {k: 0.0 for k in ["loss", "reconstruction_loss", "mask_loss"]}
        nb = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            r_batch = None if rel_t is None else rel_t[idx.numpy()].to(device)
            strong, mask = model.random_mask(x, args.mask_prob)
            out = model(strong)
            loss, parts = criterion(out, target, mask, r_batch)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss epoch {epoch}: {parts}")
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            for k in sums:
                sums[k] += parts[k]
            nb += 1
        for k in sums:
            hist[k].append(sums[k] / max(1, nb))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={hist['loss'][-1]:.4f} rec={hist['reconstruction_loss'][-1]:.4f} mask={hist['mask_loss'][-1]:.4f}")

    embedding, labels_out = extract_all(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    save_json(hist, str(save_dir / "training_history.json"))

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage})
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    baseline = BASELINES.get(dataset_name, {})
    nmi, ari = metric(eval_result, "nmi"), metric(eval_result, "ari")
    meets = bool((nmi is not None and nmi >= baseline.get("nmi", np.inf)) or (ari is not None and ari >= baseline.get("ari", np.inf)))
    save_json({"dataset": dataset_name, "method": DISPLAY_NAME, "stage": stage, "seed": int(args.seed), "n_clusters": int(n_clusters), "runtime_seconds": float(time.time() - start), "fixed_metrics": eval_result["fixed"] if eval_result else {}, "baseline": baseline, "meets_screen_baseline_any": meets, "reliability_lambda": float(args.reliability_lambda)}, str(save_dir / "summary.json"))
    print(f"[RESULT] {dataset_name} NMI={nmi} ARI={ari} meets={meets} rel_mean={rdiag.get('reliability_mean'):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
