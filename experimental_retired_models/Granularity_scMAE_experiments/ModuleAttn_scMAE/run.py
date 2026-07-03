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

from loss import ModuleAttnLoss
from model import ModuleAttnScMAE
from gene_modules import compute_gene_modules
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

METHOD_NAME = "ModuleAttn_scMAE"
DISPLAY_NAME = "scMAE + gene-module attention + DEC + variance floor"
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
    p.add_argument("--variance_weight", type=float, default=0.02)
    p.add_argument("--confidence_threshold", type=float, default=0.35)
    p.add_argument("--warmup_epochs", type=int, default=20)
    p.add_argument("--target_update_interval", type=int, default=5)
    # module-attn specific
    p.add_argument("--n_modules", type=int, default=50)
    p.add_argument("--token_dim", type=int, default=16)
    p.add_argument("--n_heads", type=int, default=4)
    p.add_argument("--n_layers", type=int, default=2)
    p.add_argument("--use_attn", type=family.str2bool, default=True, help="False = module tokens without attention (ablation)")
    return p.parse_args()


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

    print("Computing gene co-expression modules ...")
    module_of = compute_gene_modules(encoder_data, args.n_modules, seed=args.seed)
    n_modules_eff = int(module_of.max() + 1)
    np.save(save_dir / "gene_modules.npy", module_of)
    print(f"  {n_modules_eff} modules over {encoder_data.shape[1]} genes; sizes min/med/max = "
          f"{np.bincount(module_of).min()}/{int(np.median(np.bincount(module_of)))}/{np.bincount(module_of).max()}")

    dataset = ExprDataset(encoder_data, log_expr, labels)
    gen = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=gen)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = ModuleAttnScMAE(encoder_data.shape[1], n_clusters, torch.as_tensor(module_of),
                            n_modules_eff, args.token_dim, args.hidden_size,
                            args.n_heads, args.n_layers, args.dropout, args.use_attn).to(device)
    criterion = ModuleAttnLoss(args.masked_data_weight, args.mask_weight, args.cluster_weight,
                               args.variance_weight, args.confidence_threshold)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    p_targets = None
    centers_initialized = False
    hist = {k: [] for k in ["loss", "scmae_loss", "cluster_loss", "variance_loss"]}
    hist["stage"] = stage
    start = time.time()
    print(f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} clusters={n_clusters} modules={n_modules_eff} use_attn={args.use_attn}")

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch > args.warmup_epochs and ((epoch - args.warmup_epochs - 1) % max(1, args.target_update_interval) == 0 or p_targets is None):
            emb, q_full, _ = extract_all(model, full_loader, device)
            if not centers_initialized:
                km = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed).fit(emb)
                model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=device))
                emb, q_full, _ = extract_all(model, full_loader, device)
                centers_initialized = True
            p_targets = ModuleAttnScMAE.sharpen(torch.as_tensor(q_full)).numpy().astype(np.float32)

        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        model.train()
        sums = {k: 0.0 for k in ["loss", "scmae_loss", "cluster_loss", "variance_loss"]}
        nb = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            idx_np = idx.numpy()
            x = x_cpu.to(device); target = log_cpu.to(device)
            strong, mask = model.random_mask(x, args.mask_prob)
            out = model(strong)
            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)
            loss, parts = criterion(out, target, mask, p_batch, cluster_scale)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss epoch {epoch}: {parts}")
            optimizer.zero_grad(); loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0); optimizer.step()
            for k in sums:
                sums[k] += parts[k]
            nb += 1
        for k in sums:
            hist[k].append(sums[k] / max(1, nb))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={hist['loss'][-1]:.4f} scmae={hist['scmae_loss'][-1]:.4f} cluster={hist['cluster_loss'][-1]:.4f} var={hist['variance_loss'][-1]:.4f} scale={cluster_scale:.2f}")

    embedding, q_out, labels_out = extract_all(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    save_json(hist, str(save_dir / "training_history.json"))

    eval_result = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage})
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    baseline = BASELINES.get(dataset_name, {})
    nmi, ari = metric(eval_result, "nmi"), metric(eval_result, "ari")
    meets = bool((nmi is not None and nmi >= baseline.get("nmi", np.inf)) or (ari is not None and ari >= baseline.get("ari", np.inf)))
    save_json({"dataset": dataset_name, "method": DISPLAY_NAME, "stage": stage, "seed": int(args.seed), "n_clusters": int(n_clusters), "n_modules": n_modules_eff, "use_attn": bool(args.use_attn), "runtime_seconds": float(time.time() - start), "fixed_metrics": eval_result["fixed"] if eval_result else {}, "baseline": baseline, "meets_screen_baseline_any": meets}, str(save_dir / "summary.json"))
    print(f"[RESULT] {dataset_name} NMI={nmi} ARI={ari} meets={meets} modules={n_modules_eff} use_attn={args.use_attn}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
