#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VarFloor-scMAE
==============

Formal benchmark runner for zero-mask scMAE + DEC + fixed per-dimension
standard-deviation floor.

This method was promoted from the Granularity scMAE diagnostic experiments after
it emerged as the most reliable deep intervention: DEC sharpening provides the
cluster objective, while a VICReg-style std-floor prevents latent dimensional
variance collapse. Labels are used only for final benchmark evaluation.
"""

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
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
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
    def __init__(self, enc: np.ndarray, target: np.ndarray, labels: np.ndarray):
        self.enc = torch.as_tensor(enc, dtype=torch.float32)
        self.target = torch.as_tensor(target, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.enc.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.enc[idx], self.target[idx], self.labels[idx]


class ScMAEDECStdFloor(nn.Module):
    def __init__(self, num_genes: int, n_clusters: int, hidden_size: int = 128, dropout: float = 0.05):
        super().__init__()
        self.num_genes = int(num_genes)
        self.n_clusters = int(n_clusters)
        self.hidden_size = int(hidden_size)
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(self.num_genes, 256),
            nn.LayerNorm(256),
            nn.Mish(),
            nn.Linear(256, hidden_size),
            nn.LayerNorm(hidden_size),
            nn.Mish(),
            nn.Linear(hidden_size, hidden_size),
        )
        self.mask_predictor = nn.Linear(hidden_size, self.num_genes)
        self.decoder = nn.Linear(hidden_size + self.num_genes, self.num_genes)
        self.cluster_centers = nn.Parameter(torch.randn(self.n_clusters, hidden_size) * 0.02)

    def forward(self, x: torch.Tensor) -> dict[str, torch.Tensor]:
        latent = self.encoder(x)
        mask_logits = self.mask_predictor(latent)
        reconstruction = self.decoder(torch.cat([latent, mask_logits], dim=1))
        q = self.student_q(latent)
        return {"latent": latent, "mask_logits": mask_logits, "reconstruction": reconstruction, "cluster_q": q}

    def student_q(self, latent: torch.Tensor) -> torch.Tensor:
        dist = torch.cdist(latent, self.cluster_centers).pow(2)
        q = 1.0 / (1.0 + dist)
        return q / q.sum(dim=1, keepdim=True).clamp_min(1e-8)

    @torch.no_grad()
    def initialize_centers(self, centers: torch.Tensor) -> None:
        if centers.shape != self.cluster_centers.shape:
            raise ValueError(f"center shape {tuple(centers.shape)} != {tuple(self.cluster_centers.shape)}")
        self.cluster_centers.copy_(centers)

    @staticmethod
    def sharpen(q: torch.Tensor) -> torch.Tensor:
        weight = q.pow(2) / q.sum(dim=0, keepdim=True).clamp_min(1e-8)
        return weight / weight.sum(dim=1, keepdim=True).clamp_min(1e-8)

    def random_mask(self, x: torch.Tensor, mask_prob: float) -> tuple[torch.Tensor, torch.Tensor]:
        mask = (torch.rand_like(x) < float(mask_prob)).float()
        empty = mask.sum(dim=1) == 0
        if bool(empty.any()):
            cols = torch.randint(0, x.shape[1], (int(empty.sum()),), device=x.device)
            mask[empty, cols] = 1.0
        return x.masked_fill(mask.bool(), 0.0), mask


def parse_args():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default="VarFloor-scMAE")
    parser.add_argument("--variant_name", default="varfloor_scmae")
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--warmup_epochs", type=int, default=20)
    parser.add_argument("--target_update_interval", type=int, default=5)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.05)
    parser.add_argument("--mask_prob", type=float, default=0.4)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_weight", type=float, default=0.65)
    parser.add_argument("--cluster_weight", type=float, default=0.35)
    parser.add_argument("--consistency_weight", type=float, default=0.05)
    parser.add_argument("--variance_weight", type=float, default=0.02)
    parser.add_argument("--confidence_threshold", type=float, default=0.35)
    parser.add_argument("--skip_eval", type=family.str2bool, default=False)
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


def std_floor_loss(z: torch.Tensor) -> torch.Tensor:
    std = torch.sqrt(z.var(dim=0) + 1e-4)
    return F.relu(1.0 - std).mean()


def scmae_loss(out: dict[str, torch.Tensor], target: torch.Tensor, mask: torch.Tensor, masked_data_weight: float, mask_weight: float):
    w = mask * float(masked_data_weight) + (1.0 - mask) * (1.0 - float(masked_data_weight))
    rec = (w * F.smooth_l1_loss(out["reconstruction"], target, reduction="none")).mean()
    mask_loss = F.binary_cross_entropy_with_logits(out["mask_logits"], mask.float())
    return (1.0 - float(mask_weight)) * rec + float(mask_weight) * mask_loss


def dec_loss(out: dict[str, torch.Tensor], p_target: torch.Tensor | None, confidence_threshold: float) -> torch.Tensor:
    if p_target is None:
        return out["latent"].new_tensor(0.0)
    q = out["cluster_q"].clamp_min(1e-8)
    kl = (p_target * (torch.log(p_target.clamp_min(1e-8)) - torch.log(q))).sum(dim=1)
    conf = p_target.max(dim=1).values
    keep = (conf >= float(confidence_threshold)).float()
    return (kl * keep).sum() / keep.sum().clamp_min(1.0)


@torch.no_grad()
def extract_all(model: ScMAEDECStdFloor, loader: DataLoader, device: torch.device):
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
    save_json(target_bundle.preprocess_config, str(save_dir / "preprocess_config.json"))

    dataset = ExprDataset(encoder_data, log_expr, labels)
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = ScMAEDECStdFloor(encoder_data.shape[1], n_clusters, args.hidden_size, args.dropout).to(device)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    p_targets = None
    centers_initialized = False
    history = {k: [] for k in ["loss", "scmae_loss", "dec_loss", "variance_loss", "consistency_loss"]}
    start = time.time()
    print(f"Device={device} dataset={dataset_name} cells={encoder_data.shape[0]} clusters={n_clusters}")

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
            p_targets = ScMAEDECStdFloor.sharpen(torch.as_tensor(q_full)).numpy().astype(np.float32)

        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        sums = {k: 0.0 for k in history}
        batches = 0
        model.train()
        for idx, x_cpu, target_cpu, _ in train_loader:
            idx_np = idx.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)
            target = target_cpu.to(device)
            strong, mask = model.random_mask(x, args.mask_prob)
            weak, _ = model.random_mask(x, max(0.05, args.mask_prob * 0.5))
            out = model(strong)
            weak_out = model(weak)
            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)

            sc_loss = scmae_loss(out, target, mask, args.masked_data_weight, args.mask_weight)
            consistency = (F.normalize(out["latent"], dim=1) - F.normalize(weak_out["latent"], dim=1)).pow(2).sum(dim=1).mean()
            d_loss = dec_loss(out, p_batch, args.confidence_threshold)
            v_loss = std_floor_loss(out["latent"])
            loss = sc_loss + args.consistency_weight * consistency + cluster_scale * args.cluster_weight * d_loss + args.variance_weight * v_loss
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()

            sums["loss"] += float(loss.detach().cpu())
            sums["scmae_loss"] += float(sc_loss.detach().cpu())
            sums["dec_loss"] += float(d_loss.detach().cpu())
            sums["variance_loss"] += float(v_loss.detach().cpu())
            sums["consistency_loss"] += float(consistency.detach().cpu())
            batches += 1

        for key in sums:
            history[key].append(sums[key] / max(1, batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"scmae={history['scmae_loss'][-1]:.4f} dec={history['dec_loss'][-1]:.4f} "
                f"var={history['variance_loss'][-1]:.4f}"
            )

    embedding, q_out, labels_out = extract_all(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    std_profile = effective_dimensionality(embedding.std(axis=0))

    eval_result = None
    preds = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            args.method_name,
            args.seed,
            embedding,
            labels_out,
            n_clusters,
            {
                "variant": args.variant_name,
                "std_floor_weight": float(args.variance_weight),
                "cluster_method_train": "DEC_sharpened_target",
                "mask_strategy": "zero",
            },
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64) if preds is not None else np.zeros(n_clusters)
    frac = counts / max(1.0, counts.sum())
    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "method_raw": args.variant_name,
        "seed": int(args.seed),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "std_profile": std_profile,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
    }
    save_json(summary, str(save_dir / "summary.json"))
    ari = summary["fixed_metrics"].get("kmeans_known_k", {}).get("ari")
    nmi = summary["fixed_metrics"].get("kmeans_known_k", {}).get("nmi")
    print(f"[RESULT] {dataset_name} VarFloor-scMAE ARI={ari} NMI={nmi}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
