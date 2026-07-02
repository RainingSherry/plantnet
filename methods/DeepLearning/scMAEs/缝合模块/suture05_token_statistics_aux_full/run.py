#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from pathlib import Path

import h5py
import numpy as np
import torch
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import DataLoader

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in [HERE, *HERE.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import save_json
from model import TokenStatisticsAuxScMAE
from loss import TokenStatisticsAuxLoss


METHOD_DIR = "suture05_token_statistics_aux_full"
METHOD_NAME = "scMAE + gene token statistics auxiliary branch"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


def parse_args():
    p = argparse.ArgumentParser(description=METHOD_NAME)
    p.add_argument("--data_path", required=True)
    p.add_argument("--save_dir", required=True)
    p.add_argument("--dataset_name", required=True)
    p.add_argument("--label_key", default="auto")
    p.add_argument("--input_mode", default="auto")
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
    p.add_argument("--no_save_h5ad", action="store_true")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--hidden_dim", type=int, default=512)
    p.add_argument("--latent_dim", type=int, default=32)
    p.add_argument("--token_dim", type=int, default=32)
    p.add_argument("--dropout", type=float, default=0.1)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--aux_weight", type=float, default=0.15)
    return p.parse_args()


def save_h5(path, embedding, labels):
    with h5py.File(path, "w") as h:
        h.create_dataset("X", data=embedding.astype(np.float32))
        h.create_dataset("labels", data=labels.astype(np.int64))


def append_csv(row):
    root = Path(__file__).resolve().parents[2]
    for path in [root / "新模型独立快筛单次结果.csv", root / "新模型独立快筛汇总结果.csv"]:
        exists = path.exists()
        with path.open("a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(row.keys()))
            if not exists:
                writer.writeheader()
            writer.writerow(row)


def collect_importance(model, loader, device):
    model.eval()
    vals = []
    with torch.no_grad():
        for _, x, _ in loader:
            _, imp = model.token_stats(x.to(device))
            vals.append(imp.detach().cpu().numpy())
    return np.concatenate(vals, axis=0).astype(np.float32)


def diagnostics(embedding, labels, pred, importance):
    if pred is None:
        pred = np.zeros_like(labels)
    masses = np.bincount(pred.astype(np.int64), minlength=max(1, int(pred.max()) + 1)).astype(np.float64)
    masses /= max(float(masses.sum()), 1.0)
    idx = NearestNeighbors(n_neighbors=min(16, len(labels))).fit(embedding).kneighbors(embedding, return_distance=False)[:, 1:]
    entropy = -(importance * np.log(np.clip(importance, 1e-8, 1.0))).sum(axis=1) / math.log(importance.shape[1])
    emb_var = float(np.var(embedding, axis=0).mean())
    return {
        "edge_survival": 1.0,
        "neighbor_purity_proxy": float((labels[idx] == labels[:, None]).mean()) if idx.size else 1.0,
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(entropy.mean()),
        "rare_risk_fraction": float((importance.max(axis=1) > 0.10).mean()),
        "embedding_variance": emb_var,
        "cluster_mass_min": float(masses[masses > 0].min()) if np.any(masses > 0) else 0.0,
        "cluster_mass_max": float(masses.max(initial=0.0)),
        "collapse_warning": bool(emb_var < 1e-5 or masses.max(initial=0.0) > 0.90 or np.sum(masses > 0) <= 1),
        "token_entropy_mean": float(entropy.mean()),
        "token_importance_max_mean": float(importance.max(axis=1).mean()),
        "diagnostic_note": "TSSA-inspired branch summarizes gene-token statistics linearly; no pairwise attention, no image reshape, and no cell mixing.",
    }


def main():
    args = parse_args()
    family.set_seed(args.seed)
    out = Path(args.save_dir)
    out.mkdir(parents=True, exist_ok=True)
    save_json(vars(args), str(out / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    print(f"Using device: {device}")
    target = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed)
    enc = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, args.scale_input, args.label_key, args.seed)
    save_json(target.profile, str(out / "dataset_profile.json"))
    save_json(target.preprocess_config, str(out / "preprocess_config.json"))
    np.save(out / "labels.npy", target.labels.astype(np.int64))
    target_tensor = torch.as_tensor(target.data, dtype=torch.float32, device=device)
    ds = family.IndexedExpressionDataset(enc.data, target.labels)
    loader = DataLoader(ds, batch_size=args.batch_size, shuffle=True, drop_last=False)
    eval_loader = DataLoader(ds, batch_size=args.batch_size, shuffle=False, drop_last=False)
    model = TokenStatisticsAuxScMAE(enc.data.shape[1], args.hidden_dim, args.latent_dim, args.token_dim, args.dropout, args.mask_prob, args.aux_weight).to(device)
    loss_fn = TokenStatisticsAuxLoss()
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    stage = "smoke" if args.smoke else "screen"
    print(f"Dataset={args.dataset_name} cells={enc.data.shape[0]} genes={enc.data.shape[1]} clusters={args.n_clusters}")
    print(f"Method={METHOD_DIR} stage={stage} epochs={args.epochs}")
    history, start = [], time.time()
    for epoch in range(1, args.epochs + 1):
        model.train()
        rows = []
        for idx, x, _ in loader:
            x = x.to(device)
            outputs = model(x)
            loss, parts = loss_fn(outputs, target_tensor[idx.to(device)])
            if not torch.isfinite(loss):
                raise FloatingPointError("loss contains NaN or Inf")
            opt.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            rows.append(parts)
        rec = {k: float(np.mean([r[k] for r in rows])) for k in rows[0]}
        rec["epoch"] = epoch
        history.append(rec)
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={rec['total']:.4f} rec={rec['reconstruction']:.4f} mask={rec['mask_bce']:.4f} stat={rec['stat']:.4f} ent={rec['token_entropy']:.3f}")

    embedding, labels = family.extract_embedding(model, eval_loader, device)
    np.save(out / "embedding_final.npy", embedding)
    save_h5(out / "embedding.h5", embedding, labels)
    save_json({"history": history}, str(out / "training_history.json"))
    torch.save({"model_state_dict": model.state_dict(), "args": vars(args)}, out / "model_checkpoint.pth")
    importance = collect_importance(model, eval_loader, device)
    np.save(out / "token_importance.npy", importance)
    fixed, pred = {}, None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(out, args.dataset_name, METHOD_DIR, args.seed, embedding, labels, args.n_clusters, {"stage": stage})
        fixed = result["fixed"]
        pred = result["preds"]["kmeans_known_k"]
        save_json(fixed, str(out / "metrics.json"))
    diag = diagnostics(embedding, labels, pred, importance)
    save_json(diag, str(out / "diagnostics.json"))
    metric = fixed.get("kmeans_known_k", {}) if fixed else {}
    base = BASELINES.get(args.dataset_name, {})
    meets = bool(metric and ((metric.get("nmi", -1.0) >= base.get("nmi", float("inf"))) or (metric.get("ari", -1.0) >= base.get("ari", float("inf")))) and not diag["collapse_warning"])
    summary = {
        "dataset": args.dataset_name,
        "method": METHOD_NAME,
        "method_dir": METHOD_DIR,
        "stage": stage,
        "seed": args.seed,
        "n_cells": int(enc.data.shape[0]),
        "n_genes": int(enc.data.shape[1]),
        "n_clusters": int(args.n_clusters),
        "runtime_seconds": float(time.time() - start),
        "embedding_path": str((out / "embedding_final.npy").resolve()),
        "fixed_metrics": fixed,
        "diagnostics": diag,
        "baseline": base,
        "meets_screen_baseline_any": meets,
        "note": "Quick screen only; this file is not appended to 全benchmark结果.csv.",
    }
    save_json(summary, str(out / "summary.json"))
    if metric:
        append_csv({
            "method_dir": METHOD_DIR,
            "dataset": args.dataset_name,
            "stage": stage,
            "seed": args.seed,
            "nmi": metric.get("nmi", np.nan),
            "ari": metric.get("ari", np.nan),
            "acc": metric.get("acc", np.nan),
            "collapse_warning": diag["collapse_warning"],
            "meets_screen_baseline_any": meets,
            "save_dir": str(out),
        })
    print(f"Completed {METHOD_DIR}. Results saved to: {out}")


if __name__ == "__main__":
    main()
