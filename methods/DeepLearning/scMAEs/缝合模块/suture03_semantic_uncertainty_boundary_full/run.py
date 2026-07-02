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
from model import SemanticUncertaintyScMAE
from loss import SemanticUncertaintyLoss


METHOD_DIR = "suture03_semantic_uncertainty_boundary_full"
METHOD_NAME = "scMAE + semantic uncertainty boundary gate"
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
    p.add_argument("--dropout", type=float, default=0.1)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--module_weight", type=float, default=0.2)
    return p.parse_args()


def normalize01(x: np.ndarray) -> np.ndarray:
    lo, hi = np.quantile(x, [0.05, 0.95])
    return np.clip((x - lo) / (hi - lo + 1e-6), 0.0, 1.0).astype(np.float32)


def build_gate_targets(input_expr: np.ndarray, log_expr: np.ndarray, k: int = 15) -> tuple[np.ndarray, dict]:
    n = input_expr.shape[0]
    k = min(k + 1, n)
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(input_expr)
    dist, _ = nn.kneighbors(input_expr, return_distance=True)
    boundary = normalize01(dist[:, 1:].mean(axis=1))
    dropout = (log_expr <= 1e-6).mean(axis=1)
    total = np.maximum(log_expr.sum(axis=1), 1e-6)
    concentration = np.sort(log_expr, axis=1)[:, -max(1, log_expr.shape[1] // 20):].sum(axis=1) / total
    rare = normalize01(0.6 * dropout + 0.4 * concentration)
    core = np.clip(1.0 - np.maximum(boundary, rare), 0.05, 1.0)
    mat = np.stack([core, boundary + 0.05, rare + 0.05], axis=1).astype(np.float32)
    mat /= mat.sum(axis=1, keepdims=True)
    return mat, {
        "boundary_score_mean": float(boundary.mean()),
        "rare_score_mean": float(rare.mean()),
        "core_score_mean": float(core.mean()),
    }


def save_h5(path: Path, embedding: np.ndarray, labels: np.ndarray) -> None:
    with h5py.File(path, "w") as h:
        h.create_dataset("X", data=embedding.astype(np.float32))
        h.create_dataset("labels", data=labels.astype(np.int64))


def append_csv(row: dict) -> None:
    root = Path(__file__).resolve().parents[2]
    for path in [root / "新模型独立快筛单次结果.csv", root / "新模型独立快筛汇总结果.csv"]:
        exists = path.exists()
        with path.open("a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=list(row.keys()))
            if not exists:
                writer.writeheader()
            writer.writerow(row)


def diagnostics(embedding, labels, pred, gates, target_stats):
    if pred is None:
        pred = np.zeros_like(labels)
    masses = np.bincount(pred.astype(np.int64), minlength=max(1, int(pred.max()) + 1)).astype(np.float64)
    masses /= max(float(masses.sum()), 1.0)
    nn = NearestNeighbors(n_neighbors=min(16, len(labels))).fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    entropy = -(gates * np.log(np.clip(gates, 1e-6, 1.0))).sum(axis=1) / math.log(3.0)
    emb_var = float(np.var(embedding, axis=0).mean())
    return {
        "edge_survival": float(1.0 - gates[:, 1].mean()),
        "neighbor_purity_proxy": float((labels[idx] == labels[:, None]).mean()) if idx.size else 1.0,
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(entropy.mean()),
        "rare_risk_fraction": float((gates[:, 2] > 0.45).mean()),
        "embedding_variance": emb_var,
        "cluster_mass_min": float(masses[masses > 0].min()) if np.any(masses > 0) else 0.0,
        "cluster_mass_max": float(masses.max(initial=0.0)),
        "collapse_warning": bool(emb_var < 1e-5 or masses.max(initial=0.0) > 0.90 or np.sum(masses > 0) <= 1),
        "gate_core_mean": float(gates[:, 0].mean()),
        "gate_boundary_mean": float(gates[:, 1].mean()),
        "gate_rare_mean": float(gates[:, 2].mean()),
        "target_stats": target_stats,
        "diagnostic_note": "SID-inspired semantic decoupling is used as a cell-level uncertainty gate; no cell mixing is used.",
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
    gate_targets, target_stats = build_gate_targets(enc.data, target.data)
    np.save(out / "gate_targets.npy", gate_targets)
    target_tensor = torch.as_tensor(target.data, dtype=torch.float32, device=device)
    gate_tensor = torch.as_tensor(gate_targets, dtype=torch.float32, device=device)

    ds = family.IndexedExpressionDataset(enc.data, target.labels)
    loader = DataLoader(ds, batch_size=args.batch_size, shuffle=True, drop_last=False)
    eval_loader = DataLoader(ds, batch_size=args.batch_size, shuffle=False, drop_last=False)
    model = SemanticUncertaintyScMAE(enc.data.shape[1], args.hidden_dim, args.latent_dim, args.dropout, args.mask_prob, args.module_weight).to(device)
    loss_fn = SemanticUncertaintyLoss()
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    stage = "smoke" if args.smoke else "screen"
    print(f"Dataset={args.dataset_name} cells={enc.data.shape[0]} genes={enc.data.shape[1]} clusters={args.n_clusters}")
    print(f"Method={METHOD_DIR} stage={stage} epochs={args.epochs}")
    history = []
    start = time.time()
    for epoch in range(1, args.epochs + 1):
        model.train()
        rows = []
        for idx, x, _ in loader:
            x = x.to(device)
            idx = idx.to(device)
            outputs = model(x)
            loss, parts = loss_fn(outputs, target_tensor[idx], gate_tensor[idx])
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
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={rec['total']:.4f} rec={rec['reconstruction']:.4f} "
                f"mask={rec['mask_bce']:.4f} gates=({rec['gate_core']:.2f},{rec['gate_boundary']:.2f},{rec['gate_rare']:.2f})"
            )

    embedding, labels = family.extract_embedding(model, eval_loader, device)
    np.save(out / "embedding_final.npy", embedding)
    save_h5(out / "embedding.h5", embedding, labels)
    save_json({"history": history}, str(out / "training_history.json"))
    torch.save({"model_state_dict": model.state_dict(), "args": vars(args)}, out / "model_checkpoint.pth")

    with torch.no_grad():
        all_gates = []
        for _, x, _ in eval_loader:
            _, _, _, g = model.encode_with_gate(x.to(device))
            all_gates.append(g.detach().cpu().numpy())
    gates_np = np.concatenate(all_gates, axis=0).astype(np.float32)
    np.save(out / "semantic_gates.npy", gates_np)

    fixed, pred = {}, None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(out, args.dataset_name, METHOD_DIR, args.seed, embedding, labels, args.n_clusters, {"stage": stage})
        fixed = result["fixed"]
        pred = result["preds"]["kmeans_known_k"]
        save_json(fixed, str(out / "metrics.json"))
    diag = diagnostics(embedding, labels, pred, gates_np, target_stats)
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
