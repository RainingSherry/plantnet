#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import math
import os
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[_var] = "1"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import h5py
import numpy as np
import torch
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import DataLoader

HERE = Path(__file__).resolve().parent
ROOT = next(parent for parent in [HERE, *HERE.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(HERE))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import save_json
from model import AttentionMatchingTeacherScMAE
from loss import AttentionMatchingTeacherLoss


METHOD_DIR = "suture12_attention_matching_teacher_consistency_full"
METHOD_NAME = "scMAE + attention matching teacher consistency"
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
    p.add_argument("--stage", choices=["screen", "formal"], default="screen")
    p.add_argument("--hidden_dim", type=int, default=512)
    p.add_argument("--latent_dim", type=int, default=32)
    p.add_argument("--dropout", type=float, default=0.1)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--match_weight", type=float, default=0.035)
    p.add_argument("--match_loss_weight", type=float, default=0.03)
    p.add_argument("--target_gate", type=float, default=0.08)
    return p.parse_args()


def save_h5(path: Path, embedding: np.ndarray, labels: np.ndarray):
    with h5py.File(path, "w") as h:
        h.create_dataset("X", data=embedding.astype(np.float32))
        h.create_dataset("labels", data=labels.astype(np.int64))


def append_csv(row: dict):
    root = Path(__file__).resolve().parents[2]
    fieldnames = ["stage", "method", "dataset", "seed", "acc", "nmi", "ari", "f1_macro", "baseline_nmi", "baseline_ari", "meets_baseline_any", "collapse_warning", "embedding_variance", "neighbor_purity_proxy", "mixed_cell_fraction", "run_dir"]
    for path in [root / "新模型独立快筛单次结果.csv", root / "新模型独立快筛汇总结果.csv"]:
        exists = path.exists()
        with path.open("a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            if not exists:
                writer.writeheader()
            writer.writerow(row)


@torch.no_grad()
def extract_embedding(model, loader, device):
    model.eval()
    emb, labels = [], []
    for _, x, y in loader:
        emb.append(model.feature(x.to(device)).detach().cpu().numpy())
        labels.append(y.numpy())
    return np.concatenate(emb, axis=0).astype(np.float32), np.concatenate(labels, axis=0).astype(np.int64)


@torch.no_grad()
def collect_matching(model, loader, device):
    model.eval()
    gates, attention, similarity = [], [], []
    for _, x, _ in loader:
        out = model.encode_matched(x.to(device), x.to(device))
        gates.append(out["match_gate"].detach().cpu().numpy())
        attention.append(out["match_attention"].detach().cpu().numpy())
        similarity.append(out["match_similarity"].detach().cpu().numpy())
    return np.concatenate(gates, axis=0).astype(np.float32), np.concatenate(attention, axis=0).astype(np.float32), np.concatenate(similarity, axis=0).astype(np.float32)


def diagnostics(embedding, labels, pred, gate, attention, similarity):
    if pred is None:
        pred = np.zeros_like(labels)
    masses = np.bincount(pred.astype(np.int64), minlength=max(1, int(pred.max()) + 1)).astype(np.float64)
    masses /= max(float(masses.sum()), 1.0)
    k = min(16, len(labels))
    idx = NearestNeighbors(n_neighbors=k, n_jobs=1).fit(embedding).kneighbors(embedding, return_distance=False)[:, 1:] if k > 1 else np.empty((len(labels), 0), dtype=np.int64)
    emb_var = float(np.var(embedding, axis=0).mean())
    gate_entropy = -(gate * np.log(np.clip(gate, 1e-6, 1.0)) + (1.0 - gate) * np.log(np.clip(1.0 - gate, 1e-6, 1.0))).mean() / math.log(2.0)
    return {
        "edge_survival": float(1.0 - gate.mean()),
        "neighbor_purity_proxy": float((labels[idx] == labels[:, None]).mean()) if idx.size else 1.0,
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(gate_entropy),
        "rare_risk_fraction": float((gate[:, 0] > 0.65).mean()) if gate.ndim == 2 else float((gate > 0.65).mean()),
        "embedding_variance": emb_var,
        "cluster_mass_min": float(masses[masses > 0].min()) if np.any(masses > 0) else 0.0,
        "cluster_mass_max": float(masses.max(initial=0.0)),
        "collapse_warning": bool(emb_var < 1e-5 or masses.max(initial=0.0) > 0.90 or np.sum(masses > 0) <= 1),
        "match_gate_mean": float(gate.mean()),
        "match_attention_mean": float(attention.mean()),
        "match_similarity_mean": float(similarity.mean()),
        "diagnostic_note": "Attention matching is applied between masked and detached clean latent views; no graph or cell mixing is used.",
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
    model = AttentionMatchingTeacherScMAE(enc.data.shape[1], args.hidden_dim, args.latent_dim, args.dropout, args.mask_prob, args.match_weight).to(device)
    loss_fn = AttentionMatchingTeacherLoss(match_weight=args.match_loss_weight, target_gate=args.target_gate)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    stage = "smoke" if args.smoke else args.stage
    print(f"Dataset={args.dataset_name} cells={enc.data.shape[0]} genes={enc.data.shape[1]} clusters={args.n_clusters}")
    print(f"Method={METHOD_DIR} stage={stage} epochs={args.epochs}")
    history, start = [], time.time()
    for epoch in range(1, args.epochs + 1):
        model.train()
        rows = []
        for idx, x, _ in loader:
            idx = idx.to(device)
            x = x.to(device)
            outputs = model(x)
            loss, parts = loss_fn(outputs, target_tensor[idx])
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
            print(f"Epoch {epoch:03d}/{args.epochs} loss={rec['total']:.4f} rec={rec['reconstruction']:.4f} mask={rec['mask_bce']:.4f} match={rec['match_loss']:.4f} gate={rec['match_gate_mean']:.3f}")
    embedding, labels = extract_embedding(model, eval_loader, device)
    np.save(out / "embedding_final.npy", embedding)
    save_h5(out / "embedding.h5", embedding, labels)
    save_json({"history": history}, str(out / "training_history.json"))
    torch.save({"model_state_dict": model.state_dict(), "args": vars(args)}, out / "model_checkpoint.pth")
    gate, attention, similarity = collect_matching(model, eval_loader, device)
    np.save(out / "match_gate.npy", gate)
    np.save(out / "match_attention.npy", attention)
    np.save(out / "match_similarity.npy", similarity)
    fixed, pred = {}, None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(out, args.dataset_name, METHOD_DIR, args.seed, embedding, labels, args.n_clusters, {"stage": stage})
        fixed = result["fixed"]
        pred = result["preds"]["kmeans_known_k"]
        save_json(fixed, str(out / "metrics.json"))
    diag = diagnostics(embedding, labels, pred, gate, attention, similarity)
    save_json(diag, str(out / "diagnostics.json"))
    metric = fixed.get("kmeans_known_k", {}) if fixed else {}
    base = BASELINES.get(args.dataset_name, {})
    meets = bool(metric and ((metric.get("nmi", -1.0) >= base.get("nmi", float("inf"))) or (metric.get("ari", -1.0) >= base.get("ari", float("inf")))) and not diag["collapse_warning"])
    summary = {"dataset": args.dataset_name, "method": METHOD_NAME, "method_dir": METHOD_DIR, "stage": stage, "seed": args.seed, "n_cells": int(enc.data.shape[0]), "n_genes": int(enc.data.shape[1]), "n_clusters": int(args.n_clusters), "runtime_seconds": float(time.time() - start), "embedding_path": str((out / "embedding_final.npy").resolve()), "fixed_metrics": fixed, "diagnostics": diag, "baseline": base, "meets_screen_baseline_any": meets, "note": "Quick screen only; this file is not appended to 全benchmark结果.csv."}
    save_json(summary, str(out / "summary.json"))
    if metric and stage != "formal":
        append_csv({"stage": stage, "method": METHOD_DIR, "dataset": args.dataset_name, "seed": args.seed, "acc": metric.get("acc", np.nan), "nmi": metric.get("nmi", np.nan), "ari": metric.get("ari", np.nan), "f1_macro": metric.get("f1_macro", np.nan), "baseline_nmi": base.get("nmi", np.nan), "baseline_ari": base.get("ari", np.nan), "meets_baseline_any": meets, "collapse_warning": diag.get("collapse_warning", np.nan), "embedding_variance": diag.get("embedding_variance", np.nan), "neighbor_purity_proxy": diag.get("neighbor_purity_proxy", np.nan), "mixed_cell_fraction": diag.get("mixed_cell_fraction", np.nan), "run_dir": str(out.resolve())})
    print(f"Completed {METHOD_DIR}. Results saved to: {out}")


if __name__ == "__main__":
    main()
