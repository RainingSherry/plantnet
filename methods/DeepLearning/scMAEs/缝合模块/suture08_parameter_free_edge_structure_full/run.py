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
from model import ParameterFreeEdgeStructureScMAE
from loss import ParameterFreeEdgeStructureLoss


METHOD_DIR = "suture08_parameter_free_edge_structure_full"
METHOD_NAME = "scMAE + parameter-free edge/structure reliability"
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
    p.add_argument("--adapter_weight", type=float, default=0.06)
    p.add_argument("--edge_k", type=int, default=15)
    p.add_argument("--max_exact_cells", type=int, default=20000)
    p.add_argument("--anchor_count", type=int, default=1024)
    p.add_argument("--anchor_chunk", type=int, default=48)
    return p.parse_args()


def save_h5(path: Path, embedding: np.ndarray, labels: np.ndarray):
    with h5py.File(path, "w") as h:
        h.create_dataset("X", data=embedding.astype(np.float32))
        h.create_dataset("labels", data=labels.astype(np.int64))


def append_csv(row: dict):
    root = Path(__file__).resolve().parents[2]
    fieldnames = [
        "stage", "method", "dataset", "seed", "acc", "nmi", "ari", "f1_macro",
        "baseline_nmi", "baseline_ari", "meets_baseline_any", "collapse_warning",
        "embedding_variance", "neighbor_purity_proxy", "mixed_cell_fraction", "run_dir",
    ]
    for path in [root / "新模型独立快筛单次结果.csv", root / "新模型独立快筛汇总结果.csv"]:
        exists = path.exists()
        with path.open("a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            if not exists:
                writer.writeheader()
            writer.writerow(row)


def exact_graph_structure(log_expr: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = log_expr.shape[0]
    k_eff = min(k + 1, n)
    nn = NearestNeighbors(n_neighbors=k_eff, metric="euclidean", n_jobs=1).fit(log_expr)
    dist, idx = nn.kneighbors(log_expr, return_distance=True)
    neighbors = idx[:, 1:]
    d = dist[:, 1:]
    sigma = np.median(d, axis=1, keepdims=True) + 1e-6
    weights = np.exp(-(d ** 2) / (2.0 * sigma ** 2)).astype(np.float32)
    weights /= np.maximum(weights.sum(axis=1, keepdims=True), 1e-6)
    structure = np.einsum("nk,nkd->nd", weights, log_expr[neighbors]).astype(np.float32)
    return structure, neighbors.astype(np.int64), weights


def anchor_graph_structure(log_expr: np.ndarray, k: int, anchor_count: int, chunk_size: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    n, d = log_expr.shape
    n_anchor = min(max(k + 1, anchor_count), n)
    anchor_idx = rng.choice(n, size=n_anchor, replace=False)
    anchors = log_expr[anchor_idx].astype(np.float32, copy=False)
    k_eff = min(k, n_anchor)
    structure = np.zeros_like(log_expr, dtype=np.float32)
    neighbors = np.zeros((n, k_eff), dtype=np.int64)
    weights = np.zeros((n, k_eff), dtype=np.float32)
    anchor_norm = np.sum(anchors * anchors, axis=1)
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        chunk = log_expr[start:end].astype(np.float32, copy=False)
        d2 = np.sum(chunk * chunk, axis=1, keepdims=True) + anchor_norm[None, :] - 2.0 * (chunk @ anchors.T)
        d2 = np.maximum(d2, 0.0)
        part = np.argpartition(d2, kth=k_eff - 1, axis=1)[:, :k_eff]
        row = np.arange(end - start)[:, None]
        order = np.argsort(d2[row, part], axis=1)
        part = part[row, order]
        picked_d2 = d2[row, part]
        sigma = np.median(np.sqrt(picked_d2 + 1e-8), axis=1, keepdims=True) + 1e-6
        w = np.exp(-picked_d2 / (2.0 * sigma ** 2)).astype(np.float32)
        w /= np.maximum(w.sum(axis=1, keepdims=True), 1e-6)
        structure[start:end] = np.einsum("nk,nkd->nd", w, anchors[part]).astype(np.float32)
        neighbors[start:end] = anchor_idx[part]
        weights[start:end] = w
    return structure, neighbors, weights


def build_edge_structure(log_expr: np.ndarray, k: int, max_exact_cells: int, anchor_count: int, chunk_size: int, seed: int) -> dict:
    x = np.nan_to_num(log_expr.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    if x.shape[0] <= max_exact_cells:
        structure, neighbors, weights = exact_graph_structure(x, k)
        mode = "exact_cell_knn"
    else:
        structure, neighbors, weights = anchor_graph_structure(x, k, anchor_count, chunk_size, seed)
        mode = "anchor_knn_fallback"
    edge = x - structure
    x_norm = np.linalg.norm(x, axis=1)
    low_norm = np.linalg.norm(structure, axis=1)
    high_norm = np.linalg.norm(edge, axis=1)
    low_scaled = low_norm / (np.median(low_norm) + 1e-6)
    high_scaled = high_norm / (np.median(high_norm) + 1e-6)
    reliability = low_scaled / (low_scaled + high_scaled + 1e-6)
    reliability = np.clip(reliability, 0.02, 0.98).astype(np.float32)
    edge_score = np.clip(high_scaled / (low_scaled + high_scaled + 1e-6), 0.0, 1.0).astype(np.float32)
    return {
        "structure": structure.astype(np.float32),
        "edge": edge.astype(np.float32),
        "neighbors": neighbors,
        "weights": weights.astype(np.float32),
        "reliability": reliability[:, None],
        "edge_score": edge_score[:, None],
        "mode": mode,
        "x_norm_mean": float(x_norm.mean()),
    }


@torch.no_grad()
def extract_embedding(model, loader, structure_tensor, reliability_tensor, device):
    model.eval()
    emb, labels = [], []
    for idx, x, y in loader:
        idx = idx.to(device)
        z = model.feature(x.to(device), structure_tensor[idx], reliability_tensor[idx])
        emb.append(z.detach().cpu().numpy())
        labels.append(y.numpy())
    return np.concatenate(emb, axis=0).astype(np.float32), np.concatenate(labels, axis=0).astype(np.int64)


@torch.no_grad()
def collect_adapter(model, loader, structure_tensor, reliability_tensor, device):
    model.eval()
    rel, delta_norm = [], []
    for idx, x, _ in loader:
        idx = idx.to(device)
        out = model.encode_with_reliability(x.to(device), structure_tensor[idx], reliability_tensor[idx])
        rel.append(out["reliability"].detach().cpu().numpy())
        delta_norm.append(out["adapter_delta"].norm(dim=1).detach().cpu().numpy())
    return np.concatenate(rel, axis=0).astype(np.float32), np.concatenate(delta_norm, axis=0).astype(np.float32)


def diagnostics(embedding, labels, pred, reliability, edge_score, neighbors, graph_mode):
    if pred is None:
        pred = np.zeros_like(labels)
    masses = np.bincount(pred.astype(np.int64), minlength=max(1, int(pred.max()) + 1)).astype(np.float64)
    masses /= max(float(masses.sum()), 1.0)
    emb_var = float(np.var(embedding, axis=0).mean())
    purity = float((labels[neighbors] == labels[:, None]).mean()) if neighbors.size else 1.0
    boundary_prob = 1.0 - reliability
    entropy = -(boundary_prob * np.log(np.clip(boundary_prob, 1e-6, 1.0)) + reliability * np.log(np.clip(reliability, 1e-6, 1.0))).mean() / math.log(2.0)
    collapse = bool(emb_var < 1e-5 or masses.max(initial=0.0) > 0.90 or np.sum(masses > 0) <= 1)
    return {
        "edge_survival": float(1.0 - 0.06 * reliability.mean()),
        "neighbor_purity_proxy": purity,
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(entropy),
        "rare_risk_fraction": float((boundary_prob > 0.65).mean()),
        "embedding_variance": emb_var,
        "cluster_mass_min": float(masses[masses > 0].min()) if np.any(masses > 0) else 0.0,
        "cluster_mass_max": float(masses.max(initial=0.0)),
        "collapse_warning": collapse,
        "reliability_mean": float(reliability.mean()),
        "reliability_std": float(reliability.std()),
        "edge_score_mean": float(edge_score.mean()),
        "graph_mode": graph_mode,
        "diagnostic_note": "Parameter-free reliability from log-expression graph low-frequency structure and high-frequency residual; no image FFT and no cell mixing.",
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

    graph = build_edge_structure(target.data, args.edge_k, args.max_exact_cells, args.anchor_count, args.anchor_chunk, args.seed)
    np.save(out / "structure_context.npy", graph["structure"])
    np.save(out / "edge_residual.npy", graph["edge"])
    np.save(out / "edge_neighbors.npy", graph["neighbors"])
    np.save(out / "edge_weights.npy", graph["weights"])
    np.save(out / "reliability.npy", graph["reliability"])
    np.save(out / "edge_score.npy", graph["edge_score"])

    target_tensor = torch.as_tensor(target.data, dtype=torch.float32, device=device)
    structure_tensor = torch.as_tensor(graph["structure"], dtype=torch.float32, device=device)
    reliability_tensor = torch.as_tensor(graph["reliability"], dtype=torch.float32, device=device)
    ds = family.IndexedExpressionDataset(enc.data, target.labels)
    loader = DataLoader(ds, batch_size=args.batch_size, shuffle=True, drop_last=False)
    eval_loader = DataLoader(ds, batch_size=args.batch_size, shuffle=False, drop_last=False)
    model = ParameterFreeEdgeStructureScMAE(enc.data.shape[1], args.hidden_dim, args.latent_dim, args.dropout, args.mask_prob, args.adapter_weight).to(device)
    loss_fn = ParameterFreeEdgeStructureLoss()
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    stage = "smoke" if args.smoke else args.stage
    print(f"Dataset={args.dataset_name} cells={enc.data.shape[0]} genes={enc.data.shape[1]} clusters={args.n_clusters}")
    print(f"Method={METHOD_DIR} graph_mode={graph['mode']} stage={stage} epochs={args.epochs}")
    history, start = [], time.time()
    for epoch in range(1, args.epochs + 1):
        model.train()
        rows = []
        for idx, x, _ in loader:
            idx = idx.to(device)
            x = x.to(device)
            outputs = model(x, structure_tensor[idx], reliability_tensor[idx])
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
            print(f"Epoch {epoch:03d}/{args.epochs} loss={rec['total']:.4f} rec={rec['reconstruction']:.4f} mask={rec['mask_bce']:.4f} align={rec['structure_align']:.4f} rel={rec['reliability_mean']:.3f}")

    embedding, labels = extract_embedding(model, eval_loader, structure_tensor, reliability_tensor, device)
    np.save(out / "embedding_final.npy", embedding)
    save_h5(out / "embedding.h5", embedding, labels)
    save_json({"history": history}, str(out / "training_history.json"))
    torch.save({"model_state_dict": model.state_dict(), "args": vars(args)}, out / "model_checkpoint.pth")
    rel_seen, delta_norm = collect_adapter(model, eval_loader, structure_tensor, reliability_tensor, device)
    np.save(out / "adapter_reliability_seen.npy", rel_seen)
    np.save(out / "adapter_delta_norm.npy", delta_norm)

    fixed, pred = {}, None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(out, args.dataset_name, METHOD_DIR, args.seed, embedding, labels, args.n_clusters, {"stage": stage})
        fixed = result["fixed"]
        pred = result["preds"]["kmeans_known_k"]
        save_json(fixed, str(out / "metrics.json"))
    diag = diagnostics(embedding, labels, pred, graph["reliability"], graph["edge_score"], graph["neighbors"], graph["mode"])
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
    if metric and stage != "formal":
        append_csv({
            "stage": stage,
            "method": METHOD_DIR,
            "dataset": args.dataset_name,
            "seed": args.seed,
            "acc": metric.get("acc", np.nan),
            "nmi": metric.get("nmi", np.nan),
            "ari": metric.get("ari", np.nan),
            "f1_macro": metric.get("f1_macro", np.nan),
            "baseline_nmi": base.get("nmi", np.nan),
            "baseline_ari": base.get("ari", np.nan),
            "meets_baseline_any": meets,
            "collapse_warning": diag.get("collapse_warning", np.nan),
            "embedding_variance": diag.get("embedding_variance", np.nan),
            "neighbor_purity_proxy": diag.get("neighbor_purity_proxy", np.nan),
            "mixed_cell_fraction": diag.get("mixed_cell_fraction", np.nan),
            "run_dir": str(out.resolve()),
        })
    print(f"Completed {METHOD_DIR}. Results saved to: {out}")


if __name__ == "__main__":
    main()
