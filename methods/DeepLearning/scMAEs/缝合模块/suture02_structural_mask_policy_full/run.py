#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

import h5py
import numpy as np
import pandas as pd
import torch
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import DataLoader

THIS_DIR = Path(__file__).resolve().parent
ROOT = next(parent for parent in [THIS_DIR, *THIS_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import save_json
from model import StructuralMaskScMAE
from loss import StructuralMaskLoss


METHOD_DIR = "suture02_structural_mask_policy_full"
METHOD_NAME = "scMAE + structural gene mask policy"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=METHOD_NAME)
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", required=True)
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto")
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
    parser.add_argument("--no_save_h5ad", action="store_true")
    parser.add_argument("--smoke", action="store_true")
    parser.add_argument("--hidden_dim", type=int, default=512)
    parser.add_argument("--latent_dim", type=int, default=32)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--mask_prob", type=float, default=0.4)
    parser.add_argument("--policy_weight", type=float, default=1.0)
    parser.add_argument("--num_workers", type=int, default=0)
    return parser.parse_args()


def save_numpy_h5(path: Path, embedding: np.ndarray, labels: np.ndarray) -> None:
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=embedding.astype(np.float32))
        handle.create_dataset("labels", data=labels.astype(np.int64))


def finite_or_raise(value: torch.Tensor, name: str) -> None:
    if not torch.isfinite(value).all():
        raise FloatingPointError(f"{name} contains NaN or Inf")


def gene_stats_from_log(log_expr: np.ndarray, device: torch.device) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    mean = log_expr.mean(axis=0).astype(np.float32)
    var = log_expr.var(axis=0).astype(np.float32)
    dropout = (log_expr <= 1e-6).mean(axis=0).astype(np.float32)
    return (
        torch.as_tensor(mean, dtype=torch.float32, device=device),
        torch.as_tensor(var, dtype=torch.float32, device=device),
        torch.as_tensor(dropout, dtype=torch.float32, device=device),
    )


def neighbor_purity(embedding: np.ndarray, labels: np.ndarray, k: int = 15) -> float:
    if embedding.shape[0] <= 2:
        return 1.0
    k = min(k + 1, embedding.shape[0])
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean").fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    return float((labels[idx] == labels[:, None]).mean())


def diagnostics(
    embedding: np.ndarray,
    labels: np.ndarray,
    pred: np.ndarray | None,
    mask_prob: np.ndarray,
    marker_risk: np.ndarray,
) -> dict:
    emb_var = float(np.var(embedding, axis=0).mean())
    if pred is None:
        pred = np.zeros_like(labels)
    counts = np.bincount(pred.astype(np.int64), minlength=max(int(pred.max()) + 1, 1)).astype(np.float64)
    masses = counts / max(float(counts.sum()), 1.0)
    high_risk = marker_risk >= np.quantile(marker_risk, 0.75)
    low_risk = marker_risk <= np.quantile(marker_risk, 0.25)
    high_mask = float(mask_prob[high_risk].mean()) if np.any(high_risk) else float(mask_prob.mean())
    low_mask = float(mask_prob[low_risk].mean()) if np.any(low_risk) else float(mask_prob.mean())
    edge_survival = float(1.0 - high_mask)
    boundary_entropy = float(np.mean(-(mask_prob * np.log(np.clip(mask_prob, 1e-6, 1.0))
                                      + (1.0 - mask_prob) * np.log(np.clip(1.0 - mask_prob, 1e-6, 1.0)))))
    boundary_entropy /= float(math.log(2.0))
    collapse = bool(emb_var < 1e-5 or masses.max(initial=0.0) > 0.90 or np.sum(masses > 0) <= 1)
    return {
        "edge_survival": edge_survival,
        "neighbor_purity_proxy": neighbor_purity(embedding, labels),
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": boundary_entropy,
        "rare_risk_fraction": float(np.mean(high_risk)),
        "embedding_variance": emb_var,
        "cluster_mass_min": float(masses[masses > 0].min()) if np.any(masses > 0) else 0.0,
        "cluster_mass_max": float(masses.max(initial=0.0)),
        "collapse_warning": collapse,
        "mask_prob_mean": float(mask_prob.mean()),
        "mask_prob_high_risk": high_mask,
        "mask_prob_low_risk": low_mask,
        "marker_risk_mean": float(marker_risk.mean()),
        "diagnostic_note": "SMMM-inspired structural saliency is used only as a gene mask policy; no cell mixing and no image-style 2D convolution.",
    }


def append_screen_csv(row: dict) -> None:
    root = Path(__file__).resolve().parents[2]
    single = root / "新模型独立快筛单次结果.csv"
    summary = root / "新模型独立快筛汇总结果.csv"
    for path in [single, summary]:
        exists = path.exists()
        with path.open("a", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row.keys()))
            if not exists:
                writer.writeheader()
            writer.writerow(row)


def main() -> None:
    args = parse_args()
    family.set_seed(args.seed)
    output_dir = Path(args.save_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    save_json(vars(args), str(output_dir / "args.json"))

    device = family.get_device(args.gpu, args.no_cuda)
    print(f"Using device: {device}")

    target_bundle = family.load_scmae_dataset(
        args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed
    )
    if args.scale_input:
        input_bundle = family.load_scmae_dataset(
            args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed
        )
    else:
        input_bundle = target_bundle
    if input_bundle.data.shape != target_bundle.data.shape:
        raise ValueError("Encoder input and log-expression target shapes diverged after preprocessing.")

    save_json(target_bundle.profile, str(output_dir / "dataset_profile.json"))
    save_json(target_bundle.preprocess_config, str(output_dir / "preprocess_config.json"))
    np.save(output_dir / "labels.npy", target_bundle.labels.astype(np.int64))

    target_tensor = torch.as_tensor(target_bundle.data, dtype=torch.float32, device=device)
    gene_mean, gene_var, gene_dropout = gene_stats_from_log(target_bundle.data, device)
    dataset = family.IndexedExpressionDataset(input_bundle.data, target_bundle.labels)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=False, num_workers=args.num_workers)
    eval_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=False, drop_last=False, num_workers=args.num_workers)

    model = StructuralMaskScMAE(
        input_dim=input_bundle.data.shape[1],
        hidden_dim=args.hidden_dim,
        latent_dim=args.latent_dim,
        dropout=args.dropout,
        mask_prob=args.mask_prob,
        policy_weight=args.policy_weight,
    ).to(device)
    criterion = StructuralMaskLoss(target_mask_prob=args.mask_prob)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    stage = "smoke" if args.smoke else "screen"
    print(f"Dataset={args.dataset_name} cells={input_bundle.data.shape[0]} genes={input_bundle.data.shape[1]} clusters={args.n_clusters}")
    print(f"Method={METHOD_DIR} stage={stage} epochs={args.epochs}")
    history: list[dict] = []
    start = time.time()
    for epoch in range(1, args.epochs + 1):
        model.train()
        losses = []
        parts_accum: dict[str, list[float]] = {}
        for idx, x, _ in loader:
            x = x.to(device)
            y = target_tensor[idx.to(device)]
            out = model(x, gene_mean, gene_var, gene_dropout)
            loss, parts = criterion(out, y)
            finite_or_raise(loss, "loss")
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            losses.append(float(loss.detach().cpu()))
            for key, value in parts.items():
                parts_accum.setdefault(key, []).append(float(value))
        rec = {key: float(np.mean(values)) for key, values in parts_accum.items()}
        rec["epoch"] = epoch
        rec["loss"] = float(np.mean(losses))
        history.append(rec)
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={rec['loss']:.4f} "
                f"rec={rec.get('reconstruction', 0.0):.4f} mask={rec.get('mask_bce', 0.0):.4f} "
                f"p={rec.get('mask_prob_mean', 0.0):.3f} risk={rec.get('marker_risk_mean', 0.0):.3f}"
            )

    embedding, labels = family.extract_embedding(model, eval_loader, device)
    np.save(output_dir / "embedding_final.npy", embedding)
    save_numpy_h5(output_dir / "embedding.h5", embedding, labels)
    with torch.no_grad():
        mask_prob_t, marker_t = model.mask_policy(gene_mean, gene_var, gene_dropout)
    mask_prob = mask_prob_t.detach().cpu().numpy().astype(np.float32)
    marker_risk = marker_t.detach().cpu().numpy().astype(np.float32)
    np.save(output_dir / "mask_prob.npy", mask_prob)
    np.save(output_dir / "marker_risk.npy", marker_risk)
    save_json({"history": history}, str(output_dir / "training_history.json"))
    torch.save({"model_state_dict": model.state_dict(), "args": vars(args)}, output_dir / "model_checkpoint.pth")

    fixed_metrics = {}
    pred = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            output_dir,
            dataset=args.dataset_name,
            method=METHOD_DIR,
            seed=args.seed,
            embedding=embedding,
            labels=labels,
            n_clusters=args.n_clusters,
            extra={"stage": stage},
        )
        fixed_metrics = eval_result["fixed"]
        pred = eval_result["preds"]["kmeans_known_k"]
        save_json(fixed_metrics, str(output_dir / "metrics.json"))

    diag = diagnostics(embedding, labels, pred, mask_prob, marker_risk)
    save_json(diag, str(output_dir / "diagnostics.json"))
    baseline = BASELINES.get(args.dataset_name, {})
    metric = fixed_metrics.get("kmeans_known_k", {}) if fixed_metrics else {}
    meets = bool(metric and ((metric.get("nmi", -1.0) >= baseline.get("nmi", float("inf")))
                            or (metric.get("ari", -1.0) >= baseline.get("ari", float("inf"))))
                 and not diag["collapse_warning"])
    summary = {
        "dataset": args.dataset_name,
        "method": METHOD_NAME,
        "method_dir": METHOD_DIR,
        "stage": stage,
        "seed": args.seed,
        "n_cells": int(input_bundle.data.shape[0]),
        "n_genes": int(input_bundle.data.shape[1]),
        "n_clusters": int(args.n_clusters),
        "runtime_seconds": float(time.time() - start),
        "embedding_path": str((output_dir / "embedding_final.npy").resolve()),
        "fixed_metrics": fixed_metrics,
        "diagnostics": diag,
        "baseline": baseline,
        "meets_screen_baseline_any": meets,
        "note": "Quick screen only; this file is not appended to 全benchmark结果.csv.",
    }
    save_json(summary, str(output_dir / "summary.json"))
    if fixed_metrics:
        append_screen_csv({
            "method_dir": METHOD_DIR,
            "dataset": args.dataset_name,
            "stage": stage,
            "seed": args.seed,
            "nmi": metric.get("nmi", np.nan),
            "ari": metric.get("ari", np.nan),
            "acc": metric.get("acc", np.nan),
            "collapse_warning": diag["collapse_warning"],
            "meets_screen_baseline_any": meets,
            "save_dir": str(output_dir),
        })
    print(f"Completed {METHOD_DIR}. Results saved to: {output_dir}")


if __name__ == "__main__":
    main()
