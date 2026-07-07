#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import fcntl
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = next(parent for parent in [CURRENT_DIR, *CURRENT_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

from loss import TrajectoryGuidedLoss
from model import TrajectoryGuidedScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json

METHOD_NAME = "rank51_trajectory_guided_token_selection_full"
DISPLAY_NAME = "scMAE + trajectory-guided token selection"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class TrajectoryDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, labels: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.labels[idx]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Independent trajectory-guided token selection scMAE candidate.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    p.add_argument("--no_save_h5ad", action="store_true")
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--hidden_size", type=int, default=128)
    p.add_argument("--decoder_hidden", type=int, default=128)
    p.add_argument("--dropout", type=float, default=0.0)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--warmup_epochs", type=int, default=5)
    p.add_argument("--trajectory_momentum", type=float, default=0.92)
    p.add_argument("--velocity_momentum", type=float, default=0.85)
    p.add_argument("--temperature", type=float, default=0.8)
    p.add_argument("--masked_data_weight", type=float, default=0.75)
    p.add_argument("--mask_loss_weight", type=float, default=0.7)
    p.add_argument("--trajectory_weight", type=float, default=0.04)
    p.add_argument("--sampler_target_weight", type=float, default=0.01)
    p.add_argument("--sampler_entropy_weight", type=float, default=0.002)
    p.add_argument("--variance_weight", type=float, default=0.0)
    return p.parse_args()


def normalize_vec(x: torch.Tensor) -> torch.Tensor:
    x = x.float()
    return (x - x.mean()) / x.std(unbiased=False).clamp_min(1e-6)


def build_gene_features(level: torch.Tensor, velocity: torch.Tensor, variance: torch.Tensor, zero_fraction: torch.Tensor, device: torch.device) -> torch.Tensor:
    return torch.stack(
        [
            normalize_vec(level),
            normalize_vec(velocity),
            normalize_vec(variance),
            normalize_vec(zero_fraction),
        ],
        dim=1,
    ).to(device)


def trajectory_target(level: torch.Tensor, velocity: torch.Tensor, variance: torch.Tensor) -> torch.Tensor:
    score = 0.60 * normalize_vec(level) + 0.25 * normalize_vec(velocity.abs()) + 0.15 * normalize_vec(variance)
    return score


def adaptive_mask(
    model: TrajectoryGuidedScMAE,
    x: torch.Tensor,
    gene_features: torch.Tensor,
    mask_prob: float,
    epoch: int,
    warmup_epochs: int,
    temperature: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor, bool]:
    if epoch <= warmup_epochs:
        probs = torch.full((x.shape[1],), float(mask_prob), device=x.device)
        logits = torch.zeros_like(probs)
        adaptive = False
    else:
        with torch.no_grad():
            logits = model.sampler_logits(gene_features, temperature)
            weights = torch.sigmoid(logits).clamp(1e-4, 1.0)
            weights = weights / weights.mean().clamp_min(1e-6)
            probs = (float(mask_prob) * weights).clamp(0.02, 0.85)
        adaptive = True
    mask = (torch.rand_like(x) < probs.view(1, -1)).float()
    corrupted = x.masked_fill(mask.bool(), 0.0)
    return corrupted, mask, probs, logits, adaptive


@torch.no_grad()
def extract_embedding(model: TrajectoryGuidedScMAE, loader: DataLoader, device: torch.device, gene_features: torch.Tensor) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model.eval()
    emb, labels, mask_scores = [], [], []
    for _, x_cpu, _, y in loader:
        out = model(x_cpu.to(device), gene_features=gene_features)
        emb.append(out["embedding"].detach().cpu().numpy())
        mask_scores.append(torch.sigmoid(out["mask_logits"]).mean(dim=1).detach().cpu().numpy())
        labels.append(y.numpy())
    return (
        np.nan_to_num(np.concatenate(emb).astype(np.float32)),
        np.concatenate(labels).astype(np.int64),
        np.nan_to_num(np.concatenate(mask_scores).astype(np.float32)),
    )


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    idx = NearestNeighbors(n_neighbors=min(k + 1, embedding.shape[0])).fit(embedding).kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[idx] == labels[:, None]))


def diagnostics(
    embedding: np.ndarray,
    labels: np.ndarray,
    mask_scores: np.ndarray,
    n_clusters: int,
    seed: int,
    preds: np.ndarray | None,
    mask_rate: float,
    mask_probs: np.ndarray,
    trajectory_level: np.ndarray,
    sampler_entropy: float,
) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    frac = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = frac / max(1.0, frac.sum())
    probs = frac[frac > 0]
    label_counts = np.bincount(labels.astype(np.int64)).astype(np.float64)
    rare = set(np.where(label_counts <= max(5.0, 0.01 * labels.shape[0]))[0].tolist())
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    return {
        "edge_survival": float(1.0 - mask_rate),
        "neighbor_purity_proxy": neighbor_purity(labels, embedding),
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(-(probs * np.log(probs)).sum() / max(np.log(max(2, n_clusters)), 1e-8)),
        "rare_risk_fraction": float(np.mean([x in rare for x in labels])) if labels.size else 0.0,
        "embedding_variance": var,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
        "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or (frac.size and frac.min() < 0.001) or (frac.size and frac.max() > 0.95)),
        "edge_confidence_mean": float("nan"),
        "edge_confidence_std": float("nan"),
        "mask_score_mean": float(np.mean(mask_scores)) if mask_scores.size else float("nan"),
        "trajectory_level_cv": float(np.std(trajectory_level) / (np.mean(trajectory_level) + 1e-8)),
        "mask_probability_mean": float(np.mean(mask_probs)),
        "mask_probability_std": float(np.std(mask_probs)),
        "sampler_entropy": float(sampler_entropy),
        "diagnostic_note": "No NeighborMix is used. Rank51 adapts video trajectory-guided token selection to per-gene reconstruction-error trajectories.",
    }


def metric(eval_result: dict | None, name: str) -> float | None:
    if not eval_result:
        return None
    v = eval_result.get("fixed", {}).get("kmeans_known_k", {}).get(name)
    return None if v is None else float(v)


def append_row(row: dict) -> None:
    fields = ["stage", "method", "dataset", "seed", "acc", "nmi", "ari", "f1_macro", "baseline_nmi", "baseline_ari", "meets_baseline_any", "collapse_warning", "embedding_variance", "neighbor_purity_proxy", "mixed_cell_fraction", "run_dir"]
    with SINGLE_RESULT_CSV.open("a+", newline="", encoding="utf-8") as h:
        fcntl.flock(h.fileno(), fcntl.LOCK_EX)
        h.seek(0)
        exists = bool(h.read(1))
        h.seek(0, 2)
        w = csv.DictWriter(h, fieldnames=fields)
        if not exists:
            w.writeheader()
        w.writerow({k: row.get(k, "") for k in fields})
        fcntl.flock(h.fileno(), fcntl.LOCK_UN)
    df = pd.read_csv(SINGLE_RESULT_CSV)
    out = []
    for (stage, method), group in df.groupby(["stage", "method"], dropna=False):
        screen = group[group["stage"] == "screen"]
        meets = int(screen.get("meets_baseline_any", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
        noncollapse = int((~screen.get("collapse_warning", pd.Series(dtype=bool)).fillna(True).astype(bool)).sum())
        out.append({
            "stage": stage,
            "method": method,
            "n_rows": int(len(group)),
            "n_datasets": int(group["dataset"].nunique()),
            "mean_ari": float(group["ari"].dropna().mean()),
            "mean_nmi": float(group["nmi"].dropna().mean()),
            "screen_meets_baseline_count": meets,
            "screen_noncollapse_count": noncollapse,
            "effective_screen_candidate": bool(meets >= 2 and noncollapse >= 2),
            "datasets": ";".join(sorted(str(x) for x in group["dataset"].dropna().unique())),
        })
    pd.DataFrame(out).to_csv(SUMMARY_RESULT_CSV, index=False)


def main() -> int:
    args = parse_args()
    if args.gpu in {0, 7} and not args.no_cuda:
        raise ValueError("GPU 0 and GPU 7 are forbidden. Choose GPU 1-6 or --no_cuda.")
    if args.smoke:
        args.epochs = min(args.epochs, 3)
        args.warmup_epochs = min(args.warmup_epochs, 1)
    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem
    stage = "smoke" if args.smoke else "screen"

    target = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed)
    if args.scale_input:
        enc = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed)
        encoder_data = np.asarray(enc.data, dtype=np.float32)
    else:
        enc = target
        encoder_data = np.asarray(target.data, dtype=np.float32)
    log_expr = np.asarray(target.data, dtype=np.float32)
    labels = np.asarray(target.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    gene_variance = torch.as_tensor(np.var(log_expr, axis=0).astype(np.float32) + 1e-6)
    zero_fraction = torch.as_tensor((log_expr <= 1e-8).mean(axis=0).astype(np.float32))
    trajectory_level = gene_variance / gene_variance.mean().clamp_min(1e-6)
    trajectory_velocity = torch.zeros_like(trajectory_level)
    last_mask_probs = np.full(log_expr.shape[1], args.mask_prob, dtype=np.float32)
    np.save(save_dir / "gene_names.npy", target.gene_names.astype(str))
    save_json(target.profile, str(save_dir / "dataset_profile.json"))
    save_json(
        {
            **target.preprocess_config,
            "encoder_scale_input": bool(args.scale_input),
            "target_scale_input": False,
            "rank51_core": "trajectory-guided gene-token masking from EMA reconstruction difficulty and velocity",
        },
        str(save_dir / "preprocess_config.json"),
    )

    dataset = TrajectoryDataset(encoder_data, log_expr, labels)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=torch.Generator().manual_seed(args.seed))
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False)
    model = TrajectoryGuidedScMAE(encoder_data.shape[1], args.hidden_size, args.decoder_hidden, args.dropout).to(device)
    criterion = TrajectoryGuidedLoss(
        args.masked_data_weight,
        args.mask_loss_weight,
        args.trajectory_weight,
        args.sampler_target_weight,
        args.sampler_entropy_weight,
        args.variance_weight,
    )
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    keys = ["loss", "scmae_loss", "recon_loss", "mask_loss", "trajectory_loss", "sampler_target_loss", "sampler_entropy_loss", "variance_loss", "effective_mask_rate", "adaptive_mask_fraction", "sampler_entropy"]
    history = {k: [] for k in keys}
    history["stage"] = stage
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        gene_features = build_gene_features(trajectory_level, trajectory_velocity, gene_variance, zero_fraction, device)
        traj_target = trajectory_target(trajectory_level, trajectory_velocity, gene_variance).to(device)
        sums = {k: 0.0 for k in keys}
        nb = 0
        for _, x_cpu, log_cpu, _ in train_loader:
            x = x_cpu.to(device)
            y = log_cpu.to(device)
            corrupted, mask, mask_probs, sampler_logits_detached, adaptive = adaptive_mask(model, x, gene_features, args.mask_prob, epoch, args.warmup_epochs, args.temperature)
            out = model(corrupted, gene_features=gene_features, temperature=args.temperature)
            loss, parts, gene_error = criterion(out, y, mask, traj_target, out.get("sampler_logits"))
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()

            old_level = trajectory_level.clone()
            scaled_error = gene_error.cpu() / gene_error.mean().clamp_min(1e-6).cpu()
            new_level = args.trajectory_momentum * trajectory_level + (1.0 - args.trajectory_momentum) * scaled_error
            delta = new_level - old_level
            trajectory_velocity = args.velocity_momentum * trajectory_velocity + (1.0 - args.velocity_momentum) * delta
            trajectory_level = new_level.clamp(0.05, 10.0)
            last_mask_probs = mask_probs.detach().cpu().numpy().astype(np.float32)

            probs = torch.softmax(out["sampler_logits"].detach(), dim=0).clamp_min(1e-8)
            entropy = float((-(probs * probs.log()).sum() / torch.log(torch.tensor(float(probs.numel()), device=probs.device).clamp_min(2.0))).cpu())
            for k, v in parts.items():
                sums[k] += v
            sums["effective_mask_rate"] += float(mask.mean().detach().cpu())
            sums["adaptive_mask_fraction"] += float(adaptive)
            sums["sampler_entropy"] += entropy
            nb += 1
        for k in keys:
            history[k].append(sums[k] / max(1, nb))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"scmae={history['scmae_loss'][-1]:.4f} mask={history['effective_mask_rate'][-1]:.3f} "
                f"adaptive={history['adaptive_mask_fraction'][-1]:.2f} entropy={history['sampler_entropy'][-1]:.3f}"
            )

    final_gene_features = build_gene_features(trajectory_level, trajectory_velocity, gene_variance, zero_fraction, device)
    embedding, labels_out, mask_scores = extract_embedding(model, full_loader, device, final_gene_features)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "edge_confidence.npy", np.full(labels_out.shape[0], np.nan, dtype=np.float32))
    np.save(save_dir / "trajectory_level.npy", trajectory_level.numpy().astype(np.float32))
    np.save(save_dir / "trajectory_velocity.npy", trajectory_velocity.numpy().astype(np.float32))
    np.save(save_dir / "mask_probabilities.npy", last_mask_probs.astype(np.float32))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "optimizer": opt.state_dict(), "args": vars(args), "gene_names": target.gene_names.astype(str)}, save_dir / "model_checkpoint.pth")

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            DISPLAY_NAME,
            args.seed,
            embedding,
            labels_out,
            n_clusters,
            {"variant": METHOD_NAME, "stage": stage, "preprocessing": "scaled encoder + log-expression targets + trajectory-guided mask token selection"},
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    diag = diagnostics(
        embedding,
        labels_out,
        mask_scores,
        n_clusters,
        args.seed,
        preds,
        history["effective_mask_rate"][-1],
        last_mask_probs,
        trajectory_level.numpy().astype(np.float32),
        history["sampler_entropy"][-1],
    )
    save_json(diag, str(save_dir / "diagnostics.json"))
    if not args.no_save_h5ad:
        enc.adata.obsm["X_trajectory_guided_scmae"] = embedding
        enc.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(enc.adata)
        enc.adata.write_h5ad(save_dir / "adata_trajectory_guided_scmae.h5ad", compression="gzip")

    baseline = BASELINES.get(dataset_name, {})
    nmi, ari, acc, f1 = metric(eval_result, "nmi"), metric(eval_result, "ari"), metric(eval_result, "acc"), metric(eval_result, "f1_macro")
    meets = bool((nmi is not None and nmi >= baseline.get("nmi", np.inf)) or (ari is not None and ari >= baseline.get("ari", np.inf)))
    summary = {
        "dataset": dataset_name,
        "method": DISPLAY_NAME,
        "method_dir": METHOD_NAME,
        "stage": stage,
        "seed": int(args.seed),
        "n_cells": int(encoder_data.shape[0]),
        "n_genes": int(encoder_data.shape[1]),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(time.time() - start),
        "embedding_path": str((save_dir / "embedding_final.npy").resolve()),
        "fixed_metrics": eval_result["fixed"] if eval_result else {},
        "diagnostics": diag,
        "baseline": baseline,
        "meets_screen_baseline_any": meets,
        "note": "Screen result is candidate evidence only and is not appended to 全benchmark结果.csv.",
    }
    save_json(summary, str(save_dir / "summary.json"))
    if not args.skip_eval:
        append_row({
            "stage": stage,
            "method": METHOD_NAME,
            "dataset": dataset_name,
            "seed": int(args.seed),
            "acc": acc,
            "nmi": nmi,
            "ari": ari,
            "f1_macro": f1,
            "baseline_nmi": baseline.get("nmi"),
            "baseline_ari": baseline.get("ari"),
            "meets_baseline_any": meets,
            "collapse_warning": diag["collapse_warning"],
            "embedding_variance": diag["embedding_variance"],
            "neighbor_purity_proxy": diag["neighbor_purity_proxy"],
            "mixed_cell_fraction": diag["mixed_cell_fraction"],
            "run_dir": str(save_dir.resolve()),
        })
    print(f"Completed {METHOD_NAME}. Results saved to: {save_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
