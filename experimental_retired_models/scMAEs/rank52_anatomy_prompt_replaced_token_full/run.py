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
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = next(parent for parent in [CURRENT_DIR, *CURRENT_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

from loss import AnatomyPromptLoss
from model import AnatomyPromptScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json

METHOD_NAME = "rank52_anatomy_prompt_replaced_token_full"
DISPLAY_NAME = "scMAE + anatomy prompt replaced token"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class AnatomyPromptDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, token_ids: np.ndarray, labels: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.token_ids = torch.as_tensor(token_ids, dtype=torch.long)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.token_ids[idx], self.labels[idx]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Independent anatomy-prompt replaced-token scMAE candidate.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    p.add_argument("--n_modules", type=int, default=32)
    p.add_argument("--n_token_bins", type=int, default=8)
    p.add_argument("--mask_prob_start", type=float, default=0.28)
    p.add_argument("--mask_prob_end", type=float, default=0.45)
    p.add_argument("--replace_prob", type=float, default=0.5)
    p.add_argument("--difficulty_momentum", type=float, default=0.95)
    p.add_argument("--masked_data_weight", type=float, default=0.75)
    p.add_argument("--mask_loss_weight", type=float, default=0.6)
    p.add_argument("--replaced_weight", type=float, default=0.15)
    p.add_argument("--token_weight", type=float, default=0.05)
    p.add_argument("--variance_weight", type=float, default=0.0)
    return p.parse_args()


def make_gene_modules(log_expr: np.ndarray, n_modules: int, seed: int) -> np.ndarray:
    features = np.stack(
        [
            np.mean(log_expr, axis=0),
            np.var(log_expr, axis=0),
            (log_expr <= 1e-8).mean(axis=0),
        ],
        axis=1,
    ).astype(np.float32)
    features = StandardScaler().fit_transform(features)
    n_modules = int(min(max(2, n_modules), log_expr.shape[1]))
    return KMeans(n_clusters=n_modules, n_init=20, random_state=seed).fit_predict(features).astype(np.int64)


def make_token_ids(log_expr: np.ndarray, n_bins: int) -> tuple[np.ndarray, np.ndarray]:
    qs = np.quantile(log_expr, np.linspace(0.0, 1.0, n_bins + 1)[1:-1], axis=0).astype(np.float32)
    token_ids = np.zeros(log_expr.shape, dtype=np.uint8)
    for q in qs:
        token_ids += (log_expr > q[None, :]).astype(np.uint8)
    return token_ids, qs


def curriculum_mask(
    x: torch.Tensor,
    difficulty: torch.Tensor,
    epoch: int,
    epochs: int,
    mask_start: float,
    mask_end: float,
    replace_prob: float,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, float]:
    progress = (epoch - 1) / max(1, epochs - 1)
    mask_ratio = float(mask_start * (1.0 - progress) + mask_end * progress)
    weights = difficulty.to(x.device).float().clamp_min(1e-4)
    weights = weights / weights.mean().clamp_min(1e-6)
    probs = (mask_ratio * weights).clamp(0.02, 0.85)
    mask = (torch.rand_like(x) < probs.view(1, -1)).float()
    swap = (torch.rand_like(x) < float(replace_prob)) & mask.bool()
    shuffled = x[torch.randperm(x.shape[0], device=x.device)]
    corrupted = torch.where(swap, shuffled, x)
    corrupted = corrupted.masked_fill(mask.bool() & ~swap, 0.0)
    return corrupted, mask, mask, mask_ratio


@torch.no_grad()
def extract_embedding(model: AnatomyPromptScMAE, loader: DataLoader, device: torch.device, module_ids: torch.Tensor) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    model.eval()
    emb, labels, mask_scores = [], [], []
    for _, x_cpu, _, _, y in loader:
        out = model(x_cpu.to(device), module_ids)
        emb.append(out["embedding"].detach().cpu().numpy())
        mask_scores.append(torch.sigmoid(out["mask_logits"]).mean(dim=1).detach().cpu().numpy())
        labels.append(y.numpy())
    return np.nan_to_num(np.concatenate(emb).astype(np.float32)), np.concatenate(labels).astype(np.int64), np.nan_to_num(np.concatenate(mask_scores).astype(np.float32))


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    idx = NearestNeighbors(n_neighbors=min(k + 1, embedding.shape[0])).fit(embedding).kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[idx] == labels[:, None]))


def diagnostics(embedding: np.ndarray, labels: np.ndarray, mask_scores: np.ndarray, n_clusters: int, seed: int, preds: np.ndarray | None, mask_rate: float, difficulty: np.ndarray, n_modules: int) -> dict:
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
        "mask_score_mean": float(np.mean(mask_scores)) if mask_scores.size else float("nan"),
        "gene_difficulty_cv": float(np.std(difficulty) / (np.mean(difficulty) + 1e-8)),
        "n_gene_modules": int(n_modules),
        "diagnostic_note": "No NeighborMix is used. Anatomy guidance is represented by unsupervised gene modules and module prompts.",
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
        out.append({"stage": stage, "method": method, "n_rows": int(len(group)), "n_datasets": int(group["dataset"].nunique()), "mean_ari": float(group["ari"].dropna().mean()), "mean_nmi": float(group["nmi"].dropna().mean()), "screen_meets_baseline_count": meets, "screen_noncollapse_count": noncollapse, "effective_screen_candidate": bool(meets >= 2 and noncollapse >= 2), "datasets": ";".join(sorted(str(x) for x in group["dataset"].dropna().unique()))})
    pd.DataFrame(out).to_csv(SUMMARY_RESULT_CSV, index=False)


def main() -> int:
    args = parse_args()
    if args.gpu in {0, 7} and not args.no_cuda:
        raise ValueError("GPU 0 and GPU 7 are forbidden. Choose GPU 1-6 or --no_cuda.")
    if args.smoke:
        args.epochs = min(args.epochs, 3)
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
    module_ids_np = make_gene_modules(log_expr, args.n_modules, args.seed)
    n_modules = int(module_ids_np.max() + 1)
    token_ids, token_quantiles = make_token_ids(log_expr, args.n_token_bins)
    gene_difficulty = torch.as_tensor((np.var(log_expr, axis=0) + 1e-6).astype(np.float32))
    gene_difficulty = gene_difficulty / gene_difficulty.mean().clamp_min(1e-6)
    np.save(save_dir / "gene_names.npy", target.gene_names.astype(str))
    np.save(save_dir / "gene_module_ids.npy", module_ids_np.astype(np.int64))
    np.save(save_dir / "token_quantiles.npy", token_quantiles.astype(np.float32))
    save_json(target.profile, str(save_dir / "dataset_profile.json"))
    save_json({**target.preprocess_config, "encoder_scale_input": bool(args.scale_input), "target_scale_input": False, "rank52_core": "gene-module anatomy prompts with corrupted-expression detection and gene rank token targets"}, str(save_dir / "preprocess_config.json"))

    dataset = AnatomyPromptDataset(encoder_data, log_expr, token_ids, labels)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=torch.Generator().manual_seed(args.seed))
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False)
    module_ids = torch.as_tensor(module_ids_np, dtype=torch.long, device=device)
    model = AnatomyPromptScMAE(encoder_data.shape[1], args.hidden_size, args.decoder_hidden, n_modules, args.n_token_bins, args.dropout).to(device)
    criterion = AnatomyPromptLoss(args.masked_data_weight, args.mask_loss_weight, args.replaced_weight, args.token_weight, args.variance_weight)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    keys = ["loss", "scmae_loss", "recon_loss", "mask_loss", "replaced_loss", "token_loss", "variance_loss", "effective_mask_rate", "mask_ratio"]
    history = {k: [] for k in keys}
    history["stage"] = stage
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")
    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {k: 0.0 for k in keys}
        nb = 0
        for _, x_cpu, log_cpu, tok_cpu, _ in train_loader:
            x, y, tok = x_cpu.to(device), log_cpu.to(device), tok_cpu.to(device)
            corrupted, mask, replaced_label, mask_ratio = curriculum_mask(x, gene_difficulty, epoch, args.epochs, args.mask_prob_start, args.mask_prob_end, args.replace_prob)
            out = model(corrupted, module_ids)
            loss, parts, gene_error = criterion(out, y, mask, replaced_label, tok)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            scaled_error = gene_error.cpu() / gene_error.mean().clamp_min(1e-6).cpu()
            gene_difficulty.mul_(args.difficulty_momentum).add_(scaled_error, alpha=1.0 - args.difficulty_momentum).clamp_(0.05, 10.0)
            for k, v in parts.items():
                sums[k] += v
            sums["effective_mask_rate"] += float(mask.mean().detach().cpu())
            sums["mask_ratio"] += mask_ratio
            nb += 1
        for k in keys:
            history[k].append(sums[k] / max(1, nb))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} scmae={history['scmae_loss'][-1]:.4f} replaced={history['replaced_loss'][-1]:.4f} token={history['token_loss'][-1]:.4f} mask={history['effective_mask_rate'][-1]:.3f}")

    embedding, labels_out, mask_scores = extract_embedding(model, full_loader, device, module_ids)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "edge_confidence.npy", np.full(labels_out.shape[0], np.nan, dtype=np.float32))
    np.save(save_dir / "gene_difficulty.npy", gene_difficulty.numpy().astype(np.float32))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "optimizer": opt.state_dict(), "args": vars(args), "gene_names": target.gene_names.astype(str), "module_ids": module_ids_np.astype(np.int64)}, save_dir / "model_checkpoint.pth")
    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage, "preprocessing": "scaled encoder + log-expression targets + module prompts/replaced-token objectives"})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    diag = diagnostics(embedding, labels_out, mask_scores, n_clusters, args.seed, preds, history["effective_mask_rate"][-1], gene_difficulty.numpy().astype(np.float32), n_modules)
    save_json(diag, str(save_dir / "diagnostics.json"))
    if not args.no_save_h5ad:
        enc.adata.obsm["X_anatomy_prompt_scmae"] = embedding
        enc.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(enc.adata)
        enc.adata.write_h5ad(save_dir / "adata_anatomy_prompt_scmae.h5ad", compression="gzip")
    baseline = BASELINES.get(dataset_name, {})
    nmi, ari, acc, f1 = metric(eval_result, "nmi"), metric(eval_result, "ari"), metric(eval_result, "acc"), metric(eval_result, "f1_macro")
    meets = bool((nmi is not None and nmi >= baseline.get("nmi", np.inf)) or (ari is not None and ari >= baseline.get("ari", np.inf)))
    summary = {"dataset": dataset_name, "method": DISPLAY_NAME, "method_dir": METHOD_NAME, "stage": stage, "seed": int(args.seed), "n_cells": int(encoder_data.shape[0]), "n_genes": int(encoder_data.shape[1]), "n_clusters": int(n_clusters), "runtime_seconds": float(time.time() - start), "embedding_path": str((save_dir / "embedding_final.npy").resolve()), "fixed_metrics": eval_result["fixed"] if eval_result else {}, "diagnostics": diag, "baseline": baseline, "meets_screen_baseline_any": meets, "note": "Screen result is candidate evidence only and is not appended to 全benchmark结果.csv."}
    save_json(summary, str(save_dir / "summary.json"))
    if not args.skip_eval:
        append_row({"stage": stage, "method": METHOD_NAME, "dataset": dataset_name, "seed": int(args.seed), "acc": acc, "nmi": nmi, "ari": ari, "f1_macro": f1, "baseline_nmi": baseline.get("nmi"), "baseline_ari": baseline.get("ari"), "meets_baseline_any": meets, "collapse_warning": diag["collapse_warning"], "embedding_variance": diag["embedding_variance"], "neighbor_purity_proxy": diag["neighbor_purity_proxy"], "mixed_cell_fraction": diag["mixed_cell_fraction"], "run_dir": str(save_dir.resolve())})
    print(f"Completed {METHOD_NAME}. Results saved to: {save_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
