#!/usr/bin/env python3
from __future__ import annotations

import argparse
import copy
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

from loss import Data2VecLatentLoss
from model import Data2VecScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json


METHOD_NAME = "rank19_data2vec_latent_teacher_full"
DISPLAY_NAME = "scMAE + data2vec latent teacher"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class ExprDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, labels: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.labels[idx]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Independent data2vec latent teacher scMAE candidate.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
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
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.05)
    parser.add_argument("--mask_prob", type=float, default=0.4)
    parser.add_argument("--teacher_jitter", type=float, default=0.01)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_weight", type=float, default=0.65)
    parser.add_argument("--distill_weight", type=float, default=0.08)
    parser.add_argument("--variance_weight", type=float, default=0.02)
    parser.add_argument("--ema_start", type=float, default=0.99)
    parser.add_argument("--ema_end", type=float, default=0.999)
    parser.add_argument("--warmup_epochs", type=int, default=10)
    return parser.parse_args()


@torch.no_grad()
def update_ema(teacher: torch.nn.Module, student: torch.nn.Module, decay: float) -> None:
    for tp, sp in zip(teacher.parameters(), student.parameters()):
        tp.data.mul_(float(decay)).add_(sp.data, alpha=1.0 - float(decay))
    for tb, sb in zip(teacher.buffers(), student.buffers()):
        tb.copy_(sb)


@torch.no_grad()
def extract_embedding(model: Data2VecScMAE, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    emb, labels = [], []
    for _, x, _, y in loader:
        emb.append(model.feature(x.to(device)).detach().cpu().numpy())
        labels.append(y.numpy())
    return np.nan_to_num(np.concatenate(emb).astype(np.float32)), np.concatenate(labels).astype(np.int64)


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    nn = NearestNeighbors(n_neighbors=min(k + 1, embedding.shape[0]))
    nn.fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[idx] == labels[:, None]))


def diagnostics(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int, preds: np.ndarray | None, teacher_std: float) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    probs = frac[frac > 0]
    label_counts = np.bincount(labels.astype(np.int64)).astype(np.float64)
    rare_cut = max(5.0, 0.01 * labels.shape[0])
    rare = set(np.where(label_counts <= rare_cut)[0].tolist())
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    mass_min = float(frac.min()) if frac.size else 0.0
    mass_max = float(frac.max()) if frac.size else 0.0
    return {
        "edge_survival": 1.0,
        "neighbor_purity_proxy": neighbor_purity(labels, embedding),
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(-(probs * np.log(probs)).sum() / max(np.log(max(2, n_clusters)), 1e-8)),
        "rare_risk_fraction": float(np.mean([x in rare for x in labels])) if labels.size else 0.0,
        "embedding_variance": var,
        "cluster_mass_min": mass_min,
        "cluster_mass_max": mass_max,
        "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or mass_min < 0.001 or mass_max > 0.95 or teacher_std < 1e-4),
        "teacher_target_std": float(teacher_std),
        "diagnostic_note": "No NeighborMix graph is used; EMA teacher latent targets are regressed from strong masked student views.",
    }


def metric(eval_result: dict | None, name: str) -> float | None:
    if not eval_result:
        return None
    value = eval_result.get("fixed", {}).get("kmeans_known_k", {}).get(name)
    return None if value is None else float(value)


def append_row(row: dict) -> None:
    fields = ["stage", "method", "dataset", "seed", "acc", "nmi", "ari", "f1_macro", "baseline_nmi", "baseline_ari", "meets_baseline_any", "collapse_warning", "embedding_variance", "neighbor_purity_proxy", "mixed_cell_fraction", "run_dir"]
    with SINGLE_RESULT_CSV.open("a+", newline="", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        exists = bool(handle.read(1))
        handle.seek(0, 2)
        writer = csv.DictWriter(handle, fieldnames=fields)
        if not exists:
            writer.writeheader()
        writer.writerow({k: row.get(k, "") for k in fields})
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    df = pd.read_csv(SINGLE_RESULT_CSV)
    out = []
    for (stage, method_name), group in df.groupby(["stage", "method"], dropna=False):
        screen = group[group["stage"] == "screen"]
        meets = int(screen.get("meets_baseline_any", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
        noncollapse = int((~screen.get("collapse_warning", pd.Series(dtype=bool)).fillna(True).astype(bool)).sum())
        out.append({"stage": stage, "method": method_name, "n_rows": int(len(group)), "n_datasets": int(group["dataset"].nunique()), "mean_ari": float(group["ari"].dropna().mean()), "mean_nmi": float(group["nmi"].dropna().mean()), "screen_meets_baseline_count": meets, "screen_noncollapse_count": noncollapse, "effective_screen_candidate": bool(meets >= 2 and noncollapse >= 2), "datasets": ";".join(sorted(str(x) for x in group["dataset"].dropna().unique()))})
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

    target_bundle = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed)
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed)
        if not np.array_equal(encoder_bundle.gene_names.astype(str), target_bundle.gene_names.astype(str)):
            raise ValueError("Scaled encoder genes and log-expression target genes differ.")
        encoder_data = encoder_bundle.data
    else:
        encoder_bundle = target_bundle
        encoder_data = target_bundle.data
    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json({**target_bundle.preprocess_config, "encoder_scale_input": bool(args.scale_input), "target_scale_input": False, "teacher_target": "EMA teacher normalized latent target from weak/full view"}, str(save_dir / "preprocess_config.json"))
    np.save(save_dir / "gene_names.npy", target_bundle.gene_names.astype(str))

    dataset = ExprDataset(encoder_data, log_expr, labels)
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)
    student = Data2VecScMAE(encoder_data.shape[1], args.hidden_size, args.dropout).to(device)
    teacher = copy.deepcopy(student).to(device)
    teacher.eval()
    for param in teacher.parameters():
        param.requires_grad_(False)
    criterion = Data2VecLatentLoss(args.masked_data_weight, args.mask_weight, args.distill_weight, args.variance_weight)
    optimizer = torch.optim.Adam(student.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {"loss": [], "scmae_loss": [], "reconstruction_loss": [], "mask_loss": [], "distill_loss": [], "variance_loss": [], "effective_mask_rate": [], "teacher_target_std": [], "ema_decay": [], "distill_scale": [], "stage": stage}
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")

    total_epochs = max(1, args.epochs)
    for epoch in range(1, total_epochs + 1):
        student.train()
        sums = {k: 0.0 for k in ["loss", "scmae_loss", "reconstruction_loss", "mask_loss", "distill_loss", "variance_loss"]}
        mask_sum = 0.0
        target_std_sum = 0.0
        n_batches = 0
        progress = epoch / max(1, total_epochs)
        ema_decay = args.ema_end - (args.ema_end - args.ema_start) * (1.0 - progress) ** 2
        distill_scale = min(1.0, epoch / max(1, args.warmup_epochs))
        for _, x_cpu, log_cpu, _ in train_loader:
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            strong, mask = student.mask_view(x, args.mask_prob, 0.0)
            weak, _ = student.mask_view(x, max(0.02, args.mask_prob * 0.1), args.teacher_jitter)
            out = student(strong)
            with torch.no_grad():
                teacher_out = teacher(weak)
                teacher_target = teacher_out["teacher_layer"]
                teacher_std = float(teacher_target.std(dim=0).mean().detach().cpu())
            loss, parts = criterion(out, teacher_target, target, mask, distill_scale)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(student.parameters(), 5.0)
            optimizer.step()
            update_ema(teacher, student, ema_decay)
            for key in sums:
                sums[key] += parts[key]
            mask_sum += float(mask.mean().detach().cpu())
            target_std_sum += teacher_std
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        history["effective_mask_rate"].append(mask_sum / max(1, n_batches))
        history["teacher_target_std"].append(target_std_sum / max(1, n_batches))
        history["ema_decay"].append(float(ema_decay))
        history["distill_scale"].append(float(distill_scale))
        if epoch == 1 or epoch == total_epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{total_epochs} loss={history['loss'][-1]:.4f} scmae={history['scmae_loss'][-1]:.4f} distill={history['distill_loss'][-1]:.4f} var={history['variance_loss'][-1]:.4f} tstd={history['teacher_target_std'][-1]:.4f} mask_rate={history['effective_mask_rate'][-1]:.4f}")

    embedding, labels_out = extract_embedding(teacher, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"student": student.state_dict(), "teacher": teacher.state_dict(), "optimizer": optimizer.state_dict(), "args": vars(args), "gene_names": target_bundle.gene_names.astype(str)}, save_dir / "model_checkpoint.pth")

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage, "preprocessing": "scaled encoder + data2vec EMA latent target"})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    diag = diagnostics(embedding, labels_out, n_clusters, args.seed, preds, history["teacher_target_std"][-1])
    save_json(diag, str(save_dir / "diagnostics.json"))
    if not args.no_save_h5ad:
        encoder_bundle.adata.obsm["X_data2vec_scmae"] = embedding
        encoder_bundle.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(encoder_bundle.adata)
        encoder_bundle.adata.write_h5ad(save_dir / "adata_data2vec_scmae.h5ad", compression="gzip")

    baseline = BASELINES.get(dataset_name, {})
    nmi, ari, acc, f1 = metric(eval_result, "nmi"), metric(eval_result, "ari"), metric(eval_result, "acc"), metric(eval_result, "f1_macro")
    meets = bool((nmi is not None and nmi >= baseline.get("nmi", np.inf)) or (ari is not None and ari >= baseline.get("ari", np.inf)))
    summary = {"dataset": dataset_name, "method": DISPLAY_NAME, "method_dir": METHOD_NAME, "stage": stage, "seed": int(args.seed), "n_cells": int(encoder_data.shape[0]), "n_genes": int(encoder_data.shape[1]), "n_clusters": int(n_clusters), "runtime_seconds": float(time.time() - start), "embedding_path": str((save_dir / "embedding_final.npy").resolve()), "fixed_metrics": eval_result["fixed"] if eval_result is not None else {}, "diagnostics": diag, "baseline": baseline, "meets_screen_baseline_any": meets, "note": "Screen result is candidate evidence only and is not appended to 全benchmark结果.csv."}
    save_json(summary, str(save_dir / "summary.json"))
    if not args.skip_eval:
        append_row({"stage": stage, "method": METHOD_NAME, "dataset": dataset_name, "seed": int(args.seed), "acc": acc, "nmi": nmi, "ari": ari, "f1_macro": f1, "baseline_nmi": baseline.get("nmi"), "baseline_ari": baseline.get("ari"), "meets_baseline_any": meets, "collapse_warning": diag["collapse_warning"], "embedding_variance": diag["embedding_variance"], "neighbor_purity_proxy": diag["neighbor_purity_proxy"], "mixed_cell_fraction": diag["mixed_cell_fraction"], "run_dir": str(save_dir.resolve())})
    print(f"Completed {METHOD_NAME}. Results saved to: {save_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

