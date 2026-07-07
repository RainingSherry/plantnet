#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import fcntl
import json
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
ROOT = next(
    parent
    for parent in [CURRENT_DIR, *CURRENT_DIR.parents]
    if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists()
)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

from loss import ScDiVaLoss
from model import ScDiVaAbsorbingScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json


METHOD_NAME = "rank01_scdiva_absorbing_diffusion_full"
DISPLAY_NAME = "scMAE + ScDiVa absorbing token target"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class ScDiVaDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, tokens: np.ndarray, labels: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.tokens = torch.as_tensor(tokens, dtype=torch.long)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.tokens[idx], self.labels[idx]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Independent-full ScDiVa-inspired absorbing-mask scMAE candidate.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
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
    parser.add_argument("--dropout", type=float, default=0.0)
    parser.add_argument("--token_bins", type=int, default=8)
    parser.add_argument("--t_min", type=float, default=0.15)
    parser.add_argument("--t_max", type=float, default=0.65)
    parser.add_argument("--expression_weight", type=float, default=1.0)
    parser.add_argument("--mask_weight", type=float, default=0.3)
    parser.add_argument("--token_weight", type=float, default=0.5)
    parser.add_argument("--huber_beta", type=float, default=1.0)
    return parser.parse_args()


def compute_quantile_tokens(log_expr: np.ndarray, token_bins: int) -> tuple[np.ndarray, np.ndarray]:
    token_bins = int(token_bins)
    n_genes = log_expr.shape[1]
    edges = np.zeros((n_genes, max(0, token_bins - 1)), dtype=np.float32)
    tokens = np.zeros(log_expr.shape, dtype=np.int64)
    quantiles = np.linspace(0.0, 1.0, token_bins + 1, dtype=np.float32)[1:-1]
    for gene_idx in range(n_genes):
        values = np.asarray(log_expr[:, gene_idx], dtype=np.float32)
        finite = values[np.isfinite(values)]
        if finite.size == 0 or token_bins <= 1:
            continue
        gene_edges = np.quantile(finite, quantiles).astype(np.float32)
        gene_edges = np.maximum.accumulate(gene_edges)
        edges[gene_idx] = gene_edges
        tokens[:, gene_idx] = np.searchsorted(gene_edges, values, side="right")
    return tokens.astype(np.int64), edges


def extract_embedding(model: ScDiVaAbsorbingScMAE, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    embeddings = []
    labels = []
    with torch.no_grad():
        for _, x, _, _, y in loader:
            z = model.feature(x.to(device))
            embeddings.append(z.detach().cpu().numpy())
            labels.append(y.numpy())
    emb = np.nan_to_num(np.concatenate(embeddings, axis=0).astype(np.float32))
    labels_np = np.concatenate(labels, axis=0).astype(np.int64)
    return emb, labels_np


def entropy(values: np.ndarray) -> float:
    values = values.astype(np.float64)
    total = float(values.sum())
    if total <= 0.0:
        return 0.0
    probs = values / total
    probs = probs[probs > 0.0]
    return float(-(probs * np.log(probs)).sum())


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    k_eff = min(k + 1, embedding.shape[0])
    nn = NearestNeighbors(n_neighbors=k_eff, metric="euclidean")
    nn.fit(embedding)
    indices = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[indices] == labels[:, None]))


def compute_diagnostics(
    embedding: np.ndarray,
    labels: np.ndarray,
    n_clusters: int,
    seed: int,
    preds: np.ndarray | None,
) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    pred_counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    pred_frac = pred_counts / max(1.0, float(pred_counts.sum()))
    label_counts = np.bincount(labels.astype(np.int64)).astype(np.float64)
    rare_cut = max(5.0, 0.01 * float(labels.shape[0]))
    rare_labels = set(np.where(label_counts <= rare_cut)[0].tolist())
    rare_risk = float(np.mean([label in rare_labels for label in labels])) if labels.size else 0.0
    cluster_entropy = entropy(pred_counts) / max(np.log(max(2, n_clusters)), 1e-8)
    embedding_variance = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    cluster_mass_min = float(pred_frac.min()) if pred_frac.size else 0.0
    cluster_mass_max = float(pred_frac.max()) if pred_frac.size else 0.0
    collapse_warning = bool(
        (not np.isfinite(embedding_variance))
        or embedding_variance < 1e-8
        or cluster_mass_min < 0.001
        or cluster_mass_max > 0.95
    )
    return {
        "edge_survival": 1.0,
        "neighbor_purity_proxy": neighbor_purity(labels, embedding),
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(cluster_entropy),
        "rare_risk_fraction": rare_risk,
        "embedding_variance": embedding_variance,
        "cluster_mass_min": cluster_mass_min,
        "cluster_mass_max": cluster_mass_max,
        "collapse_warning": collapse_warning,
        "diagnostic_note": "No NeighborMix graph is used; edge_survival is the non-pruned diagnostic baseline.",
    }


def first_metric(eval_result: dict | None, metric: str) -> float | None:
    if not eval_result:
        return None
    fixed = eval_result.get("fixed", {}).get("kmeans_known_k", {})
    value = fixed.get(metric)
    return None if value is None else float(value)


def append_screen_row(row: dict) -> None:
    SINGLE_RESULT_CSV.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "stage",
        "method",
        "dataset",
        "seed",
        "acc",
        "nmi",
        "ari",
        "f1_macro",
        "baseline_nmi",
        "baseline_ari",
        "meets_baseline_any",
        "collapse_warning",
        "embedding_variance",
        "neighbor_purity_proxy",
        "mixed_cell_fraction",
        "run_dir",
    ]
    with SINGLE_RESULT_CSV.open("a+", newline="", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        exists = bool(handle.read(1))
        handle.seek(0, 2)
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if not exists:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in fieldnames})
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    refresh_summary_csv()


def refresh_summary_csv() -> None:
    if not SINGLE_RESULT_CSV.exists():
        return
    df = pd.read_csv(SINGLE_RESULT_CSV)
    if df.empty:
        return
    rows = []
    for (stage, method), group in df.groupby(["stage", "method"], dropna=False):
        screen = group[group["stage"] == "screen"]
        count_meets = int(screen.get("meets_baseline_any", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
        count_noncollapse = int((~screen.get("collapse_warning", pd.Series(dtype=bool)).fillna(True).astype(bool)).sum())
        rows.append(
            {
                "stage": stage,
                "method": method,
                "n_rows": int(len(group)),
                "n_datasets": int(group["dataset"].nunique()),
                "mean_ari": float(group["ari"].dropna().mean()) if "ari" in group else np.nan,
                "mean_nmi": float(group["nmi"].dropna().mean()) if "nmi" in group else np.nan,
                "screen_meets_baseline_count": count_meets,
                "screen_noncollapse_count": count_noncollapse,
                "effective_screen_candidate": bool(count_meets >= 2 and count_noncollapse >= 2),
                "datasets": ";".join(sorted(str(x) for x in group["dataset"].dropna().unique())),
            }
        )
    pd.DataFrame(rows).to_csv(SUMMARY_RESULT_CSV, index=False)


def main() -> int:
    args = parse_args()
    if args.gpu in {0, 7} and not args.no_cuda:
        raise ValueError("GPU 0 and GPU 7 are forbidden. Choose GPU 1-6 or --no_cuda.")
    if args.smoke:
        args.epochs = min(int(args.epochs), 3)

    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = family.get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem
    stage = "smoke" if args.smoke else "screen"

    target_bundle = family.load_scmae_dataset(
        file_path=args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        scale_input=False,
        label_key=args.label_key,
        seed=args.seed,
    )
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(
            file_path=args.data_path,
            input_mode=args.input_mode,
            n_top_genes=args.n_top_genes,
            target_sum=args.target_sum,
            scale_input=True,
            label_key=args.label_key,
            seed=args.seed,
        )
        encoder_data = encoder_bundle.data
        if not np.array_equal(encoder_bundle.gene_names.astype(str), target_bundle.gene_names.astype(str)):
            raise ValueError("Scaled encoder genes and log-expression target genes differ.")
    else:
        encoder_bundle = target_bundle
        encoder_data = target_bundle.data

    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    tokens, token_edges = compute_quantile_tokens(log_expr, args.token_bins)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))

    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(
        {
            **target_bundle.preprocess_config,
            "encoder_scale_input": bool(args.scale_input),
            "target_scale_input": False,
            "token_bins": int(args.token_bins),
            "token_source": "log_expr gene-specific quantile bins",
        },
        str(save_dir / "preprocess_config.json"),
    )
    np.save(save_dir / "gene_names.npy", target_bundle.gene_names.astype(str))
    np.save(save_dir / "token_quantile_edges.npy", token_edges)

    dataset = ScDiVaDataset(encoder_data, log_expr, tokens, labels)
    generator = torch.Generator()
    generator.manual_seed(args.seed)
    train_loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=True,
        generator=generator,
    )
    test_loader = DataLoader(
        dataset,
        batch_size=max(args.batch_size * 4, 512),
        shuffle=False,
        drop_last=False,
    )

    model = ScDiVaAbsorbingScMAE(
        num_genes=encoder_data.shape[1],
        hidden_size=args.hidden_size,
        token_bins=args.token_bins,
        dropout=args.dropout,
    ).to(device)
    criterion = ScDiVaLoss(
        expression_weight=args.expression_weight,
        mask_weight=args.mask_weight,
        token_weight=args.token_weight,
        huber_beta=args.huber_beta,
    )
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    history = {
        "loss": [],
        "expression_loss": [],
        "mask_loss": [],
        "token_loss": [],
        "effective_mask_rate": [],
        "t_mean": [],
        "stage": stage,
    }
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {"loss": 0.0, "expression_loss": 0.0, "mask_loss": 0.0, "token_loss": 0.0}
        mask_sum = 0.0
        t_sum = 0.0
        n_batches = 0
        for _, x_enc_cpu, x_log_cpu, token_cpu, _ in train_loader:
            x_enc = x_enc_cpu.to(device)
            x_log = x_log_cpu.to(device)
            token_target = token_cpu.to(device)
            corrupted, mask, t = model.sample_absorbing_mask(x_enc, args.t_min, args.t_max)
            outputs = model(corrupted, t)
            loss, parts = criterion(outputs, x_log, token_target, mask)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()

            for key in sums:
                sums[key] += parts[key]
            mask_sum += float(mask.mean().detach().cpu())
            t_sum += float(t.mean().detach().cpu())
            n_batches += 1

        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        history["effective_mask_rate"].append(mask_sum / max(1, n_batches))
        history["t_mean"].append(t_sum / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} "
                f"loss={history['loss'][-1]:.4f} "
                f"expr={history['expression_loss'][-1]:.4f} "
                f"mask={history['mask_loss'][-1]:.4f} "
                f"tok={history['token_loss'][-1]:.4f} "
                f"mask_rate={history['effective_mask_rate'][-1]:.4f}"
            )

    embedding, labels_out = extract_embedding(model, test_loader, device)
    if not np.isfinite(embedding).all():
        raise FloatingPointError("Embedding contains NaN or Inf after nan_to_num safeguard unexpectedly failed.")
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save(
        {
            "model": model.state_dict(),
            "optimizer": optimizer.state_dict(),
            "args": vars(args),
            "gene_names": target_bundle.gene_names.astype(str),
            "token_quantile_edges": token_edges,
        },
        save_dir / "model_checkpoint.pth",
    )

    eval_result = None
    preds = None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(
            output_dir=save_dir,
            dataset=dataset_name,
            method=DISPLAY_NAME,
            seed=args.seed,
            embedding=embedding,
            labels=labels_out,
            n_clusters=n_clusters,
            extra={
                "variant": METHOD_NAME,
                "stage": stage,
                "preprocessing": "scaled encoder + log_expr targets",
                "token_bins": int(args.token_bins),
                "t_min": float(args.t_min),
                "t_max": float(args.t_max),
            },
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    diagnostics = compute_diagnostics(embedding, labels_out, n_clusters, args.seed, preds)
    save_json(diagnostics, str(save_dir / "diagnostics.json"))

    if not args.no_save_h5ad:
        encoder_bundle.adata.obsm["X_scdiva_absorbing_scmae"] = embedding
        encoder_bundle.adata.uns["rank01_scdiva_absorbing_diffusion_full"] = {
            "method": DISPLAY_NAME,
            "variant": METHOD_NAME,
            "stage": stage,
        }
        sanitize_anndata_for_write(encoder_bundle.adata)
        encoder_bundle.adata.write_h5ad(save_dir / "adata_scdiva_absorbing_scmae.h5ad", compression="gzip")

    runtime = time.time() - start
    baseline = BASELINES.get(dataset_name, {})
    nmi = first_metric(eval_result, "nmi")
    ari = first_metric(eval_result, "ari")
    acc = first_metric(eval_result, "acc")
    f1_macro = first_metric(eval_result, "f1_macro")
    meets_baseline = bool(
        (nmi is not None and "nmi" in baseline and nmi >= baseline["nmi"])
        or (ari is not None and "ari" in baseline and ari >= baseline["ari"])
    )

    summary = {
        "dataset": dataset_name,
        "method": DISPLAY_NAME,
        "method_dir": METHOD_NAME,
        "stage": stage,
        "seed": int(args.seed),
        "n_cells": int(encoder_data.shape[0]),
        "n_genes": int(encoder_data.shape[1]),
        "n_clusters": int(n_clusters),
        "runtime_seconds": float(runtime),
        "embedding_path": str((save_dir / "embedding_final.npy").resolve()),
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "diagnostics": diagnostics,
        "baseline": baseline,
        "meets_screen_baseline_any": meets_baseline,
        "note": "Screen result is candidate evidence only and is not appended to 全benchmark结果.csv.",
    }
    save_json(summary, str(save_dir / "summary.json"))

    if not args.skip_eval:
        append_screen_row(
            {
                "stage": stage,
                "method": METHOD_NAME,
                "dataset": dataset_name,
                "seed": int(args.seed),
                "acc": acc,
                "nmi": nmi,
                "ari": ari,
                "f1_macro": f1_macro,
                "baseline_nmi": baseline.get("nmi"),
                "baseline_ari": baseline.get("ari"),
                "meets_baseline_any": meets_baseline,
                "collapse_warning": diagnostics["collapse_warning"],
                "embedding_variance": diagnostics["embedding_variance"],
                "neighbor_purity_proxy": diagnostics["neighbor_purity_proxy"],
                "mixed_cell_fraction": diagnostics["mixed_cell_fraction"],
                "run_dir": str(save_dir.resolve()),
            }
        )

    print(f"Completed {METHOD_NAME}. Results saved to: {save_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
