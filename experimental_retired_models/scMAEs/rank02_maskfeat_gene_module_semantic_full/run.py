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

from loss import MaskFeatModuleLoss
from model import MaskFeatGeneModuleScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json


METHOD_NAME = "rank02_maskfeat_gene_module_semantic_full"
DISPLAY_NAME = "scMAE + MaskFeat gene-module semantic target"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class ModuleDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, semantic: np.ndarray, labels: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.semantic = torch.as_tensor(semantic, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.semantic[idx], self.labels[idx]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Independent-full MaskFeat-style gene-module semantic scMAE candidate.",
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
    parser.add_argument("--num_modules", type=int, default=32)
    parser.add_argument("--mask_prob", type=float, default=0.25)
    parser.add_argument("--module_mask_prob", type=float, default=0.08)
    parser.add_argument("--reconstruction_weight", type=float, default=1.0)
    parser.add_argument("--mask_weight", type=float, default=0.3)
    parser.add_argument("--semantic_weight", type=float, default=0.5)
    parser.add_argument("--huber_beta", type=float, default=1.0)
    return parser.parse_args()


def build_gene_modules(log_expr: np.ndarray, num_modules: int, seed: int) -> np.ndarray:
    mean = log_expr.mean(axis=0)
    std = log_expr.std(axis=0)
    zero_rate = (log_expr <= 1e-8).mean(axis=0)
    q25, q50, q75 = np.quantile(log_expr, [0.25, 0.5, 0.75], axis=0)
    dispersion = std / (mean + 1e-4)
    features = np.stack([mean, std, zero_rate, q25, q50, q75, dispersion], axis=1).astype(np.float32)
    features = np.nan_to_num(features)
    features = (features - features.mean(axis=0, keepdims=True)) / (features.std(axis=0, keepdims=True) + 1e-6)
    n_modules = min(int(num_modules), log_expr.shape[1])
    module_ids = KMeans(n_clusters=n_modules, n_init=20, random_state=seed).fit_predict(features)
    return module_ids.astype(np.int64)


def compute_module_features(log_expr: np.ndarray, module_ids: np.ndarray, num_modules: int) -> np.ndarray:
    n_cells, n_genes = log_expr.shape
    order = np.argsort(log_expr, axis=1)
    ranks = np.empty_like(order, dtype=np.float32)
    ranks[np.arange(n_cells)[:, None], order] = np.linspace(0.0, 1.0, n_genes, dtype=np.float32)
    out = np.zeros((n_cells, num_modules, 3), dtype=np.float32)
    for module in range(num_modules):
        cols = module_ids == module
        if not np.any(cols):
            continue
        block = log_expr[:, cols]
        out[:, module, 0] = block.mean(axis=1)
        out[:, module, 1] = (block <= 1e-8).mean(axis=1)
        out[:, module, 2] = ranks[:, cols].mean(axis=1)
    return out


@torch.no_grad()
def extract_embedding(model: MaskFeatGeneModuleScMAE, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    embeddings, labels = [], []
    for _, x, _, _, y in loader:
        embeddings.append(model.feature(x.to(device)).detach().cpu().numpy())
        labels.append(y.numpy())
    return np.nan_to_num(np.concatenate(embeddings).astype(np.float32)), np.concatenate(labels).astype(np.int64)


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    k_eff = min(k + 1, embedding.shape[0])
    nn = NearestNeighbors(n_neighbors=k_eff)
    nn.fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[idx] == labels[:, None]))


def compute_diagnostics(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int, preds: np.ndarray | None) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, float(counts.sum()))
    probs = frac[frac > 0]
    boundary_entropy = float(-(probs * np.log(probs)).sum() / max(np.log(max(2, n_clusters)), 1e-8))
    label_counts = np.bincount(labels.astype(np.int64)).astype(np.float64)
    rare_cut = max(5.0, 0.01 * float(labels.shape[0]))
    rare = set(np.where(label_counts <= rare_cut)[0].tolist())
    variance = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    mass_min = float(frac.min()) if frac.size else 0.0
    mass_max = float(frac.max()) if frac.size else 0.0
    return {
        "edge_survival": 1.0,
        "neighbor_purity_proxy": neighbor_purity(labels, embedding),
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": boundary_entropy,
        "rare_risk_fraction": float(np.mean([x in rare for x in labels])) if labels.size else 0.0,
        "embedding_variance": variance,
        "cluster_mass_min": mass_min,
        "cluster_mass_max": mass_max,
        "collapse_warning": bool((not np.isfinite(variance)) or variance < 1e-8 or mass_min < 0.001 or mass_max > 0.95),
        "diagnostic_note": "No NeighborMix graph is used; MaskFeat semantic target is module-level.",
    }


def metric(eval_result: dict | None, name: str) -> float | None:
    if not eval_result:
        return None
    value = eval_result.get("fixed", {}).get("kmeans_known_k", {}).get(name)
    return None if value is None else float(value)


def append_screen_row(row: dict) -> None:
    fields = [
        "stage", "method", "dataset", "seed", "acc", "nmi", "ari", "f1_macro",
        "baseline_nmi", "baseline_ari", "meets_baseline_any", "collapse_warning",
        "embedding_variance", "neighbor_purity_proxy", "mixed_cell_fraction", "run_dir",
    ]
    SINGLE_RESULT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with SINGLE_RESULT_CSV.open("a+", newline="", encoding="utf-8") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        handle.seek(0)
        exists = bool(handle.read(1))
        handle.seek(0, 2)
        writer = csv.DictWriter(handle, fieldnames=fields)
        if not exists:
            writer.writeheader()
        writer.writerow({key: row.get(key, "") for key in fields})
        fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    refresh_summary_csv()


def refresh_summary_csv() -> None:
    if not SINGLE_RESULT_CSV.exists():
        return
    df = pd.read_csv(SINGLE_RESULT_CSV)
    rows = []
    for (stage, method_name), group in df.groupby(["stage", "method"], dropna=False):
        screen = group[group["stage"] == "screen"]
        meets = int(screen.get("meets_baseline_any", pd.Series(dtype=bool)).fillna(False).astype(bool).sum())
        noncollapse = int((~screen.get("collapse_warning", pd.Series(dtype=bool)).fillna(True).astype(bool)).sum())
        rows.append({
            "stage": stage,
            "method": method_name,
            "n_rows": int(len(group)),
            "n_datasets": int(group["dataset"].nunique()),
            "mean_ari": float(group["ari"].dropna().mean()),
            "mean_nmi": float(group["nmi"].dropna().mean()),
            "screen_meets_baseline_count": meets,
            "screen_noncollapse_count": noncollapse,
            "effective_screen_candidate": bool(meets >= 2 and noncollapse >= 2),
            "datasets": ";".join(sorted(str(x) for x in group["dataset"].dropna().unique())),
        })
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
        args.data_path, args.input_mode, args.n_top_genes, args.target_sum, False, args.label_key, args.seed
    )
    if args.scale_input:
        encoder_bundle = family.load_scmae_dataset(
            args.data_path, args.input_mode, args.n_top_genes, args.target_sum, True, args.label_key, args.seed
        )
        if not np.array_equal(encoder_bundle.gene_names.astype(str), target_bundle.gene_names.astype(str)):
            raise ValueError("Scaled encoder genes and log-expression target genes differ.")
        encoder_data = encoder_bundle.data
    else:
        encoder_bundle = target_bundle
        encoder_data = target_bundle.data

    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    module_ids = build_gene_modules(log_expr, args.num_modules, args.seed)
    num_modules = int(module_ids.max()) + 1
    semantic = compute_module_features(log_expr, module_ids, num_modules)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))

    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(
        {
            **target_bundle.preprocess_config,
            "encoder_scale_input": bool(args.scale_input),
            "target_scale_input": False,
            "semantic_target": "module mean expression, zero-rate, rank-profile",
            "num_modules": int(num_modules),
        },
        str(save_dir / "preprocess_config.json"),
    )
    np.save(save_dir / "gene_names.npy", target_bundle.gene_names.astype(str))
    np.save(save_dir / "gene_module_ids.npy", module_ids.astype(np.int64))

    dataset = ModuleDataset(encoder_data, log_expr, semantic, labels)
    generator = torch.Generator()
    generator.manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    test_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = MaskFeatGeneModuleScMAE(
        num_genes=encoder_data.shape[1],
        num_modules=num_modules,
        module_feature_dim=3,
        hidden_size=args.hidden_size,
        dropout=args.dropout,
    ).to(device)
    module_ids_t = torch.as_tensor(module_ids, dtype=torch.long, device=device)
    criterion = MaskFeatModuleLoss(args.reconstruction_weight, args.mask_weight, args.semantic_weight, args.huber_beta)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)

    history = {"loss": [], "reconstruction_loss": [], "mask_loss": [], "semantic_loss": [], "effective_mask_rate": [], "stage": stage}
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} modules={num_modules} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {"loss": 0.0, "reconstruction_loss": 0.0, "mask_loss": 0.0, "semantic_loss": 0.0}
        mask_sum = 0.0
        n_batches = 0
        for _, x_enc_cpu, x_log_cpu, sem_cpu, _ in train_loader:
            x_enc = x_enc_cpu.to(device)
            x_log = x_log_cpu.to(device)
            sem = sem_cpu.to(device)
            corrupted, gene_mask, module_mask = model.sample_module_mask(
                x_enc, module_ids_t, args.mask_prob, args.module_mask_prob
            )
            outputs = model(corrupted)
            loss, parts = criterion(outputs, x_log, sem, gene_mask, module_mask)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
            for key in sums:
                sums[key] += parts[key]
            mask_sum += float(gene_mask.mean().detach().cpu())
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        history["effective_mask_rate"].append(mask_sum / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"rec={history['reconstruction_loss'][-1]:.4f} mask={history['mask_loss'][-1]:.4f} "
                f"sem={history['semantic_loss'][-1]:.4f} mask_rate={history['effective_mask_rate'][-1]:.4f}"
            )

    embedding, labels_out = extract_embedding(model, test_loader, device)
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
            "gene_module_ids": module_ids.astype(np.int64),
        },
        save_dir / "model_checkpoint.pth",
    )

    eval_result, preds = None, None
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
                "preprocessing": "scaled encoder + log_expr module semantic targets",
                "num_modules": int(num_modules),
                "mask_prob": float(args.mask_prob),
                "module_mask_prob": float(args.module_mask_prob),
            },
        )
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))

    diagnostics = compute_diagnostics(embedding, labels_out, n_clusters, args.seed, preds)
    save_json(diagnostics, str(save_dir / "diagnostics.json"))

    if not args.no_save_h5ad:
        encoder_bundle.adata.obsm["X_maskfeat_module_scmae"] = embedding
        encoder_bundle.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(encoder_bundle.adata)
        encoder_bundle.adata.write_h5ad(save_dir / "adata_maskfeat_module_scmae.h5ad", compression="gzip")

    baseline = BASELINES.get(dataset_name, {})
    nmi, ari, acc, f1_macro = metric(eval_result, "nmi"), metric(eval_result, "ari"), metric(eval_result, "acc"), metric(eval_result, "f1_macro")
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
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "diagnostics": diagnostics,
        "baseline": baseline,
        "meets_screen_baseline_any": meets,
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
                "meets_baseline_any": meets,
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

