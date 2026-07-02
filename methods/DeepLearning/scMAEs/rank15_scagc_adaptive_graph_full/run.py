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

from loss import ScAGCGraphLoss
from model import ScAGCAdaptiveGraphScMAE, confidence_gated_mix
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json


METHOD_NAME = "rank15_scagc_adaptive_graph_full"
DISPLAY_NAME = "scMAE + scAGC adaptive graph"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class GraphDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, labels: np.ndarray, neighbor_idx: np.ndarray, edge_conf: np.ndarray, rare_risk: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)
        self.neighbor_data = torch.as_tensor(encoder_data[neighbor_idx], dtype=torch.float32)
        self.edge_conf = torch.as_tensor(edge_conf, dtype=torch.float32)
        self.rare_risk = torch.as_tensor(rare_risk, dtype=torch.float32)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.neighbor_data[idx], self.log_expr[idx], self.labels[idx], self.edge_conf[idx], self.rare_risk[idx]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Independent scAGC adaptive graph scMAE candidate.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument("--mask_prob", type=float, default=0.35)
    parser.add_argument("--neighbor_k", type=int, default=15)
    parser.add_argument("--knn_features", type=int, default=64)
    parser.add_argument("--mix_alpha", type=float, default=0.85)
    parser.add_argument("--edge_threshold", type=float, default=0.55)
    parser.add_argument("--rare_threshold", type=float, default=0.85)
    parser.add_argument("--dropedge_keep_prob", type=float, default=0.75)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_weight", type=float, default=0.65)
    parser.add_argument("--edge_weight", type=float, default=0.05)
    parser.add_argument("--graph_weight", type=float, default=0.05)
    parser.add_argument("--dropedge_weight", type=float, default=0.02)
    return parser.parse_args()


def build_graph(encoder_data: np.ndarray, k: int, max_features: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(encoder_data[:, : min(max_features, encoder_data.shape[1])], dtype=np.float32)
    x = np.nan_to_num(x)
    n_neighbors = min(max(2, k + 1), x.shape[0])
    nn = NearestNeighbors(n_neighbors=n_neighbors, algorithm="auto", metric="euclidean")
    nn.fit(x)
    dist, idx = nn.kneighbors(x, return_distance=True)
    neigh = idx[:, 1:]
    d = dist[:, 1:]
    radius = d.mean(axis=1)
    sigma = float(np.median(d[:, -1])) + 1e-6
    conf_all = np.exp(-d / sigma)
    best = np.argmax(conf_all, axis=1)
    neighbor_idx = neigh[np.arange(neigh.shape[0]), best].astype(np.int64)
    edge_conf = conf_all[np.arange(conf_all.shape[0]), best].astype(np.float32)
    lo, hi = np.percentile(radius, [50.0, 95.0])
    rare_risk = np.clip((radius - lo) / max(float(hi - lo), 1e-6), 0.0, 1.0).astype(np.float32)
    return neighbor_idx, edge_conf, rare_risk


@torch.no_grad()
def extract_embedding(model: ScAGCAdaptiveGraphScMAE, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    emb, labels = [], []
    for _, x, _, _, y, _, _ in loader:
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


def diagnostics(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int, preds: np.ndarray | None, edge_conf: np.ndarray, rare_risk: np.ndarray, mixed_fraction: float, edge_threshold: float) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    probs = frac[frac > 0]
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    mass_min = float(frac.min()) if frac.size else 0.0
    mass_max = float(frac.max()) if frac.size else 0.0
    return {
        "edge_survival": float(np.mean(edge_conf >= edge_threshold)) if edge_conf.size else 0.0,
        "neighbor_purity_proxy": neighbor_purity(labels, embedding),
        "mixed_cell_fraction": float(mixed_fraction),
        "boundary_entropy": float(-(probs * np.log(probs)).sum() / max(np.log(max(2, n_clusters)), 1e-8)),
        "rare_risk_fraction": float(np.mean(rare_risk >= 0.75)) if rare_risk.size else 0.0,
        "embedding_variance": var,
        "cluster_mass_min": mass_min,
        "cluster_mass_max": mass_max,
        "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or mass_min < 0.001 or mass_max > 0.95),
        "edge_confidence_mean": float(np.mean(edge_conf)) if edge_conf.size else 0.0,
        "diagnostic_note": "Adaptive graph reliability with optional confidence-gated local mix; boundary/rare veto can set mix to zero.",
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
    neighbor_idx, edge_conf, rare_risk = build_graph(encoder_data, args.neighbor_k, args.knn_features)
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json({**target_bundle.preprocess_config, "encoder_scale_input": bool(args.scale_input), "target_scale_input": False, "graph_source": "KNN over encoder input", "mix": "optional confidence-gated local mix with rare-risk veto"}, str(save_dir / "preprocess_config.json"))
    np.save(save_dir / "gene_names.npy", target_bundle.gene_names.astype(str))
    np.save(save_dir / "neighbor_idx.npy", neighbor_idx)
    np.save(save_dir / "edge_confidence.npy", edge_conf)
    np.save(save_dir / "rare_risk.npy", rare_risk)

    dataset = GraphDataset(encoder_data, log_expr, labels, neighbor_idx, edge_conf, rare_risk)
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)
    model = ScAGCAdaptiveGraphScMAE(encoder_data.shape[1], args.hidden_size, args.dropout).to(device)
    criterion = ScAGCGraphLoss(args.masked_data_weight, args.mask_weight, args.edge_weight, args.graph_weight, args.dropedge_weight)
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {"loss": [], "scmae_loss": [], "reconstruction_loss": [], "mask_loss": [], "edge_loss": [], "graph_loss": [], "dropedge_loss": [], "mixed_cell_fraction": [], "effective_mask_rate": [], "stage": stage}
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {k: 0.0 for k in ["loss", "scmae_loss", "reconstruction_loss", "mask_loss", "edge_loss", "graph_loss", "dropedge_loss", "mixed_cell_fraction"]}
        mask_sum = 0.0
        n_batches = 0
        for _, x_cpu, nx_cpu, log_cpu, _, conf_cpu, rare_cpu in train_loader:
            x = x_cpu.to(device)
            nx = nx_cpu.to(device)
            target = log_cpu.to(device)
            conf = conf_cpu.to(device)
            rare = rare_cpu.to(device)
            mixed, mix_gate = confidence_gated_mix(x, nx, conf, rare, args.edge_threshold, args.rare_threshold, args.mix_alpha, args.dropedge_keep_prob)
            corrupted, mask = model.random_mask(mixed, args.mask_prob)
            out = model(corrupted)
            neighbor_out = model(nx)
            dropped, _ = confidence_gated_mix(x, nx, conf, rare, args.edge_threshold, args.rare_threshold, args.mix_alpha, max(0.0, args.dropedge_keep_prob * 0.5))
            dropped, _ = model.random_mask(dropped, args.mask_prob)
            dropped_out = model(dropped)
            neg = nx[torch.randperm(nx.shape[0], device=device)]
            pos_logits = model.edge_logits(out["latent"], neighbor_out["latent"])
            with torch.no_grad():
                neg_latent = model.feature(neg)
            neg_logits = model.edge_logits(out["latent"], neg_latent)
            loss, parts = criterion(out, neighbor_out, dropped_out, target, mask, pos_logits, neg_logits, conf, mix_gate)
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            for key in sums:
                sums[key] += parts[key]
            mask_sum += float(mask.mean().detach().cpu())
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        history["effective_mask_rate"].append(mask_sum / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} scmae={history['scmae_loss'][-1]:.4f} edge={history['edge_loss'][-1]:.4f} graph={history['graph_loss'][-1]:.4f} mix={history['mixed_cell_fraction'][-1]:.4f} mask_rate={history['effective_mask_rate'][-1]:.4f}")

    embedding, labels_out = extract_embedding(model, full_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "optimizer": optimizer.state_dict(), "args": vars(args), "gene_names": target_bundle.gene_names.astype(str)}, save_dir / "model_checkpoint.pth")

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage, "preprocessing": "scaled encoder + adaptive graph reliability"})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    mixed_fraction = float(np.mean(history["mixed_cell_fraction"])) if history["mixed_cell_fraction"] else 0.0
    diag = diagnostics(embedding, labels_out, n_clusters, args.seed, preds, edge_conf, rare_risk, mixed_fraction, args.edge_threshold)
    save_json(diag, str(save_dir / "diagnostics.json"))
    if not args.no_save_h5ad:
        encoder_bundle.adata.obsm["X_scagc_graph_scmae"] = embedding
        encoder_bundle.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(encoder_bundle.adata)
        encoder_bundle.adata.write_h5ad(save_dir / "adata_scagc_graph_scmae.h5ad", compression="gzip")

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

