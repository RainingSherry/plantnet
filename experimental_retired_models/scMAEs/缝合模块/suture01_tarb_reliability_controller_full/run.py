#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import fcntl
import os
import sys
import time
from pathlib import Path

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_var, "1")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")

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

from loss import SutureTARBLoss
from model import SutureTARBScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json


METHOD_NAME = "suture01_tarb_reliability_controller_full"
DISPLAY_NAME = "scMAE + suture TARB reliability controller"
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
    parser = argparse.ArgumentParser(description="Independent suture01 TARB reliability-controller scMAE.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument("--module_weight", type=float, default=0.15)
    parser.add_argument("--masked_data_weight", type=float, default=0.75)
    parser.add_argument("--mask_weight", type=float, default=0.65)
    parser.add_argument("--cluster_weight", type=float, default=0.35)
    parser.add_argument("--consistency_weight", type=float, default=0.05)
    parser.add_argument("--confidence_threshold", type=float, default=0.35)
    parser.add_argument("--warmup_epochs", type=int, default=20)
    parser.add_argument("--target_update_interval", type=int, default=5)
    parser.add_argument("--neighbor_k", type=int, default=15)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--balance_weight", type=float, default=0.01)
    parser.add_argument("--variance_weight", type=float, default=0.02)
    parser.add_argument("--conservative_weight", type=float, default=0.01)
    parser.add_argument("--variance_floor", type=float, default=0.01)
    return parser.parse_args()


def _register_null_h5ad_reader() -> None:
    try:
        import h5py
        from anndata._io.specs.registry import _REGISTRY, IOSpec

        def _read_null(*args, **kwargs):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


def build_knn_graph(data_np: np.ndarray, k: int, pca_dim: int, seed: int) -> np.ndarray | None:
    if data_np.shape[0] <= 2:
        return None
    max_k = min(int(k), data_np.shape[0] - 1)
    if max_k <= 0:
        return None
    dim = min(int(pca_dim), min(data_np.shape) - 1)
    if dim >= 2:
        from sklearn.decomposition import PCA

        work = PCA(n_components=dim, random_state=seed).fit_transform(data_np.astype(np.float64))
    else:
        work = data_np.astype(np.float64)
    nn = NearestNeighbors(n_neighbors=max_k + 1, metric="cosine").fit(work)
    return nn.kneighbors(work, return_distance=False)[:, 1:].astype(np.int64)


def compute_reliability(embedding: np.ndarray, q: np.ndarray, neighbor_indices: np.ndarray | None, k: int = 15) -> tuple[np.ndarray, dict]:
    n_cells = int(embedding.shape[0])
    q = np.asarray(q, dtype=np.float64)
    assign = q.argmax(axis=1)
    confidence = q.max(axis=1)
    if neighbor_indices is None or neighbor_indices.shape[1] == 0:
        agree = np.ones(n_cells, dtype=np.float64)
        density = np.ones(n_cells, dtype=np.float64)
    else:
        kk = min(int(k), neighbor_indices.shape[1])
        nb = neighbor_indices[:, :kk]
        agree = (assign[nb] == assign[:, None]).mean(axis=1).astype(np.float64)
        mean_dist = np.empty(n_cells, dtype=np.float64)
        step = 4096
        for start in range(0, n_cells, step):
            end = min(start + step, n_cells)
            diff = embedding[start:end, None, :] - embedding[nb[start:end]]
            mean_dist[start:end] = np.linalg.norm(diff, axis=2).mean(axis=1)
        lo, hi = np.percentile(mean_dist, [5.0, 95.0])
        if hi - lo < 1e-8:
            density = np.ones(n_cells, dtype=np.float64)
        else:
            density = 1.0 - np.clip((mean_dist - lo) / (hi - lo), 0.0, 1.0)
    r = np.clip(agree * density, 0.0, 1.0).astype(np.float32)
    diag = {
        "reliability_mean": float(r.mean()),
        "reliability_min": float(r.min()),
        "reliability_max": float(r.max()),
        "core_fraction": float((r >= 0.5).mean()),
        "agree_mean": float(agree.mean()),
        "density_mean": float(density.mean()),
        "confidence_mean": float(confidence.mean()),
    }
    return r, diag


@torch.no_grad()
def extract_all(model: SutureTARBScMAE, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    model.eval()
    emb, labels, q_all, op_all = [], [], [], []
    for _, x, _, y in loader:
        out = model(x.to(device), reliability=None)
        emb.append(out["latent"].detach().cpu().numpy())
        q_all.append(out["cluster_q"].detach().cpu().numpy())
        op_all.append(out["operation_weights"].detach().cpu().numpy())
        labels.append(y.numpy())
    return (
        np.nan_to_num(np.concatenate(emb).astype(np.float32)),
        np.concatenate(labels).astype(np.int64),
        np.concatenate(q_all).astype(np.float32),
        np.concatenate(op_all).astype(np.float32),
    )


def initialize_centers(model: SutureTARBScMAE, embedding: np.ndarray, n_clusters: int, seed: int, device: torch.device) -> None:
    km = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed)
    km.fit(embedding)
    model.initialize_centers(torch.as_tensor(km.cluster_centers_, dtype=torch.float32, device=device))


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    nn = NearestNeighbors(n_neighbors=min(k + 1, embedding.shape[0])).fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[idx] == labels[:, None]))


def diagnostics(embedding: np.ndarray, labels: np.ndarray, q: np.ndarray, opw: np.ndarray, n_clusters: int, seed: int, preds: np.ndarray | None, rel_diag: dict) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    probs = frac[frac > 0]
    label_counts = np.bincount(labels.astype(np.int64)).astype(np.float64)
    rare_cut = max(5.0, 0.01 * labels.shape[0])
    rare = set(np.where(label_counts <= rare_cut)[0].tolist())
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    q = np.asarray(q, dtype=np.float64)
    q_entropy = float(np.mean(-(q * np.log(np.clip(q, 1e-8, 1.0))).sum(axis=1) / max(np.log(max(2, n_clusters)), 1e-8))) if q.size else 0.0
    mass_min = float(frac.min()) if frac.size else 0.0
    mass_max = float(frac.max()) if frac.size else 0.0
    op_mean = np.asarray(opw, dtype=np.float64).mean(axis=0).tolist() if opw.size else []
    op_entropy = float(np.mean(-(opw * np.log(np.clip(opw, 1e-8, 1.0))).sum(axis=1) / np.log(opw.shape[1]))) if opw.size else 0.0
    d = {
        "edge_survival": 1.0,
        "neighbor_purity_proxy": neighbor_purity(labels, embedding),
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": q_entropy,
        "rare_risk_fraction": float(np.mean([x in rare for x in labels])) if labels.size else 0.0,
        "embedding_variance": var,
        "cluster_mass_min": mass_min,
        "cluster_mass_max": mass_max,
        "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or mass_min < 0.001 or mass_max > 0.95),
        "module_gate_entropy": op_entropy,
        "operation_weight_mean": op_mean,
        "cluster_confidence_mean": float(np.max(q, axis=1).mean()) if q.size else 0.0,
        "diagnostic_note": "TARB-style latent controller; no NeighborMix cell mixing; reliability gates only risky auxiliary operations, not DEC.",
    }
    d.update(rel_diag)
    return d


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
        out.append(
            {
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
            }
        )
    pd.DataFrame(out).to_csv(SUMMARY_RESULT_CSV, index=False)


def main() -> int:
    _register_null_h5ad_reader()
    args = parse_args()
    if args.gpu in {0, 7} and not args.no_cuda:
        raise ValueError("GPU 0 and GPU 7 are forbidden. Choose GPU 1-6 or --no_cuda.")
    if args.smoke:
        args.epochs = min(args.epochs, 3)
        args.warmup_epochs = min(args.warmup_epochs, args.epochs + 1)
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
    save_json({**target_bundle.preprocess_config, "encoder_scale_input": bool(args.scale_input), "target_scale_input": False, "module": "TARB latent operation controller"}, str(save_dir / "preprocess_config.json"))

    nb_indices = build_knn_graph(encoder_data, args.neighbor_k, args.knn_pca_dim, args.seed)
    dataset = ExprDataset(encoder_data, log_expr, labels)
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)
    model = SutureTARBScMAE(encoder_data.shape[1], n_clusters, args.hidden_size, args.dropout, args.module_weight).to(device)
    criterion = SutureTARBLoss(
        args.masked_data_weight,
        args.mask_weight,
        args.cluster_weight,
        args.consistency_weight,
        args.confidence_threshold,
        args.balance_weight,
        args.variance_weight,
        args.conservative_weight,
        args.variance_floor,
    )
    optimizer = torch.optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    p_targets: np.ndarray | None = None
    cluster_initialized = False
    reliability = np.ones(encoder_data.shape[0], dtype=np.float32)
    rel_diag = {"reliability_mean": 1.0, "reliability_min": 1.0, "reliability_max": 1.0, "core_fraction": 1.0, "agree_mean": 1.0, "density_mean": 1.0, "confidence_mean": 0.0}
    history = {k: [] for k in ["loss", "scmae_loss", "reconstruction_loss", "mask_loss", "cluster_loss", "consistency_loss", "operation_balance_loss", "variance_loss", "conservative_loss", "confidence_fraction", "latent_std", "effective_mask_rate", "cluster_scale", "reliability_mean"]}
    history["stage"] = stage
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch > args.warmup_epochs and ((epoch - args.warmup_epochs - 1) % max(1, args.target_update_interval) == 0 or p_targets is None):
            emb, _, q_full, _ = extract_all(model, full_loader, device)
            if not cluster_initialized:
                initialize_centers(model, emb, n_clusters, args.seed, device)
                emb, _, q_full, _ = extract_all(model, full_loader, device)
                cluster_initialized = True
            p_targets = SutureTARBScMAE.target_distribution(torch.as_tensor(q_full, dtype=torch.float32)).numpy().astype(np.float32)
            reliability, rel_diag = compute_reliability(emb, q_full, nb_indices, k=args.neighbor_k)
        cluster_scale = 0.0 if epoch <= args.warmup_epochs else min(1.0, (epoch - args.warmup_epochs) / max(1, args.warmup_epochs))
        model.train()
        sums = {k: 0.0 for k in ["loss", "scmae_loss", "reconstruction_loss", "mask_loss", "cluster_loss", "consistency_loss", "operation_balance_loss", "variance_loss", "conservative_loss", "confidence_fraction", "latent_std"]}
        mask_sum = 0.0
        n_batches = 0
        for idx, x_cpu, log_cpu, _ in train_loader:
            idx_np = idx.numpy()
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            r_batch = torch.as_tensor(reliability[idx_np], dtype=torch.float32, device=device)
            strong, mask = model.random_mask(x, args.mask_prob)
            weak, _ = model.random_mask(x, max(0.05, args.mask_prob * 0.5))
            out = model(strong, reliability=r_batch)
            weak_out = model(weak, reliability=r_batch)
            p_batch = None if p_targets is None else torch.as_tensor(p_targets[idx_np], dtype=torch.float32, device=device)
            loss, parts = criterion(out, weak_out, target, mask, p_batch, cluster_scale, r_batch)
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
        history["cluster_scale"].append(cluster_scale)
        history["reliability_mean"].append(float(rel_diag["reliability_mean"]))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} scmae={history['scmae_loss'][-1]:.4f} cluster={history['cluster_loss'][-1]:.4f} opbal={history['operation_balance_loss'][-1]:.4f} r={history['reliability_mean'][-1]:.3f} conf={history['confidence_fraction'][-1]:.3f}")

    embedding, labels_out, q_out, opw_out = extract_all(model, full_loader, device)
    reliability, rel_diag = compute_reliability(embedding, q_out, nb_indices, k=args.neighbor_k)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "cluster_q.npy", q_out.astype(np.float32))
    np.save(save_dir / "operation_weights.npy", opw_out.astype(np.float32))
    np.save(save_dir / "reliability.npy", reliability.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "optimizer": optimizer.state_dict(), "args": vars(args), "gene_names": target_bundle.gene_names.astype(str)}, save_dir / "model_checkpoint.pth")

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage, "preprocessing": "scaled encoder + log target + TARB latent controller"})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    diag = diagnostics(embedding, labels_out, q_out, opw_out, n_clusters, args.seed, preds, rel_diag)
    save_json(diag, str(save_dir / "diagnostics.json"))
    if not args.no_save_h5ad:
        encoder_bundle.adata.obsm["X_suture01_tarb"] = embedding
        encoder_bundle.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(encoder_bundle.adata)
        encoder_bundle.adata.write_h5ad(save_dir / "adata_suture01_tarb.h5ad", compression="gzip")

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
        "fixed_metrics": eval_result["fixed"] if eval_result is not None else {},
        "diagnostics": diag,
        "baseline": baseline,
        "meets_screen_baseline_any": meets,
        "note": "Screen result is candidate evidence only and is not appended to 全benchmark结果.csv.",
    }
    save_json(summary, str(save_dir / "summary.json"))
    if not args.skip_eval:
        append_row({"stage": stage, "method": METHOD_NAME, "dataset": dataset_name, "seed": int(args.seed), "acc": acc, "nmi": nmi, "ari": ari, "f1_macro": f1, "baseline_nmi": baseline.get("nmi"), "baseline_ari": baseline.get("ari"), "meets_baseline_any": meets, "collapse_warning": diag["collapse_warning"], "embedding_variance": diag["embedding_variance"], "neighbor_purity_proxy": diag["neighbor_purity_proxy"], "mixed_cell_fraction": diag["mixed_cell_fraction"], "run_dir": str(save_dir.resolve())})
    print(f"Completed {METHOD_NAME}. Results saved to: {save_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
