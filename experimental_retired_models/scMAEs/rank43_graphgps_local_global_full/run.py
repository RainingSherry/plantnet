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
from sklearn.decomposition import TruncatedSVD
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

from loss import GraphGPSLoss
from model import GraphGPSScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json


METHOD_NAME = "rank43_graphgps_local_global_full"
DISPLAY_NAME = "scMAE + GraphGPS local-global adapter"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {
    "Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029},
    "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275},
    "Macosko": {"nmi": 0.657465, "ari": 0.494268},
}


class GraphGPSDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, labels: np.ndarray, knn_index: np.ndarray, pe: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)
        self.knn_index = torch.as_tensor(knn_index, dtype=torch.long)
        self.pe = torch.as_tensor(pe, dtype=torch.float32)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.labels[idx], self.knn_index[idx], self.pe[idx]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Independent GraphGPS local-global scMAE candidate.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    parser.add_argument("--decoder_hidden", type=int, default=128)
    parser.add_argument("--num_heads", type=int, default=4)
    parser.add_argument("--dropout", type=float, default=0.05)
    parser.add_argument("--neighbor_k", type=int, default=12)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--pe_dim", type=int, default=8)
    parser.add_argument("--mask_prob", type=float, default=0.4)
    parser.add_argument("--dropedge_prob", type=float, default=0.35)
    parser.add_argument("--recon_weight", type=float, default=1.0)
    parser.add_argument("--mask_weight", type=float, default=0.05)
    parser.add_argument("--pe_weight", type=float, default=0.05)
    parser.add_argument("--edge_weight", type=float, default=0.02)
    parser.add_argument("--consistency_weight", type=float, default=0.05)
    parser.add_argument("--variance_weight", type=float, default=0.01)
    return parser.parse_args()


def build_graph_context(data: np.ndarray, k: int, pca_dim: int, pe_dim: int, seed: int) -> tuple[np.ndarray, np.ndarray]:
    x = np.asarray(data, dtype=np.float32)
    dim = int(min(max(2, pca_dim), x.shape[1] - 1, x.shape[0] - 1))
    rep = TruncatedSVD(n_components=dim, random_state=seed).fit_transform(x) if dim < x.shape[1] else x
    nn = NearestNeighbors(n_neighbors=min(k + 1, x.shape[0]), metric="cosine")
    nn.fit(rep)
    dist, idx = nn.kneighbors(rep, return_distance=True)
    idx, dist = idx[:, 1:], dist[:, 1:]
    if idx.shape[1] < k:
        pad = k - idx.shape[1]
        idx = np.concatenate([idx, np.repeat(idx[:, -1:], pad, axis=1)], axis=1)
        dist = np.concatenate([dist, np.repeat(dist[:, -1:], pad, axis=1)], axis=1)

    pe = np.zeros((x.shape[0], pe_dim), dtype=np.float32)
    use = min(max(0, pe_dim - 2), rep.shape[1])
    if use:
        coords = rep[:, :use]
        coords = (coords - coords.mean(axis=0, keepdims=True)) / (coords.std(axis=0, keepdims=True) + 1e-6)
        pe[:, :use] = coords.astype(np.float32)
    density = 1.0 - np.clip(dist.mean(axis=1), 0.0, 2.0) / 2.0
    indeg = np.bincount(idx.reshape(-1), minlength=x.shape[0]).astype(np.float32)
    indeg = (indeg - indeg.mean()) / (indeg.std() + 1e-6)
    pe[:, -2] = density.astype(np.float32)
    pe[:, -1] = indeg.astype(np.float32)
    return idx.astype(np.int64), np.nan_to_num(pe.astype(np.float32))


def gather_neighbors(all_x: torch.Tensor, all_pe: torch.Tensor, neighbor_ids: torch.Tensor, device: torch.device, drop_prob: float) -> tuple[torch.Tensor, torch.Tensor, float]:
    neighbor_ids = neighbor_ids.to(all_x.device)
    neigh_x = all_x[neighbor_ids.reshape(-1)].view(neighbor_ids.shape[0], neighbor_ids.shape[1], all_x.shape[1])
    neigh_pe = all_pe[neighbor_ids.reshape(-1)].view(neighbor_ids.shape[0], neighbor_ids.shape[1], all_pe.shape[1])
    keep = torch.rand(neighbor_ids.shape, device=neighbor_ids.device) >= float(drop_prob)
    no_keep = ~keep.any(dim=1)
    if bool(no_keep.any()):
        cols = torch.randint(0, neighbor_ids.shape[1], (int(no_keep.sum()),), device=neighbor_ids.device)
        keep[no_keep, cols] = True
    scale = neighbor_ids.shape[1] / keep.sum(dim=1).clamp_min(1).view(-1, 1, 1)
    keep3 = keep.unsqueeze(-1)
    visible_x = neigh_x.masked_fill(~keep3, 0.0) * scale
    visible_pe = neigh_pe.masked_fill(~keep3, 0.0) * scale
    return visible_x.to(device), visible_pe.to(device), float(keep.float().mean().detach().cpu())


def node_base_embedding(model: GraphGPSScMAE, x: torch.Tensor, pe: torch.Tensor) -> torch.Tensor:
    return model.encode_base(x) + model.gps.pe_proj(pe)


@torch.no_grad()
def extract_embedding(model: GraphGPSScMAE, loader: DataLoader, all_x: torch.Tensor, all_pe: torch.Tensor, device: torch.device) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    model.eval()
    emb, labels, edge_conf, gate_vals = [], [], [], []
    for _, x_cpu, _, y, neigh_ids, pe_cpu in loader:
        x = x_cpu.to(device)
        pe = pe_cpu.to(device)
        neigh_x, neigh_pe, _ = gather_neighbors(all_x, all_pe, neigh_ids, device, 0.0)
        out = model(x, neigh_x, pe, neigh_pe)
        pos_ids = neigh_ids[:, 0].to(all_x.device)
        neg_ids = torch.randint(0, all_x.shape[0], (x.shape[0],), device=all_x.device)
        pos_z = node_base_embedding(model, all_x[pos_ids].to(device), all_pe[pos_ids].to(device))
        neg_z = node_base_embedding(model, all_x[neg_ids].to(device), all_pe[neg_ids].to(device))
        pos_logit = model.edge_logits(out["embedding"], pos_z)
        neg_logit = model.edge_logits(out["embedding"], neg_z)
        emb.append(out["embedding"].detach().cpu().numpy())
        edge_conf.append((torch.sigmoid(pos_logit) - torch.sigmoid(neg_logit)).detach().cpu().numpy())
        gate_vals.append(out["local_gate"].mean(dim=1).detach().cpu().numpy())
        labels.append(y.numpy())
    return (
        np.nan_to_num(np.concatenate(emb).astype(np.float32)),
        np.concatenate(labels).astype(np.int64),
        np.nan_to_num(np.concatenate(edge_conf).astype(np.float32)),
        np.nan_to_num(np.concatenate(gate_vals).astype(np.float32)),
    )


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    nn = NearestNeighbors(n_neighbors=min(k + 1, embedding.shape[0]))
    nn.fit(embedding)
    idx = nn.kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[idx] == labels[:, None]))


def diagnostics(embedding: np.ndarray, labels: np.ndarray, edge_conf: np.ndarray, gate_vals: np.ndarray, n_clusters: int, seed: int, preds: np.ndarray | None, edge_survival: float) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    counts = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = counts / max(1.0, counts.sum())
    probs = frac[frac > 0]
    label_counts = np.bincount(labels.astype(np.int64)).astype(np.float64)
    rare_cut = max(5.0, 0.01 * labels.shape[0])
    rare = set(np.where(label_counts <= rare_cut)[0].tolist())
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    finite_conf = edge_conf[np.isfinite(edge_conf)]
    return {
        "edge_survival": float(edge_survival),
        "neighbor_purity_proxy": neighbor_purity(labels, embedding),
        "mixed_cell_fraction": 0.0,
        "boundary_entropy": float(-(probs * np.log(probs)).sum() / max(np.log(max(2, n_clusters)), 1e-8)),
        "rare_risk_fraction": float(np.mean([x in rare for x in labels])) if labels.size else 0.0,
        "embedding_variance": var,
        "cluster_mass_min": float(frac.min()) if frac.size else 0.0,
        "cluster_mass_max": float(frac.max()) if frac.size else 0.0,
        "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or (frac.size and frac.min() < 0.001) or (frac.size and frac.max() > 0.95)),
        "edge_confidence_mean": float(finite_conf.mean()) if finite_conf.size else float("nan"),
        "edge_confidence_std": float(finite_conf.std()) if finite_conf.size else float("nan"),
        "local_gate_mean": float(np.mean(gate_vals)) if gate_vals.size else float("nan"),
        "diagnostic_note": "No NeighborMix is used. GraphGPS local-global adapter combines structural encodings, local KNN messages, and batch global attention.",
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
        encoder_data = np.asarray(encoder_bundle.data, dtype=np.float32)
    else:
        encoder_bundle = target_bundle
        encoder_data = np.asarray(target_bundle.data, dtype=np.float32)
    log_expr = np.asarray(target_bundle.data, dtype=np.float32)
    labels = np.asarray(target_bundle.labels, dtype=np.int64)
    n_clusters = int(args.n_clusters if args.n_clusters > 0 else len(np.unique(labels)))
    knn_index, pe = build_graph_context(encoder_data, args.neighbor_k, args.knn_pca_dim, args.pe_dim, args.seed)
    save_json(target_bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json({**target_bundle.preprocess_config, "encoder_scale_input": bool(args.scale_input), "target_scale_input": False, "graph_knn": f"TruncatedSVD({args.knn_pca_dim}) + cosine KNN(k={args.neighbor_k})", "graphgps_pe": "SVD coordinates + local density + KNN in-degree proxy", "graphgps_core": "shallow local KNN message plus batch global multi-head attention"}, str(save_dir / "preprocess_config.json"))
    np.save(save_dir / "gene_names.npy", target_bundle.gene_names.astype(str))
    np.save(save_dir / "knn_index.npy", knn_index)
    np.save(save_dir / "graph_pe.npy", pe)

    dataset = GraphGPSDataset(encoder_data, log_expr, labels, knn_index, pe)
    all_encoder = dataset.encoder_data
    all_pe = dataset.pe
    generator = torch.Generator().manual_seed(args.seed)
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)
    model = GraphGPSScMAE(encoder_data.shape[1], args.pe_dim, args.hidden_size, args.decoder_hidden, args.num_heads, args.dropout).to(device)
    criterion = GraphGPSLoss(args.recon_weight, args.mask_weight, args.pe_weight, args.edge_weight, args.consistency_weight, args.variance_weight)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    keys = ["loss", "recon_loss", "mask_loss", "pe_loss", "edge_loss", "consistency_loss", "variance_loss", "edge_confidence", "edge_negative_confidence", "edge_proxy_accuracy", "edge_survival", "effective_mask_rate", "local_gate"]
    history = {k: [] for k in keys}
    history["stage"] = stage
    start = time.time()
    print(f"Using device: {device}")
    print(f"Dataset={dataset_name} cells={encoder_data.shape[0]} genes={encoder_data.shape[1]} clusters={n_clusters}")
    print(f"Method={METHOD_NAME} stage={stage} epochs={args.epochs}")

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {k: 0.0 for k in keys}
        n_batches = 0
        for _, x_cpu, log_cpu, _, neigh_ids, pe_cpu in train_loader:
            x = x_cpu.to(device)
            target = log_cpu.to(device)
            pe_batch = pe_cpu.to(device)
            neigh_x, neigh_pe, edge_survival = gather_neighbors(all_encoder, all_pe, neigh_ids, device, args.dropedge_prob)
            neigh_x2, neigh_pe2, _ = gather_neighbors(all_encoder, all_pe, neigh_ids, device, min(0.75, args.dropedge_prob + 0.25))
            corrupted, mask = model.mask_view(x, args.mask_prob)
            out = model(corrupted, neigh_x, pe_batch, neigh_pe)
            with torch.no_grad():
                out2 = model(corrupted, neigh_x2, pe_batch, neigh_pe2)
            pos_ids = neigh_ids[:, 0].to(all_encoder.device)
            neg_ids = torch.randint(0, len(dataset), (x.shape[0],), device=all_encoder.device)
            pos_z = node_base_embedding(model, all_encoder[pos_ids].to(device), all_pe[pos_ids].to(device))
            neg_z = node_base_embedding(model, all_encoder[neg_ids].to(device), all_pe[neg_ids].to(device))
            pos_logits = model.edge_logits(out["embedding"], pos_z)
            neg_logits = model.edge_logits(out["embedding"], neg_z)
            loss, parts = criterion(out, target, mask, pe_batch, pos_logits, neg_logits, out2["embedding"])
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            optimizer.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            for key, value in parts.items():
                sums[key] += value
            sums["edge_survival"] += edge_survival
            sums["effective_mask_rate"] += float(mask.mean().detach().cpu())
            sums["local_gate"] += float(out["local_gate"].mean().detach().cpu())
            n_batches += 1
        for key in keys:
            history[key].append(sums[key] / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} recon={history['recon_loss'][-1]:.4f} pe={history['pe_loss'][-1]:.4f} edge={history['edge_loss'][-1]:.4f} gate={history['local_gate'][-1]:.4f}")

    embedding, labels_out, edge_conf, gate_vals = extract_embedding(model, full_loader, all_encoder, all_pe, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "edge_confidence.npy", edge_conf.astype(np.float32))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "optimizer": optimizer.state_dict(), "args": vars(args), "gene_names": target_bundle.gene_names.astype(str)}, save_dir / "model_checkpoint.pth")

    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage, "preprocessing": "scaled encoder + log-expression targets + GraphGPS local-global adapter"})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    diag = diagnostics(embedding, labels_out, edge_conf, gate_vals, n_clusters, args.seed, preds, history["edge_survival"][-1])
    save_json(diag, str(save_dir / "diagnostics.json"))
    if not args.no_save_h5ad:
        encoder_bundle.adata.obsm["X_graphgps_scmae"] = embedding
        encoder_bundle.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(encoder_bundle.adata)
        encoder_bundle.adata.write_h5ad(save_dir / "adata_graphgps_scmae.h5ad", compression="gzip")

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
