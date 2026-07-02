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
from sklearn.cluster import KMeans, MiniBatchKMeans
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

from loss import MoleBERTContextLoss
from model import MoleBERTContextScMAE
from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, sanitize_anndata_for_write, save_json

METHOD_NAME = "rank44_molebert_context_token_full"
DISPLAY_NAME = "scMAE + Mole-BERT context token modeling"
RESULT_ROOT = ROOT / "methods" / "DeepLearning" / "scMAEs"
SINGLE_RESULT_CSV = RESULT_ROOT / "新模型独立快筛单次结果.csv"
SUMMARY_RESULT_CSV = RESULT_ROOT / "新模型独立快筛汇总结果.csv"
BASELINES = {"Melanoma_5K": {"nmi": 0.735414, "ari": 0.668029}, "Quake_10x_Spleen": {"nmi": 0.851730, "ari": 0.922275}, "Macosko": {"nmi": 0.657465, "ari": 0.494268}}


class ContextTokenDataset(Dataset):
    def __init__(self, encoder_data: np.ndarray, log_expr: np.ndarray, labels: np.ndarray, knn_index: np.ndarray, token_labels: np.ndarray):
        self.encoder_data = torch.as_tensor(encoder_data, dtype=torch.float32)
        self.log_expr = torch.as_tensor(log_expr, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)
        self.knn_index = torch.as_tensor(knn_index, dtype=torch.long)
        self.token_labels = torch.as_tensor(token_labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.encoder_data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.encoder_data[idx], self.log_expr[idx], self.labels[idx], self.knn_index[idx], self.token_labels[idx]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Independent Mole-BERT context-token scMAE candidate.", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
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
    p.add_argument("--dropout", type=float, default=0.05)
    p.add_argument("--neighbor_k", type=int, default=12)
    p.add_argument("--knn_pca_dim", type=int, default=50)
    p.add_argument("--num_tokens", type=int, default=128)
    p.add_argument("--mask_prob", type=float, default=0.4)
    p.add_argument("--dropedge_prob", type=float, default=0.25)
    p.add_argument("--recon_weight", type=float, default=1.0)
    p.add_argument("--mask_weight", type=float, default=0.05)
    p.add_argument("--token_weight", type=float, default=0.12)
    p.add_argument("--edge_weight", type=float, default=0.02)
    p.add_argument("--variance_weight", type=float, default=0.01)
    return p.parse_args()


def build_knn_and_tokens(data: np.ndarray, k: int, pca_dim: int, num_tokens: int, seed: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(data, dtype=np.float32)
    dim = int(min(max(2, pca_dim), x.shape[1] - 1, x.shape[0] - 1))
    rep = TruncatedSVD(n_components=dim, random_state=seed).fit_transform(x) if dim < x.shape[1] else x
    nn = NearestNeighbors(n_neighbors=min(k + 1, x.shape[0]), metric="cosine").fit(rep)
    idx = nn.kneighbors(rep, return_distance=False)[:, 1:]
    if idx.shape[1] < k:
        idx = np.concatenate([idx, np.repeat(idx[:, -1:], k - idx.shape[1], axis=1)], axis=1)
    ctx = rep[idx].mean(axis=1)
    tok_input = np.concatenate([rep, ctx, np.abs(rep - ctx)], axis=1).astype(np.float32)
    n_tok = int(min(max(8, num_tokens), max(8, x.shape[0] // 4)))
    km = MiniBatchKMeans(n_clusters=n_tok, random_state=seed, batch_size=min(4096, max(256, x.shape[0])), n_init=5)
    token_labels = km.fit_predict(tok_input).astype(np.int64)
    return idx.astype(np.int64), token_labels, tok_input.astype(np.float32)


def gather_neighbors(all_x: torch.Tensor, neighbor_ids: torch.Tensor, device: torch.device, drop_prob: float) -> tuple[torch.Tensor, float]:
    neighbor_ids = neighbor_ids.to(all_x.device)
    neigh = all_x[neighbor_ids.reshape(-1)].view(neighbor_ids.shape[0], neighbor_ids.shape[1], all_x.shape[1])
    keep = torch.rand(neighbor_ids.shape, device=neighbor_ids.device) >= float(drop_prob)
    none = ~keep.any(dim=1)
    if bool(none.any()):
        keep[none, torch.randint(0, neighbor_ids.shape[1], (int(none.sum()),), device=neighbor_ids.device)] = True
    visible = neigh.masked_fill(~keep.unsqueeze(-1), 0.0) * (neighbor_ids.shape[1] / keep.sum(dim=1).clamp_min(1).view(-1, 1, 1))
    return visible.to(device), float(keep.float().mean().detach().cpu())


@torch.no_grad()
def extract_embedding(model: MoleBERTContextScMAE, loader: DataLoader, all_x: torch.Tensor, device: torch.device) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    model.eval()
    emb, labels, edge_conf, gates = [], [], [], []
    for _, x_cpu, _, y, neigh_ids, _ in loader:
        x = x_cpu.to(device)
        neigh, _ = gather_neighbors(all_x, neigh_ids, device, 0.0)
        out = model(x, neigh)
        pos_ids = neigh_ids[:, 0].to(all_x.device)
        neg_ids = torch.randint(0, all_x.shape[0], (x.shape[0],), device=all_x.device)
        pos_z = model.encode_base(all_x[pos_ids].to(device))
        neg_z = model.encode_base(all_x[neg_ids].to(device))
        edge_conf.append((torch.sigmoid(model.edge_logits(out["embedding"], pos_z)) - torch.sigmoid(model.edge_logits(out["embedding"], neg_z))).detach().cpu().numpy())
        emb.append(out["embedding"].detach().cpu().numpy())
        gates.append(out["context_gate"].mean(dim=1).detach().cpu().numpy())
        labels.append(y.numpy())
    return np.nan_to_num(np.concatenate(emb).astype(np.float32)), np.concatenate(labels).astype(np.int64), np.nan_to_num(np.concatenate(edge_conf).astype(np.float32)), np.nan_to_num(np.concatenate(gates).astype(np.float32))


def neighbor_purity(labels: np.ndarray, embedding: np.ndarray, k: int = 10) -> float:
    if embedding.shape[0] <= 2:
        return float("nan")
    idx = NearestNeighbors(n_neighbors=min(k + 1, embedding.shape[0])).fit(embedding).kneighbors(embedding, return_distance=False)[:, 1:]
    return float(np.mean(labels[idx] == labels[:, None]))


def diagnostics(embedding: np.ndarray, labels: np.ndarray, edge_conf: np.ndarray, gates: np.ndarray, n_clusters: int, seed: int, preds: np.ndarray | None, edge_survival: float, token_labels: np.ndarray) -> dict:
    if preds is None:
        preds = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(embedding)
    frac = np.bincount(preds.astype(np.int64), minlength=n_clusters).astype(np.float64)
    frac = frac / max(1.0, frac.sum())
    probs = frac[frac > 0]
    label_counts = np.bincount(labels.astype(np.int64)).astype(np.float64)
    rare = set(np.where(label_counts <= max(5.0, 0.01 * labels.shape[0]))[0].tolist())
    var = float(np.var(embedding, axis=0).mean()) if embedding.size else 0.0
    return {"edge_survival": float(edge_survival), "neighbor_purity_proxy": neighbor_purity(labels, embedding), "mixed_cell_fraction": 0.0, "boundary_entropy": float(-(probs * np.log(probs)).sum() / max(np.log(max(2, n_clusters)), 1e-8)), "rare_risk_fraction": float(np.mean([x in rare for x in labels])) if labels.size else 0.0, "embedding_variance": var, "cluster_mass_min": float(frac.min()) if frac.size else 0.0, "cluster_mass_max": float(frac.max()) if frac.size else 0.0, "collapse_warning": bool((not np.isfinite(var)) or var < 1e-8 or (frac.size and frac.min() < 0.001) or (frac.size and frac.max() > 0.95)), "edge_confidence_mean": float(np.nanmean(edge_conf)), "edge_confidence_std": float(np.nanstd(edge_conf)), "context_gate_mean": float(np.mean(gates)) if gates.size else float("nan"), "token_count_used": int(np.unique(token_labels).size), "diagnostic_note": "No NeighborMix is used. Mole-BERT adaptation uses context-aware discrete cell tokens and masked expression reconstruction."}


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
    knn_index, token_labels, token_features = build_knn_and_tokens(encoder_data, args.neighbor_k, args.knn_pca_dim, args.num_tokens, args.seed)
    np.save(save_dir / "gene_names.npy", target.gene_names.astype(str))
    np.save(save_dir / "knn_index.npy", knn_index)
    np.save(save_dir / "context_token_labels.npy", token_labels)
    save_json(target.profile, str(save_dir / "dataset_profile.json"))
    save_json({**target.preprocess_config, "encoder_scale_input": bool(args.scale_input), "target_scale_input": False, "context_tokenizer": f"MiniBatchKMeans over SVD cell + KNN context features, tokens={np.unique(token_labels).size}"}, str(save_dir / "preprocess_config.json"))
    dataset = ContextTokenDataset(encoder_data, log_expr, labels, knn_index, token_labels)
    all_encoder = dataset.encoder_data
    train_loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=torch.Generator().manual_seed(args.seed))
    full_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False)
    model = MoleBERTContextScMAE(encoder_data.shape[1], int(np.unique(token_labels).size), args.hidden_size, args.decoder_hidden, args.dropout).to(device)
    criterion = MoleBERTContextLoss(args.recon_weight, args.mask_weight, args.token_weight, args.edge_weight, args.variance_weight)
    opt = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    keys = ["loss", "recon_loss", "mask_loss", "token_loss", "edge_loss", "variance_loss", "token_accuracy_proxy", "edge_confidence", "edge_negative_confidence", "edge_proxy_accuracy", "edge_survival", "effective_mask_rate", "context_gate"]
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
        for _, x_cpu, log_cpu, _, neigh_ids, tok_cpu in train_loader:
            x, y, toks = x_cpu.to(device), log_cpu.to(device), tok_cpu.to(device)
            neigh, edge_survival = gather_neighbors(all_encoder, neigh_ids, device, args.dropedge_prob)
            corrupted, mask = model.mask_view(x, args.mask_prob)
            out = model(corrupted, neigh)
            pos_ids = neigh_ids[:, 0].to(all_encoder.device)
            neg_ids = torch.randint(0, len(dataset), (x.shape[0],), device=all_encoder.device)
            pos_z = model.encode_base(all_encoder[pos_ids].to(device))
            neg_z = model.encode_base(all_encoder[neg_ids].to(device))
            loss, parts = criterion(out, y, mask, toks, model.edge_logits(out["embedding"], pos_z), model.edge_logits(out["embedding"], neg_z))
            if not torch.isfinite(loss):
                raise FloatingPointError(f"Non-finite loss at epoch {epoch}: {parts}")
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            opt.step()
            for k, v in parts.items():
                sums[k] += v
            sums["edge_survival"] += edge_survival
            sums["effective_mask_rate"] += float(mask.mean().detach().cpu())
            sums["context_gate"] += float(out["context_gate"].mean().detach().cpu())
            nb += 1
        for k in keys:
            history[k].append(sums[k] / max(1, nb))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} recon={history['recon_loss'][-1]:.4f} token={history['token_loss'][-1]:.4f} tok_acc={history['token_accuracy_proxy'][-1]:.3f} gate={history['context_gate'][-1]:.4f}")
    embedding, labels_out, edge_conf, gates = extract_embedding(model, full_loader, all_encoder, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "edge_confidence.npy", edge_conf.astype(np.float32))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "optimizer": opt.state_dict(), "args": vars(args), "gene_names": target.gene_names.astype(str)}, save_dir / "model_checkpoint.pth")
    eval_result, preds = None, None
    if not args.skip_eval:
        eval_result = family.write_kmeans_known_k_outputs(save_dir, dataset_name, DISPLAY_NAME, args.seed, embedding, labels_out, n_clusters, {"variant": METHOD_NAME, "stage": stage, "preprocessing": "scaled encoder + log-expression targets + Mole-BERT context tokens"})
        preds = eval_result["preds"]["kmeans_known_k"]
        save_json(eval_result["fixed"], str(save_dir / "metrics.json"))
    diag = diagnostics(embedding, labels_out, edge_conf, gates, n_clusters, args.seed, preds, history["edge_survival"][-1], token_labels)
    save_json(diag, str(save_dir / "diagnostics.json"))
    if not args.no_save_h5ad:
        enc.adata.obsm["X_molebert_context_scmae"] = embedding
        enc.adata.uns[METHOD_NAME] = {"method": DISPLAY_NAME, "variant": METHOD_NAME, "stage": stage}
        sanitize_anndata_for_write(enc.adata)
        enc.adata.write_h5ad(save_dir / "adata_molebert_context_scmae.h5ad", compression="gzip")
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
