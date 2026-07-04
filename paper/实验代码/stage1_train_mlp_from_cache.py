"""Train a minimal scMAE-compatible MLP from cached HVG arrays.

Input cache is produced by stage1_scmae_compatible.py --save-cache.
This script avoids h5py/scipy and is intended for Torch environments where
the sparse scientific stack may be fragile.
"""

from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.metrics import adjusted_rand_score, f1_score, normalized_mutual_info_score

try:
    from scipy.optimize import linear_sum_assignment
except Exception:  # pragma: no cover - optional dependency guard for fragile local envs.
    linear_sum_assignment = None


def ensure_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)


def write_json(path: Path, obj: Dict) -> None:
    with path.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2, ensure_ascii=False)


def load_labels(path: Path) -> Optional[np.ndarray]:
    if not path.exists():
        return None
    return np.load(path, allow_pickle=True)


def encode_labels(labels: np.ndarray) -> Tuple[np.ndarray, Dict[str, int]]:
    labels_str = np.asarray([str(x) for x in labels], dtype=object)
    unique = sorted(set(labels_str.tolist()))
    mapping = {label: i for i, label in enumerate(unique)}
    y = np.asarray([mapping[x] for x in labels_str], dtype=np.int64)
    return y, mapping


def simple_kmeans(
    x: np.ndarray,
    n_clusters: int,
    seed: int,
    max_iter: int = 100,
) -> np.ndarray:
    """Small NumPy KMeans to avoid fragile threadpoolctl paths on Windows."""
    rng = np.random.default_rng(seed)
    n_samples = x.shape[0]
    if n_clusters <= 0 or n_clusters > n_samples:
        raise ValueError(f"Invalid n_clusters={n_clusters} for n_samples={n_samples}")
    centers = x[rng.choice(n_samples, size=n_clusters, replace=False)].astype(np.float32)
    labels = np.zeros(n_samples, dtype=np.int64)
    for _ in range(max_iter):
        distances = ((x[:, None, :] - centers[None, :, :]) ** 2).sum(axis=2)
        new_labels = distances.argmin(axis=1).astype(np.int64)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for k in range(n_clusters):
            members = x[labels == k]
            if len(members) == 0:
                centers[k] = x[rng.integers(0, n_samples)]
            else:
                centers[k] = members.mean(axis=0)
    return labels


def greedy_assignment(score: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Fallback assignment when SciPy is unavailable; exact Hungarian is preferred."""
    remaining_rows = set(range(score.shape[0]))
    remaining_cols = set(range(score.shape[1]))
    rows = []
    cols = []
    while remaining_rows and remaining_cols:
        best = None
        for r in remaining_rows:
            for c in remaining_cols:
                item = (score[r, c], r, c)
                if best is None or item > best:
                    best = item
        _, row, col = best
        rows.append(row)
        cols.append(col)
        remaining_rows.remove(row)
        remaining_cols.remove(col)
    return np.asarray(rows, dtype=np.int64), np.asarray(cols, dtype=np.int64)


def best_map_cluster_labels(y_true: np.ndarray, y_pred: np.ndarray, n_classes: int) -> Tuple[np.ndarray, Dict]:
    contingency = np.zeros((n_classes, n_classes), dtype=np.int64)
    for true_label, pred_label in zip(y_true, y_pred):
        if 0 <= true_label < n_classes and 0 <= pred_label < n_classes:
            contingency[int(true_label), int(pred_label)] += 1

    if linear_sum_assignment is not None:
        true_ind, pred_ind = linear_sum_assignment(contingency.max() - contingency)
        method = "hungarian"
    else:
        true_ind, pred_ind = greedy_assignment(contingency)
        method = "greedy_fallback"

    pred_to_true = {int(pred): int(true) for true, pred in zip(true_ind, pred_ind)}
    mapped = np.asarray([pred_to_true.get(int(label), -1) for label in y_pred], dtype=np.int64)
    mapping = {
        "method": method,
        "pred_to_true": {str(k): int(v) for k, v in pred_to_true.items()},
        "contingency": contingency.tolist(),
    }
    return mapped, mapping


class ScMAECompatibleMLP(nn.Module):
    def __init__(self, n_genes: int, latent_dim: int = 64, hidden_dim: int = 512):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(n_genes, hidden_dim),
            nn.BatchNorm1d(hidden_dim),
            nn.ReLU(),
            nn.Dropout(0.1),
            nn.Linear(hidden_dim, latent_dim),
        )
        self.decoder = nn.Linear(latent_dim, n_genes)
        self.mask_head = nn.Linear(latent_dim, n_genes)

    def forward(self, x_tilde: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        z = self.encoder(x_tilde)
        x_hat = self.decoder(z)
        mask_logits = self.mask_head(z)
        return z, x_hat, mask_logits


def gene_wise_shuffle(x: torch.Tensor) -> torch.Tensor:
    columns = []
    batch_size, n_genes = x.shape
    for j in range(n_genes):
        perm = torch.randperm(batch_size, device=x.device)
        columns.append(x[perm, j])
    return torch.stack(columns, dim=1)


def corrupt_batch(x: torch.Tensor, mask_ratio: float) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    mask = (torch.rand_like(x) < mask_ratio).float()
    x_prime = gene_wise_shuffle(x)
    x_tilde = x * (1.0 - mask) + x_prime * mask
    return x_tilde, mask, x_prime


def batch_mask_diagnostics(x: torch.Tensor, mask: torch.Tensor, x_prime: torch.Tensor) -> Dict[str, float]:
    with torch.no_grad():
        masked = mask > 0
        masked_count = masked.sum().clamp_min(1)
        observed = x > 0
        changed = masked & (x_prime != x)
        zero_to_zero = masked & (x == 0) & (x_prime == 0)
        return {
            "masked_count": float(masked.sum().item()),
            "masked_observed_count": float((masked & observed).sum().item()),
            "zero_to_zero_count": float(zero_to_zero.sum().item()),
            "effective_changed_count": float(changed.sum().item()),
            "masked_zero_count": float((masked & (x == 0)).sum().item()),
            "masked_nonzero_to_zero_count": float((masked & (x > 0) & (x_prime == 0)).sum().item()),
            "masked_zero_to_nonzero_count": float((masked & (x == 0) & (x_prime > 0)).sum().item()),
            "total_count": float(mask.numel()),
            "observed_count": float(observed.sum().item()),
            "mask_gene_sum": mask.sum(dim=0).detach().cpu().numpy(),
        }


def merge_diag(acc: Dict, item: Dict) -> None:
    for key, value in item.items():
        if key == "mask_gene_sum":
            if key not in acc:
                acc[key] = value.copy()
            else:
                acc[key] += value
        else:
            acc[key] = acc.get(key, 0.0) + float(value)


def finalize_diag(acc: Dict, n_genes: int) -> Dict:
    masked_count = max(acc.get("masked_count", 0.0), 1.0)
    total_count = max(acc.get("total_count", 0.0), 1.0)
    observed_count = max(acc.get("observed_count", 0.0), 1.0)
    gene_frequency = acc.get("mask_gene_sum", np.zeros(n_genes)) / max(total_count / n_genes, 1.0)
    probs = gene_frequency / max(float(gene_frequency.sum()), 1e-12)
    probs = probs[probs > 0]
    entropy = float(-(probs * np.log(probs + 1e-12)).sum()) if probs.size else 0.0
    return {
        "actual_mask_ratio_global": float(acc.get("masked_count", 0.0) / total_count),
        "actual_mask_ratio_observed": float(acc.get("masked_observed_count", 0.0) / observed_count),
        "masked_observed_fraction_among_masked": float(acc.get("masked_observed_count", 0.0) / masked_count),
        "zero_to_zero_fraction_among_masked": float(acc.get("zero_to_zero_count", 0.0) / masked_count),
        "effective_changed_fraction_among_masked": float(
            acc.get("effective_changed_count", 0.0) / masked_count
        ),
        "masked_zero_fraction_among_masked": float(acc.get("masked_zero_count", 0.0) / masked_count),
        "masked_nonzero_to_zero_fraction_among_masked": float(
            acc.get("masked_nonzero_to_zero_count", 0.0) / masked_count
        ),
        "masked_zero_to_nonzero_fraction_among_masked": float(
            acc.get("masked_zero_to_nonzero_count", 0.0) / masked_count
        ),
        "mask_entropy": entropy,
        "mask_entropy_normalized": float(entropy / np.log(n_genes)),
        "gene_mask_frequency_mean": float(gene_frequency.mean()),
        "gene_mask_frequency_std": float(gene_frequency.std()),
        "gene_mask_frequency_max": float(gene_frequency.max()),
    }


def clustering_metrics(
    embedding: np.ndarray,
    labels: Optional[np.ndarray],
    seed: int,
    name: str,
    save_dir: Path,
) -> Dict:
    if labels is None:
        return {}
    y, mapping = encode_labels(labels)
    k = len(mapping)
    pred = simple_kmeans(embedding.astype(np.float32), n_clusters=k, seed=seed)
    mapped_pred, cluster_mapping = best_map_cluster_labels(y, pred, k)

    np.save(save_dir / f"cluster_labels_{name}_kmeans.npy", pred.astype(np.int64))
    np.save(save_dir / f"cluster_labels_{name}_kmeans_mapped.npy", mapped_pred.astype(np.int64))
    write_json(
        save_dir / f"cluster_mapping_{name}_kmeans.json",
        {
            "label_to_id": mapping,
            **cluster_mapping,
        },
    )

    return {
        f"{name}_n_classes": int(k),
        f"{name}_ari": float(adjusted_rand_score(y, pred)),
        f"{name}_nmi": float(normalized_mutual_info_score(y, pred)),
        f"{name}_macro_f1_unmapped": float(f1_score(y, pred, average="macro")),
        f"{name}_acc_mapped": float((mapped_pred == y).mean()),
        f"{name}_macro_f1_mapped": float(
            f1_score(y, mapped_pred, labels=np.arange(k), average="macro", zero_division=0)
        ),
        f"{name}_cluster_mapping_method": cluster_mapping["method"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, required=True)
    parser.add_argument("--save-dir", type=Path, required=True)
    parser.add_argument("--epochs", type=int, default=2)
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--mask-ratio", type=float, default=0.3)
    parser.add_argument("--latent-dim", type=int, default=64)
    parser.add_argument("--hidden-dim", type=int, default=512)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--lambda-masked", type=float, default=4.0)
    parser.add_argument("--gamma-mask", type=float, default=1.0)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--device", type=str, default="auto")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    start = time.time()
    if args.save_dir.exists() and any(args.save_dir.iterdir()) and not args.force:
        raise FileExistsError(f"{args.save_dir} is not empty. Use --force to overwrite/add files.")
    ensure_dir(args.save_dir)

    torch.manual_seed(args.seed)
    np.random.seed(args.seed)
    if args.device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    else:
        device = args.device

    x = np.load(args.cache_dir / "x_hvg_log1p.npy").astype(np.float32)
    n_cells, n_genes = x.shape
    x_tensor = torch.from_numpy(x)

    model = ScMAECompatibleMLP(n_genes=n_genes, latent_dim=args.latent_dim, hidden_dim=args.hidden_dim).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=1e-4)

    history = []
    diag_acc = {}
    indices = np.arange(n_cells)
    for epoch in range(args.epochs):
        np.random.shuffle(indices)
        model.train()
        epoch_loss = 0.0
        epoch_rec = 0.0
        epoch_mask = 0.0
        n_seen = 0
        for start_idx in range(0, n_cells, args.batch_size):
            batch_idx = indices[start_idx : start_idx + args.batch_size]
            xb = x_tensor[batch_idx].to(device)
            x_tilde, mask, x_prime = corrupt_batch(xb, args.mask_ratio)
            merge_diag(diag_acc, batch_mask_diagnostics(xb, mask, x_prime))

            _, x_hat, mask_logits = model(x_tilde)
            rec_weight = 1.0 + args.lambda_masked * mask
            loss_rec = (rec_weight * (x_hat - xb).pow(2)).mean()
            loss_mask = F.binary_cross_entropy_with_logits(mask_logits, mask)
            loss = loss_rec + args.gamma_mask * loss_mask

            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

            batch_n = len(batch_idx)
            epoch_loss += float(loss.detach().cpu()) * batch_n
            epoch_rec += float(loss_rec.detach().cpu()) * batch_n
            epoch_mask += float(loss_mask.detach().cpu()) * batch_n
            n_seen += batch_n

        history.append(
            {
                "epoch": epoch + 1,
                "loss": epoch_loss / n_seen,
                "loss_rec": epoch_rec / n_seen,
                "loss_mask": epoch_mask / n_seen,
            }
        )

    model.eval()
    embeddings = []
    with torch.no_grad():
        for start_idx in range(0, n_cells, args.batch_size):
            xb = x_tensor[start_idx : start_idx + args.batch_size].to(device)
            z = model.encoder(xb)
            embeddings.append(z.detach().cpu().numpy())
    embedding = np.concatenate(embeddings, axis=0)

    labels_celltype = load_labels(args.cache_dir / "labels_Celltype.npy")
    labels_seurat = load_labels(args.cache_dir / "labels_Seurat_clusters.npy")
    metrics = {}
    metrics.update(clustering_metrics(embedding, labels_celltype, args.seed, "celltype", args.save_dir))
    metrics.update(clustering_metrics(embedding, labels_seurat, args.seed, "seurat_clusters", args.save_dir))
    metrics["embedding_shape"] = [int(x) for x in embedding.shape]

    param_count = sum(p.numel() for p in model.parameters())
    trainable_param_count = sum(p.numel() for p in model.parameters() if p.requires_grad)
    mask_diagnostics = finalize_diag(diag_acc, n_genes)

    np.save(args.save_dir / "embedding_final.npy", embedding.astype(np.float32))
    write_json(args.save_dir / "metrics.json", metrics)
    write_json(args.save_dir / "training_history.json", {"history": history})
    write_json(args.save_dir / "mask_diagnostics.json", mask_diagnostics)
    write_json(
        args.save_dir / "runtime.json",
        {
            "elapsed_seconds": float(time.time() - start),
            "device": str(device),
            "cuda_available": bool(torch.cuda.is_available()),
            "training_executed": True,
        },
    )
    write_json(
        args.save_dir / "param_count.json",
        {
            "param_count": int(param_count),
            "trainable_param_count": int(trainable_param_count),
        },
    )
    write_json(
        args.save_dir / "config.json",
        {
            "cache_dir": str(args.cache_dir),
            "save_dir": str(args.save_dir),
            "epochs": int(args.epochs),
            "batch_size": int(args.batch_size),
            "mask_ratio": float(args.mask_ratio),
            "latent_dim": int(args.latent_dim),
            "hidden_dim": int(args.hidden_dim),
            "lr": float(args.lr),
            "lambda_masked": float(args.lambda_masked),
            "gamma_mask": float(args.gamma_mask),
            "seed": int(args.seed),
            "device": str(device),
        },
    )

    print(json.dumps({"history": history, "metrics": metrics, "mask_diagnostics": mask_diagnostics}, indent=2))


if __name__ == "__main__":
    main()
