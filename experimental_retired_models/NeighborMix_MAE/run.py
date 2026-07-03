#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import random
import sys
from pathlib import Path
from typing import Optional

import h5py
import numpy as np
import scipy.sparse as sp
import torch
import torch.nn as nn
import torch.nn.functional as F
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler, normalize
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent


def _find_project_root(start: Path) -> Path:
    for parent in [start, *start.parents]:
        if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists():
            return parent
    raise RuntimeError("Could not locate project root containing methods/DeepLearning/scMAE_family.py")


ROOT = _find_project_root(CURRENT_DIR)
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experimental_retired_models.PlantSPADE_LGCL.data import load_lgcl_dataset, write_dataset_artifacts
from experimental_retired_models.PlantSPADE_LGCL.eval import write_evaluation_outputs
from experimental_retired_models.PlantSPADE_LGCL.utils import ensure_dir, save_json


class IndexedExpressionDataset(Dataset):
    def __init__(self, data: np.ndarray, labels: Optional[np.ndarray]):
        self.data = torch.as_tensor(data, dtype=torch.float32)
        if labels is None:
            labels = np.zeros(data.shape[0], dtype=np.int64)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.data.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.data[idx], self.labels[idx]


class NeighborMixMAE(nn.Module):
    def __init__(
        self,
        input_dim: int,
        latent_dim: int = 32,
        hidden_dim: int = 512,
        bottleneck_dim: int = 128,
        dropout: float = 0.1,
    ):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Dropout(float(dropout)),
            nn.Linear(input_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, bottleneck_dim),
            nn.LayerNorm(bottleneck_dim),
            nn.GELU(),
            nn.Linear(bottleneck_dim, latent_dim),
        )
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, bottleneck_dim),
            nn.LayerNorm(bottleneck_dim),
            nn.GELU(),
            nn.Linear(bottleneck_dim, hidden_dim),
            nn.LayerNorm(hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, input_dim),
        )

    def forward(self, x: torch.Tensor):
        z = self.encoder(x)
        recon = self.decoder(z)
        return z, recon

    @torch.no_grad()
    def feature(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)


def str2bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "t", "yes", "y"}:
        return True
    if value in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got {value!r}")


def parse_float_list(value: str):
    if value is None or str(value).strip() == "":
        return []
    return [float(item) for item in str(value).split(",") if item.strip()]


def parse_args():
    parser = argparse.ArgumentParser(description="NeighborMix-MAE fixed-protocol runner")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default="neighbormix_mae")
    parser.add_argument("--variant_name", default="neighbormix_mae")
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=2000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--svd_dim", type=int, default=32)
    parser.add_argument("--svd_iter", type=int, default=7)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--latent_dim", type=int, default=32)
    parser.add_argument("--hidden_dim", type=int, default=512)
    parser.add_argument("--bottleneck_dim", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-5)
    parser.add_argument("--grad_clip", type=float, default=5.0)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
    parser.add_argument("--mask_strategy", default="zero", choices=["zero", "random_swap"])
    parser.add_argument("--masked_weight", type=float, default=1.0)
    parser.add_argument("--unmasked_weight", type=float, default=0.2)
    parser.add_argument("--self_weight", type=float, default=1.0)
    parser.add_argument("--denoise_weight", type=float, default=1.0)
    parser.add_argument("--consistency_weight", type=float, default=0.0)
    parser.add_argument("--mix_mode", default="neighbor", choices=["none", "random", "neighbor"])
    parser.add_argument("--target_mode", default="original", choices=["original", "mixed"])
    parser.add_argument("--alpha", type=float, default=0.8)
    parser.add_argument("--neighbor_k", type=int, default=15)
    parser.add_argument("--mix_neighbors", type=int, default=4)
    parser.add_argument("--tau", type=float, default=0.2)
    parser.add_argument("--knn_pca_dim", type=int, default=50)
    parser.add_argument("--eval_neighbors", type=int, default=15)
    parser.add_argument("--leiden_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--louvain_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--leiden_resolutions", default="0.2,0.4,0.6,0.8,1.0,1.2")
    parser.add_argument("--include_louvain", type=str2bool, default=False)
    parser.add_argument("--run_oracle_sweep", type=str2bool, default=False)
    parser.add_argument("--sweep_max_cells", type=int, default=10000)
    parser.add_argument("--silhouette_sample_size", type=int, default=3000)
    parser.add_argument("--skip_eval", type=str2bool, default=False)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--no_save_h5ad", action="store_true")
    return parser.parse_args()


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible:
        visible_ids = [item.strip() for item in visible.split(",") if item.strip()]
        if set(visible_ids).intersection({"0", "7"}):
            raise ValueError("CUDA_VISIBLE_DEVICES includes forbidden physical GPU 0 or 7.")
        if len(visible_ids) == 1:
            return torch.device("cuda:0")
        if str(gpu) in visible_ids:
            return torch.device(f"cuda:{visible_ids.index(str(gpu))}")
        if 0 <= gpu < len(visible_ids):
            return torch.device(f"cuda:{gpu}")
        raise ValueError(f"--gpu {gpu} is outside isolated CUDA_VISIBLE_DEVICES={visible!r}.")
    if gpu in {0, 7}:
        raise ValueError("Physical GPU 0 and GPU 7 are forbidden. Use 1,2,3,4,5,6 or --no_cuda.")
    return torch.device(f"cuda:{gpu}")


def matrix_to_dense_float32(matrix) -> np.ndarray:
    if sp.issparse(matrix):
        arr = matrix.toarray()
    else:
        arr = np.asarray(matrix)
    arr = np.asarray(arr, dtype=np.float32)
    arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    arr[arr < 0.0] = 0.0
    return arr


def build_knn_distribution(
    data: np.ndarray,
    k: int,
    pca_dim: int,
    tau: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray, dict]:
    n_cells, n_genes = data.shape
    k = max(1, min(int(k), n_cells - 1))
    n_components = max(2, min(int(pca_dim), n_cells - 1, n_genes - 1))
    if n_cells < 3:
        indices = np.zeros((n_cells, 1), dtype=np.int64)
        probs = np.ones((n_cells, 1), dtype=np.float32)
        return indices, probs, {"neighbor_k": 1, "pca_dim": 0, "note": "too_few_cells"}

    if n_genes > n_components:
        reducer = PCA(n_components=n_components, random_state=seed) if not sp.issparse(data) else TruncatedSVD(n_components=n_components, random_state=seed)
        emb = reducer.fit_transform(data).astype(np.float32)
    else:
        emb = data.astype(np.float32, copy=False)
    if emb.shape[1] > 1:
        emb = StandardScaler().fit_transform(emb).astype(np.float32)
    emb = normalize(emb, norm="l2", axis=1, copy=False).astype(np.float32)

    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine")
    nn.fit(emb)
    distances, neighbors = nn.kneighbors(emb, return_distance=True)
    neighbors = neighbors[:, 1:].astype(np.int64, copy=False)
    distances = distances[:, 1:].astype(np.float32, copy=False)
    sim = np.clip(1.0 - distances, a_min=0.0, a_max=None)
    logits = sim / max(float(tau), 1e-6)
    logits = logits - logits.max(axis=1, keepdims=True)
    probs = np.exp(logits).astype(np.float64)
    probs = probs / probs.sum(axis=1, keepdims=True).clip(min=1e-12)
    return neighbors, probs.astype(np.float32), {
        "neighbor_k": int(k),
        "pca_dim": int(n_components),
        "tau": float(tau),
        "mean_max_neighbor_prob": float(probs.max(axis=1).mean()),
        "mean_neighbor_similarity": float(sim.mean()),
    }


def apply_mask(x: torch.Tensor, mask_ratio: float, strategy: str) -> tuple[torch.Tensor, torch.Tensor]:
    mask = (torch.rand_like(x) < float(mask_ratio)).float()
    if strategy == "zero":
        corrupted = x * (1.0 - mask)
    elif strategy == "random_swap":
        if x.shape[0] <= 1:
            replacement = torch.zeros_like(x)
        else:
            replacement = x[torch.randperm(x.shape[0], device=x.device)]
        corrupted = torch.where(mask.bool(), replacement, x)
    else:
        raise ValueError(f"Unknown mask_strategy: {strategy}")
    return corrupted, mask


def reconstruction_loss(
    recon: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
    masked_weight: float,
    unmasked_weight: float,
) -> torch.Tensor:
    weight = mask * float(masked_weight) + (1.0 - mask) * float(unmasked_weight)
    return (weight * (recon - target).pow(2)).mean()


def sample_mix(
    data_np: np.ndarray,
    batch_indices: np.ndarray,
    batch_x: torch.Tensor,
    mode: str,
    alpha: float,
    mix_neighbors: int,
    rng: np.random.Generator,
    neighbor_indices: Optional[np.ndarray] = None,
    neighbor_probs: Optional[np.ndarray] = None,
) -> torch.Tensor:
    if mode == "none":
        return batch_x
    n_cells = data_np.shape[0]
    bsz = int(batch_indices.shape[0])
    mix_neighbors = max(1, int(mix_neighbors))
    if mode == "random":
        sampled = rng.integers(0, n_cells, size=(bsz, mix_neighbors), dtype=np.int64)
        weights = np.full((bsz, mix_neighbors), 1.0 / mix_neighbors, dtype=np.float32)
    elif mode == "neighbor":
        if neighbor_indices is None or neighbor_probs is None:
            raise ValueError("NeighborMix requires neighbor_indices and neighbor_probs.")
        sampled = np.empty((bsz, mix_neighbors), dtype=np.int64)
        weights = np.empty((bsz, mix_neighbors), dtype=np.float32)
        for pos, cell in enumerate(batch_indices):
            probs = neighbor_probs[cell]
            choices = rng.choice(neighbor_indices.shape[1], size=mix_neighbors, replace=True, p=probs)
            sampled[pos] = neighbor_indices[cell, choices]
            picked = probs[choices].astype(np.float32, copy=False)
            weights[pos] = picked / max(float(picked.sum()), 1e-12)
    else:
        raise ValueError(f"Unknown mix mode: {mode}")

    neighbor_expr = data_np[sampled]
    neighbor_mean = np.sum(neighbor_expr * weights[:, :, None], axis=1).astype(np.float32)
    neighbor_t = torch.as_tensor(neighbor_mean, dtype=batch_x.dtype, device=batch_x.device)
    alpha = float(alpha)
    return alpha * batch_x + (1.0 - alpha) * neighbor_t


@torch.no_grad()
def extract_embedding(model: NeighborMixMAE, loader: DataLoader, device: torch.device) -> tuple[np.ndarray, np.ndarray]:
    model.eval()
    embeddings = []
    labels = []
    for _, x, y in loader:
        x = x.to(device)
        z = model.feature(x)
        embeddings.append(z.detach().cpu().numpy())
        labels.append(y.numpy())
    emb = np.concatenate(embeddings, axis=0).astype(np.float32)
    labels_np = np.concatenate(labels, axis=0).astype(np.int64)
    emb = np.nan_to_num(emb, nan=0.0, posinf=0.0, neginf=0.0)
    emb = normalize(emb, norm="l2", axis=1, copy=False).astype(np.float32)
    return emb, labels_np


def save_embedding_h5(path: Path, embedding: np.ndarray, labels: Optional[np.ndarray]) -> None:
    with h5py.File(path, "w") as handle:
        handle.create_dataset("X", data=embedding.astype(np.float32))
        if labels is not None:
            handle.create_dataset("labels", data=labels.astype(np.int64))


def main():
    args = parse_args()
    set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = get_device(args.gpu, args.no_cuda)
    print(f"Using device: {device}")

    bundle = load_lgcl_dataset(
        args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        svd_dim=args.svd_dim,
        svd_iter=args.svd_iter,
        seed=args.seed,
        label_key=args.label_key,
    )
    write_dataset_artifacts(bundle, str(save_dir))
    if bundle.labels is None:
        raise ValueError("Labels are required for fixed benchmark evaluation.")

    data_np = matrix_to_dense_float32(bundle.amplitude)
    labels = bundle.labels.astype(np.int64)
    n_clusters = int(args.n_clusters) if args.n_clusters and args.n_clusters > 0 else int(len(np.unique(labels)))
    dataset_name = args.dataset_name or Path(args.data_path).stem
    print(f"Cells={data_np.shape[0]} genes={data_np.shape[1]} clusters={n_clusters} variant={args.variant_name}")

    neighbor_indices = None
    neighbor_probs = None
    neighbor_profile = {"enabled": False}
    if args.mix_mode == "neighbor":
        neighbor_indices, neighbor_probs, neighbor_profile = build_knn_distribution(
            data_np,
            k=args.neighbor_k,
            pca_dim=args.knn_pca_dim,
            tau=args.tau,
            seed=args.seed,
        )
        neighbor_profile["enabled"] = True
    save_json(neighbor_profile, str(save_dir / "neighbor_graph_profile.json"))

    dataset = IndexedExpressionDataset(data_np, labels)
    generator = torch.Generator()
    generator.manual_seed(args.seed)
    train_loader = DataLoader(
        dataset,
        batch_size=args.batch_size,
        shuffle=True,
        drop_last=False,
        generator=generator,
    )
    eval_loader = DataLoader(dataset, batch_size=max(args.batch_size * 4, 512), shuffle=False, drop_last=False)

    model = NeighborMixMAE(
        input_dim=data_np.shape[1],
        latent_dim=args.latent_dim,
        hidden_dim=args.hidden_dim,
        bottleneck_dim=args.bottleneck_dim,
        dropout=args.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=max(1, args.epochs))
    rng = np.random.default_rng(args.seed + 2027)
    history = {
        "loss": [],
        "self_loss": [],
        "denoise_loss": [],
        "consistency_loss": [],
        "lr": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        totals = {key: 0.0 for key in ["loss", "self_loss", "denoise_loss", "consistency_loss"]}
        n_batches = 0
        for idx_t, x_cpu, _ in train_loader:
            idx_np = idx_t.numpy().astype(np.int64, copy=False)
            x = x_cpu.to(device)
            x_self_corrupt, self_mask = apply_mask(x, args.mask_ratio, args.mask_strategy)
            z_self, recon_self = model(x_self_corrupt)
            self_loss = reconstruction_loss(
                recon_self,
                x,
                self_mask,
                masked_weight=args.masked_weight,
                unmasked_weight=args.unmasked_weight,
            )

            denoise_loss = torch.zeros((), dtype=torch.float32, device=device)
            consistency_loss = torch.zeros((), dtype=torch.float32, device=device)
            if args.mix_mode != "none" and args.denoise_weight > 0.0:
                x_mix = sample_mix(
                    data_np=data_np,
                    batch_indices=idx_np,
                    batch_x=x,
                    mode=args.mix_mode,
                    alpha=args.alpha,
                    mix_neighbors=args.mix_neighbors,
                    rng=rng,
                    neighbor_indices=neighbor_indices,
                    neighbor_probs=neighbor_probs,
                )
                x_mix_corrupt, mix_mask = apply_mask(x_mix, args.mask_ratio, args.mask_strategy)
                z_mix, recon_mix = model(x_mix_corrupt)
                target = x if args.target_mode == "original" else x_mix
                denoise_loss = reconstruction_loss(
                    recon_mix,
                    target,
                    mix_mask,
                    masked_weight=args.masked_weight,
                    unmasked_weight=args.unmasked_weight,
                )
                if args.consistency_weight > 0.0:
                    consistency_loss = (
                        1.0
                        - F.cosine_similarity(F.normalize(z_mix, dim=1), F.normalize(z_self.detach(), dim=1), dim=1)
                    ).mean()

            loss = float(args.self_weight) * self_loss
            loss = loss + float(args.denoise_weight) * denoise_loss
            loss = loss + float(args.consistency_weight) * consistency_loss
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            if args.grad_clip and args.grad_clip > 0:
                torch.nn.utils.clip_grad_norm_(model.parameters(), float(args.grad_clip))
            optimizer.step()

            totals["loss"] += float(loss.detach().cpu())
            totals["self_loss"] += float(self_loss.detach().cpu())
            totals["denoise_loss"] += float(denoise_loss.detach().cpu())
            totals["consistency_loss"] += float(consistency_loss.detach().cpu())
            n_batches += 1
        scheduler.step()
        for key in totals:
            history[key].append(totals[key] / max(1, n_batches))
        history["lr"].append(float(scheduler.get_last_lr()[0]))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"Epoch {epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"self={history['self_loss'][-1]:.4f} denoise={history['denoise_loss'][-1]:.4f} "
                f"cons={history['consistency_loss'][-1]:.4f}"
            )

    embedding, labels_out = extract_embedding(model, eval_loader, device)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels_out.astype(np.int64))
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    save_json(history, str(save_dir / "training_history.json"))
    torch.save(
        {
            "model_state": model.state_dict(),
            "args": vars(args),
            "gene_names": bundle.gene_names.astype(str),
            "neighbor_profile": neighbor_profile,
        },
        save_dir / "model.pt",
    )

    result = None
    if not args.skip_eval:
        result = write_evaluation_outputs(
            output_dir=str(save_dir),
            dataset=dataset_name,
            method=args.method_name,
            seed=args.seed,
            embedding=embedding,
            labels=labels_out,
            n_clusters=n_clusters,
            n_neighbors=args.eval_neighbors,
            leiden_fixed_resolution=args.leiden_fixed_resolution,
            louvain_fixed_resolution=args.louvain_fixed_resolution,
            leiden_sweep_resolutions=parse_float_list(args.leiden_resolutions),
            include_louvain=args.include_louvain,
            run_oracle_sweep=args.run_oracle_sweep,
            sweep_max_cells=args.sweep_max_cells,
            silhouette_sample_size=args.silhouette_sample_size,
            prefix="eval",
            extra={
                "variant": args.variant_name,
                "mix_mode": args.mix_mode,
                "target_mode": args.target_mode,
                "alpha": float(args.alpha),
                "neighbor_k": int(args.neighbor_k),
                "mask_ratio": float(args.mask_ratio),
                "consistency_weight": float(args.consistency_weight),
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
        save_embedding_h5(save_dir / "embedding.h5", embedding, labels_out)

    if not args.no_save_h5ad:
        bundle.adata.obsm["X_neighbormix_mae"] = embedding
        bundle.adata.uns["neighbormix_mae"] = {
            "method": args.method_name,
            "variant": args.variant_name,
            "mix_mode": args.mix_mode,
            "target_mode": args.target_mode,
            "alpha": float(args.alpha),
            "neighbor_k": int(args.neighbor_k),
            "mask_ratio": float(args.mask_ratio),
            "consistency_weight": float(args.consistency_weight),
        }
        bundle.adata.write_h5ad(save_dir / "adata_neighbormix_mae.h5ad", compression="gzip")

    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "variant": args.variant_name,
        "seed": int(args.seed),
        "n_cells": int(data_np.shape[0]),
        "n_genes": int(data_np.shape[1]),
        "n_clusters": int(n_clusters),
        "embedding_path": str((save_dir / "embedding_final.npy").resolve()),
        "fixed_metrics": result["fixed"] if result is not None else {},
        "note": "Enhanced cells are used only as training corruptions; evaluation uses original cells only.",
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(f"Results saved to: {save_dir}")


if __name__ == "__main__":
    main()
