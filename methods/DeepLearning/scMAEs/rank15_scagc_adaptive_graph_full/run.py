#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
os.environ.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
os.makedirs(os.environ["MPLCONFIGDIR"], exist_ok=True)
os.makedirs(os.environ["NUMBA_CACHE_DIR"], exist_ok=True)

import numpy as np
import torch
from sklearn.cluster import KMeans
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json
from loss import apply_mask_corruption, build_knn_adjacency, scagc_loss
from model import ScAGCAdaptiveGraphScMAE


class IndexedDualTargetDataset(Dataset):
    def __init__(self, x_scaled: np.ndarray, x_zinb: np.ndarray, labels: np.ndarray) -> None:
        if x_scaled.shape != x_zinb.shape:
            raise ValueError(f"x_scaled and x_zinb shapes differ: {x_scaled.shape} vs {x_zinb.shape}")
        self.x_scaled = torch.as_tensor(x_scaled, dtype=torch.float32)
        self.x_zinb = torch.as_tensor(x_zinb, dtype=torch.float32)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.x_scaled.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.x_scaled[idx], self.x_zinb[idx], self.labels[idx]


def parse_args():
    parser = argparse.ArgumentParser("rank15 scAGC adaptive graph scMAE")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--hidden_size", type=int, default=128)
    parser.add_argument("--latent_size", type=int, default=64)
    parser.add_argument("--hop_order", type=int, default=2)
    parser.add_argument("--graph_top_k", type=int, default=15)
    parser.add_argument("--graph_sigma", type=float, default=1.0)
    parser.add_argument("--graph_temperature", type=float, default=0.7)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
    parser.add_argument("--graph_weight", type=float, default=0.1)
    parser.add_argument("--zinb_weight", type=float, default=0.05)
    parser.add_argument("--contrastive_weight", type=float, default=0.1)
    parser.add_argument("--clustering_weight", type=float, default=0.05)
    parser.add_argument("--reconstruction_weight", type=float, default=0.2)
    parser.add_argument("--mask_weight", type=float, default=0.05)
    parser.add_argument("--ridge_lambda", type=float, default=0.0)
    parser.add_argument("--contrastive_temperature", type=float, default=0.7)
    parser.add_argument("--prototype_update_interval", type=int, default=5)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--weight_decay", type=float, default=1e-4)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--skip_eval", type=family.str2bool, default=False)
    return parser.parse_args()


def resolve_input_path(data_path: str, dataset_name: str) -> str:
    path = Path(data_path).resolve()
    if path.suffix.lower() != ".h5":
        return str(path)
    out_dir = CURRENT_DIR.parent / "benchmark_data"
    out_dir.mkdir(parents=True, exist_ok=True)
    converted = out_dir / f"{dataset_name}.h5ad"
    if converted.exists():
        return str(converted)
    cmd = [
        sys.executable,
        str(ROOT / "scripts" / "prepare_dataset.py"),
        "--input_path",
        str(path),
        "--dataset_name",
        dataset_name,
        "--output_dir",
        str(out_dir),
        "--force",
    ]
    result = subprocess.run(cmd, cwd=str(ROOT), text=True, capture_output=True, timeout=3600)
    if result.returncode != 0:
        raise RuntimeError(
            "Failed to convert input .h5 to .h5ad\n"
            f"Command: {' '.join(cmd)}\nSTDOUT:\n{result.stdout}\nSTDERR:\n{result.stderr}"
        )
    return str(converted)


def load_bundles(data_path: str, args: argparse.Namespace):
    scaled = family.load_scmae_dataset(
        data_path,
        args.input_mode,
        args.n_top_genes,
        args.target_sum,
        args.scale_input,
        args.label_key,
        args.seed,
    )
    non_scaled = family.load_scmae_dataset(
        data_path,
        args.input_mode,
        args.n_top_genes,
        args.target_sum,
        False,
        args.label_key,
        args.seed,
    )
    if scaled.data.shape != non_scaled.data.shape:
        raise RuntimeError("scaled and non-scaled bundles produced different shapes")
    if not np.array_equal(scaled.gene_names.astype(str), non_scaled.gene_names.astype(str)):
        raise RuntimeError("scaled and non-scaled bundles selected different HVG genes")
    return scaled, non_scaled


@torch.no_grad()
def extract_features(model: ScAGCAdaptiveGraphScMAE, x: np.ndarray, device: torch.device, batch_size: int, graph_top_k: int) -> np.ndarray:
    model.eval()
    outputs = []
    for start in range(0, x.shape[0], int(batch_size)):
        xb = torch.as_tensor(x[start:start + int(batch_size)], dtype=torch.float32, device=device)
        adjacency = build_knn_adjacency(xb, graph_top_k)
        outputs.append(model.feature(xb, adjacency).detach().cpu().numpy())
    return np.concatenate(outputs, axis=0).astype(np.float32)


def update_cluster_centers(
    model: ScAGCAdaptiveGraphScMAE,
    x: np.ndarray,
    n_clusters: int,
    device: torch.device,
    batch_size: int,
    graph_top_k: int,
    seed: int,
) -> torch.Tensor:
    emb = extract_features(model, x, device, batch_size, graph_top_k)
    labels = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(emb)
    centers = np.zeros((n_clusters, emb.shape[1]), dtype=np.float32)
    for cluster_id in range(n_clusters):
        members = emb[labels == cluster_id]
        if members.size:
            centers[cluster_id] = members.mean(axis=0)
    return torch.as_tensor(centers, dtype=torch.float32, device=device)


def main():
    args = parse_args()
    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    (save_dir / "source_manifest.json").write_text(
        (CURRENT_DIR / "source_manifest.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    device = family.get_device(args.gpu, args.no_cuda)
    dataset_name = args.dataset_name or Path(args.data_path).stem
    data_path = resolve_input_path(args.data_path, dataset_name)
    bundle, zinb_bundle = load_bundles(data_path, args)
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    preprocess_config = dict(bundle.preprocess_config)
    preprocess_config["zinb_target"] = "same HVG genes, scale_input=False, nonnegative log-normalized expression"
    preprocess_config["adaptive_graph"] = "mini-batch KNN initialization followed by RBF Gumbel-TopK adaptive graph"
    save_json(preprocess_config, str(save_dir / "preprocess_config.json"))

    x = bundle.data.astype(np.float32, copy=False)
    x_zinb = np.clip(zinb_bundle.data.astype(np.float32, copy=False), 0.0, None)
    labels_true = bundle.labels.astype(np.int64)
    n_clusters = int(args.n_clusters) if args.n_clusters > 0 else int(len(np.unique(labels_true)))

    dataset = IndexedDualTargetDataset(x, x_zinb, labels_true)
    generator = torch.Generator().manual_seed(args.seed)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    model = ScAGCAdaptiveGraphScMAE(
        num_genes=x.shape[1],
        hidden_size=args.hidden_size,
        latent_size=args.latent_size,
        hop_order=args.hop_order,
        graph_top_k=args.graph_top_k,
        graph_sigma=args.graph_sigma,
        graph_temperature=args.graph_temperature,
        dropout=args.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    centers = update_cluster_centers(
        model, x, n_clusters, device, max(512, args.batch_size * 4), args.graph_top_k, args.seed
    )
    history = {
        "loss": [],
        "graph_loss": [],
        "zinb_loss": [],
        "contrastive_loss": [],
        "clustering_loss": [],
        "reconstruction_loss": [],
        "mask_loss": [],
        "mask_rate": [],
        "adaptive_density": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch == 1 or epoch % max(1, args.prototype_update_interval) == 0:
            centers = update_cluster_centers(
                model, x, n_clusters, device, max(512, args.batch_size * 4), args.graph_top_k, args.seed + epoch
            )
        model.train()
        sums = {key: 0.0 for key in history}
        n_batches = 0
        for _, xb_cpu, xz_cpu, _ in loader:
            clean = xb_cpu.to(device)
            zinb_target = xz_cpu.to(device)
            corrupted, mask = apply_mask_corruption(clean, args.mask_ratio)
            initial_adjacency = build_knn_adjacency(clean, args.graph_top_k)
            outputs = model(corrupted, initial_adjacency)
            parts = scagc_loss(
                outputs,
                clean,
                zinb_target,
                mask,
                centers,
                graph_weight=args.graph_weight,
                zinb_weight=args.zinb_weight,
                contrastive_weight=args.contrastive_weight,
                clustering_weight=args.clustering_weight,
                reconstruction_weight=args.reconstruction_weight,
                mask_weight=args.mask_weight,
                ridge_lambda=args.ridge_lambda,
                contrastive_temperature=args.contrastive_temperature,
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            sums["loss"] += float(parts.total.detach().cpu())
            sums["graph_loss"] += float(parts.graph.detach().cpu())
            sums["zinb_loss"] += float(parts.zinb.detach().cpu())
            sums["contrastive_loss"] += float(parts.contrastive.detach().cpu())
            sums["clustering_loss"] += float(parts.clustering.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            sums["adaptive_density"] += float(parts.adaptive_density.detach().cpu())
            n_batches += 1
        for key in history:
            history[key].append(sums[key] / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank15 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"zinb={history['zinb_loss'][-1]:.4f} graph={history['graph_loss'][-1]:.4f} "
                f"cg={history['contrastive_loss'][-1]:.4f}",
                flush=True,
            )

    emb = extract_features(model, x, device, max(512, args.batch_size * 4), args.graph_top_k)
    np.save(save_dir / "embedding_final.npy", emb)
    np.save(save_dir / "labels.npy", labels_true)
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    family.save_embedding_h5(save_dir / "embedding.h5", emb, labels_true)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "args": vars(args)}, save_dir / "model_checkpoint.pth")

    result = None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            "rank15_scagc_adaptive_graph_full",
            args.seed,
            emb,
            labels_true,
            n_clusters,
            {
                "rank": 15,
                "source_paper": "scAGC",
                "mask_semantics": "1 = expression gene replaced for the auxiliary scMAE reconstruction branch",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 15,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
