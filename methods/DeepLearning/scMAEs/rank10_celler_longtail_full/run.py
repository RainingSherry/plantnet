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
from torch.utils.data import DataLoader

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json
from loss import apply_mask_corruption, celler_longtail_loss
from model import CellerLongTailScMAE


def parse_args():
    parser = argparse.ArgumentParser("rank10 Celler long-tail scMAE")
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
    parser.add_argument("--depth", type=int, default=2)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
    parser.add_argument("--prototype_update_interval", type=int, default=5)
    parser.add_argument("--ginf_sigma", type=float, default=0.1)
    parser.add_argument("--rare_strength", type=float, default=0.5)
    parser.add_argument("--hard_strength", type=float, default=0.25)
    parser.add_argument("--mask_weight", type=float, default=0.1)
    parser.add_argument("--ginf_weight", type=float, default=0.15)
    parser.add_argument("--compactness_weight", type=float, default=0.1)
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


@torch.no_grad()
def extract_features(model: CellerLongTailScMAE, x: np.ndarray, device: torch.device, batch_size: int) -> np.ndarray:
    model.eval()
    outputs = []
    for start in range(0, x.shape[0], int(batch_size)):
        xb = torch.as_tensor(x[start:start + int(batch_size)], dtype=torch.float32, device=device)
        outputs.append(model.feature(xb).detach().cpu().numpy())
    return np.concatenate(outputs, axis=0).astype(np.float32)


def update_pseudo_state(
    model: CellerLongTailScMAE,
    x: np.ndarray,
    n_clusters: int,
    device: torch.device,
    batch_size: int,
    seed: int,
) -> tuple[np.ndarray, torch.Tensor, torch.Tensor]:
    emb = extract_features(model, x, device, batch_size)
    labels = KMeans(n_clusters=n_clusters, n_init=20, random_state=seed).fit_predict(emb).astype(np.int64)
    counts_np = np.bincount(labels, minlength=n_clusters).astype(np.float32)
    centers = np.zeros((n_clusters, emb.shape[1]), dtype=np.float32)
    for cluster_id in range(n_clusters):
        members = emb[labels == cluster_id]
        if members.size:
            centers[cluster_id] = members.mean(axis=0)
    return labels, torch.as_tensor(counts_np, dtype=torch.float32, device=device), torch.as_tensor(centers, dtype=torch.float32, device=device)


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
    bundle = family.load_scmae_dataset(
        data_path,
        args.input_mode,
        args.n_top_genes,
        args.target_sum,
        args.scale_input,
        args.label_key,
        args.seed,
    )
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    x = bundle.data.astype(np.float32, copy=False)
    labels_true = bundle.labels.astype(np.int64)
    n_clusters = int(args.n_clusters) if args.n_clusters > 0 else int(len(np.unique(labels_true)))

    dataset = family.IndexedExpressionDataset(x, labels_true)
    generator = torch.Generator().manual_seed(args.seed)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    model = CellerLongTailScMAE(
        num_genes=x.shape[1],
        n_prototypes=n_clusters,
        hidden_size=args.hidden_size,
        depth=args.depth,
        dropout=args.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    pseudo_labels, class_counts, prototype_centers = update_pseudo_state(
        model, x, n_clusters, device, max(512, args.batch_size * 4), args.seed
    )
    history = {
        "loss": [],
        "reconstruction_loss": [],
        "mask_loss": [],
        "ginf_loss": [],
        "compactness_loss": [],
        "hard_weight_mean": [],
        "mask_rate": [],
        "pseudo_min_count": [],
        "pseudo_max_count": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        if epoch == 1 or epoch % max(1, args.prototype_update_interval) == 0:
            pseudo_labels, class_counts, prototype_centers = update_pseudo_state(
                model, x, n_clusters, device, max(512, args.batch_size * 4), args.seed + epoch
            )
        model.train()
        sums = {key: 0.0 for key in history if not key.startswith("pseudo_")}
        n_batches = 0
        for idx_cpu, xb_cpu, _ in loader:
            idx = idx_cpu.numpy()
            clean = xb_cpu.to(device)
            corrupted, mask = apply_mask_corruption(clean, args.mask_ratio)
            outputs = model(corrupted)
            batch_pseudo = torch.as_tensor(pseudo_labels[idx], dtype=torch.long, device=device)
            parts = celler_longtail_loss(
                outputs,
                clean,
                mask,
                batch_pseudo,
                class_counts,
                prototype_centers,
                ginf_sigma=args.ginf_sigma,
                rare_strength=args.rare_strength,
                hard_strength=args.hard_strength,
                mask_weight=args.mask_weight,
                ginf_weight=args.ginf_weight,
                compactness_weight=args.compactness_weight,
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            sums["loss"] += float(parts.total.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["ginf_loss"] += float(parts.ginf.detach().cpu())
            sums["compactness_loss"] += float(parts.compactness.detach().cpu())
            sums["hard_weight_mean"] += float(parts.hard_weight_mean.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        counts_cpu = class_counts.detach().cpu().numpy()
        history["pseudo_min_count"].append(float(counts_cpu.min()))
        history["pseudo_max_count"].append(float(counts_cpu.max()))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank10 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"recon={history['reconstruction_loss'][-1]:.4f} ginf={history['ginf_loss'][-1]:.4f} "
                f"count_range={counts_cpu.min():.0f}-{counts_cpu.max():.0f}",
                flush=True,
            )

    emb = extract_features(model, x, device, max(512, args.batch_size * 4))
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
            "rank10_celler_longtail_full",
            args.seed,
            emb,
            labels_true,
            n_clusters,
            {
                "rank": 10,
                "source_paper": "Celler",
                "mask_semantics": "1 = expression gene replaced for masked reconstruction branch",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 10,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
