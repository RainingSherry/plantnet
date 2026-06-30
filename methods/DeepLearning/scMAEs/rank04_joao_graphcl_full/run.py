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
from torch.utils.data import DataLoader

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json
from loss import (
    apply_graph_augmentation,
    joao_scmae_loss,
    joao_update_probabilities,
    nt_xent_loss,
    sample_augmentation_pair,
)
from model import JOAOScMAEGraphEncoder, build_knn_adjacency


def parse_args():
    parser = argparse.ArgumentParser("rank04 JOAO GraphCL scMAE")
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
    parser.add_argument("--depth", type=int, default=3)
    parser.add_argument("--projection_size", type=int, default=128)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--knn_k", type=int, default=15)
    parser.add_argument("--aug_ratio", type=float, default=0.2)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
    parser.add_argument("--temperature", type=float, default=0.2)
    parser.add_argument("--reconstruction_weight", type=float, default=0.4)
    parser.add_argument("--mask_weight", type=float, default=0.1)
    parser.add_argument("--joao_beta", type=float, default=0.1)
    parser.add_argument("--joao_gamma", type=float, default=0.1)
    parser.add_argument("--joao_update_interval", type=int, default=5)
    parser.add_argument("--joao_eval_batches", type=int, default=8)
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


def random_expression_mask(x: torch.Tensor, mask_ratio: float) -> tuple[torch.Tensor, torch.Tensor]:
    mask = (torch.rand_like(x) < float(mask_ratio)).float()
    return x * (1.0 - mask), mask


@torch.no_grad()
def evaluate_augmentation_losses(
    model: JOAOScMAEGraphEncoder,
    loader: DataLoader,
    device: torch.device,
    args: argparse.Namespace,
) -> np.ndarray:
    model.eval()
    losses = np.zeros(5, dtype=np.float64)
    counts = np.zeros(5, dtype=np.float64)
    for batch_index, (_, xb_cpu, _) in enumerate(loader):
        if batch_index >= args.joao_eval_batches:
            break
        xb = xb_cpu.to(device)
        adj = build_knn_adjacency(xb, args.knn_k)
        for aug_id in range(5):
            x1, a1, _ = apply_graph_augmentation(xb, adj, aug_id, args.aug_ratio)
            x2, a2, _ = apply_graph_augmentation(xb, adj, aug_id, args.aug_ratio)
            z1 = model(x1, a1)["projection"]
            z2 = model(x2, a2)["projection"]
            losses[aug_id] += float(nt_xent_loss(z1, z2, args.temperature).detach().cpu())
            counts[aug_id] += 1.0
    if np.any(counts == 0):
        raise RuntimeError("JOAO augmentation evaluation received no batches")
    return losses / counts


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
    labels = bundle.labels.astype(np.int64)
    n_clusters = int(args.n_clusters) if args.n_clusters > 0 else int(len(np.unique(labels)))

    dataset = family.IndexedExpressionDataset(x, labels)
    generator = torch.Generator().manual_seed(args.seed)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    model = JOAOScMAEGraphEncoder(
        num_genes=x.shape[1],
        hidden_size=args.hidden_size,
        depth=args.depth,
        projection_size=args.projection_size,
        dropout=args.dropout,
        knn_k=args.knn_k,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    rng = np.random.default_rng(args.seed)
    aug_prob = np.ones(5, dtype=np.float64) / 5.0
    history = {
        "loss": [],
        "contrastive_loss": [],
        "reconstruction_loss": [],
        "mask_loss": [],
        "mask_rate": [],
        "aug_prob": [],
        "aug_eval_loss": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {"loss": 0.0, "contrastive_loss": 0.0, "reconstruction_loss": 0.0, "mask_loss": 0.0, "mask_rate": 0.0}
        n_batches = 0
        for _, xb_cpu, _ in loader:
            xb = xb_cpu.to(device)
            adj = build_knn_adjacency(xb, args.knn_k)
            aug1, aug2 = sample_augmentation_pair(aug_prob, rng)
            x1, a1, _ = apply_graph_augmentation(xb, adj, aug1, args.aug_ratio)
            x2, a2, _ = apply_graph_augmentation(xb, adj, aug2, args.aug_ratio)
            x_rec, rec_mask = random_expression_mask(xb, args.mask_ratio)
            out1 = model(x1, a1)
            out2 = model(x2, a2)
            out_rec = model(x_rec, adj)
            parts = joao_scmae_loss(
                out1,
                out2,
                out_rec,
                xb,
                rec_mask,
                temperature=args.temperature,
                reconstruction_weight=args.reconstruction_weight,
                mask_weight=args.mask_weight,
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            sums["loss"] += float(parts.total.detach().cpu())
            sums["contrastive_loss"] += float(parts.contrastive.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))

        aug_eval = np.full(5, np.nan, dtype=np.float64)
        if epoch == 1 or epoch % max(1, args.joao_update_interval) == 0:
            aug_eval = evaluate_augmentation_losses(model, loader, device, args)
            aug_prob = joao_update_probabilities(aug_prob, aug_eval, args.joao_beta, args.joao_gamma)
        history["aug_prob"].append(aug_prob.tolist())
        history["aug_eval_loss"].append(aug_eval.tolist())
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank04 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"contrast={history['contrastive_loss'][-1]:.4f} recon={history['reconstruction_loss'][-1]:.4f} "
                f"augP={','.join(f'{p:.2f}' for p in aug_prob)}",
                flush=True,
            )

    model.eval()
    embeddings = []
    eval_batch = max(512, args.batch_size * 4)
    with torch.no_grad():
        for start in range(0, x.shape[0], eval_batch):
            xb = torch.as_tensor(x[start:start + eval_batch], dtype=torch.float32, device=device)
            embeddings.append(model.feature(xb).detach().cpu().numpy())
    emb = np.concatenate(embeddings, axis=0).astype(np.float32)
    np.save(save_dir / "embedding_final.npy", emb)
    np.save(save_dir / "labels.npy", labels)
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    family.save_embedding_h5(save_dir / "embedding.h5", emb, labels)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "args": vars(args), "aug_prob": aug_prob.tolist()}, save_dir / "model_checkpoint.pth")

    result = None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            "rank04_joao_graphcl_full",
            args.seed,
            emb,
            labels,
            n_clusters,
            {
                "rank": 4,
                "source_paper": "JOAO",
                "mask_semantics": "1 = expression gene masked in scMAE reconstruction branch",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 4,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
            "final_aug_prob": aug_prob.tolist(),
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
