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
    contrastive_sequence_loss,
    gather_patch_targets,
    make_fixed_count_mask_indices,
    mask_sc_loss,
)
from model import MaskSCClusteringModel, SequenceLevelEncoder


def parse_args():
    parser = argparse.ArgumentParser("rank13 mask-sc single-cell clustering")
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
    parser.add_argument("--decoder_size", type=int, default=128)
    parser.add_argument("--target_dim", type=int, default=64)
    parser.add_argument("--encoder_depth", type=int, default=2)
    parser.add_argument("--decoder_depth", type=int, default=2)
    parser.add_argument("--num_heads", type=int, default=4)
    parser.add_argument("--matrix_patch_size", type=int, default=4)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--mask_ratio", type=float, default=0.75)
    parser.add_argument("--variance_weight", type=float, default=0.02)
    parser.add_argument("--target_pretrain_epochs", type=int, default=5)
    parser.add_argument("--target_temperature", type=float, default=0.2)
    parser.add_argument("--target_noise_std", type=float, default=0.05)
    parser.add_argument("--target_drop_prob", type=float, default=0.15)
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=128)
    parser.add_argument("--lr", type=float, default=1e-3)
    parser.add_argument("--target_lr", type=float, default=1e-3)
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


def augment_patch_vectors(patches: torch.Tensor, noise_std: float, drop_prob: float) -> torch.Tensor:
    if patches.ndim != 2:
        raise ValueError(f"patches must be [batch, raw_patch_dim], got {tuple(patches.shape)}")
    keep = (torch.rand_like(patches) > float(drop_prob)).to(dtype=patches.dtype)
    noise = torch.randn_like(patches) * float(noise_std)
    return patches * keep + noise


def train_sequence_target_encoder(
    model: MaskSCClusteringModel,
    target_encoder: SequenceLevelEncoder,
    loader: DataLoader,
    device: torch.device,
    args: argparse.Namespace,
) -> list[float]:
    optimizer = torch.optim.AdamW(target_encoder.parameters(), lr=args.target_lr, weight_decay=args.weight_decay)
    history = []
    for epoch in range(1, max(0, args.target_pretrain_epochs) + 1):
        target_encoder.train()
        total = 0.0
        n_batches = 0
        for _, xb_cpu, _ in loader:
            clean = xb_cpu.to(device)
            raw = model.patcher.raw_patches(clean)
            pick = torch.randint(0, raw.shape[1], (raw.shape[0],), device=device)
            selected = raw[torch.arange(raw.shape[0], device=device), pick]
            view1 = augment_patch_vectors(selected, args.target_noise_std, args.target_drop_prob).unsqueeze(1)
            view2 = augment_patch_vectors(selected, args.target_noise_std, args.target_drop_prob).unsqueeze(1)
            z1 = target_encoder(view1).squeeze(1)
            z2 = target_encoder(view2).squeeze(1)
            loss = contrastive_sequence_loss(z1, z2, args.target_temperature)
            optimizer.zero_grad(set_to_none=True)
            loss.backward()
            torch.nn.utils.clip_grad_norm_(target_encoder.parameters(), 5.0)
            optimizer.step()
            total += float(loss.detach().cpu())
            n_batches += 1
        history.append(total / max(1, n_batches))
        if epoch == 1 or epoch == args.target_pretrain_epochs:
            print(f"rank13 target_pretrain={epoch:03d}/{args.target_pretrain_epochs} loss={history[-1]:.4f}", flush=True)
    target_encoder.eval()
    for param in target_encoder.parameters():
        param.requires_grad_(False)
    return history


@torch.no_grad()
def extract_features(model: MaskSCClusteringModel, x: np.ndarray, device: torch.device, batch_size: int) -> np.ndarray:
    model.eval()
    outputs = []
    for start in range(0, x.shape[0], int(batch_size)):
        xb = torch.as_tensor(x[start:start + int(batch_size)], dtype=torch.float32, device=device)
        outputs.append(model.feature(xb).detach().cpu().numpy())
    return np.concatenate(outputs, axis=0).astype(np.float32)


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
    model = MaskSCClusteringModel(
        num_genes=x.shape[1],
        hidden_size=args.hidden_size,
        decoder_size=args.decoder_size,
        target_dim=args.target_dim,
        encoder_depth=args.encoder_depth,
        decoder_depth=args.decoder_depth,
        num_heads=args.num_heads,
        matrix_patch_size=args.matrix_patch_size,
        dropout=args.dropout,
    ).to(device)
    target_encoder = SequenceLevelEncoder(
        raw_patch_dim=model.patcher.raw_patch_dim,
        target_dim=args.target_dim,
        hidden_size=args.hidden_size,
        dropout=args.dropout,
    ).to(device)
    target_history = train_sequence_target_encoder(model, target_encoder, loader, device, args)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {
        "loss": [],
        "reconstruction_loss": [],
        "token_variance_loss": [],
        "mask_rate": [],
        "target_pretrain_loss": target_history,
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {"loss": 0.0, "reconstruction_loss": 0.0, "token_variance_loss": 0.0, "mask_rate": 0.0}
        n_batches = 0
        for _, xb_cpu, _ in loader:
            clean = xb_cpu.to(device)
            visible_indices, masked_indices, patch_mask = make_fixed_count_mask_indices(
                clean.shape[0], model.num_patches, args.mask_ratio, clean.device
            )
            with torch.no_grad():
                all_target_features = target_encoder(model.patcher.raw_patches(clean))
                masked_target_features = gather_patch_targets(all_target_features, masked_indices)
            outputs = model(clean, visible_indices, masked_indices)
            parts = mask_sc_loss(outputs, masked_target_features, patch_mask, variance_weight=args.variance_weight)
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            sums["loss"] += float(parts.total.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["token_variance_loss"] += float(parts.token_variance.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank13 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"recon={history['reconstruction_loss'][-1]:.4f} var={history['token_variance_loss'][-1]:.4f}",
                flush=True,
            )

    emb = extract_features(model, x, device, max(512, args.batch_size * 4))
    np.save(save_dir / "embedding_final.npy", emb)
    np.save(save_dir / "labels.npy", labels_true)
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))
    family.save_embedding_h5(save_dir / "embedding.h5", emb, labels_true)
    save_json(history, str(save_dir / "training_history.json"))
    torch.save({"model": model.state_dict(), "target_encoder": target_encoder.state_dict(), "args": vars(args)}, save_dir / "model_checkpoint.pth")

    result = None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            "rank13_masked_sc_clustering_full",
            args.seed,
            emb,
            labels_true,
            n_clusters,
            {
                "rank": 13,
                "source_paper": "Masked Modeling for Single-cell Clustering of scRNA-seq Data",
                "mask_semantics": "1 = expression matrix patch removed from encoder and predicted by sequence-guided decoder",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 13,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
