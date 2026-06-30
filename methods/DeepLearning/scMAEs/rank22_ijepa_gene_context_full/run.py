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
from loss import ijepa_loss, make_context_target_masks, patch_mask_to_gene_mask
from model import IJEPAGeneContextScMAE


def parse_args():
    parser = argparse.ArgumentParser("rank22 I-JEPA gene-context scMAE")
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
    parser.add_argument("--num_heads", type=int, default=4)
    parser.add_argument("--patch_size", type=int, default=20)
    parser.add_argument("--predictor_size", type=int, default=64)
    parser.add_argument("--predictor_depth", type=int, default=2)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--target_ratio", type=float, default=0.35)
    parser.add_argument("--context_ratio", type=float, default=0.55)
    parser.add_argument("--loss_beta", type=float, default=1.0)
    parser.add_argument("--ema_decay", type=float, default=0.996)
    parser.add_argument("--ema_end_decay", type=float, default=0.999)
    parser.add_argument("--ema_anneal_epochs", type=int, default=20)
    parser.add_argument("--reconstruction_weight", type=float, default=0.2)
    parser.add_argument("--mask_weight", type=float, default=0.05)
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


def ema_momentum_for_epoch(args: argparse.Namespace, epoch: int) -> float:
    if args.ema_anneal_epochs <= 0 or epoch >= args.ema_anneal_epochs:
        return float(args.ema_end_decay)
    pct = float(epoch) / float(args.ema_anneal_epochs)
    return float(args.ema_decay) + (float(args.ema_end_decay) - float(args.ema_decay)) * pct


@torch.no_grad()
def extract_features(model: IJEPAGeneContextScMAE, x: np.ndarray, device: torch.device, batch_size: int) -> np.ndarray:
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
    model = IJEPAGeneContextScMAE(
        num_genes=x.shape[1],
        patch_size=args.patch_size,
        hidden_size=args.hidden_size,
        depth=args.depth,
        num_heads=args.num_heads,
        predictor_size=args.predictor_size,
        predictor_depth=args.predictor_depth,
        dropout=args.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {
        "loss": [],
        "latent_loss": [],
        "reconstruction_loss": [],
        "mask_loss": [],
        "target_patch_rate": [],
        "ema_momentum": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {key: 0.0 for key in history if key != "ema_momentum"}
        n_batches = 0
        momentum = ema_momentum_for_epoch(args, epoch)
        for _, xb_cpu, _ in loader:
            clean = xb_cpu.to(device)
            context_mask, target_mask = make_context_target_masks(
                clean.shape[0],
                model.num_patches,
                args.target_ratio,
                args.context_ratio,
                clean.device,
            )
            gene_target_mask = patch_mask_to_gene_mask(target_mask, args.patch_size, model.num_genes)
            target_features = model.target_features(clean)
            outputs = model(clean, context_mask, target_mask)
            parts = ijepa_loss(
                outputs,
                clean,
                target_features,
                target_mask,
                gene_target_mask,
                beta=args.loss_beta,
                reconstruction_weight=args.reconstruction_weight,
                mask_weight=args.mask_weight,
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            model.update_target(momentum)
            sums["loss"] += float(parts.total.detach().cpu())
            sums["latent_loss"] += float(parts.latent.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["target_patch_rate"] += float(parts.target_patch_rate.detach().cpu())
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        history["ema_momentum"].append(momentum)
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank22 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"latent={history['latent_loss'][-1]:.4f} ema={momentum:.5f}",
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
            "rank22_ijepa_gene_context_full",
            args.seed,
            emb,
            labels_true,
            n_clusters,
            {
                "rank": 22,
                "source_paper": "I-JEPA",
                "mask_semantics": "1 = target gene patch predicted by predictor from non-overlapping context patches",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 22,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
