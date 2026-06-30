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
from loss import ibot_distillation_loss, make_patch_mask, patch_mask_to_gene_mask
from model import IBOTOnlineTokenizerScMAE


def parse_args():
    parser = argparse.ArgumentParser("rank33 iBOT online tokenizer scMAE")
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
    parser.add_argument("--out_dim", type=int, default=256)
    parser.add_argument("--bottleneck_size", type=int, default=64)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
    parser.add_argument("--view_noise_std", type=float, default=0.03)
    parser.add_argument("--view_dropout", type=float, default=0.05)
    parser.add_argument("--student_temp", type=float, default=0.1)
    parser.add_argument("--teacher_temp", type=float, default=0.04)
    parser.add_argument("--teacher_patch_temp", type=float, default=0.07)
    parser.add_argument("--center_momentum", type=float, default=0.9)
    parser.add_argument("--ema_decay", type=float, default=0.996)
    parser.add_argument("--ema_end_decay", type=float, default=0.9995)
    parser.add_argument("--ema_anneal_epochs", type=int, default=20)
    parser.add_argument("--cls_weight", type=float, default=1.0)
    parser.add_argument("--patch_weight", type=float, default=1.0)
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


def ema_decay_for_epoch(args: argparse.Namespace, epoch: int) -> float:
    if args.ema_anneal_epochs <= 0 or epoch >= args.ema_anneal_epochs:
        return float(args.ema_end_decay)
    pct = float(epoch) / float(args.ema_anneal_epochs)
    return float(args.ema_decay) + (float(args.ema_end_decay) - float(args.ema_decay)) * pct


def augment_expression(x: torch.Tensor, noise_std: float, dropout_prob: float) -> torch.Tensor:
    out = x
    if float(dropout_prob) > 0.0:
        keep = (torch.rand_like(out) >= float(dropout_prob)).to(dtype=out.dtype)
        out = out * keep
    if float(noise_std) > 0.0:
        out = out + torch.randn_like(out) * float(noise_std)
    return out


@torch.no_grad()
def extract_features(model: IBOTOnlineTokenizerScMAE, x: np.ndarray, device: torch.device, batch_size: int) -> np.ndarray:
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
    model = IBOTOnlineTokenizerScMAE(
        num_genes=x.shape[1],
        patch_size=args.patch_size,
        hidden_size=args.hidden_size,
        depth=args.depth,
        num_heads=args.num_heads,
        dropout=args.dropout,
        out_dim=args.out_dim,
        bottleneck_size=args.bottleneck_size,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {
        "loss": [],
        "cls_loss": [],
        "patch_loss": [],
        "reconstruction_loss": [],
        "mask_loss": [],
        "mask_rate": [],
        "ema_decay": [],
        "cls_center_norm": [],
        "patch_center_norm": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {key: 0.0 for key in history if key not in {"ema_decay", "cls_center_norm", "patch_center_norm"}}
        n_batches = 0
        decay = ema_decay_for_epoch(args, epoch)
        for _, xb_cpu, _ in loader:
            clean = xb_cpu.to(device)
            view1 = augment_expression(clean, args.view_noise_std, args.view_dropout)
            view2 = augment_expression(clean, args.view_noise_std, args.view_dropout)
            mask1 = make_patch_mask(clean.shape[0], model.num_patches, args.mask_ratio, clean.device)
            mask2 = make_patch_mask(clean.shape[0], model.num_patches, args.mask_ratio, clean.device)
            gene_mask1 = patch_mask_to_gene_mask(mask1, args.patch_size, model.num_genes)
            gene_mask2 = patch_mask_to_gene_mask(mask2, args.patch_size, model.num_genes)
            outputs = model(view1, mask1, view2, mask2)
            parts = ibot_distillation_loss(
                outputs,
                clean,
                mask1,
                mask2,
                gene_mask1,
                gene_mask2,
                model.cls_center.detach(),
                model.patch_center.detach(),
                student_temp=args.student_temp,
                teacher_temp=args.teacher_temp,
                teacher_patch_temp=args.teacher_patch_temp,
                cls_weight=args.cls_weight,
                patch_weight=args.patch_weight,
                reconstruction_weight=args.reconstruction_weight,
                mask_weight=args.mask_weight,
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            model.update_teacher(decay)
            model.update_centers(outputs["teacher1"], outputs["teacher2"], args.center_momentum)
            sums["loss"] += float(parts.total.detach().cpu())
            sums["cls_loss"] += float(parts.cls.detach().cpu())
            sums["patch_loss"] += float(parts.patch.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            n_batches += 1
        for key in sums:
            history[key].append(sums[key] / max(1, n_batches))
        history["ema_decay"].append(decay)
        history["cls_center_norm"].append(float(model.cls_center.norm().detach().cpu()))
        history["patch_center_norm"].append(float(model.patch_center.norm().detach().cpu()))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank33 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"cls={history['cls_loss'][-1]:.4f} patch={history['patch_loss'][-1]:.4f} ema={decay:.5f}",
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
            "rank33_ibot_online_tokenizer_full",
            args.seed,
            emb,
            labels_true,
            n_clusters,
            {
                "rank": 33,
                "source_paper": "iBOT: Image BERT Pre-Training with Online Tokenizer",
                "mask_semantics": "1 = expression patch replaced by learned iBOT mask token",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 33,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
