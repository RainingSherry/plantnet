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
from loss import gene_patch_feature_target, maskfeat_loss
from model import GenePatchMaskFeatMAE


def parse_args():
    p = argparse.ArgumentParser("rank02 MaskFeat gene patch prediction")
    p.add_argument("--data_path", required=True)
    p.add_argument("--save_dir", required=True)
    p.add_argument("--dataset_name", default=None)
    p.add_argument("--label_key", default="auto")
    p.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    p.add_argument("--n_top_genes", type=int, default=1000)
    p.add_argument("--target_sum", type=float, default=10000.0)
    p.add_argument("--scale_input", type=family.str2bool, default=True)
    p.add_argument("--n_clusters", type=int, default=0)
    p.add_argument("--hidden_size", type=int, default=128)
    p.add_argument("--depth", type=int, default=4)
    p.add_argument("--heads", type=int, default=4)
    p.add_argument("--patch_size", type=int, default=16)
    p.add_argument("--mask_ratio", type=float, default=0.75)
    p.add_argument("--epochs", type=int, default=80)
    p.add_argument("--batch_size", type=int, default=128)
    p.add_argument("--lr", type=float, default=1e-3)
    p.add_argument("--weight_decay", type=float, default=1e-4)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--gpu", type=int, default=1)
    p.add_argument("--no_cuda", action="store_true")
    p.add_argument("--skip_eval", type=family.str2bool, default=False)
    return p.parse_args()


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
    model = GenePatchMaskFeatMAE(
        num_genes=x.shape[1],
        patch_size=args.patch_size,
        hidden_size=args.hidden_size,
        depth=args.depth,
        num_heads=args.heads,
        dropout=0.1,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {"loss": [], "feature_loss": [], "reconstruction_loss": [], "mask_loss": [], "mask_rate": []}

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {key: 0.0 for key in history}
        n_batches = 0
        for _, xb_cpu, _ in loader:
            xb = xb_cpu.to(device)
            out = model(xb, mask_ratio=args.mask_ratio)
            target = gene_patch_feature_target(out["patches"])
            parts = maskfeat_loss(
                out["feature_pred"],
                target,
                out["reconstruction_pred"],
                out["patches"],
                out["mask_logits"],
                out["patch_mask"],
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            sums["loss"] += float(parts.total.detach().cpu())
            sums["feature_loss"] += float(parts.feature.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            n_batches += 1
        for key in history:
            history[key].append(sums[key] / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank02 epoch={epoch:03d}/{args.epochs} "
                f"loss={history['loss'][-1]:.4f} feature={history['feature_loss'][-1]:.4f} "
                f"mask={history['mask_rate'][-1]:.4f}",
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
    torch.save({"model": model.state_dict(), "args": vars(args)}, save_dir / "model_checkpoint.pth")

    result = None
    if not args.skip_eval:
        result = family.write_kmeans_known_k_outputs(
            save_dir,
            dataset_name,
            "rank02_maskfeat_gene_features_full",
            args.seed,
            emb,
            labels,
            n_clusters,
            {
                "rank": 2,
                "source_paper": "MaskFeat",
                "mask_semantics": "1 = masked gene patch predicted through feature target",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 2,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
