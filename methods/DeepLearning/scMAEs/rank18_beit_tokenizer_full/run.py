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
from torch.utils.data import DataLoader, Dataset

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json
from loss import beit_expression_loss, make_patch_mask, patch_mask_to_gene_mask
from model import ExpressionBEiTScMAE


class IndexedTokenDataset(Dataset):
    def __init__(self, x_scaled: np.ndarray, patch_tokens: np.ndarray, labels: np.ndarray) -> None:
        if x_scaled.shape[0] != patch_tokens.shape[0]:
            raise ValueError("x_scaled and patch_tokens must have the same number of cells")
        self.x_scaled = torch.as_tensor(x_scaled, dtype=torch.float32)
        self.patch_tokens = torch.as_tensor(patch_tokens, dtype=torch.long)
        self.labels = torch.as_tensor(labels, dtype=torch.long)

    def __len__(self) -> int:
        return int(self.x_scaled.shape[0])

    def __getitem__(self, idx: int):
        return int(idx), self.x_scaled[idx], self.patch_tokens[idx], self.labels[idx]


def parse_args():
    parser = argparse.ArgumentParser("rank18 BEiT expression tokenizer scMAE")
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
    parser.add_argument("--num_heads", type=int, default=4)
    parser.add_argument("--patch_size", type=int, default=20)
    parser.add_argument("--vocab_size", type=int, default=32)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--mask_ratio", type=float, default=0.4)
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


def patch_means(x: np.ndarray, patch_size: int) -> np.ndarray:
    cells, genes = x.shape
    num_patches = int(np.ceil(genes / patch_size))
    pad = num_patches * patch_size - genes
    padded = np.pad(x, ((0, 0), (0, pad)), mode="constant") if pad else x
    return padded.reshape(cells, num_patches, patch_size).mean(axis=2)


def build_token_edges(values: np.ndarray, vocab_size: int) -> np.ndarray:
    quantiles = np.linspace(0.0, 1.0, int(vocab_size) + 1, dtype=np.float64)[1:-1]
    edges = np.quantile(values.reshape(-1), quantiles).astype(np.float32)
    return np.maximum.accumulate(edges)


def tokenize_patches(values: np.ndarray, edges: np.ndarray) -> np.ndarray:
    return np.searchsorted(edges, values, side="right").astype(np.int64)


@torch.no_grad()
def extract_features(model: ExpressionBEiTScMAE, x: np.ndarray, device: torch.device, batch_size: int) -> np.ndarray:
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
    bundle, token_bundle = load_bundles(data_path, args)
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    preprocess_config = dict(bundle.preprocess_config)
    preprocess_config["beit_tokenizer"] = "same HVG genes, scale_input=False, patch mean quantile tokens"
    save_json(preprocess_config, str(save_dir / "preprocess_config.json"))

    x = bundle.data.astype(np.float32, copy=False)
    token_source = np.clip(token_bundle.data.astype(np.float32, copy=False), 0.0, None)
    token_values = patch_means(token_source, args.patch_size)
    token_edges = build_token_edges(token_values, args.vocab_size)
    patch_tokens = tokenize_patches(token_values, token_edges)
    np.save(save_dir / "token_edges.npy", token_edges)
    labels_true = bundle.labels.astype(np.int64)
    n_clusters = int(args.n_clusters) if args.n_clusters > 0 else int(len(np.unique(labels_true)))

    dataset = IndexedTokenDataset(x, patch_tokens, labels_true)
    generator = torch.Generator().manual_seed(args.seed)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    model = ExpressionBEiTScMAE(
        num_genes=x.shape[1],
        vocab_size=args.vocab_size,
        patch_size=args.patch_size,
        hidden_size=args.hidden_size,
        depth=args.depth,
        num_heads=args.num_heads,
        dropout=args.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    history = {
        "loss": [],
        "token_loss": [],
        "reconstruction_loss": [],
        "mask_loss": [],
        "mask_rate": [],
        "token_accuracy": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {key: 0.0 for key in history}
        n_batches = 0
        for _, xb_cpu, token_cpu, _ in loader:
            clean = xb_cpu.to(device)
            patch_targets = token_cpu.to(device)
            patch_mask = make_patch_mask(clean.shape[0], model.num_patches, args.mask_ratio, clean.device)
            gene_mask = patch_mask_to_gene_mask(patch_mask, args.patch_size, model.num_genes)
            outputs = model(clean, patch_mask)
            parts = beit_expression_loss(
                outputs,
                clean,
                patch_targets,
                patch_mask,
                gene_mask,
                reconstruction_weight=args.reconstruction_weight,
                mask_weight=args.mask_weight,
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            sums["loss"] += float(parts.total.detach().cpu())
            sums["token_loss"] += float(parts.token.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            sums["token_accuracy"] += float(parts.token_accuracy.detach().cpu())
            n_batches += 1
        for key in history:
            history[key].append(sums[key] / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank18 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"token={history['token_loss'][-1]:.4f} acc={history['token_accuracy'][-1]:.4f}",
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
            "rank18_beit_tokenizer_full",
            args.seed,
            emb,
            labels_true,
            n_clusters,
            {
                "rank": 18,
                "source_paper": "BEiT",
                "mask_semantics": "1 = expression patch replaced by learned BEiT mask token",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 18,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
