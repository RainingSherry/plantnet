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
from loss import apply_feature_mask, build_knn_adjacency, drop_edges, graph_barlow_loss
from model import GraphBarlowTwinsScMAE


def parse_args():
    parser = argparse.ArgumentParser("rank27 Graph Barlow Twins scMAE")
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
    parser.add_argument("--projection_size", type=int, default=64)
    parser.add_argument("--projector_hidden", type=int, default=256)
    parser.add_argument("--dropout", type=float, default=0.1)
    parser.add_argument("--graph_top_k", type=int, default=15)
    parser.add_argument("--edge_drop", type=float, default=0.2)
    parser.add_argument("--feature_drop_prob", type=float, default=0.25)
    parser.add_argument("--element_mask_ratio", type=float, default=0.1)
    parser.add_argument("--lambda_offdiag", type=float, default=-1.0)
    parser.add_argument("--barlow_eps", type=float, default=1e-5)
    parser.add_argument("--barlow_weight", type=float, default=0.02)
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


@torch.no_grad()
def extract_features(model: GraphBarlowTwinsScMAE, x: np.ndarray, device: torch.device, batch_size: int, graph_top_k: int) -> np.ndarray:
    model.eval()
    outputs = []
    for start in range(0, x.shape[0], int(batch_size)):
        xb = torch.as_tensor(x[start:start + int(batch_size)], dtype=torch.float32, device=device)
        adjacency = build_knn_adjacency(xb, graph_top_k)
        outputs.append(model.feature(xb, adjacency).detach().cpu().numpy())
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
    preprocess_config = dict(bundle.preprocess_config)
    preprocess_config["graph_barlow"] = "mini-batch KNN graph with two feature/edge dropout views"
    preprocess_config["mask_semantics"] = "1 = expression feature zeroed by Graph Barlow feature corruption"
    save_json(preprocess_config, str(save_dir / "preprocess_config.json"))
    x = bundle.data.astype(np.float32, copy=False)
    labels_true = bundle.labels.astype(np.int64)
    n_clusters = int(args.n_clusters) if args.n_clusters > 0 else int(len(np.unique(labels_true)))

    dataset = family.IndexedExpressionDataset(x, labels_true)
    generator = torch.Generator().manual_seed(args.seed)
    loader = DataLoader(dataset, batch_size=args.batch_size, shuffle=True, drop_last=True, generator=generator)
    model = GraphBarlowTwinsScMAE(
        num_genes=x.shape[1],
        hidden_size=args.hidden_size,
        latent_size=args.latent_size,
        projection_size=args.projection_size,
        projector_hidden=args.projector_hidden,
        dropout=args.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    lambda_offdiag = None if args.lambda_offdiag < 0.0 else float(args.lambda_offdiag)
    history = {
        "loss": [],
        "barlow_loss": [],
        "on_diag_loss": [],
        "off_diag_loss": [],
        "reconstruction_loss": [],
        "mask_loss": [],
        "mask_rate": [],
        "edge_density1": [],
        "edge_density2": [],
    }

    for epoch in range(1, max(1, args.epochs) + 1):
        model.train()
        sums = {key: 0.0 for key in history}
        n_batches = 0
        for _, xb_cpu, _ in loader:
            clean = xb_cpu.to(device)
            view1, mask1 = apply_feature_mask(clean, args.feature_drop_prob, args.element_mask_ratio)
            view2, mask2 = apply_feature_mask(clean, args.feature_drop_prob, args.element_mask_ratio)
            adjacency1 = drop_edges(build_knn_adjacency(view1, args.graph_top_k), args.edge_drop)
            adjacency2 = drop_edges(build_knn_adjacency(view2, args.graph_top_k), args.edge_drop)
            outputs = model(view1, adjacency1, view2, adjacency2)
            parts = graph_barlow_loss(
                outputs,
                clean,
                mask1,
                mask2,
                adjacency1,
                adjacency2,
                lambda_offdiag=lambda_offdiag,
                eps=args.barlow_eps,
                barlow_weight=args.barlow_weight,
                reconstruction_weight=args.reconstruction_weight,
                mask_weight=args.mask_weight,
            )
            optimizer.zero_grad(set_to_none=True)
            parts.total.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 5.0)
            optimizer.step()
            sums["loss"] += float(parts.total.detach().cpu())
            sums["barlow_loss"] += float(parts.barlow.detach().cpu())
            sums["on_diag_loss"] += float(parts.on_diag.detach().cpu())
            sums["off_diag_loss"] += float(parts.off_diag.detach().cpu())
            sums["reconstruction_loss"] += float(parts.reconstruction.detach().cpu())
            sums["mask_loss"] += float(parts.mask.detach().cpu())
            sums["mask_rate"] += float(parts.mask_rate.detach().cpu())
            sums["edge_density1"] += float(parts.edge_density1.detach().cpu())
            sums["edge_density2"] += float(parts.edge_density2.detach().cpu())
            n_batches += 1
        for key in history:
            history[key].append(sums[key] / max(1, n_batches))
        if epoch == 1 or epoch == args.epochs or epoch % 10 == 0:
            print(
                f"rank27 epoch={epoch:03d}/{args.epochs} loss={history['loss'][-1]:.4f} "
                f"barlow={history['barlow_loss'][-1]:.4f} mask={history['mask_rate'][-1]:.3f}",
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
            "rank27_graph_barlow_twins_full",
            args.seed,
            emb,
            labels_true,
            n_clusters,
            {
                "rank": 27,
                "source_paper": "Graph Barlow Twins",
                "mask_semantics": "1 = expression feature zeroed by graph-view corruption",
            },
        )
        save_json(result["fixed"], str(save_dir / "metrics.json"))
    save_json(
        {
            "dataset": dataset_name,
            "rank": 27,
            "source_manifest": str((CURRENT_DIR / "source_manifest.json").resolve()),
            "fixed_metrics": result["fixed"] if result else {},
        },
        str(save_dir / "summary.json"),
    )


if __name__ == "__main__":
    main()
