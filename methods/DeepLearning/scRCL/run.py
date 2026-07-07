#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PlantNet benchmark wrapper for the full upstream scRCL implementation.

The upstream THPengL/scRCL tree is vendored under ``scRCL_upstream``.  This
file only adapts PlantNet's benchmark CLI/data/output contract to that
upstream implementation; model, graph construction, loss functions and
clustering are imported from the vendored source.
"""

from __future__ import annotations

import argparse
import json
import os
import random
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np
import torch
from sklearn.preprocessing import LabelEncoder
from tqdm import tqdm


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[2]
METHODS_DIR = PROJECT_ROOT / "methods"
UPSTREAM_DIR = SCRIPT_DIR / "scRCL_upstream"

for path in (str(METHODS_DIR), str(PROJECT_ROOT)):
    if path not in sys.path:
        sys.path.insert(0, path)

from methods.utils import save as benchmark_save  # noqa: E402

if str(UPSTREAM_DIR) not in sys.path:
    sys.path.insert(0, str(UPSTREAM_DIR))

from clustering import clustering as scrcl_clustering  # noqa: E402
from dataset import load_h5_data  # noqa: E402
from model import Model  # noqa: E402
from utils import get_fusion, set_random_seed  # noqa: E402


LABEL_CANDIDATES = (
    "resolved_label",
    "maintype",
    "celltype",
    "cell_type",
    "Celltype",
    "cell_label",
    "label",
    "Y",
)


def register_null_h5ad_reader() -> None:
    """Allow AnnData files with null-encoded uns entries to be read."""
    try:
        import h5py
        from anndata._io.specs.registry import IOSpec, _REGISTRY

        def _read_null(*args, **kwargs):
            return None

        for typ in (h5py.Dataset, h5py.Group):
            try:
                _REGISTRY.register_read(typ, IOSpec("null", "0.1.0"))(_read_null)
            except Exception:
                pass
    except Exception:
        pass


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False
    set_random_seed(seed)


def parse_tms(value: str) -> list[int]:
    parts = [part.strip() for part in str(value).split(",") if part.strip()]
    if not parts:
        return [0, 0]
    parsed = [int(part) for part in parts]
    if len(parsed) == 1:
        parsed.append(0)
    return parsed


def infer_label_column(adata: Any, explicit: str = "auto") -> str:
    if explicit and explicit != "auto":
        if explicit not in adata.obs.columns:
            raise KeyError(f"label_key={explicit!r} not found; obs columns={list(adata.obs.columns)}")
        return explicit
    for candidate in LABEL_CANDIDATES:
        if candidate in adata.obs.columns:
            return candidate
    raise KeyError(f"No label column found; obs columns={list(adata.obs.columns)}")


def build_scrcl_input(data_path: Path, tmp_dir: Path, label_key: str) -> tuple[Path, np.ndarray, dict[str, Any]]:
    import scanpy as sc

    register_null_h5ad_reader()
    adata = sc.read_h5ad(data_path)
    resolved_label = infer_label_column(adata, label_key)
    raw_labels = adata.obs[resolved_label].astype(str).to_numpy()

    work = adata.copy()
    # Upstream load_h5_data uses cell_type1.values.codes and densifies feature
    # only on this branch, so provide its expected categorical column.
    work.obs["cell_type1"] = work.obs[resolved_label].astype(str).astype("category")
    work.uns["plantnet_label_key"] = resolved_label
    out_path = tmp_dir / "scrcl_input.h5ad"
    work.write_h5ad(out_path)
    meta = {
        "source_data_path": str(data_path),
        "label_key": resolved_label,
        "n_cells": int(work.n_obs),
        "n_genes": int(work.n_vars),
    }
    return out_path, raw_labels, meta


def make_config(args: argparse.Namespace, n_clusters: int) -> dict[str, Any]:
    out_dim = args.out_dim if args.out_dim > 0 else n_clusters
    return {
        "dataset": args.dataset_name or Path(args.data_path).stem,
        "task": "SC",
        "lr": float(args.lr),
        "epochs": int(args.epochs),
        "tms": parse_tms(args.tms),
        "dropout": float(args.dropout),
        "weight_decay": float(args.weight_decay),
        "hvg": int(args.hvg),
        "lambda1": float(args.lambda1),
        "lambda2": float(args.lambda2),
        "k": int(args.n_neighbors),
        "t": float(args.temperature),
        "hid_dim": int(args.hid_dim),
        "out_dim": int(out_dim),
        "loss_chunk_size": int(args.loss_chunk_size),
    }


def train_scrcl(
    data_path: Path,
    config: dict[str, Any],
    n_clusters: int,
    seed: int,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, list[dict[str, float]]]:
    data, _adata = load_h5_data(
        str(data_path),
        config["dataset"],
        config["hvg"],
        config["k"],
        config["tms"],
    )

    features = data["features"]
    config["n_samples"] = int(features[0].shape[0])
    config["n_classes"] = int(n_clusters)
    config["in_dim"] = int(features[0].shape[-1])
    config["gene_dim"] = int(features[0].shape[0])
    if config["out_dim"] <= 0:
        config["out_dim"] = int(n_clusters)

    adj = torch.tensor(data["adj"], dtype=torch.float32).fill_diagonal_(0).to(device)
    edge_index = torch.tensor(data["edge_index"], dtype=torch.int64).to(device)
    edge_index_g = torch.tensor(data["edge_index_g"], dtype=torch.int64).to(device)
    in_features = [torch.FloatTensor(feat).to(device) for feat in features]

    set_seed(seed)
    model = Model(
        in_dim=config["in_dim"],
        hid_dim=config["hid_dim"],
        out_dim=config["out_dim"],
        gene_dim=config["gene_dim"],
        dropout=config["dropout"],
        device=device,
    ).to(device)
    model.loss_chunk_size = int(config.get("loss_chunk_size", 0) or 0)
    optimizer = torch.optim.Adam(
        params=model.parameters(),
        lr=float(config["lr"]),
        weight_decay=float(config["weight_decay"]),
    )

    history: list[dict[str, float]] = []
    final_embedding: np.ndarray | None = None
    for epoch in tqdm(range(int(config["epochs"])), desc="scRCL", leave=False):
        model.train()
        emb_1, emb_2, z_1, z_2 = model(in_features, edge_index, edge_index_g)
        loss_hea = model.embedding_distribution_alignment(x_1=emb_1, x_2=emb_2, t=config["t"])
        loss_ndc = model.neighborhood_contrastive_alignment(x_1=emb_1, x_2=emb_2, adj=adj, t=config["t"])
        loss_cvc = model.cross_view_consistency(Z1=z_1, Z2=z_2, adj=adj)
        loss = loss_hea + config["lambda1"] * loss_ndc + config["lambda2"] * loss_cvc

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        history.append(
            {
                "epoch": float(epoch + 1),
                "loss": float(loss.detach().cpu().item()),
                "loss_hea": float(loss_hea.detach().cpu().item()),
                "loss_ndc": float(loss_ndc.detach().cpu().item()),
                "loss_cvc": float(loss_cvc.detach().cpu().item()),
            }
        )

    model.eval()
    with torch.no_grad():
        _emb_1, _emb_2, z_1, z_2 = model(in_features, edge_index, edge_index_g)
        emb_fusion = get_fusion(z1=z_1.detach(), z2=z_2.detach())
        final_embedding = emb_fusion.detach().cpu().numpy().astype(np.float32)
        label_pred, _, _ = scrcl_clustering(
            feature=emb_fusion,
            cluster_num=int(n_clusters),
            seed=int(seed),
            device=device,
        )

    return np.asarray(label_pred, dtype=np.int64), final_embedding, history


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="scRCL wrapper for the PlantNet benchmark",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, required=True)
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--epochs", type=int, default=200)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--dataset_name", type=str, default="")
    parser.add_argument("--label_key", type=str, default="auto")

    parser.add_argument("--hvg", type=int, default=2000)
    parser.add_argument("--n_neighbors", type=int, default=15)
    parser.add_argument("--tms", type=str, default="0,0")
    parser.add_argument("--temperature", type=float, default=0.1)
    parser.add_argument("--lambda1", type=float, default=100.0)
    parser.add_argument("--lambda2", type=float, default=1.0)
    parser.add_argument("--hid_dim", type=int, default=1500)
    parser.add_argument("--out_dim", type=int, default=0)
    parser.add_argument("--loss_chunk_size", type=int, default=4096)
    parser.add_argument("--dropout", type=float, default=0.3)
    parser.add_argument("--lr", type=float, default=1e-4)
    parser.add_argument("--weight_decay", type=float, default=0.1)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    save_dir = Path(args.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)

    use_cuda = not args.no_cuda and torch.cuda.is_available()
    device = torch.device(f"cuda:{args.gpu}" if use_cuda else "cpu")
    if use_cuda:
        torch.cuda.set_device(args.gpu)
    print(f"Using device: {device}")

    config = make_config(args, args.n_clusters)
    with tempfile.TemporaryDirectory(prefix="plantnet_scrcl_") as tmp:
        tmp_dir = Path(tmp)
        scrcl_input, raw_labels, input_meta = build_scrcl_input(Path(args.data_path), tmp_dir, args.label_key)
        config["hvg"] = max(1, min(int(config["hvg"]), int(input_meta["n_genes"])))
        config["k"] = max(1, min(int(config["k"]), int(input_meta["n_cells"]) - 1))

        pred, embedding, history = train_scrcl(
            data_path=scrcl_input,
            config=config,
            n_clusters=args.n_clusters,
            seed=args.seed,
            device=device,
        )

    encoder = LabelEncoder()
    y_true = encoder.fit_transform(raw_labels).astype(np.int64)
    benchmark_save(
        str(save_dir),
        y_true,
        pred,
        int(args.epochs),
        embedding,
        args=vars(args),
        preprocess_config={
            "adapter": "methods/DeepLearning/scRCL/run.py",
            "upstream": "THPengL/scRCL",
            "upstream_commit": (SCRIPT_DIR / "scRCL_upstream_commit.txt").read_text(encoding="utf-8").strip(),
            "upstream_dir": str(UPSTREAM_DIR.relative_to(PROJECT_ROOT)),
            "input": input_meta,
            "config": config,
            "checkpoint_policy": "final_epoch_no_label_selection",
        },
    )
    with (save_dir / "training_history.json").open("w", encoding="utf-8") as handle:
        json.dump(history, handle, indent=2)
    print(f"scRCL completed. Results saved to: {save_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
