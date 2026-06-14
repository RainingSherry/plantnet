#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import torch
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import normalize

CURRENT_DIR = Path(__file__).resolve().parent
ROOT = CURRENT_DIR.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning import scMAE_family as family
from methods.shared_utils import ensure_dir, save_json


SCGPT_DIR = CURRENT_DIR
SCGPT_MODEL_FILE = SCGPT_DIR / "best_model.pt"
SCGPT_VOCAB_FILE = SCGPT_DIR / "vocab.json"
SCGPT_ARGS_FILE = SCGPT_DIR / "args.json"
FORBIDDEN_GPUS = {0, 7}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PlantNet-safe scGPT embedding runner.")
    parser.add_argument("--data_path", required=True)
    parser.add_argument("--save_dir", required=True)
    parser.add_argument("--dataset_name", default=None)
    parser.add_argument("--method_name", default="scGPT")
    parser.add_argument("--label_key", default="auto")
    parser.add_argument("--input_mode", default="auto", choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--scale_input", type=family.str2bool, default=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--pca_dim", type=int, default=50)
    parser.add_argument("--batch_size", type=int, default=64)
    parser.add_argument("--max_length", type=int, default=1200)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    return parser.parse_args()


def get_device(gpu: int, no_cuda: bool) -> torch.device:
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu")
    if int(gpu) in FORBIDDEN_GPUS:
        raise ValueError("GPU 0 and 7 are forbidden for this experiment.")
    visible = [x.strip() for x in str(__import__("os").environ.get("CUDA_VISIBLE_DEVICES", "")).split(",") if x.strip()]
    if visible:
        if any(int(x) in FORBIDDEN_GPUS for x in visible if x.lstrip("-").isdigit()):
            raise ValueError("CUDA_VISIBLE_DEVICES contains forbidden GPU 0 or 7.")
        return torch.device("cuda:0" if str(gpu) in visible and len(visible) == 1 else f"cuda:{visible.index(str(gpu))}")
    return torch.device(f"cuda:{int(gpu)}")


def load_vocab() -> dict[str, int]:
    if not SCGPT_VOCAB_FILE.exists():
        return {}
    with SCGPT_VOCAB_FILE.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def try_scgpt_embedding(data_np: np.ndarray, gene_names: np.ndarray, args: argparse.Namespace, device: torch.device):
    """Best-effort scGPT embedding.

    In this plant project the bundled scGPT folder contains the human vocabulary
    but not the checkpoint. Most datasets use plant gene identifiers, so true
    transformer inference is expected to be unavailable. We still keep this path
    explicit so future checkpoints can be dropped into the same folder.
    """
    vocab = load_vocab()
    matched = [str(g) for g in gene_names.astype(str) if str(g) in vocab]
    profile = {
        "scgpt_package_available": False,
        "checkpoint_exists": SCGPT_MODEL_FILE.exists(),
        "vocab_exists": SCGPT_VOCAB_FILE.exists(),
        "args_exists": SCGPT_ARGS_FILE.exists(),
        "n_input_genes": int(len(gene_names)),
        "n_vocab_genes": int(len(vocab)),
        "n_matched_genes": int(len(matched)),
        "matched_gene_fraction": float(len(matched) / max(1, len(gene_names))),
        "used_scgpt_transformer": False,
        "fallback_reason": "",
    }
    try:
        import scgpt  # noqa: F401
        profile["scgpt_package_available"] = True
    except Exception as exc:
        profile["fallback_reason"] = f"scgpt package unavailable: {exc}"
        return None, profile

    if not SCGPT_MODEL_FILE.exists():
        profile["fallback_reason"] = f"missing checkpoint: {SCGPT_MODEL_FILE}"
        return None, profile
    if len(matched) == 0:
        profile["fallback_reason"] = "no input genes matched bundled scGPT vocabulary"
        return None, profile

    profile["fallback_reason"] = "transformer path not enabled in this safe runner"
    return None, profile


def pca_fallback_embedding(data_np: np.ndarray, pca_dim: int, seed: int) -> np.ndarray:
    dim = max(1, min(int(pca_dim), data_np.shape[1], data_np.shape[0] - 1))
    emb = PCA(n_components=dim, random_state=seed).fit_transform(data_np)
    return normalize(np.nan_to_num(emb, nan=0.0, posinf=0.0, neginf=0.0), axis=1).astype(np.float32)


def main() -> int:
    args = parse_args()
    family.set_seed(args.seed)
    save_dir = Path(ensure_dir(args.save_dir))
    save_json(vars(args), str(save_dir / "args.json"))
    device = get_device(args.gpu, args.no_cuda)
    print(f"Using device: {device}", flush=True)

    bundle = family.load_scmae_dataset(
        file_path=args.data_path,
        input_mode=args.input_mode,
        n_top_genes=args.n_top_genes,
        target_sum=args.target_sum,
        scale_input=args.scale_input,
        label_key=args.label_key,
        seed=args.seed,
    )
    data_np = bundle.data.astype(np.float32)
    labels = bundle.labels.astype(np.int64)
    dataset_name = args.dataset_name or Path(args.data_path).stem
    n_clusters = int(args.n_clusters) if args.n_clusters and args.n_clusters > 0 else int(len(np.unique(labels)))
    save_json(bundle.profile, str(save_dir / "dataset_profile.json"))
    save_json(bundle.preprocess_config, str(save_dir / "preprocess_config.json"))
    np.save(save_dir / "gene_names.npy", bundle.gene_names.astype(str))

    embedding, scgpt_profile = try_scgpt_embedding(data_np, bundle.gene_names, args, device)
    if embedding is None:
        print(f"[scGPT fallback] {scgpt_profile['fallback_reason']}", flush=True)
        embedding = pca_fallback_embedding(data_np, args.pca_dim, args.seed)
        embedding_source = "pca_fallback"
    else:
        embedding_source = "scgpt_transformer"

    kmeans = KMeans(n_clusters=n_clusters, n_init=20, random_state=args.seed)
    pred = kmeans.fit_predict(embedding)
    np.save(save_dir / "embedding_final.npy", embedding.astype(np.float32))
    np.save(save_dir / "embeddings_base.npy", embedding.astype(np.float32))
    np.save(save_dir / "labels.npy", labels.astype(np.int64))
    np.save(save_dir / "eval_kmeans_known_k.npy", pred.astype(np.int64))
    family.save_embedding_h5(save_dir / "embedding.h5", embedding, labels)
    result = family.write_kmeans_known_k_outputs(
        output_dir=save_dir,
        dataset=dataset_name,
        method=args.method_name,
        seed=args.seed,
        embedding=embedding,
        labels=labels,
        n_clusters=n_clusters,
        extra={
            "variant": embedding_source,
            "used_scgpt_transformer": bool(scgpt_profile["used_scgpt_transformer"]),
            "n_matched_genes": int(scgpt_profile["n_matched_genes"]),
            "matched_gene_fraction": float(scgpt_profile["matched_gene_fraction"]),
            "fallback_reason": scgpt_profile["fallback_reason"],
        },
    )
    save_json(result["fixed"], str(save_dir / "metrics.json"))
    summary = {
        "dataset": dataset_name,
        "method": args.method_name,
        "variant": embedding_source,
        "seed": int(args.seed),
        "n_cells": int(data_np.shape[0]),
        "n_genes": int(data_np.shape[1]),
        "n_clusters": int(n_clusters),
        "scgpt_profile": scgpt_profile,
        "fixed_metrics": result["fixed"],
        "label_leakage": False,
    }
    save_json(summary, str(save_dir / "summary.json"))
    print(f"Results saved to: {save_dir}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
