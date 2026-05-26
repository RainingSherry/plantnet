#!/usr/bin/env python
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")

import numpy as np
import scipy.sparse as sp
import torch
import torch.nn.functional as F

SCRIPT_DIR = Path(__file__).resolve().parent
PKG_DIR = SCRIPT_DIR.parent
ROOT = SCRIPT_DIR.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning.PlantSPADE_LGCL.eval import write_evaluation_outputs
from methods.DeepLearning.PlantSPADE_LGCL.support_gene_attention import SupportGeneAttention
from methods.DeepLearning.PlantSPADE_LGCL.utils import save_json


def str2bool(value):
    if isinstance(value, bool):
        return value
    value = str(value).strip().lower()
    if value in {"1", "true", "t", "yes", "y"}:
        return True
    if value in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected true/false, got {value!r}")


def parse_float_list(value: str):
    if value is None or str(value).strip() == "":
        return []
    return [float(item) for item in str(value).split(",") if item.strip()]


def parse_args():
    parser = argparse.ArgumentParser(description="CPU-only evaluation recovery from saved embeddings.")
    parser.add_argument("--run_dir", required=True)
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--method_name", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--n_clusters", type=int, default=0)
    parser.add_argument("--embedding_path", default=None)
    parser.add_argument("--labels_path", default=None)
    parser.add_argument("--variant_name", default="baseline")
    parser.add_argument("--prefix", default=None)
    parser.add_argument("--use_support_attention", type=str2bool, default=False)
    parser.add_argument("--support_path", default=None)
    parser.add_argument("--amplitude_path", default=None)
    parser.add_argument("--gene_embedding_path", default=None)
    parser.add_argument("--attention_topk_genes", type=int, default=128)
    parser.add_argument("--attention_beta", type=float, default=0.1)
    parser.add_argument("--attention_gamma", type=float, default=0.1)
    parser.add_argument("--attention_eta", type=float, default=0.5)
    parser.add_argument("--attention_dropout", type=float, default=0.0)
    parser.add_argument("--normalize_embedding", type=str2bool, default=True)
    parser.add_argument("--eval_neighbors", type=int, default=15)
    parser.add_argument("--leiden_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--louvain_fixed_resolution", type=float, default=1.0)
    parser.add_argument("--leiden_resolutions", default="0.2,0.4,0.6,0.8,1.0,1.2")
    parser.add_argument("--sweep_max_cells", type=int, default=10000)
    parser.add_argument("--include_louvain", type=str2bool, default=False)
    parser.add_argument("--run_oracle_sweep", type=str2bool, default=False)
    parser.add_argument("--silhouette_sample_size", type=int, default=3000)
    return parser.parse_args()


def load_array(path: Path, name: str) -> np.ndarray:
    if not path.exists():
        raise FileNotFoundError(f"Missing {name}: {path}")
    return np.load(path)


def compute_gene_idf(support: sp.csr_matrix) -> np.ndarray:
    df = np.diff(support.tocsc().indptr).astype(np.float32)
    n_cells = float(support.shape[0])
    return np.log1p(n_cells / (1.0 + df)).astype(np.float32)


@torch.no_grad()
def support_attention_embedding(args, run_dir: Path, base_embedding: np.ndarray) -> np.ndarray:
    support_path = Path(args.support_path) if args.support_path else run_dir / "support_matrix.npz"
    amplitude_path = Path(args.amplitude_path) if args.amplitude_path else run_dir / "amplitude_matrix.npz"
    gene_embedding_path = Path(args.gene_embedding_path) if args.gene_embedding_path else run_dir / "gene_embedding.npy"
    support = sp.load_npz(support_path).tocsr()
    amplitude = sp.load_npz(amplitude_path).tocsr()
    gene_embedding = load_array(gene_embedding_path, "gene embedding").astype(np.float32)
    gene_idf = compute_gene_idf(support)

    device = torch.device("cpu")
    attention = SupportGeneAttention(
        support=support,
        amplitude=amplitude,
        gene_idf=torch.as_tensor(gene_idf, dtype=torch.float32, device=device),
        top_k_genes=args.attention_topk_genes,
        beta=args.attention_beta,
        gamma=args.attention_gamma,
        eta=args.attention_eta,
        dropout=args.attention_dropout,
    ).to(device)
    attention.eval()
    cell_t = torch.as_tensor(base_embedding, dtype=torch.float32, device=device)
    gene_t = torch.as_tensor(gene_embedding, dtype=torch.float32, device=device)
    refined = attention(cell_t, gene_t, return_attention=False)
    if args.normalize_embedding:
        refined = F.normalize(refined, dim=1)
    return refined.detach().cpu().numpy().astype(np.float32)


def main():
    args = parse_args()
    run_dir = Path(args.run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    embedding_path = Path(args.embedding_path) if args.embedding_path else run_dir / "embedding_final.npy"
    labels_path = Path(args.labels_path) if args.labels_path else run_dir / "labels.npy"
    base_embedding = load_array(embedding_path, "embedding").astype(np.float32)
    labels = load_array(labels_path, "labels").astype(np.int64)
    n_clusters = int(args.n_clusters) if args.n_clusters and args.n_clusters > 0 else int(len(np.unique(labels)))
    resolutions = parse_float_list(args.leiden_resolutions)
    variant_name = str(args.variant_name)
    prefix = args.prefix or f"eval_{variant_name}"

    if args.use_support_attention:
        embedding = support_attention_embedding(args, run_dir, base_embedding)
        np.save(run_dir / f"embedding_{variant_name}.npy", embedding.astype(np.float32))
        np.save(run_dir / "embedding_final.npy", embedding.astype(np.float32))
    else:
        embedding = base_embedding
        np.save(run_dir / f"embedding_{variant_name}.npy", embedding.astype(np.float32))

    try:
        result = write_evaluation_outputs(
            output_dir=str(run_dir),
            dataset=args.dataset,
            method=args.method_name,
            seed=args.seed,
            embedding=embedding,
            labels=labels,
            n_clusters=n_clusters,
            n_neighbors=args.eval_neighbors,
            leiden_fixed_resolution=args.leiden_fixed_resolution,
            louvain_fixed_resolution=args.louvain_fixed_resolution,
            leiden_sweep_resolutions=resolutions,
            include_louvain=args.include_louvain,
            run_oracle_sweep=args.run_oracle_sweep,
            sweep_max_cells=args.sweep_max_cells,
            silhouette_sample_size=args.silhouette_sample_size,
            prefix=prefix,
            extra={"variant": variant_name},
        )
    except Exception as exc:
        failure = {
            "dataset": args.dataset,
            "method": args.method_name,
            "seed": int(args.seed),
            "variant": variant_name,
            "error": repr(exc),
        }
        save_json(failure, str(run_dir / f"{prefix}_failure.json"))
        raise

    payload = {
        "dataset": args.dataset,
        "method": args.method_name,
        "seed": int(args.seed),
        "variant": variant_name,
        "embedding_path": str(embedding_path),
        "fixed": result["fixed"],
        "oracle": result["oracle"],
        "run_oracle_sweep": bool(args.run_oracle_sweep),
        "include_louvain": bool(args.include_louvain),
        "note": "Main results use fixed protocol; oracle sweep is supplementary only when explicitly enabled.",
    }
    save_json(payload, str(run_dir / f"{args.method_name}.json"))
    save_json(result["fixed"], str(run_dir / "metrics.json"))
    print(json.dumps({"status": "ok", "run_dir": str(run_dir), "prefix": prefix}, indent=2))


if __name__ == "__main__":
    main()
