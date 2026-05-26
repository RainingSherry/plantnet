#!/usr/bin/env python
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import StandardScaler, normalize

SCRIPT_DIR = Path(__file__).resolve().parent
PKG_DIR = SCRIPT_DIR.parent
ROOT = SCRIPT_DIR.parents[3]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from methods.DeepLearning.PlantSPADE_LGCL.data import load_lgcl_dataset, write_dataset_artifacts
from methods.DeepLearning.PlantSPADE_LGCL.eval import write_evaluation_outputs
from methods.DeepLearning.PlantSPADE_LGCL.utils import ensure_dir, save_json


PLANTSPADE_METHODS = {
    "plantspade_lgcl_baseline",
    "plantspade_lgcl_support_attention",
    "plantspade_lgcl_attention_no_idf",
    "plantspade_lgcl_attention_no_amplitude",
    "plantspade_lgcl_attention_topk_64",
    "plantspade_lgcl_attention_topk_128",
    "plantspade_lgcl_attention_topk_256",
    "plantspade_lgcl_neg_random_zero",
    "plantspade_lgcl_neg_idf_weighted_zero",
    "plantspade_lgcl_neg_neighbor_conflict_zero",
}

EXTERNAL_METHODS = {"phytocluster", "scvi", "scmae"}
TRADITIONAL_METHODS = {"sc3"}


def parse_args():
    parser = argparse.ArgumentParser(description="Run one method/dataset/seed under the PlantSPADE-LGCL protocol.")
    parser.add_argument("--dataset", required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--seed", type=int, required=True)
    parser.add_argument("--datasets_config", default=str(PKG_DIR / "configs" / "datasets_8plant.yaml"))
    parser.add_argument("--main_config", default=str(PKG_DIR / "configs" / "main_lgcl.yaml"))
    parser.add_argument("--baselines_config", default=str(PKG_DIR / "configs" / "baselines.yaml"))
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--save_h5ad", action="store_true")
    parser.add_argument("--dry_run", action="store_true")
    return parser.parse_args()


def load_yaml(path: str) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def find_dataset(cfg: dict, name: str) -> dict:
    for entry in cfg.get("datasets", []):
        if entry.get("dataset_name") == name:
            return entry
    raise KeyError(f"Dataset {name!r} not found in {cfg}")


def infer_n_clusters(entry: dict, bundle) -> int:
    expected = entry.get("expected_n_clusters")
    if expected:
        return int(expected)
    if bundle.labels is None:
        raise ValueError(f"Dataset {entry['dataset_name']} has no labels and expected_n_clusters is null.")
    return int(len(np.unique(bundle.labels)))


def guard_gpu(gpu: int | None):
    if gpu in {0, 7}:
        raise ValueError("GPU 0 and GPU 7 are forbidden. Use 1-6 or --no_cuda.")


def list_to_csv(values) -> str:
    if isinstance(values, str):
        return values
    return ",".join(str(v) for v in values)


def run_command(cmd: list[str], dry_run: bool, cwd: Path = ROOT):
    print(" ".join(cmd))
    if dry_run:
        return
    max_retries = 3
    for attempt in range(max_retries):
        try:
            subprocess.run(cmd, cwd=str(cwd), check=True)
            return
        except subprocess.CalledProcessError as e:
            if e.returncode == -11 and attempt < max_retries - 1:
                print(f"  [SIGSEGV] retry {attempt + 1}/{max_retries - 1}")
                continue
            raise


def canonical_bundle(entry: dict, main_cfg: dict, seed: int, output_dir: Path):
    prep = main_cfg.get("preprocessing", {})
    bundle = load_lgcl_dataset(
        entry["file_path"],
        input_mode=prep.get("input_mode", "auto"),
        n_top_genes=int(prep.get("n_top_genes", 2000)),
        target_sum=float(prep.get("target_sum", 10000.0)),
        svd_dim=int(prep.get("svd_dim", main_cfg.get("training", {}).get("latent_dim", 32))),
        svd_iter=int(prep.get("svd_iter", 7)),
        seed=seed,
        label_key=entry.get("label_key", "auto"),
    )
    write_dataset_artifacts(bundle, str(output_dir))
    return bundle


def write_canonical_h5ad(bundle, output_dir: Path) -> Path:
    path = output_dir / "canonical_hvg.h5ad"
    if not path.exists():
        adata = bundle.adata.copy()
        adata.X = bundle.amplitude.copy()
        adata.uns["canonical_preprocessing"] = bundle.preprocess_config
        for frame in (adata.obs, adata.var):
            for reserved in ("_index", "reserved_index"):
                if reserved in frame.columns:
                    replacement = reserved + "_renamed"
                    suffix = 1
                    while replacement in frame.columns:
                        replacement = f"{reserved}_renamed_{suffix}"
                        suffix += 1
                    frame.rename(columns={reserved: replacement}, inplace=True)
            if frame.index.name == "_index":
                frame.index.name = "cell_name" if frame is adata.obs else "gene_name"
        adata.write_h5ad(path, compression="gzip")
    return path


def pca_embedding(matrix, n_components: int, seed: int) -> np.ndarray:
    max_components = max(2, min(n_components, matrix.shape[0] - 1, matrix.shape[1] - 1))
    if sparse.issparse(matrix):
        emb = TruncatedSVD(n_components=max_components, random_state=seed).fit_transform(matrix)
    else:
        emb = PCA(n_components=max_components, random_state=seed).fit_transform(np.asarray(matrix))
    emb = np.nan_to_num(emb.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)
    if emb.shape[1] > 1:
        emb = StandardScaler().fit_transform(emb).astype(np.float32)
    return normalize(emb, norm="l2", axis=1, copy=False).astype(np.float32)


def run_traditional_pca(entry: dict, main_cfg: dict, seed: int, output_dir: Path):
    bundle = canonical_bundle(entry, main_cfg, seed, output_dir)
    n_clusters = infer_n_clusters(entry, bundle)
    baseline_cfg = load_yaml(str(PKG_DIR / "configs" / "baselines.yaml"))
    n_components = int(baseline_cfg.get("traditional", {}).get("pca_components", 50))
    embedding = pca_embedding(bundle.amplitude, n_components=n_components, seed=seed)
    np.save(output_dir / "embedding_pca.npy", embedding)
    eval_cfg = main_cfg.get("evaluation", {})
    result = write_evaluation_outputs(
        output_dir=str(output_dir),
        dataset=entry["dataset_name"],
        method="traditional_pca",
        seed=seed,
        embedding=embedding,
        labels=bundle.labels,
        n_clusters=n_clusters,
        n_neighbors=int(eval_cfg.get("n_neighbors", 15)),
        leiden_fixed_resolution=float(eval_cfg.get("leiden_fixed_resolution", 1.0)),
        louvain_fixed_resolution=float(eval_cfg.get("louvain_fixed_resolution", 1.0)),
        leiden_sweep_resolutions=eval_cfg.get("leiden_sweep_resolutions", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]),
        sweep_max_cells=int(eval_cfg.get("sweep_max_cells", 10000)),
        include_louvain=True,
        prefix="pca",
    )
    save_json({"method": "traditional_pca", "fixed": result["fixed"], "oracle": result["oracle"]}, str(output_dir / "summary.json"))


def plantspade_method_args(method: str) -> tuple[str, list[str]]:
    if method == "plantspade_lgcl_baseline":
        return "plantspade_lgcl", ["--use_support_attention", "false"]
    if method == "plantspade_lgcl_support_attention":
        return "plantspade_lgcl_sga", ["--use_support_attention", "true"]
    if method == "plantspade_lgcl_attention_no_idf":
        return "plantspade_lgcl_attention_no_idf", ["--use_support_attention", "true", "--attention_gamma", "0.0"]
    if method == "plantspade_lgcl_attention_no_amplitude":
        return "plantspade_lgcl_attention_no_amplitude", ["--use_support_attention", "true", "--attention_beta", "0.0"]
    if method.startswith("plantspade_lgcl_attention_topk_"):
        topk = method.rsplit("_", 1)[-1]
        return method, ["--use_support_attention", "true", "--attention_topk_genes", topk]
    if method == "plantspade_lgcl_neg_random_zero":
        return method, ["--negative_sampler", "random_zero"]
    if method == "plantspade_lgcl_neg_idf_weighted_zero":
        return method, ["--negative_sampler", "idf_weighted_zero"]
    if method == "plantspade_lgcl_neg_neighbor_conflict_zero":
        return method, ["--negative_sampler", "neighbor_conflict_zero"]
    raise ValueError(f"Unsupported PlantSPADE-LGCL method: {method}")


def run_plantspade(entry: dict, main_cfg: dict, method: str, seed: int, output_dir: Path, gpu: int | None, no_cuda: bool, save_h5ad: bool, dry_run: bool):
    prep = main_cfg.get("preprocessing", {})
    train = main_cfg.get("training", {})
    eval_cfg = main_cfg.get("evaluation", {})
    attention = main_cfg.get("attention", {})
    method_name, extra_args = plantspade_method_args(method)
    cmd = [
        sys.executable,
        str(PKG_DIR / "run_plantspade.py"),
        "--data_path",
        entry["file_path"],
        "--save_dir",
        str(output_dir),
        "--dataset_name",
        entry["dataset_name"],
        "--label_key",
        str(entry.get("label_key", "auto")),
        "--method_name",
        method_name,
        "--seed",
        str(seed),
        "--input_mode",
        str(prep.get("input_mode", "auto")),
        "--n_top_genes",
        str(prep.get("n_top_genes", 2000)),
        "--target_sum",
        str(prep.get("target_sum", 10000.0)),
        "--svd_dim",
        str(prep.get("svd_dim", 32)),
        "--svd_iter",
        str(prep.get("svd_iter", 7)),
        "--latent_dim",
        str(train.get("latent_dim", 32)),
        "--layers",
        str(train.get("layers", 2)),
        "--epochs",
        str(train.get("epochs", 80)),
        "--pairs_per_epoch",
        str(train.get("pairs_per_epoch", 262144)),
        "--contrastive_batch_size",
        str(train.get("contrastive_batch_size", 2048)),
        "--lr",
        str(train.get("lr", 0.001)),
        "--weight_decay",
        str(train.get("weight_decay", 0.00001)),
        "--edge_dropout",
        str(train.get("edge_dropout", 0.1)),
        "--temperature",
        str(train.get("temperature", 0.2)),
        "--contrastive_weight",
        str(train.get("contrastive_weight", 0.05)),
        "--module_weight",
        str(train.get("module_weight", 0.001)),
        "--num_modules",
        str(train.get("num_modules", 16)),
        "--module_top_k",
        str(train.get("module_top_k", 30)),
        "--negative_sampler",
        str(train.get("negative_sampler", "random_zero")),
        "--eval_neighbors",
        str(eval_cfg.get("n_neighbors", 15)),
        "--leiden_fixed_resolution",
        str(eval_cfg.get("leiden_fixed_resolution", 1.0)),
        "--louvain_fixed_resolution",
        str(eval_cfg.get("louvain_fixed_resolution", 1.0)),
        "--leiden_resolutions",
        list_to_csv(eval_cfg.get("leiden_sweep_resolutions", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2])),
        "--sweep_max_cells",
        str(eval_cfg.get("sweep_max_cells", 10000)),
        "--attention_topk_genes",
        str(attention.get("attention_topk_genes", 128)),
        "--attention_beta",
        str(attention.get("attention_beta", 0.1)),
        "--attention_gamma",
        str(attention.get("attention_gamma", 0.1)),
        "--attention_eta",
        str(attention.get("attention_eta", 0.5)),
        "--attention_dropout",
        str(attention.get("attention_dropout", 0.1)),
    ]
    if entry.get("expected_n_clusters"):
        cmd.extend(["--n_clusters", str(entry["expected_n_clusters"])])
    if gpu is not None:
        cmd.extend(["--gpu", str(gpu)])
    if no_cuda:
        cmd.append("--no_cuda")
    if not save_h5ad:
        cmd.append("--no_save_h5ad")
    cmd.extend(extra_args)
    run_command(cmd, dry_run=dry_run)


def latest_embedding(output_dir: Path) -> Path:
    candidates = [Path(p) for p in glob.glob(str(output_dir / "**" / "embedding_*.npy"), recursive=True)]
    candidates.extend(Path(p) for p in glob.glob(str(output_dir / "**" / "embedding_final.npy"), recursive=True))
    if not candidates:
        raise FileNotFoundError(f"No embedding_*.npy found under {output_dir}")
    return max(candidates, key=lambda p: p.stat().st_mtime)


def run_external(entry: dict, main_cfg: dict, baselines_cfg: dict, method: str, seed: int, output_dir: Path, gpu: int | None, no_cuda: bool, dry_run: bool):
    bundle = canonical_bundle(entry, main_cfg, seed, output_dir)
    n_clusters = infer_n_clusters(entry, bundle)
    canonical_path = write_canonical_h5ad(bundle, output_dir)
    method_cfg = baselines_cfg.get("deep", {}).get(method, {})
    runner = ROOT / method_cfg["runner"]
    data_path = entry["file_path"] if method == "scvi" else str(canonical_path)
    cmd = [
        sys.executable,
        str(runner),
        "--data_path",
        str(data_path),
        "--save_dir",
        str(output_dir),
        "--n_clusters",
        str(n_clusters),
        "--seed",
        str(seed),
    ]
    for key, value in method_cfg.get("args", {}).items():
        cmd.extend([f"--{key}", str(value)])
    if gpu is not None:
        cmd.extend(["--gpu", str(gpu)])
    if no_cuda and method in {"phytocluster", "scmae"}:
        cmd.append("--no_cuda")
    run_command(cmd, dry_run=dry_run)
    if dry_run:
        return

    embedding_path = latest_embedding(output_dir)
    embedding = np.load(embedding_path)
    eval_cfg = main_cfg.get("evaluation", {})
    result = write_evaluation_outputs(
        output_dir=str(output_dir),
        dataset=entry["dataset_name"],
        method=method,
        seed=seed,
        embedding=embedding,
        labels=bundle.labels,
        n_clusters=n_clusters,
        n_neighbors=int(eval_cfg.get("n_neighbors", 15)),
        leiden_fixed_resolution=float(eval_cfg.get("leiden_fixed_resolution", 1.0)),
        louvain_fixed_resolution=float(eval_cfg.get("louvain_fixed_resolution", 1.0)),
        leiden_sweep_resolutions=eval_cfg.get("leiden_sweep_resolutions", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]),
        sweep_max_cells=int(eval_cfg.get("sweep_max_cells", 10000)),
        include_louvain=True,
        prefix="external_eval",
    )
    save_json(
        {"method": method, "embedding_path": str(embedding_path), "fixed": result["fixed"], "oracle": result["oracle"]},
        str(output_dir / "summary_unified_eval.json"),
    )


def run_traditional_external(entry: dict, main_cfg: dict, baselines_cfg: dict, method: str, seed: int, output_dir: Path, dry_run: bool):
    bundle = canonical_bundle(entry, main_cfg, seed, output_dir)
    n_clusters = infer_n_clusters(entry, bundle)
    method_cfg = baselines_cfg.get("sc3", {})
    runner = ROOT / method_cfg.get("runner", "methods/Traditional/sc3/run.py")
    cmd = [
        sys.executable,
        str(runner),
        "--data_path",
        entry["file_path"],
        "--save_dir",
        str(output_dir),
        "--n_clusters",
        str(n_clusters),
        "--seed",
        str(seed),
    ]
    run_command(cmd, dry_run=dry_run)
    if dry_run:
        return
    embedding_path = latest_embedding(output_dir)
    embedding = np.load(embedding_path)
    eval_cfg = main_cfg.get("evaluation", {})
    result = write_evaluation_outputs(
        output_dir=str(output_dir),
        dataset=entry["dataset_name"],
        method=method,
        seed=seed,
        embedding=embedding,
        labels=bundle.labels,
        n_clusters=n_clusters,
        n_neighbors=int(eval_cfg.get("n_neighbors", 15)),
        leiden_fixed_resolution=float(eval_cfg.get("leiden_fixed_resolution", 1.0)),
        louvain_fixed_resolution=float(eval_cfg.get("louvain_fixed_resolution", 1.0)),
        leiden_sweep_resolutions=eval_cfg.get("leiden_sweep_resolutions", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2]),
        sweep_max_cells=int(eval_cfg.get("sweep_max_cells", 10000)),
        include_louvain=True,
        prefix="sc3_eval",
    )
    save_json(
        {"method": method, "embedding_path": str(embedding_path), "fixed": result["fixed"], "oracle": result["oracle"]},
        str(output_dir / "summary_sc3.json"),
    )


def main():
    args = parse_args()
    datasets_cfg = load_yaml(args.datasets_config)
    main_cfg = load_yaml(args.main_config)
    baselines_cfg = load_yaml(args.baselines_config)
    entry = find_dataset(datasets_cfg, args.dataset)
    gpu = args.gpu if args.gpu is not None else main_cfg.get("gpu", 1)
    if args.no_cuda:
        gpu = None
    guard_gpu(gpu)

    base_output = Path(args.output_dir or main_cfg.get("output_dir", ROOT / "results" / "PlantSPADE_LGCL_protocol"))
    output_dir = Path(ensure_dir(base_output / entry["dataset_name"] / args.method / f"seed_{args.seed}"))
    save_json({"dataset": entry, "method": args.method, "seed": args.seed}, str(output_dir / "run_single_config.json"))

    if not Path(entry["file_path"]).exists() and not args.dry_run:
        raise FileNotFoundError(f"Dataset file is missing: {entry['file_path']}")

    if args.method == "traditional_pca":
        run_traditional_pca(entry, main_cfg, args.seed, output_dir)
    elif args.method in PLANTSPADE_METHODS:
        run_plantspade(entry, main_cfg, args.method, args.seed, output_dir, gpu, args.no_cuda, args.save_h5ad, args.dry_run)
    elif args.method in EXTERNAL_METHODS:
        run_external(entry, main_cfg, baselines_cfg, args.method, args.seed, output_dir, gpu, args.no_cuda, args.dry_run)
    elif args.method in TRADITIONAL_METHODS:
        run_traditional_external(entry, main_cfg, baselines_cfg, args.method, args.seed, output_dir, args.dry_run)
    else:
        raise ValueError(f"Unsupported method {args.method!r}")


if __name__ == "__main__":
    main()
