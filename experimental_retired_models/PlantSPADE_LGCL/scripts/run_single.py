#!/usr/bin/env python
from __future__ import annotations

import argparse
import glob
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import yaml
from scipy import sparse
from sklearn.decomposition import PCA, TruncatedSVD
from sklearn.preprocessing import StandardScaler, normalize

SCRIPT_DIR = Path(__file__).resolve().parent
PKG_DIR = SCRIPT_DIR.parent
ROOT = next(parent for parent in [SCRIPT_DIR, *SCRIPT_DIR.parents] if (parent / "methods" / "DeepLearning" / "scMAE_family.py").exists())
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experimental_retired_models.PlantSPADE_LGCL.data import load_lgcl_dataset, write_dataset_artifacts
from experimental_retired_models.PlantSPADE_LGCL.utils import ensure_dir, save_json


PLANTSPADE_METHODS = {
    "plantspade_lgcl_baseline",
    "plantspade_lgcl_support_attention",
    "plantspade_lgcl_gated_fusion",
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
TRADITIONAL_EMBEDDING_METHODS = {"traditional_pca", "traditional_leiden", "traditional_louvain"}

THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMBA_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
}


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
    visible = os.environ.get("CUDA_VISIBLE_DEVICES")
    if visible:
        visible_ids = [item.strip() for item in visible.split(",") if item.strip()]
        if set(visible_ids).intersection({"0", "7"}):
            raise ValueError("CUDA_VISIBLE_DEVICES includes forbidden physical GPU 0 or 7.")
        if gpu is not None and (gpu < 0 or gpu >= len(visible_ids)):
            raise ValueError(f"--gpu {gpu} is outside isolated CUDA_VISIBLE_DEVICES={visible!r}.")
        return
    if gpu in {0, 7}:
        raise ValueError("GPU 0 and GPU 7 are forbidden. Use 1-6 or --no_cuda.")


def list_to_csv(values) -> str:
    if isinstance(values, str):
        return values
    return ",".join(str(v) for v in values)


def format_command(cmd: list[str], env: dict | None = None) -> str:
    if not env:
        return " ".join(cmd)
    shown = []
    for key in ["CUDA_VISIBLE_DEVICES", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMBA_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        if key in env:
            shown.append(f"{key}={env[key]}")
    return " ".join(shown + cmd)


def run_command(cmd: list[str], dry_run: bool, cwd: Path = ROOT, env: dict | None = None):
    display_env = dict(THREAD_ENV)
    if env:
        display_env.update(env)
    print(format_command(cmd, display_env))
    if dry_run:
        return
    merged_env = os.environ.copy()
    merged_env.update(THREAD_ENV)
    if env:
        merged_env.update(env)
    subprocess.run(cmd, cwd=str(cwd), check=True, env=merged_env)


def is_cell_embedding_path(path: Path, labels: np.ndarray | None = None) -> bool:
    name = path.name
    if name in {"gene_embedding.npy", "global_embedding_svd_projected.npy"}:
        return False
    if name.endswith("_mapped.npy"):
        return False
    if name.startswith(("eval_", "external_eval_", "pca_", "sc3_eval_")):
        return False
    if not (name.startswith("embedding_") or name == "embeddings_base.npy"):
        return False
    if labels is None:
        return True
    try:
        arr = np.load(path, mmap_mode="r")
        return bool(arr.ndim == 2 and arr.shape[0] == labels.shape[0])
    except Exception:
        return False


def find_existing_cell_embedding(output_dir: Path, labels: np.ndarray | None = None) -> Path | None:
    preferred = [
        "embedding_baseline.npy",
        "embeddings_base.npy",
        "embedding_primary.npy",
        "embedding_final.npy",
    ]
    for name in preferred:
        path = output_dir / name
        if path.exists() and is_cell_embedding_path(path, labels=labels):
            return path
    candidates = [Path(p) for p in glob.glob(str(output_dir / "**" / "embedding_*.npy"), recursive=True)]
    candidates.extend(Path(p) for p in glob.glob(str(output_dir / "**" / "embeddings_base.npy"), recursive=True))
    candidates = [path for path in candidates if is_cell_embedding_path(path, labels=labels)]
    if not candidates:
        return None
    return max(candidates, key=lambda p: p.stat().st_mtime)


def has_complete_plantspade_artifacts(output_dir: Path) -> bool:
    required = [
        "embedding_baseline.npy",
        "training_history.json",
        "labels.npy",
        "support_matrix.npz",
        "amplitude_matrix.npz",
        "gene_embedding.npy",
        "global_embedding_svd_projected.npy",
    ]
    return all((output_dir / name).exists() for name in required)


def canonical_bundle(entry: dict, main_cfg: dict, seed: int, output_dir: Path):
    prep = main_cfg.get("preprocessing", {})
    subsample = main_cfg.get("subsample", {})
    bundle = load_lgcl_dataset(
        entry["file_path"],
        input_mode=prep.get("input_mode", "auto"),
        n_top_genes=int(prep.get("n_top_genes", 2000)),
        target_sum=float(prep.get("target_sum", 10000.0)),
        svd_dim=int(prep.get("svd_dim", main_cfg.get("training", {}).get("latent_dim", 32))),
        svd_iter=int(prep.get("svd_iter", 7)),
        seed=seed,
        label_key=entry.get("label_key", "auto"),
        subsample_per_class_max=int(subsample.get("per_class_max_cells", 0) or 0),
        subsample_fallback_max=int(subsample.get("fallback_max_cells", 0) or 0),
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


def run_traditional_embedding(entry: dict, main_cfg: dict, method: str, seed: int, output_dir: Path, dry_run: bool = False):
    bundle = canonical_bundle(entry, main_cfg, seed, output_dir)
    n_clusters = infer_n_clusters(entry, bundle)
    baseline_cfg = load_yaml(str(PKG_DIR / "configs" / "baselines.yaml"))
    n_components = int(baseline_cfg.get("traditional", {}).get("pca_components", 50))
    embedding_path = output_dir / "embedding_final.npy"
    if not embedding_path.exists():
        embedding = pca_embedding(bundle.amplitude, n_components=n_components, seed=seed)
        np.save(output_dir / "embedding_pca.npy", embedding)
        np.save(embedding_path, embedding)
    save_labels(output_dir, bundle.labels)
    if method == "traditional_pca":
        prefix = "pca"
        variant_name = "pca"
        eval_overrides = None
    elif method == "traditional_leiden":
        prefix = "eval_traditional_leiden"
        variant_name = "leiden"
        eval_overrides = None
    elif method == "traditional_louvain":
        prefix = "eval_traditional_louvain"
        variant_name = "louvain"
        eval_overrides = {"include_louvain": "true"}
    else:
        raise ValueError(f"Unsupported traditional embedding method: {method}")
    run_eval_from_embedding(
        entry=entry,
        main_cfg=main_cfg,
        method=method,
        seed=seed,
        output_dir=output_dir,
        n_clusters=n_clusters,
        embedding_path=embedding_path,
        prefix=prefix,
        variant_name=variant_name,
        dry_run=dry_run,
        attention_overrides=eval_overrides,
    )


def plantspade_method_args(method: str) -> tuple[str, list[str]]:
    if method == "plantspade_lgcl_baseline":
        return "plantspade_lgcl", ["--use_support_attention", "false"]
    if method == "plantspade_lgcl_support_attention":
        return "plantspade_lgcl_sga", ["--use_support_attention", "true"]
    if method == "plantspade_lgcl_gated_fusion":
        return "plantspade_lgcl_gated_fusion", ["--use_gated_fusion", "true", "--use_support_attention", "true"]
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
    baseline_embedding_path = output_dir / "embedding_baseline.npy"
    history_path = output_dir / "training_history.json"
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
        "--subsample_per_class_max",
        str(main_cfg.get("subsample", {}).get("per_class_max_cells", 0) or 0),
        "--subsample_fallback_max",
        str(main_cfg.get("subsample", {}).get("fallback_max_cells", 0) or 0),
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
        "--include_louvain",
        str(bool(eval_cfg.get("include_louvain", False))).lower(),
        "--run_oracle_sweep",
        str(bool(eval_cfg.get("run_oracle_sweep", False))).lower(),
        "--silhouette_sample_size",
        str(eval_cfg.get("silhouette_sample_size", 3000)),
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
        "--train_only",
        "true",
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
    if has_complete_plantspade_artifacts(output_dir):
        print(f"[skip train] existing PlantSPADE artifacts found in {output_dir}")
    else:
        try:
            run_command(cmd, dry_run=dry_run)
        except subprocess.CalledProcessError:
            if has_complete_plantspade_artifacts(output_dir):
                print(f"[recover] training command failed after embeddings were saved; continuing to CPU eval: {output_dir}")
            else:
                raise
    if dry_run:
        return
    if not has_complete_plantspade_artifacts(output_dir):
        raise FileNotFoundError(f"Missing PlantSPADE training artifacts in {output_dir}")

    labels = np.load(output_dir / "labels.npy")
    n_clusters = int(entry.get("expected_n_clusters") or len(np.unique(labels)))
    prefix, variant_name, use_attention, attention_overrides = plantspade_eval_plan(method, main_cfg)
    eval_embedding = output_dir / "embedding_gated_fusion.npy" if method == "plantspade_lgcl_gated_fusion" else baseline_embedding_path
    run_eval_from_embedding(
        entry=entry,
        main_cfg=main_cfg,
        method=method,
        seed=seed,
        output_dir=output_dir,
        n_clusters=n_clusters,
        embedding_path=eval_embedding,
        prefix=prefix,
        variant_name=variant_name,
        dry_run=dry_run,
        use_support_attention=use_attention,
        attention_overrides=attention_overrides,
    )


def latest_embedding(output_dir: Path) -> Path:
    candidate = find_existing_cell_embedding(output_dir)
    if candidate is None:
        raise FileNotFoundError(f"No embedding_*.npy found under {output_dir}")
    return candidate


def save_labels(output_dir: Path, labels: np.ndarray | None) -> Path:
    if labels is None:
        raise ValueError("Labels are required for fixed benchmark evaluation.")
    labels_path = output_dir / "labels.npy"
    np.save(labels_path, labels.astype(np.int64))
    return labels_path


def evaluation_done(output_dir: Path, prefix: str) -> bool:
    path = output_dir / f"{prefix}_fixed.csv"
    return path.exists() and path.stat().st_size > 0


def eval_config_args(main_cfg: dict) -> list[str]:
    eval_cfg = main_cfg.get("evaluation", {})
    args = [
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
        "--include_louvain",
        str(bool(eval_cfg.get("include_louvain", False))).lower(),
        "--run_oracle_sweep",
        str(bool(eval_cfg.get("run_oracle_sweep", False))).lower(),
        "--silhouette_sample_size",
        str(eval_cfg.get("silhouette_sample_size", 3000)),
    ]
    return args


def run_eval_from_embedding(
    entry: dict,
    main_cfg: dict,
    method: str,
    seed: int,
    output_dir: Path,
    n_clusters: int,
    embedding_path: Path,
    prefix: str,
    variant_name: str,
    dry_run: bool,
    use_support_attention: bool = False,
    attention_overrides: dict | None = None,
) -> None:
    if evaluation_done(output_dir, prefix):
        print(f"[skip eval] {prefix}_fixed.csv already exists")
        return
    labels_path = output_dir / "labels.npy"
    if not labels_path.exists():
        raise FileNotFoundError(f"Missing labels for evaluation recovery: {labels_path}")
    cmd = [
        sys.executable,
        str(SCRIPT_DIR / "eval_from_embedding.py"),
        "--run_dir",
        str(output_dir),
        "--dataset",
        entry["dataset_name"],
        "--method_name",
        method,
        "--seed",
        str(seed),
        "--n_clusters",
        str(n_clusters),
        "--embedding_path",
        str(embedding_path),
        "--labels_path",
        str(labels_path),
        "--prefix",
        prefix,
        "--variant_name",
        variant_name,
        "--use_support_attention",
        str(bool(use_support_attention)).lower(),
    ]
    cmd.extend(eval_config_args(main_cfg))
    if attention_overrides:
        for key, value in attention_overrides.items():
            cmd.extend([f"--{key}", str(value)])
    try:
        run_command(cmd, dry_run=dry_run, env={"CUDA_VISIBLE_DEVICES": ""})
    except subprocess.CalledProcessError as exc:
        save_json(
            {
                "dataset": entry["dataset_name"],
                "method": method,
                "seed": int(seed),
                "prefix": prefix,
                "embedding_path": str(embedding_path),
                "returncode": int(exc.returncode),
                "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "note": "Evaluation failed in CPU recovery. Training outputs are left intact and will not be retrained automatically.",
            },
            str(output_dir / f"{prefix}_recovery_failure.json"),
        )
        raise


def plantspade_eval_plan(method: str, main_cfg: dict) -> tuple[str, str, bool, dict]:
    attention = main_cfg.get("attention", {})
    common = {
        "attention_topk_genes": attention.get("attention_topk_genes", 128),
        "attention_beta": attention.get("attention_beta", 0.1),
        "attention_gamma": attention.get("attention_gamma", 0.1),
        "attention_eta": attention.get("attention_eta", 0.5),
        "attention_dropout": 0.0,
    }
    if method == "plantspade_lgcl_baseline":
        return "eval_baseline", "baseline", False, {}
    if method == "plantspade_lgcl_support_attention":
        return "eval_support_attention", "support_attention", True, common
    if method == "plantspade_lgcl_gated_fusion":
        return "eval_gated_fusion", "gated_fusion", False, {}
    if method == "plantspade_lgcl_attention_no_idf":
        cfg = dict(common)
        cfg["attention_gamma"] = 0.0
        return "eval_attention_no_idf", "attention_no_idf", True, cfg
    if method == "plantspade_lgcl_attention_no_amplitude":
        cfg = dict(common)
        cfg["attention_beta"] = 0.0
        return "eval_attention_no_amplitude", "attention_no_amplitude", True, cfg
    if method.startswith("plantspade_lgcl_attention_topk_"):
        cfg = dict(common)
        cfg["attention_topk_genes"] = int(method.rsplit("_", 1)[-1])
        return f"eval_attention_topk_{cfg['attention_topk_genes']}", f"attention_topk_{cfg['attention_topk_genes']}", True, cfg
    return "eval_baseline", "baseline", False, {}


def run_external(entry: dict, main_cfg: dict, baselines_cfg: dict, method: str, seed: int, output_dir: Path, gpu: int | None, no_cuda: bool, dry_run: bool):
    bundle = canonical_bundle(entry, main_cfg, seed, output_dir)
    n_clusters = infer_n_clusters(entry, bundle)
    save_labels(output_dir, bundle.labels)
    canonical_path = write_canonical_h5ad(bundle, output_dir)
    method_cfg = baselines_cfg.get("deep", {}).get(method, {})
    runner = ROOT / method_cfg["runner"]
    data_path = entry["file_path"] if method == "scvi" else str(canonical_path)
    final_embedding = output_dir / "embedding_final.npy"
    existing_embedding = find_existing_cell_embedding(output_dir, labels=bundle.labels)
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
    if existing_embedding is not None:
        print(f"[skip train] existing cell embedding found for {method}: {existing_embedding}")
    else:
        try:
            run_command(cmd, dry_run=dry_run)
        except subprocess.CalledProcessError:
            existing_embedding = find_existing_cell_embedding(output_dir, labels=bundle.labels)
            if existing_embedding is not None:
                print(f"[recover] external runner failed after embedding was saved; continuing to CPU eval: {existing_embedding}")
            else:
                raise
    if dry_run:
        return

    embedding_path = find_existing_cell_embedding(output_dir, labels=bundle.labels)
    if embedding_path is None:
        raise FileNotFoundError(f"No valid cell embedding found under {output_dir}")
    if embedding_path != final_embedding:
        np.save(final_embedding, np.load(embedding_path).astype(np.float32))
        embedding_path = final_embedding
    run_eval_from_embedding(
        entry=entry,
        main_cfg=main_cfg,
        method=method,
        seed=seed,
        output_dir=output_dir,
        n_clusters=n_clusters,
        embedding_path=embedding_path,
        prefix="external_eval",
        variant_name="embedding",
        dry_run=dry_run,
    )


def run_traditional_external(entry: dict, main_cfg: dict, baselines_cfg: dict, method: str, seed: int, output_dir: Path, dry_run: bool):
    bundle = canonical_bundle(entry, main_cfg, seed, output_dir)
    n_clusters = infer_n_clusters(entry, bundle)
    save_labels(output_dir, bundle.labels)
    method_cfg = baselines_cfg.get("sc3", {})
    runner = ROOT / method_cfg.get("runner", "methods/Traditional/sc3/run.py")
    final_embedding = output_dir / "embedding_final.npy"
    existing_embedding = find_existing_cell_embedding(output_dir, labels=bundle.labels)
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
    if existing_embedding is not None:
        print(f"[skip train] existing cell embedding found for {method}: {existing_embedding}")
    else:
        run_command(cmd, dry_run=dry_run)
    if dry_run:
        return
    embedding_path = find_existing_cell_embedding(output_dir, labels=bundle.labels)
    if embedding_path is None:
        raise FileNotFoundError(f"No valid cell embedding found under {output_dir}")
    if embedding_path != final_embedding:
        np.save(final_embedding, np.load(embedding_path).astype(np.float32))
        embedding_path = final_embedding
    run_eval_from_embedding(
        entry=entry,
        main_cfg=main_cfg,
        method=method,
        seed=seed,
        output_dir=output_dir,
        n_clusters=n_clusters,
        embedding_path=embedding_path,
        prefix="sc3_eval",
        variant_name="embedding",
        dry_run=dry_run,
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

    if args.method in TRADITIONAL_EMBEDDING_METHODS:
        run_traditional_embedding(entry, main_cfg, args.method, args.seed, output_dir, dry_run=args.dry_run)
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
