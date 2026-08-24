#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import random
import sys
import traceback
from pathlib import Path
from typing import Any

import numpy as np

CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from methods.DeepLearning.APA_scMAE.config import build_arg_parser, config_hash, resolve_config, write_yaml
from methods.DeepLearning.APA_scMAE.data import APAExpressionDataset, load_apa_data, save_text_lines
from methods.DeepLearning.APA_scMAE.model import APAModel
from methods.DeepLearning.APA_scMAE.trainer import APATrainer, save_json


def set_seed(seed: int, deterministic: bool = True) -> None:
    import torch

    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    if deterministic:
        torch.backends.cudnn.deterministic = True
        torch.backends.cudnn.benchmark = False


def resolve_device(config: dict[str, Any]):
    import torch

    runtime = config["runtime"]
    if bool(runtime.get("no_cuda", False)) or not torch.cuda.is_available():
        return torch.device("cpu"), {"physical_gpu": None, "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES", ""), "logical_device": "cpu"}
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    forbidden = {str(x) for x in runtime.get("forbidden_gpus", [0, 7])}
    if visible:
        visible_ids = [x.strip() for x in visible.split(",") if x.strip()]
        bad = set(visible_ids).intersection(forbidden)
        if bad:
            raise ValueError(f"CUDA_VISIBLE_DEVICES={visible!r} contains forbidden GPU(s): {sorted(bad)}")
        gpu_explicit = bool(runtime.get("gpu_explicit", False))
        if gpu_explicit:
            requested_gpu = str(int(runtime["gpu"]))
            if requested_gpu not in visible_ids:
                raise ValueError(f"Requested --gpu {requested_gpu} is not in CUDA_VISIBLE_DEVICES={visible!r}.")
            logical_index = visible_ids.index(requested_gpu)
            physical = int(requested_gpu) if requested_gpu.isdigit() else None
        else:
            logical_index = 0
            physical = int(visible_ids[0]) if visible_ids[0].isdigit() else None
        logical_device = f"cuda:{logical_index}"
        return torch.device(logical_device), {
            "physical_gpu": physical,
            "cuda_visible_devices": visible,
            "logical_device": logical_device,
            "gpu_explicit": gpu_explicit,
        }
    gpu = int(runtime.get("gpu", 1))
    allowed = {int(x) for x in runtime.get("allowed_gpus", [1, 2, 3, 4, 5, 6])}
    if gpu not in allowed or str(gpu) in forbidden:
        raise ValueError(f"Physical GPU {gpu} is not allowed; use 1-6 or --no_cuda.")
    return torch.device(f"cuda:{gpu}"), {"physical_gpu": gpu, "cuda_visible_devices": "", "logical_device": f"cuda:{gpu}"}


def evaluate_embeddings(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int) -> dict[str, Any]:
    from sklearn.cluster import KMeans

    validate_embedding_inputs(embedding, labels, n_clusters)
    pred = KMeans(n_clusters=int(n_clusters), random_state=int(seed), n_init=20).fit_predict(embedding)
    metrics = mapped_clustering_metrics(labels, pred)
    metrics.update({"uses_known_k": True, "oracle-K": True, "cluster_method": "kmeans_known_k"})
    return {"metrics": metrics, "pred": pred.astype(np.int64)}


def mapped_clustering_metrics(labels: np.ndarray, pred: np.ndarray) -> dict[str, float]:
    from sklearn.metrics import (
        accuracy_score,
        adjusted_rand_score,
        completeness_score,
        f1_score,
        fowlkes_mallows_score,
        homogeneity_score,
        v_measure_score,
    )
    from sklearn.metrics.cluster import normalized_mutual_info_score

    mapped = hungarian_map(labels, pred)
    label_values = np.unique(labels)
    return {
        "acc": float(accuracy_score(labels, mapped)),
        "nmi": float(normalized_mutual_info_score(labels, pred)),
        "ari": float(adjusted_rand_score(labels, pred)),
        "f1_macro": float(f1_score(labels, mapped, average="macro", labels=label_values, zero_division=0)),
        "fmi": float(fowlkes_mallows_score(labels, pred)),
        "v_measure": float(v_measure_score(labels, pred)),
        "homogeneity": float(homogeneity_score(labels, pred)),
        "completeness": float(completeness_score(labels, pred)),
    }


def hungarian_map(y_true: np.ndarray, y_pred: np.ndarray) -> np.ndarray:
    from scipy.optimize import linear_sum_assignment

    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    true_values = np.unique(y_true)
    pred_values = np.unique(y_pred)
    n = max(len(true_values), len(pred_values))
    counts = np.zeros((n, n), dtype=np.int64)
    for i, true_value in enumerate(true_values):
        for j, pred_value in enumerate(pred_values):
            counts[i, j] = int(np.sum((y_true == true_value) & (y_pred == pred_value)))
    rows, cols = linear_sum_assignment(-counts)
    mapped = np.full_like(y_pred, fill_value=-1, dtype=np.int64)
    for row, col in zip(rows, cols):
        if row < len(true_values) and col < len(pred_values):
            mapped[y_pred == pred_values[col]] = true_values[row]
    return mapped


def validate_embedding_inputs(embedding: np.ndarray, labels: np.ndarray | None, n_clusters: int) -> None:
    emb = np.asarray(embedding)
    if emb.ndim != 2:
        raise ValueError(f"embedding must be 2D, got shape {emb.shape}")
    if not np.isfinite(emb).all():
        raise ValueError("embedding contains NaN or Inf values")
    if int(n_clusters) < 1:
        raise ValueError(f"n_clusters must be >= 1, got {n_clusters}")
    if int(n_clusters) > int(emb.shape[0]):
        raise ValueError(f"n_clusters={n_clusters} cannot exceed n_cells={emb.shape[0]}")
    if labels is not None and len(labels) != int(emb.shape[0]):
        raise ValueError(f"labels length {len(labels)} does not match embedding rows {emb.shape[0]}")


def predict_clusters(embedding: np.ndarray, n_clusters: int, seed: int) -> np.ndarray:
    from sklearn.cluster import KMeans

    validate_embedding_inputs(embedding, None, n_clusters)
    return KMeans(n_clusters=int(n_clusters), random_state=int(seed), n_init=20).fit_predict(embedding).astype(np.int64)


def try_leiden_fixed(embedding: np.ndarray, labels: np.ndarray, config: dict[str, Any]) -> dict[str, Any]:
    try:
        from methods.DeepLearning.CAAM_scMAE.evaluation.local_metrics import leiden_fixed

        _legacy_metrics, pred = leiden_fixed(
            embedding,
            labels,
            resolution=float(config["evaluation"]["leiden_fixed_resolution"]),
            n_neighbors=int(config["evaluation"]["n_neighbors"]),
            seed=int(config["seed"]),
        )
        pred = np.asarray(pred, dtype=np.int64)
        validate_embedding_inputs(embedding, labels, max(1, int(np.unique(pred).size)))
        metrics = mapped_clustering_metrics(labels, pred)
        metrics.update(
            {
                "uses_known_k": False,
                "oracle-K": False,
                "cluster_method": "leiden_fixed",
            }
        )
        return {"status": "success", "metrics": metrics, "pred": pred.astype(np.int64)}
    except (ImportError, ModuleNotFoundError) as exc:
        return {
            "status": "skipped",
            "reason": str(exc),
            "metrics": {
                "uses_known_k": False,
                "oracle-K": False,
                "cluster_method": "leiden_fixed",
                "status": "skipped",
            },
            "pred": None,
        }


def artifact_manifest(config: dict[str, Any], embedding_shape: tuple[int, int], prediction_status: str) -> dict[str, Any]:
    required = [
        "args.json",
        "runtime.json",
        "embedding_final.npy",
        "metrics.json",
        "training_history.json",
        "corruption_stats.json",
        "mask_stats.json",
        "gradient_stats.json",
        "gene_stats.json",
        "gene_stats.npy",
        "prototypes.npy",
        "selected_gene_indices.npy",
        "selected_genes.txt",
        "resolved_config.yaml",
        "model_checkpoint.pth",
        "artifact_manifest.json",
    ]
    if prediction_status != "skipped":
        required.insert(required.index("metrics.json"), "pred_labels.npy")
    if not bool(config.get("skip_eval", False)):
        required.insert(required.index("pred_labels.npy") if "pred_labels.npy" in required else required.index("metrics.json"), "labels.npy")
    return {
        "status": "complete",
        "method": config.get("method_name", "apa_scmae"),
        "dataset": config.get("dataset_name"),
        "seed": int(config["seed"]),
        "config_hash": config_hash(config),
        "embedding_shape": [int(embedding_shape[0]), int(embedding_shape[1])],
        "input_mode_resolved": config.get("preprocessing", {}).get("input_mode_resolved"),
        "feature_space_source": config.get("preprocessing", {}).get("feature_space_source"),
        "label_key": config.get("preprocessing", {}).get("label_key"),
        "labels_available": bool(config.get("preprocessing", {}).get("labels_available", False)),
        "skip_eval": bool(config.get("skip_eval", False)),
        "prediction_status": prediction_status,
        "training_status": "complete",
        "required_files": required,
    }


def validate_required_files(save_dir: Path, manifest: dict[str, Any]) -> None:
    missing = [name for name in manifest["required_files"] if name != "artifact_manifest.json" and not (save_dir / name).exists()]
    if missing:
        raise RuntimeError(f"Cannot write complete artifact manifest; missing required files: {missing}")


def attention_guard(config: dict[str, Any], n_genes: int) -> dict[str, Any]:
    batch_size = int(config["training"]["batch_size"])
    heads = int(config["model"]["attention_heads"])
    elements = int(batch_size * heads * (int(n_genes) + 1) ** 2)
    max_elements = int(config["runtime"].get("max_attention_elements", 300_000_000))
    info = {
        "attention_elements": elements,
        "max_attention_elements": max_elements,
        "force_large_attention": bool(config["runtime"].get("force_large_attention", False)),
        "warning": None,
    }
    if elements > max_elements:
        message = (
            f"APA-scMAE gene-axis attention would allocate about {elements} attention elements "
            f"(batch_size={batch_size}, heads={heads}, n_genes={n_genes}), exceeding limit {max_elements}. "
            "Lower --batch_size or --n_top_genes, or pass --force_large_attention true."
        )
        if not bool(config["runtime"].get("force_large_attention", False)):
            raise ValueError(message)
        info["warning"] = message
    return info


def embedding_extraction_batch_size(config: dict[str, Any]) -> int:
    return int(config["training"]["batch_size"])


def validate_runtime_config(config: dict[str, Any]) -> None:
    if bool(config.get("skip_eval", False)) and int(config.get("n_clusters", 0)) <= 0:
        raise ValueError("skip_eval=true requires n_clusters > 0 because prediction-only mode cannot infer K from labels.")


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    save_dir = Path(args.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    save_json(save_dir / "args.json", vars(args))
    try:
        import torch

        config = resolve_config(args)
        config["dataset_name"] = config.get("dataset_name") or Path(config["data_path"]).stem
        validate_runtime_config(config)
        set_seed(int(config["seed"]), deterministic=bool(config["runtime"]["deterministic"]))
        device, runtime_info = resolve_device(config)

        bundle = load_apa_data(
            config["data_path"],
            input_mode=str(config["preprocessing"]["input_mode"]),
            target_sum=float(config["preprocessing"]["target_sum"]),
            n_top_genes=int(config["preprocessing"]["n_top_genes"]),
            scale_input=bool(config["preprocessing"]["scale_input"]),
            n_prototypes=int(config["prototype"]["n_prototypes"]),
            pca_dim=int(config["prototype"]["pca_dim"]),
            seed=int(config["seed"]),
            require_labels=not bool(config.get("skip_eval", False)),
            label_key=config.get("evaluation", {}).get("label_key"),
        )
        config["preprocessing"].update(bundle.preprocess_config)
        runtime_info.update(attention_guard(config, int(bundle.x.shape[1])))
        save_json(save_dir / "runtime.json", runtime_info)
        write_yaml(save_dir / "resolved_config.yaml", config)
        np.save(save_dir / "gene_stats.npy", bundle.gene_stats)
        save_json(save_dir / "gene_stats.json", {"shape": list(bundle.gene_stats.shape), "columns": ["mean", "variance", "zero_rate", "hvg_rank"]})
        np.save(save_dir / "prototypes.npy", bundle.prototypes)
        np.save(save_dir / "selected_gene_indices.npy", bundle.selected_gene_indices.astype(np.int64, copy=False))
        save_text_lines(save_dir / "selected_genes.txt", bundle.gene_names)

        dataset = APAExpressionDataset(bundle.x)
        model = APAModel(
            n_genes=int(bundle.x.shape[1]),
            token_dim=int(config["model"]["token_dim"]),
            cell_dim=int(config["model"]["cell_dim"]),
            proto_dim=int(bundle.prototypes.shape[1]),
            attention_heads=int(config["model"]["attention_heads"]),
            dropout=float(config["model"]["attention_dropout"]),
            decoder_mode=str(config["model"]["decoder_mode"]),
        )
        trainer = APATrainer(
            config=config,
            model=model,
            train_dataset=dataset,
            full_x=torch.as_tensor(bundle.x, dtype=torch.float32),
            gene_stats=torch.as_tensor(bundle.gene_stats, dtype=torch.float32, device=device),
            prototypes=torch.as_tensor(bundle.prototypes, dtype=torch.float32, device=device),
            device=device,
            save_dir=save_dir,
        )
        trainer.train()
        trainer.save_diagnostics()
        embedding = trainer.extract_embeddings(batch_size=embedding_extraction_batch_size(config))
        labels = bundle.labels.astype(np.int64) if bundle.labels is not None else None
        np.save(save_dir / "embedding_final.npy", embedding)
        if labels is not None:
            np.save(save_dir / "labels.npy", labels)
        metrics: dict[str, Any] = {"diagnostics": {"embedding_shape": list(embedding.shape)}}
        prediction_status = "skipped"
        if not bool(config.get("skip_eval", False)):
            if labels is None:
                raise ValueError("skip_eval=false requires labels, but no label column was loaded.")
            configured_n_clusters = int(config.get("n_clusters", 0))
            n_clusters = configured_n_clusters if configured_n_clusters > 0 else int(len(np.unique(labels)))
            eval_result = evaluate_embeddings(embedding, labels, n_clusters, int(config["seed"]))
            metrics["kmeans_known_k"] = eval_result["metrics"]
            pred = eval_result["pred"]
            prediction_status = "kmeans_known_k"
            leiden_result = try_leiden_fixed(embedding, labels, config)
            metrics["leiden_fixed"] = leiden_result["metrics"]
            if leiden_result["pred"] is not None:
                np.save(save_dir / "eval_leiden_fixed.npy", leiden_result["pred"])
            else:
                metrics["leiden_fixed"]["skip_reason"] = leiden_result["reason"]
            np.save(save_dir / "pred_labels.npy", pred)
        else:
            n_clusters = int(config.get("n_clusters", 0))
            pred = predict_clusters(embedding, n_clusters, int(config["seed"]))
            np.save(save_dir / "pred_labels.npy", pred)
            prediction_status = "prediction_only_no_labels"
            metrics.update(
                {
                    "evaluation_status": "prediction_only_no_labels",
                    "supervised_metrics": "skipped_no_labels",
                    "prediction": {
                        "cluster_method": "kmeans",
                        "uses_known_k": True,
                        "oracle-K": False,
                        "n_clusters": int(n_clusters),
                    },
                }
            )
        save_json(save_dir / "metrics.json", metrics)
        torch.save(
            {
                "model": trainer.model.state_dict(),
                "teacher": trainer.teacher.state_dict() if trainer.teacher is not None else None,
                "embedding_prototypes": trainer.embedding_prototypes.detach().cpu() if trainer.embedding_prototypes is not None else None,
                "student_optimizer": trainer.student_optimizer.state_dict(),
                "generator_optimizer": trainer.generator_optimizer.state_dict(),
                "resolved_config": config,
            },
            save_dir / "model_checkpoint.pth",
        )
        manifest = artifact_manifest(config, embedding.shape, prediction_status)
        validate_required_files(save_dir, manifest)
        save_json(save_dir / "artifact_manifest.json", manifest)
        save_json(save_dir / "run_manifest.json", {"status": "complete", "method": config["method_name"]})
        return 0
    except Exception as exc:
        save_json(save_dir / "run_manifest.json", {"status": "failed", "error": str(exc)})
        traceback.print_exc()
        raise


if __name__ == "__main__":
    raise SystemExit(main())
