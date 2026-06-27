#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import random
import sys
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
        physical = int(visible_ids[0]) if len(visible_ids) == 1 and visible_ids[0].isdigit() else None
        return torch.device("cuda:0"), {"physical_gpu": physical, "cuda_visible_devices": visible, "logical_device": "cuda:0"}
    gpu = int(runtime.get("gpu", 1))
    allowed = {int(x) for x in runtime.get("allowed_gpus", [1, 2, 3, 4, 5, 6])}
    if gpu not in allowed or str(gpu) in forbidden:
        raise ValueError(f"Physical GPU {gpu} is not allowed; use 1-6 or --no_cuda.")
    return torch.device(f"cuda:{gpu}"), {"physical_gpu": gpu, "cuda_visible_devices": "", "logical_device": f"cuda:{gpu}"}


def evaluate_embeddings(embedding: np.ndarray, labels: np.ndarray, n_clusters: int, seed: int) -> dict[str, Any]:
    from sklearn.cluster import KMeans
    from sklearn.metrics import adjusted_rand_score, completeness_score, f1_score, fowlkes_mallows_score, homogeneity_score
    from sklearn.metrics.cluster import normalized_mutual_info_score

    from methods.evaluation import cluster_acc

    pred = KMeans(n_clusters=int(n_clusters), random_state=int(seed), n_init=20).fit_predict(embedding)
    acc, f1 = cluster_acc(labels, pred)
    metrics = {
        "acc": float(acc),
        "nmi": float(normalized_mutual_info_score(labels, pred)),
        "ari": float(adjusted_rand_score(labels, pred)),
        "f1_macro": float(f1_score(labels, pred, average="macro")),
        "f1_macro_hungarian": float(f1),
        "fmi": float(fowlkes_mallows_score(labels, pred)),
        "homogeneity": float(homogeneity_score(labels, pred)),
        "completeness": float(completeness_score(labels, pred)),
    }
    return {"metrics": metrics, "pred": pred.astype(np.int64)}


def try_leiden_fixed(embedding: np.ndarray, labels: np.ndarray, config: dict[str, Any]) -> dict[str, Any]:
    try:
        from methods.DeepLearning.CAAM_scMAE.evaluation.local_metrics import leiden_fixed

        metrics, pred = leiden_fixed(
            embedding,
            labels,
            resolution=float(config["evaluation"]["leiden_fixed_resolution"]),
            n_neighbors=int(config["evaluation"]["n_neighbors"]),
            seed=int(config["seed"]),
        )
        return {"status": "success", "metrics": metrics, "pred": pred.astype(np.int64)}
    except Exception as exc:
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


def artifact_manifest(config: dict[str, Any], embedding_shape: tuple[int, int]) -> dict[str, Any]:
    return {
        "status": "complete",
        "method": config.get("method_name", "apa_scmae"),
        "dataset": config.get("dataset_name"),
        "seed": int(config["seed"]),
        "config_hash": config_hash(config),
        "embedding_shape": [int(embedding_shape[0]), int(embedding_shape[1])],
        "required_files": [
            "embedding_final.npy",
            "labels.npy",
            "pred_labels.npy",
            "metrics.json",
            "training_history.json",
            "corruption_stats.json",
            "mask_stats.json",
            "gene_stats.json",
            "prototypes.npy",
            "selected_genes.txt",
            "resolved_config.yaml",
            "model_checkpoint.pth",
            "artifact_manifest.json",
        ],
    }


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
        set_seed(int(config["seed"]), deterministic=bool(config["runtime"]["deterministic"]))
        device, runtime_info = resolve_device(config)
        save_json(save_dir / "runtime.json", runtime_info)

        bundle = load_apa_data(
            config["data_path"],
            input_mode=str(config["preprocessing"]["input_mode"]),
            target_sum=float(config["preprocessing"]["target_sum"]),
            n_top_genes=int(config["preprocessing"]["n_top_genes"]),
            scale_input=bool(config["preprocessing"]["scale_input"]),
            n_prototypes=int(config["prototype"]["n_prototypes"]),
            pca_dim=int(config["prototype"]["pca_dim"]),
            seed=int(config["seed"]),
        )
        config["preprocessing"].update(bundle.preprocess_config)
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
        )
        trainer = APATrainer(
            config=config,
            model=model,
            train_dataset=dataset,
            full_x=torch.as_tensor(bundle.x, dtype=torch.float32, device=device),
            gene_stats=torch.as_tensor(bundle.gene_stats, dtype=torch.float32, device=device),
            prototypes=torch.as_tensor(bundle.prototypes, dtype=torch.float32, device=device),
            device=device,
            save_dir=save_dir,
        )
        trainer.train()
        trainer.save_diagnostics()
        embedding = trainer.extract_embeddings(batch_size=max(512, int(config["training"]["batch_size"]) * 2))
        labels = bundle.labels.astype(np.int64)
        np.save(save_dir / "embedding_final.npy", embedding)
        np.save(save_dir / "labels.npy", labels)
        metrics: dict[str, Any] = {"diagnostics": {"embedding_shape": list(embedding.shape)}}
        pred = np.zeros(labels.shape[0], dtype=np.int64)
        if not bool(config.get("skip_eval", False)):
            n_clusters = int(config["n_clusters"]) if int(config["n_clusters"]) > 0 else int(len(np.unique(labels)))
            eval_result = evaluate_embeddings(embedding, labels, n_clusters, int(config["seed"]))
            metrics["kmeans_known_k"] = eval_result["metrics"]
            pred = eval_result["pred"]
            leiden_result = try_leiden_fixed(embedding, labels, config)
            metrics["leiden_fixed"] = leiden_result["metrics"]
            if leiden_result["pred"] is not None:
                np.save(save_dir / "eval_leiden_fixed.npy", leiden_result["pred"])
            else:
                metrics["leiden_fixed"]["skip_reason"] = leiden_result["reason"]
        np.save(save_dir / "pred_labels.npy", pred)
        save_json(save_dir / "metrics.json", metrics)
        torch.save(
            {
                "model": trainer.model.state_dict(),
                "student_optimizer": trainer.student_optimizer.state_dict(),
                "generator_optimizer": trainer.generator_optimizer.state_dict(),
                "resolved_config": config,
            },
            save_dir / "model_checkpoint.pth",
        )
        save_json(save_dir / "artifact_manifest.json", artifact_manifest(config, embedding.shape))
        save_json(save_dir / "run_manifest.json", {"status": "complete", "method": config["method_name"]})
        return 0
    except Exception as exc:
        save_json(save_dir / "run_manifest.json", {"status": "failed", "error": str(exc)})
        print(f"APA-scMAE failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
