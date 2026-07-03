#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import random
import sys
from pathlib import Path
from typing import Any

CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parents[2]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from methods.DeepLearning.CAAM_scMAE.registry import build_arg_parser, config_hash, resolve_config, write_yaml


def save_json(path: Path, data: Any) -> None:
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, default=str)


def set_seed(seed: int, deterministic: bool = True) -> None:
    import numpy as np
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
    no_cuda = bool(runtime.get("no_cuda", False))
    visible = os.environ.get("CUDA_VISIBLE_DEVICES", "")
    gpu = int(runtime.get("gpu", 1))
    allowed = {int(x) for x in runtime.get("allowed_gpus", [1, 2, 3, 4, 5, 6])}
    forbidden = {str(x) for x in runtime.get("forbidden_gpus", [0, 7])}
    if no_cuda or not torch.cuda.is_available():
        return torch.device("cpu"), {"physical_gpu": None, "cuda_visible_devices": visible, "logical_device": "cpu"}
    if visible:
        visible_ids = [x.strip() for x in visible.split(",") if x.strip()]
        bad = set(visible_ids).intersection(forbidden)
        if bad:
            raise ValueError(f"CUDA_VISIBLE_DEVICES={visible!r} contains forbidden GPU(s): {sorted(bad)}")
        physical = int(visible_ids[0]) if len(visible_ids) == 1 and visible_ids[0].isdigit() else None
        return torch.device("cuda:0"), {"physical_gpu": physical, "cuda_visible_devices": visible, "logical_device": "cuda:0"}
    if gpu not in allowed or str(gpu) in forbidden:
        raise ValueError(f"Physical GPU {gpu} is not allowed; use 1-6 or --no_cuda.")
    return torch.device(f"cuda:{gpu}"), {"physical_gpu": gpu, "cuda_visible_devices": "", "logical_device": f"cuda:{gpu}"}


def write_environment(path: Path) -> None:
    try:
        import torch
        torch_version = torch.__version__
        cuda_available = torch.cuda.is_available()
    except Exception as exc:
        torch_version = f"unavailable: {exc}"
        cuda_available = False
    path.write_text(
        "\n".join([
            f"python={sys.version.split()[0]}",
            f"executable={sys.executable}",
            f"torch={torch_version}",
            f"cuda_available={cuda_available}",
        ]) + "\n",
        encoding="utf-8",
    )


def artifact_manifest(config: dict[str, Any], embedding_shape: tuple[int, int]) -> dict[str, Any]:
    required = [
        "metrics.json",
        "embedding_final.npy",
        "labels.npy",
        "args.json",
        "artifact_manifest.json",
        "preprocess_config.json",
        "selected_gene_indices.npy",
        "selected_genes.txt",
        "corruption_stats.json",
    ]
    return {
        "status": "complete",
        "dataset": config.get("dataset_name"),
        "method": config.get("method_name", "caam_scmae"),
        "variant": config.get("variant"),
        "seed": int(config["seed"]),
        "config_hash": config_hash(config),
        "git_commit": "unknown",
        "embedding_shape": [int(embedding_shape[0]), int(embedding_shape[1])],
        "required_files": required,
        "data_path": config.get("data_path"),
    }


def main() -> int:
    parser = build_arg_parser()
    args = parser.parse_args()
    save_dir = Path(args.save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    save_json(save_dir / "args.json", vars(args))

    try:
        import numpy as np
        import torch

        from methods.DeepLearning.CAAM_scMAE.data.context_selection import select_context_indices
        from methods.DeepLearning.CAAM_scMAE.data.dataset import CAAMExpressionDataset
        from methods.DeepLearning.CAAM_scMAE.data.donor_candidates import DonorCandidateProvider
        from methods.DeepLearning.CAAM_scMAE.data.gene_modules import build_gene_modules, normalized_assignment_dense
        from methods.DeepLearning.CAAM_scMAE.data.preprocessing import load_caam_data
        from methods.DeepLearning.CAAM_scMAE.diagnostics.embedding_stats import embedding_stats
        from methods.DeepLearning.CAAM_scMAE.evaluation.embedding import extract_embeddings
        from methods.DeepLearning.CAAM_scMAE.evaluation.local_metrics import kmeans_known_k, leiden_fixed
        from methods.DeepLearning.CAAM_scMAE.models.caam_model import build_student
        from methods.DeepLearning.CAAM_scMAE.models.common import trainable_parameter_count
        from methods.DeepLearning.CAAM_scMAE.trainers.common import CAAMTrainer

        config = resolve_config(args)
        config["dataset_name"] = config.get("dataset_name") or Path(config["data_path"]).stem
        set_seed(int(config["seed"]), deterministic=bool(config["runtime"]["deterministic"]))
        device, runtime_info = resolve_device(config)
        runtime_info.update({"amp": bool(config["runtime"]["amp"]), "num_workers": int(config["runtime"]["num_workers"])})
        save_json(save_dir / "runtime.json", runtime_info)
        write_environment(save_dir / "environment.txt")

        bundle = load_caam_data(
            config["data_path"],
            input_mode=config["preprocessing"]["input_mode"],
            target_sum=float(config["preprocessing"]["target_sum"]),
            n_top_genes=int(config["preprocessing"]["n_top_genes"]),
            scale_input=bool(config["preprocessing"]["scale_input"]),
            benchmark_mode=bool(config.get("benchmark_mode", False)),
            seed=int(config["seed"]),
        )
        config["preprocessing"].update(
            {
                "feature_space_source": bundle.preprocess_config["feature_space_source"],
                "actual_n_genes_after_selection": int(bundle.x.shape[1]),
                "selected_gene_indices_path": "selected_gene_indices.npy",
                "selected_gene_names_path": "selected_genes.txt",
            }
        )
        write_yaml(save_dir / "resolved_config.yaml", config)
        save_json(save_dir / "dataset_profile.json", bundle.profile)
        save_json(save_dir / "preprocess_config.json", bundle.preprocess_config)
        np.save(save_dir / "selected_gene_indices.npy", bundle.selected_gene_indices.astype(np.int64, copy=False))
        with open(save_dir / "selected_genes.txt", "w", encoding="utf-8") as handle:
            handle.write("\n".join(map(str, bundle.gene_names)) + "\n")

        train_dataset = CAAMExpressionDataset(bundle.x, bundle.batch_code, bundle.library_size, bundle.zero_ratio)
        donor = DonorCandidateProvider(
            bundle.x,
            bundle.batch_code,
            bundle.library_size,
            bundle.zero_ratio,
            candidate_pool_size=int(config["corruption"]["candidate_pool_size"]),
            library_size_bins=int(config["corruption"]["library_size_bins"]),
            zero_ratio_bins=int(config["corruption"]["zero_ratio_bins"]),
            atol=float(config["mask"]["changed_tolerance_abs"]),
            rtol=float(config["mask"]["changed_tolerance_rel"]),
            seed=int(config["seed"]),
        )
        donor.save(save_dir)

        assignment = None
        context_indices = None
        if config["model"]["encoder_type"] == "axial":
            _, sparse_assignment = build_gene_modules(
                bundle.x,
                int(config["axial"]["n_gene_modules"]),
                int(config["axial"]["module_svd_dim"]),
                int(config["axial"]["module_seed"]),
                save_dir,
            )
            assignment = normalized_assignment_dense(sparse_assignment)
            context_indices = select_context_indices(
                bundle.x,
                int(config["axial"]["context_size"]),
                int(config["axial"]["context_pca_dim"]),
                int(config["axial"]["context_seed"]),
                save_dir,
            )

        student = build_student(n_genes=int(bundle.x.shape[1]), config=config, assignment=assignment)
        full_x = torch.as_tensor(bundle.x, dtype=torch.float32, device=device)
        trainer = CAAMTrainer(
            config=config,
            student=student,
            train_dataset=train_dataset,
            donor_provider=donor,
            full_x=full_x,
            device=device,
            save_dir=save_dir,
            context_indices=context_indices,
        )
        trainer.train()
        trainer.save_diagnostics()

        context_tensor = None
        context_idx_tensor = None
        if context_indices is not None:
            context_idx_tensor = torch.as_tensor(context_indices, dtype=torch.long, device=device)
            context_tensor = full_x[context_idx_tensor]
        embedding = extract_embeddings(
            trainer.student,
            train_dataset,
            device,
            batch_size=max(512, int(config["training"]["batch_size"]) * 2),
            context_x=context_tensor,
            context_indices=context_idx_tensor,
        )
        embedding = embedding.astype(np.float32)
        labels = bundle.evaluation.labels.astype(np.int64)
        np.save(save_dir / "embedding_final.npy", embedding)
        np.save(save_dir / "embeddings_base.npy", embedding)
        np.save(save_dir / "labels.npy", labels)

        n_clusters = int(config["n_clusters"])
        if n_clusters <= 0:
            n_clusters = int(len(np.unique(labels)))
        metrics_kmeans, pred_kmeans = kmeans_known_k(embedding, labels, n_clusters, int(config["seed"]))
        metrics_leiden, pred_leiden = leiden_fixed(
            embedding,
            labels,
            resolution=float(config["evaluation"]["leiden_fixed_resolution"]),
            n_neighbors=int(config["evaluation"]["n_neighbors"]),
            seed=int(config["seed"]),
        )
        np.save(save_dir / "eval_kmeans_known_k.npy", pred_kmeans)
        np.save(save_dir / "eval_leiden_fixed.npy", pred_leiden)
        emb_stats = embedding_stats(embedding)
        save_json(save_dir / "metrics.json", {"kmeans_known_k": metrics_kmeans, "leiden_fixed": metrics_leiden, "diagnostics": {"embedding": emb_stats}})
        save_json(save_dir / "embedding_stats.json", emb_stats)
        torch.save(
            {
                "student": trainer.student.state_dict(),
                "generator": trainer.generator.state_dict() if trainer.generator is not None else None,
                "student_optimizer": trainer.student_optimizer.state_dict(),
                "generator_optimizer": trainer.generator_optimizer.state_dict() if trainer.generator_optimizer is not None else None,
                "resolved_config": config,
                "rng_state_torch": torch.get_rng_state(),
                "rng_state_numpy": np.random.get_state(),
                "rng_state_python": random.getstate(),
                "context_indices": context_indices,
                "gene_module_ids": np.load(save_dir / "gene_module_ids.npy") if (save_dir / "gene_module_ids.npy").exists() else None,
            },
            save_dir / "model_checkpoint_last.pt",
        )
        save_json(save_dir / "artifact_manifest.json", artifact_manifest(config, embedding.shape))
        save_json(
            save_dir / "run_manifest.json",
            {
                "status": "complete",
                "variant": config["variant"],
                "student_trainable_params": trainable_parameter_count(trainer.student),
                "generator_trainable_params": trainable_parameter_count(trainer.generator) if trainer.generator is not None else 0,
            },
        )
        return 0
    except Exception as exc:
        save_json(save_dir / "run_manifest.json", {"status": "failed", "error": str(exc)})
        print(f"CAAM-scMAE failed: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
