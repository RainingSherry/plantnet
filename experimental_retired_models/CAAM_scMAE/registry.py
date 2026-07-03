from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any


VARIANTS = ("control", "axial", "advmask", "full")
CORRUPTION_TYPES = ("scmae_shuffle", "matched_donor", "nonzero_aware_donor")


DEFAULT_CONFIG: dict[str, Any] = {
    "runtime": {
        "no_cuda": False,
        "allowed_gpus": [1, 2, 3, 4, 5, 6],
        "forbidden_gpus": [0, 7],
        "amp": False,
        "deterministic": True,
        "num_workers": 0,
        "fail_fast": True,
    },
    "preprocessing": {
        "input_mode": "log1p",
        "target_sum": 10000.0,
        "n_top_genes": 2000,
        "scale_input": False,
        "canonical_input": False,
    },
    "training": {
        "epochs": 100,
        "student_warmup_epochs": 20,
        "batch_size": 256,
        "lr_student": 0.001,
        "lr_generator": 0.0001,
        "weight_decay": 0.00001,
        "student_grad_clip": 5.0,
        "generator_grad_clip": 1.0,
        "generator_update_interval": 5,
        "fixed_epochs": True,
        "label_based_early_stop": False,
    },
    "model": {
        "latent_dim": 128,
        "encoder_type": "mlp",
        "mask_selector": "random",
        "decoder_mask_conditioning": "pred_detached",
        "mlp_hidden_dim": 256,
        "dropout": 0.0,
    },
    "mask": {
        "ratio": 0.4,
        "fixed_budget": True,
        "changed_tolerance_abs": 1.0e-6,
        "changed_tolerance_rel": 1.0e-5,
    },
    "corruption": {
        "type": "matched_donor",
        "candidate_pool_size": 32,
        "library_size_bins": 10,
        "zero_ratio_bins": 10,
        "donor_resample_attempts": 8,
        "max_budget_deficit_fraction": 0.01,
        "strict_effective_budget": False,
    },
    "axial": {
        "n_gene_modules": 64,
        "module_svd_dim": 32,
        "module_seed": 0,
        "token_dim": 128,
        "gene_attention_layers": 2,
        "gene_attention_heads": 4,
        "attention_dropout": 0.1,
        "context_size": 256,
        "context_pca_dim": 32,
        "context_seed": 0,
        "context_refresh": "per_epoch",
        "detach_context_kv": True,
    },
    "generator": {
        "hidden_dim": 64,
        "temperature_start": 1.0,
        "temperature_end": 0.3,
        "beta_mask_loss": 0.5,
        "distortion_min": 0.25,
        "distortion_max": 3.0,
    },
    "loss": {
        "lambda_visible": 0.1,
        "lambda_mask": 0.5,
        "lambda_coverage": 0.1,
        "lambda_distortion": 0.1,
        "lambda_entropy": 0.01,
    },
    "evaluation": {
        "n_neighbors": 15,
        "leiden_fixed_resolution": 1.0,
        "louvain_fixed_resolution": 1.0,
        "include_louvain": False,
        "run_oracle_sweep": False,
        "silhouette_sample_size": 3000,
    },
}


def str2bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"1", "true", "t", "yes", "y"}:
        return True
    if text in {"0", "false", "f", "no", "n"}:
        return False
    raise argparse.ArgumentTypeError(f"Expected boolean, got {value!r}")


def deep_update(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    out = copy.deepcopy(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = deep_update(out[key], value)
        else:
            out[key] = value
    return out


def load_yaml(path: str | None) -> dict[str, Any]:
    if not path:
        return {}
    import yaml

    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def apply_variant(config: dict[str, Any], variant: str) -> dict[str, Any]:
    if variant not in VARIANTS:
        raise ValueError(f"Unknown CAAM variant {variant!r}; expected one of {VARIANTS}")
    cfg = copy.deepcopy(config)
    cfg["variant"] = variant
    if variant == "control":
        cfg["model"]["encoder_type"] = "mlp"
        cfg["model"]["mask_selector"] = "random"
    elif variant == "axial":
        cfg["model"]["encoder_type"] = "axial"
        cfg["model"]["mask_selector"] = "random"
    elif variant == "advmask":
        cfg["model"]["encoder_type"] = "mlp"
        cfg["model"]["mask_selector"] = "adversarial"
    elif variant == "full":
        cfg["model"]["encoder_type"] = "axial"
        cfg["model"]["mask_selector"] = "adversarial"
    return cfg


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="CAAM-scMAE runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, required=True)
    parser.add_argument("--dataset_name", type=str, default=None)
    parser.add_argument("--method_name", type=str, default="caam_scmae")
    parser.add_argument("--variant", type=str, default="full", choices=VARIANTS)
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--gpu", type=int, default=1)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--benchmark_mode", type=str2bool, default=False)
    parser.add_argument("--input_mode", type=str, default=None, choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=None)
    parser.add_argument("--target_sum", type=float, default=None)
    parser.add_argument("--scale_input", type=str2bool, default=None)
    parser.add_argument("--skip_eval", type=str2bool, default=False)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--batch_size", type=int, default=None)
    parser.add_argument("--lr_student", type=float, default=None)
    parser.add_argument("--lr_generator", type=float, default=None)
    parser.add_argument("--mask_ratio", type=float, default=None)
    parser.add_argument("--latent_dim", type=int, default=None)
    parser.add_argument("--mlp_hidden_dim", type=int, default=None)
    parser.add_argument("--corruption_type", type=str, default=None, choices=CORRUPTION_TYPES)
    parser.add_argument("--strict_effective_budget", type=str2bool, default=None)
    return parser


def resolve_config(args: argparse.Namespace) -> dict[str, Any]:
    cfg = deep_update(DEFAULT_CONFIG, load_yaml(args.config))
    cli: dict[str, Any] = {
        "seed": int(args.seed),
        "data_path": str(args.data_path),
        "save_dir": str(args.save_dir),
        "dataset_name": args.dataset_name,
        "method_name": args.method_name,
        "n_clusters": int(args.n_clusters),
        "benchmark_mode": bool(args.benchmark_mode),
        "skip_eval": bool(args.skip_eval),
        "resume": bool(args.resume),
        "overwrite": bool(args.overwrite),
    }
    if args.no_cuda:
        cli.setdefault("runtime", {})["no_cuda"] = True
    cli.setdefault("runtime", {})["gpu"] = int(args.gpu)
    prep = cli.setdefault("preprocessing", {})
    if args.input_mode is not None:
        prep["input_mode"] = args.input_mode
    if args.n_top_genes is not None:
        prep["n_top_genes"] = int(args.n_top_genes)
    if args.target_sum is not None:
        prep["target_sum"] = float(args.target_sum)
    if args.scale_input is not None:
        prep["scale_input"] = bool(args.scale_input)
    train = cli.setdefault("training", {})
    if args.epochs is not None:
        train["epochs"] = int(args.epochs)
    if args.batch_size is not None:
        train["batch_size"] = int(args.batch_size)
    if args.lr_student is not None:
        train["lr_student"] = float(args.lr_student)
    if args.lr_generator is not None:
        train["lr_generator"] = float(args.lr_generator)
    if args.mask_ratio is not None:
        cli.setdefault("mask", {})["ratio"] = float(args.mask_ratio)
    if args.latent_dim is not None:
        cli.setdefault("model", {})["latent_dim"] = int(args.latent_dim)
    if args.mlp_hidden_dim is not None:
        cli.setdefault("model", {})["mlp_hidden_dim"] = int(args.mlp_hidden_dim)
    if args.corruption_type is not None:
        cli.setdefault("corruption", {})["type"] = str(args.corruption_type)
    if args.strict_effective_budget is not None:
        cli.setdefault("corruption", {})["strict_effective_budget"] = bool(args.strict_effective_budget)

    cfg = deep_update(cfg, cli)
    if cfg.get("benchmark_mode"):
        cfg["preprocessing"]["input_mode"] = "log1p"
        cfg["preprocessing"]["scale_input"] = False
        cfg["preprocessing"]["canonical_input"] = True
    return apply_variant(cfg, args.variant)


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


def write_yaml(path: Path, data: dict[str, Any]) -> None:
    import yaml

    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)
