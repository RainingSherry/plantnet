from __future__ import annotations

import argparse
import copy
import hashlib
import json
from pathlib import Path
from typing import Any


DEFAULT_CONFIG: dict[str, Any] = {
    "method_name": "apa_scmae",
    "seed": 42,
    "skip_eval": False,
    "runtime": {
        "no_cuda": False,
        "gpu": 1,
        "allowed_gpus": [1, 2, 3, 4, 5, 6],
        "forbidden_gpus": [0, 7],
        "deterministic": True,
        "num_workers": 0,
        "fail_fast": True,
        "max_attention_elements": 300_000_000,
        "force_large_attention": False,
    },
    "preprocessing": {
        "input_mode": "auto",
        "target_sum": 10000.0,
        "n_top_genes": 2000,
        "scale_input": False,
    },
    "prototype": {
        "n_prototypes": 16,
        "pca_dim": 32,
    },
    "model": {
        "token_dim": 64,
        "cell_dim": 64,
        "attention_heads": 4,
        "attention_dropout": 0.1,
        "dropout": 0.0,
    },
    "corruption": {
        "type": "scmae_shuffle",
        "changed_tolerance_abs": 1.0e-6,
        "changed_tolerance_rel": 1.0e-5,
    },
    "mask": {
        "ratio": 0.4,
        "temperature": 1.0,
        "masked_data_weight": 0.75,
        "generator_topk_only_effective": True,
    },
    "training": {
        "epochs": 100,
        "batch_size": 16,
        "lr_student": 1.0e-3,
        "lr_generator": 1.0e-4,
        "weight_decay": 1.0e-5,
        "student_grad_clip": 5.0,
        "generator_grad_clip": 1.0,
        "generator_update_interval": 1,
        "gamma": 0.7,
    },
    "generator_loss": {
        "lambda_entropy": 0.01,
        "lambda_balance": 0.01,
        "lambda_distortion": 0.01,
        "lambda_coverage": 0.01,
    },
    "evaluation": {
        "label_key": None,
        "n_neighbors": 15,
        "leiden_fixed_resolution": 1.0,
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


def write_yaml(path: Path, data: dict[str, Any]) -> None:
    import yaml

    with open(path, "w", encoding="utf-8") as handle:
        yaml.safe_dump(data, handle, sort_keys=False, allow_unicode=True)


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="APA-scMAE runner",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--config", type=str, default=None)
    parser.add_argument("--data_path", type=str, required=True)
    parser.add_argument("--save_dir", type=str, required=True)
    parser.add_argument("--dataset_name", type=str, default=None)
    parser.add_argument("--method_name", type=str, default=None)
    parser.add_argument("--n_clusters", type=int, required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--gpu", type=int, default=None)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--input_mode", type=str, default=None, choices=["auto", "raw", "log1p"])
    parser.add_argument("--n_top_genes", type=int, default=None)
    parser.add_argument("--target_sum", type=float, default=None)
    parser.add_argument("--scale_input", type=str2bool, default=None)
    parser.add_argument("--n_prototypes", type=int, default=None)
    parser.add_argument("--pca_dim", type=int, default=None)
    parser.add_argument("--token_dim", type=int, default=None)
    parser.add_argument("--cell_dim", type=int, default=None)
    parser.add_argument("--attention_heads", type=int, default=None)
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--batch_size", type=int, default=None)
    parser.add_argument("--lr_student", type=float, default=None)
    parser.add_argument("--lr_generator", type=float, default=None)
    parser.add_argument("--mask_ratio", type=float, default=None)
    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--corruption_type", type=str, default=None, choices=["scmae_shuffle"])
    parser.add_argument("--skip_eval", type=str2bool, default=None)
    parser.add_argument("--label_key", type=str, default=None)
    parser.add_argument("--max_attention_elements", type=int, default=None)
    parser.add_argument("--force_large_attention", type=str2bool, default=None)
    parser.add_argument("--generator_topk_only_effective", type=str2bool, default=None)
    return parser


def resolve_config(args: argparse.Namespace) -> dict[str, Any]:
    cfg = deep_update(DEFAULT_CONFIG, load_yaml(args.config))
    cli: dict[str, Any] = {
        "data_path": str(args.data_path),
        "save_dir": str(args.save_dir),
        "n_clusters": int(args.n_clusters),
    }
    if args.seed is not None:
        cli["seed"] = int(args.seed)
    if args.dataset_name is not None:
        cli["dataset_name"] = args.dataset_name
    if args.method_name is not None:
        cli["method_name"] = args.method_name
    if args.skip_eval is not None:
        cli["skip_eval"] = bool(args.skip_eval)
    if args.gpu is not None:
        cli.setdefault("runtime", {})["gpu"] = int(args.gpu)
    if args.no_cuda:
        cli.setdefault("runtime", {})["no_cuda"] = True
    if args.max_attention_elements is not None:
        cli.setdefault("runtime", {})["max_attention_elements"] = int(args.max_attention_elements)
    if args.force_large_attention is not None:
        cli.setdefault("runtime", {})["force_large_attention"] = bool(args.force_large_attention)
    prep = cli.setdefault("preprocessing", {})
    if args.input_mode is not None:
        prep["input_mode"] = args.input_mode
    if args.n_top_genes is not None:
        prep["n_top_genes"] = int(args.n_top_genes)
    if args.target_sum is not None:
        prep["target_sum"] = float(args.target_sum)
    if args.scale_input is not None:
        prep["scale_input"] = bool(args.scale_input)
    proto = cli.setdefault("prototype", {})
    if args.n_prototypes is not None:
        proto["n_prototypes"] = int(args.n_prototypes)
    if args.pca_dim is not None:
        proto["pca_dim"] = int(args.pca_dim)
    model = cli.setdefault("model", {})
    if args.token_dim is not None:
        model["token_dim"] = int(args.token_dim)
    if args.cell_dim is not None:
        model["cell_dim"] = int(args.cell_dim)
    if args.attention_heads is not None:
        model["attention_heads"] = int(args.attention_heads)
    train = cli.setdefault("training", {})
    if args.epochs is not None:
        train["epochs"] = int(args.epochs)
    if args.batch_size is not None:
        train["batch_size"] = int(args.batch_size)
    if args.lr_student is not None:
        train["lr_student"] = float(args.lr_student)
    if args.lr_generator is not None:
        train["lr_generator"] = float(args.lr_generator)
    mask = cli.setdefault("mask", {})
    if args.mask_ratio is not None:
        mask["ratio"] = float(args.mask_ratio)
    if args.temperature is not None:
        mask["temperature"] = float(args.temperature)
    if args.generator_topk_only_effective is not None:
        mask["generator_topk_only_effective"] = bool(args.generator_topk_only_effective)
    if args.corruption_type is not None:
        cli.setdefault("corruption", {})["type"] = str(args.corruption_type)
    if args.label_key is not None:
        cli.setdefault("evaluation", {})["label_key"] = str(args.label_key)
    return deep_update(cfg, cli)


def config_hash(config: dict[str, Any]) -> str:
    payload = json.dumps(config, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]
