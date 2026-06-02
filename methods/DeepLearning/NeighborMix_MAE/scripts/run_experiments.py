#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from pathlib import Path

import yaml

SCRIPT_DIR = Path(__file__).resolve().parent
METHOD_DIR = SCRIPT_DIR.parent
ROOT = METHOD_DIR.parents[2]
DATASETS_CONFIG = ROOT / "methods" / "DeepLearning" / "PlantSPADE_LGCL" / "configs" / "datasets_8plant.yaml"


THREAD_ENV = {
    "OMP_NUM_THREADS": "1",
    "MKL_NUM_THREADS": "1",
    "OPENBLAS_NUM_THREADS": "1",
    "NUMBA_NUM_THREADS": "1",
    "VECLIB_MAXIMUM_THREADS": "1",
    "PYTHONUNBUFFERED": "1",
}


def parse_args():
    parser = argparse.ArgumentParser(description="Run NeighborMix-MAE phase1/phase2 experiments.")
    parser.add_argument("--phase", choices=["phase1", "phase2"], default="phase1")
    parser.add_argument("--config", default=None)
    parser.add_argument("--datasets_config", default=str(DATASETS_CONFIG))
    parser.add_argument("--datasets", default=None, help="Comma-separated dataset override.")
    parser.add_argument("--variants", default=None, help="Comma-separated variant override.")
    parser.add_argument("--seeds", default=None, help="Comma-separated seed override.")
    parser.add_argument("--output_dir", default=None)
    parser.add_argument("--gpus", default=None, help="Allowed physical GPUs, e.g. 1,2,3,4,5,6.")
    parser.add_argument("--jobs", type=int, default=None)
    parser.add_argument("--no_cuda", action="store_true")
    parser.add_argument("--skip_eval", action="store_true")
    parser.add_argument("--dry_run", action="store_true")
    return parser.parse_args()


def split_csv(value):
    if value is None or str(value).strip() == "":
        return None
    return [item.strip() for item in str(value).split(",") if item.strip()]


def parse_int_csv(value):
    parsed = split_csv(value)
    if parsed is None:
        return None
    return [int(item) for item in parsed]


def load_yaml(path: Path) -> dict:
    with open(path, "r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def default_config_path(phase: str) -> Path:
    name = "phase1_quick.yaml" if phase == "phase1" else "phase2_ablation.yaml"
    return METHOD_DIR / "configs" / name


def dataset_map(datasets_config: Path) -> dict:
    cfg = load_yaml(datasets_config)
    return {entry["dataset_name"]: entry for entry in cfg.get("datasets", [])}


def resolve_gpus(args, cfg: dict):
    if args.no_cuda:
        return []
    if args.gpus:
        gpus = parse_int_csv(args.gpus)
    else:
        gpus = [int(gpu) for gpu in cfg.get("gpus", [1, 2, 3, 4, 5, 6])]
    if not gpus:
        raise ValueError("No GPUs configured. Use --gpus 1,2,3,4,5,6 or --no_cuda.")
    forbidden = sorted(set(gpus).intersection({0, 7}))
    if forbidden:
        raise ValueError(f"Physical GPUs {forbidden} are forbidden. Use only 1,2,3,4,5,6.")
    return gpus


def list_to_csv(values) -> str:
    if isinstance(values, str):
        return values
    return ",".join(str(v) for v in values)


def add_section_args(cmd: list[str], section: dict, mapping: dict[str, str]) -> None:
    for cfg_key, arg_name in mapping.items():
        if cfg_key in section:
            cmd.extend([f"--{arg_name}", str(section[cfg_key])])


def build_command(
    cfg: dict,
    entry: dict,
    variant: dict,
    seed: int,
    output_dir: Path,
    no_cuda: bool,
    skip_eval: bool,
    gpu: int | None,
):
    prep = cfg.get("preprocessing", {})
    model = cfg.get("model", {})
    train = cfg.get("training", {})
    mix = cfg.get("neighbor_mix", {})
    eval_cfg = cfg.get("evaluation", {})
    method_name = variant.get("method_name", variant["name"])
    run_dir = output_dir / entry["dataset_name"] / variant["name"] / f"seed_{seed}"
    cmd = [
        sys.executable,
        str(METHOD_DIR / "run.py"),
        "--data_path",
        entry["file_path"],
        "--save_dir",
        str(run_dir),
        "--dataset_name",
        entry["dataset_name"],
        "--label_key",
        str(entry.get("label_key", "auto")),
        "--method_name",
        str(method_name),
        "--variant_name",
        str(variant["name"]),
        "--seed",
        str(seed),
        "--n_clusters",
        str(entry.get("expected_n_clusters", 0) or 0),
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
        str(model.get("latent_dim", 32)),
        "--hidden_dim",
        str(model.get("hidden_dim", 512)),
        "--bottleneck_dim",
        str(model.get("bottleneck_dim", 128)),
        "--dropout",
        str(model.get("dropout", 0.1)),
        "--epochs",
        str(train.get("epochs", 80)),
        "--batch_size",
        str(train.get("batch_size", 256)),
        "--lr",
        str(train.get("lr", 0.001)),
        "--weight_decay",
        str(train.get("weight_decay", 1e-5)),
        "--mask_ratio",
        str(variant.get("mask_ratio", train.get("mask_ratio", 0.4))),
        "--mask_strategy",
        str(variant.get("mask_strategy", train.get("mask_strategy", "zero"))),
        "--masked_weight",
        str(train.get("masked_weight", 1.0)),
        "--unmasked_weight",
        str(train.get("unmasked_weight", 0.2)),
        "--mix_mode",
        str(variant.get("mix_mode", "neighbor")),
        "--target_mode",
        str(variant.get("target_mode", "original")),
        "--alpha",
        str(variant.get("alpha", mix.get("alpha", 0.8))),
        "--neighbor_k",
        str(variant.get("neighbor_k", mix.get("neighbor_k", 15))),
        "--mix_neighbors",
        str(variant.get("mix_neighbors", mix.get("mix_neighbors", 4))),
        "--tau",
        str(variant.get("tau", mix.get("tau", 0.2))),
        "--knn_pca_dim",
        str(variant.get("knn_pca_dim", mix.get("knn_pca_dim", 50))),
        "--denoise_weight",
        str(variant.get("denoise_weight", 1.0)),
        "--consistency_weight",
        str(variant.get("consistency_weight", 0.0)),
        "--eval_neighbors",
        str(eval_cfg.get("n_neighbors", 15)),
        "--leiden_fixed_resolution",
        str(eval_cfg.get("leiden_fixed_resolution", 1.0)),
        "--louvain_fixed_resolution",
        str(eval_cfg.get("louvain_fixed_resolution", 1.0)),
        "--leiden_resolutions",
        list_to_csv(eval_cfg.get("leiden_sweep_resolutions", [0.2, 0.4, 0.6, 0.8, 1.0, 1.2])),
        "--include_louvain",
        str(bool(eval_cfg.get("include_louvain", False))).lower(),
        "--run_oracle_sweep",
        str(bool(eval_cfg.get("run_oracle_sweep", False))).lower(),
        "--sweep_max_cells",
        str(eval_cfg.get("sweep_max_cells", 10000)),
        "--silhouette_sample_size",
        str(eval_cfg.get("silhouette_sample_size", 3000)),
        "--no_save_h5ad",
    ]
    if no_cuda:
        cmd.append("--no_cuda")
    else:
        cmd.extend(["--gpu", str(gpu)])
    if skip_eval:
        cmd.extend(["--skip_eval", "true"])
    return cmd, run_dir


def command_for_log(cmd: list[str], env: dict) -> str:
    shown = []
    for key in ["CUDA_VISIBLE_DEVICES", *THREAD_ENV.keys()]:
        if key in env:
            shown.append(f"{key}={env[key]}")
    return " ".join(shown + cmd)


def main():
    args = parse_args()
    config_path = Path(args.config) if args.config else default_config_path(args.phase)
    cfg = load_yaml(config_path)
    datasets = split_csv(args.datasets) or list(cfg.get("datasets", []))
    wanted_variants = set(split_csv(args.variants) or [])
    variants = cfg.get("variants", [])
    if wanted_variants:
        variants = [variant for variant in variants if variant.get("name") in wanted_variants]
    seeds = parse_int_csv(args.seeds) or [int(seed) for seed in cfg.get("seeds", [1, 2, 3])]
    output_dir = Path(args.output_dir or cfg.get("output_dir", ROOT / "results" / "NeighborMix_MAE"))
    output_dir.mkdir(parents=True, exist_ok=True)
    ds_map = dataset_map(Path(args.datasets_config))
    gpus = resolve_gpus(args, cfg)
    jobs = int(args.jobs if args.jobs is not None else cfg.get("jobs", 1))
    if jobs < 1:
        raise ValueError("--jobs must be >= 1")

    tasks = []
    for dataset_name in datasets:
        if dataset_name not in ds_map:
            raise KeyError(f"Dataset {dataset_name!r} not found in {args.datasets_config}")
        for variant in variants:
            for seed in seeds:
                tasks.append((ds_map[dataset_name], variant, int(seed)))
    print(f"Scheduled {len(tasks)} runs from {config_path}")

    running = []
    pending = list(tasks)
    next_gpu = 0
    while pending or running:
        while pending and len(running) < jobs:
            entry, variant, seed = pending.pop(0)
            gpu = None if args.no_cuda else gpus[next_gpu % len(gpus)]
            next_gpu += 1
            cmd, run_dir = build_command(cfg, entry, variant, seed, output_dir, args.no_cuda, args.skip_eval, gpu)
            run_dir.mkdir(parents=True, exist_ok=True)
            log_dir = run_dir / "logs"
            log_dir.mkdir(parents=True, exist_ok=True)
            log_path = log_dir / "run.log"
            env = os.environ.copy()
            env.update(THREAD_ENV)
            if gpu is not None:
                env["CUDA_VISIBLE_DEVICES"] = str(gpu)
            print(f"[launch] gpu={gpu if gpu is not None else 'cpu'} {entry['dataset_name']} {variant['name']} seed={seed}")
            print(command_for_log(cmd, env))
            if args.dry_run:
                continue
            handle = open(log_path, "a", encoding="utf-8")
            handle.write(f"\n\n===== launch {time.strftime('%Y-%m-%d %H:%M:%S')} =====\n")
            handle.write(command_for_log(cmd, env) + "\n")
            handle.flush()
            proc = subprocess.Popen(cmd, cwd=str(ROOT), env=env, stdout=handle, stderr=subprocess.STDOUT)
            running.append({"proc": proc, "handle": handle, "log": log_path, "task": (entry, variant, seed), "gpu": gpu})

        if args.dry_run:
            break
        time.sleep(2.0)
        still = []
        for item in running:
            rc = item["proc"].poll()
            if rc is None:
                still.append(item)
                continue
            item["handle"].close()
            entry, variant, seed = item["task"]
            status = "ok" if rc == 0 else f"failed rc={rc}"
            print(f"[done] {status} gpu={item['gpu']} {entry['dataset_name']} {variant['name']} seed={seed} log={item['log']}")
            if rc != 0:
                raise SystemExit(f"Stopping after failed run: {item['log']}")
        running = still


if __name__ == "__main__":
    main()
