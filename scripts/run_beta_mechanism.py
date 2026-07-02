#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any


EXPERIMENT_ID = "neighbormix_beta_mechanism_20260617"
CORE_DATASETS = ["Tosches", "Macosko", "worm_neuron_cell", "Melanoma_5K", "Shekhar", "Guo"]
FAILURE_DATASETS = ["Tosches", "Macosko", "worm_neuron_cell"]
OPTIONAL_DATASETS = ["Wang", "Pollen"]


STAGE_VARIANTS: dict[str, list[dict[str, Any]]] = {
    "stage1": [
        {"variant": "nm_scmae_nomix", "ablation_method": "nm_scmae_nomix"},
        {"variant": "fixed_beta_0.025", "beta_mode": "fixed", "beta_fixed": 0.025, "beta_max": 0.025},
        {"variant": "random_beta_uniform_0.05", "beta_mode": "uniform", "beta_max": 0.05, "beta_mean": 0.025},
        {"variant": "fixed_beta_0.05", "beta_mode": "fixed", "beta_fixed": 0.05, "beta_max": 0.05},
        {"variant": "random_beta_uniform_0.1", "beta_mode": "uniform", "beta_max": 0.1, "beta_mean": 0.05},
        {"variant": "fixed_beta_0.1", "beta_mode": "fixed", "beta_fixed": 0.1, "beta_max": 0.1},
        {"variant": "random_beta_uniform_0.2", "beta_mode": "uniform", "beta_max": 0.2, "beta_mean": 0.1},
        {"variant": "fixed_beta_0.2", "beta_mode": "fixed", "beta_fixed": 0.2, "beta_max": 0.2},
    ],
    "stage2": [
        {"variant": "fixed_beta_0.05", "beta_mode": "fixed", "beta_fixed": 0.05, "beta_max": 0.05},
        {"variant": "uniform_beta_0.1", "beta_mode": "uniform", "beta_max": 0.1, "beta_mean": 0.05},
        {"variant": "bernoulli_beta_0_or_0.1_p0.5", "beta_mode": "bernoulli", "beta_max": 0.1, "beta_p": 0.5, "beta_mean": 0.05},
        {
            "variant": "truncated_normal_beta_mean0.05_std0.02",
            "beta_mode": "truncated_normal",
            "beta_mean": 0.05,
            "beta_std": 0.02,
            "beta_max": 0.1,
        },
        {
            "variant": "beta_distribution_mean0.05_lowvar",
            "beta_mode": "beta_distribution",
            "beta_mean": 0.05,
            "beta_max": 0.1,
            "beta_concentration": 32.0,
        },
    ],
    "stage3": [
        {
            "variant": "anchor_target_local_mix",
            "beta_mode": "truncated_normal",
            "beta_max": 0.1,
            "beta_mean": 0.05,
            "beta_std": 0.02,
            "target_mode": "anchor",
            "noise_mode": "local_mix",
        },
        {
            "variant": "mixed_target_local_mix",
            "beta_mode": "truncated_normal",
            "beta_max": 0.1,
            "beta_mean": 0.05,
            "beta_std": 0.02,
            "target_mode": "mixed",
            "noise_mode": "local_mix",
        },
        {
            "variant": "global_random_mix_anchor_target",
            "beta_mode": "truncated_normal",
            "beta_max": 0.1,
            "beta_mean": 0.05,
            "beta_std": 0.02,
            "target_mode": "anchor",
            "noise_mode": "global_mix",
        },
        {
            "variant": "gaussian_noise_matched_anchor_target",
            "beta_mode": "truncated_normal",
            "beta_max": 0.1,
            "beta_mean": 0.05,
            "beta_std": 0.02,
            "target_mode": "anchor",
            "noise_mode": "gaussian_matched",
        },
    ],
    "full": [
        {"variant": "nm_scmae_nomix", "ablation_method": "nm_scmae_nomix"},
        {"variant": "fixed_beta_0.1", "beta_mode": "fixed", "beta_fixed": 0.1, "beta_max": 0.1},
        {"variant": "fixed_beta_0.05", "beta_mode": "fixed", "beta_fixed": 0.05, "beta_max": 0.05},
        {"variant": "random_beta_uniform_0.1", "beta_mode": "uniform", "beta_max": 0.1, "beta_mean": 0.05},
        {
            "variant": "truncated_normal_beta_mean0.05_std0.02",
            "beta_mode": "truncated_normal",
            "beta_mean": 0.05,
            "beta_std": 0.02,
            "beta_max": 0.1,
        },
        {"variant": "global_random_neighbor_control", "ablation_method": "global_random_neighbor_control"},
        {"variant": "snn_neighbormix", "ablation_method": "snn_neighbormix"},
        {
            "variant": "gaussian_noise_matched_anchor_target",
            "beta_mode": "truncated_normal",
            "beta_max": 0.1,
            "beta_mean": 0.05,
            "beta_std": 0.02,
            "target_mode": "anchor",
            "noise_mode": "gaussian_matched",
        },
    ],
}

STAGE4_BASE_VARIANTS = [
    {"variant": "fixed_beta_0.1", "beta_mode": "fixed", "beta_fixed": 0.1, "beta_max": 0.1},
    {"variant": "fixed_beta_0.05", "beta_mode": "fixed", "beta_fixed": 0.05, "beta_max": 0.05},
    {"variant": "random_beta_uniform_0.1", "beta_mode": "uniform", "beta_max": 0.1, "beta_mean": 0.05},
]
BAD_EDGE_RATIOS = [0.0, 0.1, 0.2, 0.4]


@dataclass
class Job:
    stage: str
    dataset: str
    variant: dict[str, Any]
    seed: int
    gpu: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run staged NeighborMix beta mechanism experiments")
    parser.add_argument("--root", default=f"results/experiments/{EXPERIMENT_ID}")
    parser.add_argument("--stage", choices=["stage1", "stage2", "stage3", "stage4", "full"], default="stage1")
    parser.add_argument("--datasets", default="")
    parser.add_argument("--include_optional", action="store_true")
    parser.add_argument("--seeds", default="0,1,2")
    parser.add_argument("--variants", default="", help="Comma-separated variant filter for the selected stage")
    parser.add_argument("--gpus", default="1,2,3,4,5,6")
    parser.add_argument("--epochs", type=int, default=80)
    parser.add_argument("--batch_size", type=int, default=256)
    parser.add_argument("--neighbor_k", type=int, default=5)
    parser.add_argument("--mix_neighbors", type=int, default=4)
    parser.add_argument("--n_top_genes", type=int, default=1000)
    parser.add_argument("--target_sum", type=float, default=10000.0)
    parser.add_argument("--poll_interval", type=float, default=10.0)
    parser.add_argument("--dry_run", action="store_true")
    parser.add_argument("--rerun_failed", action="store_true")
    return parser.parse_args()


def split_csv(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(",") if item.strip()]


def dataset_path(dataset: str) -> Path:
    candidates = [
        Path("data/processed_scmae") / f"{dataset}.h5ad",
        Path("data/processed") / f"{dataset}.h5ad",
        Path("data") / f"{dataset}.h5ad",
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError(f"No h5ad found for dataset {dataset!r}; tried: {candidates}")


def stage_datasets(args: argparse.Namespace) -> list[str]:
    if args.datasets:
        datasets = split_csv(args.datasets)
    elif args.stage == "stage4":
        datasets = list(FAILURE_DATASETS)
    elif args.stage == "full":
        datasets = list(CORE_DATASETS) + list(OPTIONAL_DATASETS)
    else:
        datasets = list(CORE_DATASETS)
    if args.include_optional:
        datasets += [d for d in OPTIONAL_DATASETS if d not in datasets]
    return datasets


def stage_variants(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.stage == "stage4":
        variants = []
        for base in STAGE4_BASE_VARIANTS:
            for ratio in BAD_EDGE_RATIOS:
                item = dict(base)
                item["bad_edge_ratio"] = ratio
                item["variant"] = f"{base['variant']}_bad{ratio:g}"
                variants.append(item)
    else:
        variants = [dict(item) for item in STAGE_VARIANTS[args.stage]]
    wanted = set(split_csv(args.variants))
    if wanted:
        variants = [item for item in variants if item["variant"] in wanted]
    if not variants:
        raise ValueError(f"No variants selected for {args.stage}; filter={args.variants!r}")
    return variants


def completed(run_dir: Path) -> bool:
    return (run_dir / "summary.json").exists() and (run_dir / "eval_fixed.csv").exists()


def failed(run_dir: Path) -> bool:
    return (run_dir / "FAILED").exists()


def add_opt(cmd: list[str], name: str, value: Any) -> None:
    if value is None:
        return
    cmd.extend([name, str(value)])


def command_for(job: Job, args: argparse.Namespace, run_dir: Path) -> list[str]:
    variant = job.variant
    method = str(variant.get("ablation_method", "beta_control"))
    cmd = [
        sys.executable,
        "experimental_retired_models/NeighborMix_scMAE/run_beta_mechanism.py",
        "--data_path",
        str(dataset_path(job.dataset)),
        "--save_dir",
        str(run_dir),
        "--dataset_name",
        job.dataset,
        "--ablation_method",
        method,
        "--method_name",
        str(variant["variant"]),
        "--variant_name",
        str(variant["variant"]),
        "--seed",
        str(job.seed),
        "--gpu",
        str(job.gpu),
        "--epochs",
        str(args.epochs),
        "--batch_size",
        str(args.batch_size),
        "--neighbor_k",
        str(args.neighbor_k),
        "--mix_neighbors",
        str(args.mix_neighbors),
        "--n_top_genes",
        str(args.n_top_genes),
        "--target_sum",
        str(args.target_sum),
    ]
    defaults = {
        "beta_mode": "fixed",
        "target_mode": "anchor",
        "noise_mode": "local_mix",
        "oracle_neighbor": "none",
        "bad_edge_ratio": 0.0,
    }
    for key, default in defaults.items():
        add_opt(cmd, f"--{key}", variant.get(key, default))
    for key in [
        "beta_fixed",
        "beta_max",
        "beta_mean",
        "beta_std",
        "beta_p",
        "beta_alpha",
        "beta_beta",
        "beta_concentration",
    ]:
        add_opt(cmd, f"--{key}", variant.get(key))
    return cmd


def launch(job: Job, args: argparse.Namespace, run_dir: Path) -> subprocess.Popen:
    run_dir.mkdir(parents=True, exist_ok=True)
    for marker in ["FAILED", "RUNNING"]:
        marker_path = run_dir / marker
        if marker_path.exists():
            marker_path.unlink()
    cmd = command_for(job, args, run_dir)
    (run_dir / "RUNNING").write_text(
        f"stage={job.stage} dataset={job.dataset} variant={job.variant['variant']} seed={job.seed} gpu={job.gpu}\n",
        encoding="utf-8",
    )
    (run_dir / "command.txt").write_text(" ".join(cmd) + "\n", encoding="utf-8")
    stdout = open(run_dir / "stdout.log", "a", encoding="utf-8")
    stderr = open(run_dir / "stderr.log", "a", encoding="utf-8")
    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"] = str(job.gpu)
    env.setdefault("OMP_NUM_THREADS", "4")
    env.setdefault("MKL_NUM_THREADS", "4")
    env.setdefault("OPENBLAS_NUM_THREADS", "4")
    env.setdefault("NUMBA_NUM_THREADS", "4")
    env.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
    env.setdefault("NUMBA_CACHE_DIR", "/tmp/numba-cache")
    return subprocess.Popen(cmd, stdout=stdout, stderr=stderr, env=env)


def summarize(root: Path) -> int:
    cmd = [sys.executable, "scripts/summarize_beta_mechanism.py", "--root", str(root)]
    return subprocess.run(cmd).returncode


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    datasets = stage_datasets(args)
    seeds = [int(item) for item in split_csv(args.seeds)]
    variants = stage_variants(args)
    gpus = [int(item) for item in split_csv(args.gpus)]
    forbidden = {0, 7}
    bad = sorted(set(gpus).intersection(forbidden))
    if bad:
        raise ValueError(f"GPU 0 and 7 are forbidden; got {bad}")

    jobs = [Job(stage=args.stage, dataset=d, variant=v, seed=s) for d in datasets for s in seeds for v in variants]
    pending: list[tuple[Job, Path]] = []
    skipped = 0
    for job in jobs:
        run_dir = root / job.stage / job.dataset / str(job.variant["variant"]) / f"seed{job.seed}"
        if completed(run_dir):
            skipped += 1
            continue
        if failed(run_dir) and not args.rerun_failed:
            skipped += 1
            continue
        pending.append((job, run_dir))

    print(f"root={root}")
    print(f"stage={args.stage}")
    print(f"datasets={datasets}")
    print(f"seeds={seeds}")
    print(f"variants={[v['variant'] for v in variants]}")
    print(f"gpus={gpus}")
    print(f"jobs_total={len(jobs)} skipped={skipped} pending={len(pending)}")
    if args.dry_run:
        for job, run_dir in pending[:30]:
            print("DRY", job, run_dir)
        return 0

    active: dict[int, tuple[subprocess.Popen, Job, Path]] = {}
    queue = list(pending)
    completed_count = 0
    failed_count = 0
    while queue or active:
        free_gpus = [gpu for gpu in gpus if gpu not in active]
        while free_gpus and queue:
            gpu = free_gpus.pop(0)
            job, run_dir = queue.pop(0)
            job.gpu = gpu
            print(
                f"[launch] gpu={gpu} stage={job.stage} dataset={job.dataset} "
                f"seed={job.seed} variant={job.variant['variant']}",
                flush=True,
            )
            active[gpu] = (launch(job, args, run_dir), job, run_dir)

        time.sleep(float(args.poll_interval))
        for gpu, (proc, job, run_dir) in list(active.items()):
            code = proc.poll()
            if code is None:
                continue
            running = run_dir / "RUNNING"
            if running.exists():
                running.unlink()
            if code == 0 and completed(run_dir):
                completed_count += 1
                print(
                    f"[done] gpu={gpu} stage={job.stage} dataset={job.dataset} "
                    f"seed={job.seed} variant={job.variant['variant']}",
                    flush=True,
                )
            else:
                failed_count += 1
                (run_dir / "FAILED").write_text(f"returncode={code}\n", encoding="utf-8")
                print(
                    f"[fail] gpu={gpu} code={code} stage={job.stage} dataset={job.dataset} "
                    f"seed={job.seed} variant={job.variant['variant']}",
                    flush=True,
                )
            del active[gpu]

    print(f"completed_now={completed_count} failed_now={failed_count}", flush=True)
    summarize_code = summarize(root)
    if summarize_code != 0:
        return summarize_code
    return 1 if failed_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
