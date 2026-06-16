#!/usr/bin/env python
from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path


CORE_DATASETS = ["Tosches", "Macosko", "worm_neuron_cell", "Melanoma_5K", "Shekhar", "Guo"]
OPTIONAL_DATASETS = ["Wang", "Pollen"]
METHODS = [
    "nm_scmae_nomix",
    "neighbormix_scmae",
    "random_pseudo_gate_p0.5",
    "random_edge_dropout_keep0.5",
    "random_beta_uniform_0.1",
    "mutual_knn_neighbormix",
    "snn_neighbormix",
    "consensus_neighbormix_threshold0.4",
    "global_random_neighbor_control",
]


@dataclass
class Job:
    dataset: str
    method: str
    seed: int
    gpu: int | None = None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run NeighborMix stochastic ablation matrix")
    parser.add_argument("--root", default="results/experiments/neighbormix_stochastic_regularization_20260616")
    parser.add_argument("--datasets", default=",".join(CORE_DATASETS))
    parser.add_argument("--include_optional", action="store_true")
    parser.add_argument("--seeds", default="0,1,2")
    parser.add_argument("--methods", default=",".join(METHODS))
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


def completed(run_dir: Path) -> bool:
    return (run_dir / "summary.json").exists() and (run_dir / "eval_fixed.csv").exists()


def failed(run_dir: Path) -> bool:
    return (run_dir / "FAILED").exists()


def command_for(job: Job, args: argparse.Namespace, run_dir: Path) -> list[str]:
    cmd = [
        sys.executable,
        "methods/DeepLearning/NeighborMix_scMAE/run_stochastic_ablation.py",
        "--data_path",
        str(dataset_path(job.dataset)),
        "--save_dir",
        str(run_dir),
        "--dataset_name",
        job.dataset,
        "--ablation_method",
        job.method,
        "--method_name",
        job.method,
        "--variant_name",
        job.method,
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
    return cmd


def launch(job: Job, args: argparse.Namespace, run_dir: Path) -> subprocess.Popen:
    run_dir.mkdir(parents=True, exist_ok=True)
    for marker in ["FAILED", "RUNNING"]:
        marker_path = run_dir / marker
        if marker_path.exists():
            marker_path.unlink()
    (run_dir / "RUNNING").write_text(f"dataset={job.dataset} method={job.method} seed={job.seed} gpu={job.gpu}\n", encoding="utf-8")
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
    cmd = command_for(job, args, run_dir)
    return subprocess.Popen(cmd, stdout=stdout, stderr=stderr, env=env)


def summarize(root: Path) -> int:
    cmd = [
        sys.executable,
        "scripts/summarize_neighbormix_stochastic_ablation.py",
        "--root",
        str(root),
    ]
    return subprocess.run(cmd).returncode


def main() -> int:
    args = parse_args()
    root = Path(args.root)
    datasets = split_csv(args.datasets)
    if args.include_optional:
        datasets += [d for d in OPTIONAL_DATASETS if d not in datasets]
    seeds = [int(item) for item in split_csv(args.seeds)]
    methods = split_csv(args.methods)
    gpus = [int(item) for item in split_csv(args.gpus)]
    forbidden = {0, 7}
    bad = sorted(set(gpus).intersection(forbidden))
    if bad:
        raise ValueError(f"GPU 0 and 7 are forbidden; got {bad}")

    jobs = [Job(dataset=d, method=m, seed=s) for d in datasets for s in seeds for m in methods]
    pending: list[tuple[Job, Path]] = []
    skipped = 0
    for job in jobs:
        run_dir = root / job.dataset / job.method / f"seed{job.seed}"
        if completed(run_dir):
            skipped += 1
            continue
        if failed(run_dir) and not args.rerun_failed:
            skipped += 1
            continue
        pending.append((job, run_dir))

    print(f"root={root}")
    print(f"datasets={datasets}")
    print(f"seeds={seeds}")
    print(f"methods={methods}")
    print(f"gpus={gpus}")
    print(f"jobs_total={len(jobs)} skipped={skipped} pending={len(pending)}")
    if args.dry_run:
        for job, run_dir in pending[:20]:
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
            print(f"[launch] gpu={gpu} dataset={job.dataset} seed={job.seed} method={job.method}", flush=True)
            proc = launch(job, args, run_dir)
            active[gpu] = (proc, job, run_dir)

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
                print(f"[done] gpu={gpu} dataset={job.dataset} seed={job.seed} method={job.method}", flush=True)
            else:
                failed_count += 1
                (run_dir / "FAILED").write_text(f"returncode={code}\n", encoding="utf-8")
                print(f"[fail] gpu={gpu} code={code} dataset={job.dataset} seed={job.seed} method={job.method}", flush=True)
            del active[gpu]

    print(f"completed_now={completed_count} failed_now={failed_count}")
    summarize_code = summarize(root)
    if summarize_code != 0:
        return summarize_code
    return 1 if failed_count else 0


if __name__ == "__main__":
    raise SystemExit(main())
