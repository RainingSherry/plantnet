#!/usr/bin/env python3
"""并行调度器：conditional-shuffle 消融矩阵分发到 GPU 1-6（禁用 0/7）。

矩阵 = 5 corruption arm (zero / swap_global=S0 / swap_lib=S1 / swap_ndet=S2 / swap_zerolib=S3)
     × 3 数据集 (Macosko / Melanoma_5K / Quake_10x_Spleen)
     × 3 种子 (42/43/44) = 45 runs。
zero = 复现赢家 (zero-mask + DEC + std-floor) 同-harness 对照。
已完成(summary.json)的自动跳过，可断点续跑。
"""
from __future__ import annotations
import subprocess, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY = "/data/luolie/conda/envs/scssl_bench_py310/bin/python"
RUN = str(Path(__file__).resolve().parent / "run_ablation.py")
OUT = Path(__file__).resolve().parent / "runs"

DATASETS = {
    "Macosko": ("methods/DeepLearning/scMAEs/benchmark_data/Macosko.h5ad", 12),
    "Melanoma_5K": ("methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad", 9),
    "Quake_10x_Spleen": ("methods/DeepLearning/scMAEs/benchmark_data/Quake_10x_Spleen.h5ad", 5),
}
ARMS = ["zero", "swap_global", "swap_lib", "swap_ndet", "swap_zerolib"]
SEEDS = [42, 43, 44]
GPUS = [1, 2, 3, 4, 5, 6]


def build_jobs():
    jobs = []
    for ds, (path, k) in DATASETS.items():
        for arm in ARMS:
            for seed in SEEDS:
                name = f"{ds}__{arm}__seed{seed}"
                save_dir = OUT / name
                if (save_dir / "summary.json").exists():
                    continue
                cmd = [PY, RUN,
                       "--data_path", str(ROOT / path), "--save_dir", str(save_dir),
                       "--dataset_name", ds, "--n_clusters", str(k), "--label_key", "resolved_label",
                       "--epochs", "80", "--variance_weight", "0.02", "--force_gate", "1.0",
                       "--corruption", arm, "--seed", str(seed)]
                jobs.append((name, cmd))
    return jobs


def main():
    jobs = build_jobs()
    if not jobs:
        print("All runs already complete.")
        return 0
    print(f"Dispatching {len(jobs)} run(s) across GPUs {GPUS}.")
    running = {}
    free = list(GPUS)
    queue = list(jobs)
    logdir = OUT / "_logs"
    logdir.mkdir(parents=True, exist_ok=True)
    while queue or running:
        while free and queue:
            gpu = free.pop(0)
            name, cmd = queue.pop(0)
            cmd = cmd + ["--gpu", str(gpu)]
            log = open(logdir / f"{name}.log", "w")
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
            running[gpu] = (name, p, log)
            print(f"[launch] GPU{gpu} <- {name}", flush=True)
        time.sleep(10)
        done = []
        for gpu, (name, p, log) in running.items():
            if p.poll() is not None:
                log.close()
                print(f"[done]   GPU{gpu} -> {name} {'OK' if p.returncode==0 else f'FAIL(rc={p.returncode})'}", flush=True)
                done.append(gpu)
        for gpu in done:
            running.pop(gpu); free.append(gpu)
    print("ALL DISPATCHED RUNS FINISHED", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
