#!/usr/bin/env python3
"""并行调度器：把 decoder-bypass 消融的完整 run 矩阵分发到多张 GPU。

矩阵 = 3 decoder_mode (concat=D0 原始 / none=D1 / lowrank=D2)
     × 3 数据集 (Macosko / Melanoma_5K / Quake_10x_Spleen)
     × 3 种子 (42/43/44)
     = 27 runs。concat 即复现赢家 scMAE+DEC+std-floor（对照）。

每张 GPU 同时最多跑 1 个 run；GPU 池默认 1-6（禁用 0/7）。已完成(summary.json 存在)的
run 自动跳过，可断点续跑。
"""
from __future__ import annotations

import subprocess
import sys
import time
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
DECODER_MODES = ["concat", "none", "lowrank"]
SEEDS = [42, 43, 44]
GPUS = [1, 2, 3, 4, 5, 6]


def build_jobs():
    jobs = []
    for ds, (path, k) in DATASETS.items():
        for mode in DECODER_MODES:
            for seed in SEEDS:
                name = f"{ds}__{mode}__seed{seed}"
                save_dir = OUT / name
                if (save_dir / "summary.json").exists():
                    continue
                cmd = [
                    PY, RUN,
                    "--data_path", str(ROOT / path),
                    "--save_dir", str(save_dir),
                    "--dataset_name", ds,
                    "--n_clusters", str(k),
                    "--label_key", "resolved_label",
                    "--epochs", "80",
                    "--variance_weight", "0.02",
                    "--force_gate", "1.0",
                    "--decoder_mode", mode,
                    "--seed", str(seed),
                ]
                jobs.append((name, cmd))
    return jobs


def main():
    jobs = build_jobs()
    if not jobs:
        print("All runs already complete.")
        return 0
    print(f"Dispatching {len(jobs)} run(s) across GPUs {GPUS}.")
    running = {}  # gpu -> (name, Popen, logfile handle)
    free = list(GPUS)
    queue = list(jobs)
    logdir = OUT / "_logs"
    logdir.mkdir(parents=True, exist_ok=True)

    while queue or running:
        # launch on free GPUs
        while free and queue:
            gpu = free.pop(0)
            name, cmd = queue.pop(0)
            cmd = cmd + ["--gpu", str(gpu)]
            log = open(logdir / f"{name}.log", "w")
            p = subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT)
            running[gpu] = (name, p, log)
            print(f"[launch] GPU{gpu} <- {name}")
        # poll
        time.sleep(10)
        done = []
        for gpu, (name, p, log) in running.items():
            if p.poll() is not None:
                rc = p.returncode
                log.close()
                tag = "OK" if rc == 0 else f"FAIL(rc={rc})"
                print(f"[done]   GPU{gpu} -> {name} {tag}")
                done.append(gpu)
        for gpu in done:
            running.pop(gpu)
            free.append(gpu)

    print("ALL DISPATCHED RUNS FINISHED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
