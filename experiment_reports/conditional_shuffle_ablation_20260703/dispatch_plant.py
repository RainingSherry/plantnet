#!/usr/bin/env python3
"""植物数据关键检验：nuisance-matched swap 在项目真实目标(植物 scRNA)上像 Macosko(受益)
还是像 Melanoma(持平)？

矩阵 = 4 arm (zero / swap_global / swap_lib / swap_ndet) × 植物数据集 × 3 seed。
跳过 swap_zerolib(动物实验里双峰/无信息)。输出到 runs_plant/。GPU 1-6，禁用 0/7。
"""
from __future__ import annotations
import subprocess, time
from pathlib import Path

BASE = Path(__file__).resolve().parent
PY = "/data/luolie/conda/envs/scssl_bench_py310/bin/python"
RUN = str(BASE / "run_ablation.py")
OUT = BASE / "runs_plant"
DATA_ROOT = "/data/luolie/biopipeline/scCluBench/data"

DATASETS = {  # name -> (file, n_clusters, label_key)
    "SRP182008": (f"{DATA_ROOT}/SRP182008.h5ad", 15, "Celltype"),
    "CRA002977_1": (f"{DATA_ROOT}/CRA002977_1.h5ad", 7, "Celltype"),
}
ARMS = ["zero", "swap_global", "swap_lib", "swap_ndet"]
SEEDS = [42, 43, 44]
GPUS = [1, 2, 3, 4, 5, 6]


def build_jobs():
    jobs = []
    for ds, (path, k, lk) in DATASETS.items():
        for arm in ARMS:
            for seed in SEEDS:
                name = f"{ds}__{arm}__seed{seed}"
                if (OUT / name / "summary.json").exists():
                    continue
                jobs.append((name, [
                    PY, RUN, "--data_path", path, "--save_dir", str(OUT / name),
                    "--dataset_name", ds, "--n_clusters", str(k), "--label_key", lk,
                    "--epochs", "80", "--variance_weight", "0.02", "--force_gate", "1.0",
                    "--corruption", arm, "--seed", str(seed)]))
    return jobs


def main():
    jobs = build_jobs()
    if not jobs:
        print("All plant runs complete."); return 0
    print(f"Dispatching {len(jobs)} plant run(s) across GPUs {GPUS}.", flush=True)
    running, free, queue = {}, list(GPUS), list(jobs)
    logdir = OUT / "_logs"; logdir.mkdir(parents=True, exist_ok=True)
    while queue or running:
        while free and queue:
            gpu = free.pop(0); name, cmd = queue.pop(0)
            log = open(logdir / f"{name}.log", "w")
            p = subprocess.Popen(cmd + ["--gpu", str(gpu)], stdout=log, stderr=subprocess.STDOUT)
            running[gpu] = (name, p, log); print(f"[launch] GPU{gpu} <- {name}", flush=True)
        time.sleep(10)
        done = []
        for gpu, (name, p, log) in running.items():
            if p.poll() is not None:
                log.close(); print(f"[done] GPU{gpu} -> {name} {'OK' if p.returncode==0 else 'FAIL'}", flush=True); done.append(gpu)
        for gpu in done:
            running.pop(gpu); free.append(gpu)
    print("PLANT DONE", flush=True); return 0


if __name__ == "__main__":
    raise SystemExit(main())
