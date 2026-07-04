#!/usr/bin/env python3
"""Phase 1b：追 Macosko 上 lowrank/none 的 0.872 高峰。

Phase 1 发现去掉/压缩 mask 旁路后 Macosko 出现双峰（lowrank seed43=0.872 >> 赢家 0.70，
但 seed42/44 崩到 ~0.50）。这里加种子刻画分布 + 扫 mask_rank，判断高峰是否可稳定命中。

矩阵:
  - lowrank × mask_rank{4,8,16,32} × seed{42..49}   (主力：mask_rank 是否影响命中率)
  - none    × seed{45..49}                            (补 D1 的种子，与 Phase1 的 42/43/44 合并)
GPU 1-6，禁用 0/7。已完成(summary.json)自动跳过。
"""
from __future__ import annotations
import subprocess, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PY = "/data/luolie/conda/envs/scssl_bench_py310/bin/python"
RUN = str(Path(__file__).resolve().parent / "run_ablation.py")
OUT = Path(__file__).resolve().parent / "runs_phase1b"
DATA = "methods/DeepLearning/scMAEs/benchmark_data/Macosko.h5ad"
GPUS = [1, 2, 3, 4, 5, 6]


def build_jobs():
    jobs = []
    # lowrank × mask_rank × seeds 42-49
    for mr in [4, 8, 16, 32]:
        for seed in range(42, 50):
            name = f"lowrank_mr{mr}_seed{seed}"
            if (OUT / name / "summary.json").exists():
                continue
            jobs.append((name, [
                PY, RUN, "--data_path", str(ROOT / DATA), "--save_dir", str(OUT / name),
                "--dataset_name", "Macosko", "--n_clusters", "12", "--label_key", "resolved_label",
                "--epochs", "80", "--variance_weight", "0.02", "--force_gate", "1.0",
                "--decoder_mode", "lowrank", "--mask_rank", str(mr), "--seed", str(seed)]))
    # none × seeds 45-49 (42/43/44 already in Phase 1 runs/)
    for seed in range(45, 50):
        name = f"none_seed{seed}"
        if (OUT / name / "summary.json").exists():
            continue
        jobs.append((name, [
            PY, RUN, "--data_path", str(ROOT / DATA), "--save_dir", str(OUT / name),
            "--dataset_name", "Macosko", "--n_clusters", "12", "--label_key", "resolved_label",
            "--epochs", "80", "--variance_weight", "0.02", "--force_gate", "1.0",
            "--decoder_mode", "none", "--seed", str(seed)]))
    return jobs


def main():
    jobs = build_jobs()
    if not jobs:
        print("All Phase 1b runs already complete.")
        return 0
    print(f"Dispatching {len(jobs)} Phase-1b run(s) across GPUs {GPUS}.", flush=True)
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
    print("PHASE1B DONE", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
