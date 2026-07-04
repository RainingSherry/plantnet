#!/usr/bin/env python3
"""诊断: VarFloor 的 DEC 中心初始化来源(embedding vs PCA-KMeans) × 特征数据集 × 3 种子。
配合事后 fusion_eval.py(Step A: PCA⊕VF 融合)。GPU 1-6, 禁 0/7。"""
from __future__ import annotations
import subprocess, time
from pathlib import Path

BASE = Path(__file__).resolve().parent
PY = "/data/luolie/conda/envs/scssl_bench_py310/bin/python"
RUN = str(BASE / "run_variant.py")
OUT = BASE / "runs"

# 按 PCA-vs-VF 两极客观挑: VF占优(Pollen/Quake_10x) + PCA占优(Macosko/Tosches)
DATASETS = {
    "Pollen":           ("/data/luolie/biopipeline/scCluBench/data/processed/Pollen.h5ad", 11, "VF>>PCA"),
    "Quake_10x_Spleen": ("/data/luolie/biopipeline/scCluBench/data/processed_scmae/Quake_10x_Spleen.h5ad", 5, "VF>>PCA"),
    "Macosko":          ("/data/luolie/biopipeline/scCluBench/data/processed_scmae/Macosko.h5ad", 12, "PCA>>VF"),
    "Tosches":          ("/data/luolie/biopipeline/scCluBench/data/processed_scmae/Tosches.h5ad", 15, "PCA>>VF"),
}
INITS = ["embedding", "pca_kmeans"]
SEEDS = [42, 2024, 3407]   # 与全benchmark同种子
GPUS = [1, 2, 3, 4, 5, 6]


def build():
    jobs = []
    for ds, (path, k, _) in DATASETS.items():
        for ci in INITS:
            for s in SEEDS:
                name = f"{ds}__{ci}__seed{s}"
                if (OUT / name / "summary.json").exists():
                    continue
                jobs.append((name, [PY, RUN, "--data_path", path, "--save_dir", str(OUT / name),
                    "--dataset_name", ds, "--n_clusters", str(k), "--label_key", "resolved_label",
                    "--center_init", ci, "--epochs", "80", "--warmup_epochs", "20", "--seed", str(s)]))
    return jobs


def main():
    jobs = build()
    if not jobs:
        print("all done"); return 0
    print(f"Dispatching {len(jobs)} runs across {GPUS}", flush=True)
    running, free, q = {}, list(GPUS), list(jobs)
    logd = OUT / "_logs"; logd.mkdir(parents=True, exist_ok=True)
    while q or running:
        while free and q:
            g = free.pop(0); name, cmd = q.pop(0)
            log = open(logd / f"{name}.log", "w")
            p = subprocess.Popen(cmd + ["--gpu", str(g)], stdout=log, stderr=subprocess.STDOUT)
            running[g] = (name, p, log); print(f"[launch] GPU{g} <- {name}", flush=True)
        time.sleep(10)
        done = []
        for g, (name, p, log) in running.items():
            if p.poll() is not None:
                log.close(); print(f"[done] GPU{g} {name} {'OK' if p.returncode==0 else 'FAIL'}", flush=True); done.append(g)
        for g in done:
            running.pop(g); free.append(g)
    print("ALL DONE", flush=True); return 0


if __name__ == "__main__":
    raise SystemExit(main())
