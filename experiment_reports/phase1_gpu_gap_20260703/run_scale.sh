#!/usr/bin/env bash
# Phase 1 可扩展性曲线：经典步骤(PCA+KNN/KMeans/eval)随细胞数的超线性增长。
# epochs=2 最小化训练(训练非重点); 每档超时 900s, 超时=该规模 CPU 经典步骤不可行(本身是发现)。
set -uo pipefail
cd "$(dirname "$0")/../.."
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python
M=experiment_reports/phase1_gpu_gap_20260703/measure.py
OUT=experiment_reports/phase1_gpu_gap_20260703/runs
S=/home/luolie/.claude/jobs/23ea9258/tmp/scale_bench/hvg1000
GPU="${1:-2}"

run(){ # name file
  local dir="$OUT/$1"
  [ -f "$dir/gap.json" ] && { echo "skip $1"; return; }
  [ -f "$2" ] || { echo "missing $2"; return; }
  echo "=== $1 $(date +%H:%M:%S) ==="
  timeout 900 $PY $M --data_path "$2" --save_dir "$dir" --dataset_name "$1" \
    --n_clusters 15 --label_key Celltype --epochs 2 --warmup_epochs 1 \
    --device cuda --gpu "$GPU" --cpu_threads 48 || echo "$1 TIMEOUT/FAIL(=CPU不可行)"
}

run scale_50k   $S/real_50k_hvg1000.h5ad
run scale_86k   $S/real_86k_hvg1000.h5ad
run scale_200k  $S/real_200k_hvg1000_xspecies.h5ad
run scale_500k  $S/synth_500k_hvg1000.h5ad
run scale_1M    $S/synth_1M_hvg1000.h5ad
echo "SCALE DONE $(date +%H:%M:%S)"
