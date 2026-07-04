#!/usr/bin/env bash
# Phase 1 缺口实测：顺序跑(避免 CPU 经典步骤争抢污染 wall-time)。
# 单 GPU，线程绑到 48 核。cuda 档跑全部；cpu 档只跑小数据(拿 train 的 CPU/GPU 比)。
set -uo pipefail
cd "$(dirname "$0")/../.."
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python
M=experiment_reports/phase1_gpu_gap_20260703/measure.py
OUT=experiment_reports/phase1_gpu_gap_20260703/runs
BD=methods/DeepLearning/scMAEs/benchmark_data
GPU="${1:-1}"

run(){ # name data k label device threads
  local dir="$OUT/$1"
  [ -f "$dir/gap.json" ] && { echo "skip $1 (done)"; return; }
  echo "=== $1 (device=$5) $(date +%H:%M:%S) ==="
  $PY $M --data_path "$2" --save_dir "$dir" --dataset_name "$1" \
    --n_clusters "$3" --label_key "$4" --epochs 30 --warmup_epochs 10 \
    --device "$5" --gpu "$GPU" --cpu_threads "$6"
}

# cuda 档(全部)
run Melanoma_cuda   $BD/Melanoma_5K.h5ad       9  resolved_label cuda 48
run Quake_cuda      $BD/Quake_10x_Spleen.h5ad  5  resolved_label cuda 48
run Macosko_cuda    $BD/Macosko.h5ad           12 resolved_label cuda 48
# cpu 档(仅小数据, 拿 train CPU/GPU 比; Macosko CPU 训练太慢, 跳过)
run Melanoma_cpu    $BD/Melanoma_5K.h5ad       9  resolved_label cpu  48
run Quake_cpu       $BD/Quake_10x_Spleen.h5ad  5  resolved_label cpu  48

echo "LADDER DONE $(date +%H:%M:%S)"
