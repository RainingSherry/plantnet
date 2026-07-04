#!/usr/bin/env bash
# Phase 2 —— 条件化 / nuisance 匹配 swap-noise 消融 的完整 run 矩阵生成器。
#
# 矩阵 = 5 corruption arm x 3 数据集 x 3 种子 = 45 条 run。
#   corruption: zero(复现赢家) / swap_global(S0) / swap_lib(S1) / swap_ndet(S2) / swap_zerolib(S3)
#   数据集    : Macosko(k=12) / Melanoma_5K(k=9) / Quake_10x_Spleen(k=5)
#   种子      : 42 43 44
#
# 默认 DRY-RUN：只把 45 条命令逐行打印到 stdout（每条一行、可被外部调度器分发到不同 GPU）。
# 真正执行请显式加 `go`。GPU 号可参数化（默认 1；禁止用 GPU 0/7，runner 内已有保护）。
#
# 用法：
#   bash run_all.sh                 # dry-run，打印 45 条命令（GPU=1）
#   bash run_all.sh 3               # dry-run，打印命令（GPU=3）
#   bash run_all.sh 3 go            # 在 GPU 3 上顺序执行（跳过已完成的）
#   bash run_all.sh | parallel -j6  # 交给外部调度器并行分发（自行覆盖 --gpu）
set -euo pipefail

cd "$(dirname "$0")/../.."
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python
RUN=experiment_reports/conditional_shuffle_ablation_20260703/run_ablation.py
OUT=experiment_reports/conditional_shuffle_ablation_20260703/runs

GPU="${1:-1}"
MODE="${2:-dry}"   # dry | go

# 数据集: name|path|n_clusters
DATASETS=(
  "Macosko|methods/DeepLearning/scMAEs/benchmark_data/Macosko.h5ad|12"
  "Melanoma_5K|methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad|9"
  "Quake_10x_Spleen|methods/DeepLearning/scMAEs/benchmark_data/Quake_10x_Spleen.h5ad|5"
)
CORRUPTIONS=(zero swap_global swap_lib swap_ndet swap_zerolib)
SEEDS=(42 43 44)

count=0
for entry in "${DATASETS[@]}"; do
  IFS="|" read -r name path k <<< "$entry"
  for corr in "${CORRUPTIONS[@]}"; do
    for seed in "${SEEDS[@]}"; do
      dir="$OUT/${name}_${corr}_seed${seed}"
      cmd="$PY $RUN --data_path $path --dataset_name $name --n_clusters $k --label_key resolved_label \
--corruption $corr --epochs 80 --variance_weight 0.02 --force_gate 1.0 \
--n_nuisance_bins 10 --n_joint_bins 5 --gpu $GPU --seed $seed --save_dir $dir"
      # 折行归一成单行，便于外部调度器逐条分发
      cmd="$(echo "$cmd" | tr -s ' \\\n' ' ')"
      count=$((count + 1))
      if [[ "$MODE" == "go" ]]; then
        if [[ -f "$dir/summary.json" ]]; then echo "skip $dir (done)"; continue; fi
        echo "=== [$count/45] $name $corr seed$seed (GPU $GPU) ==="
        eval "$cmd"
      else
        echo "$cmd"
      fi
    done
  done
done

if [[ "$MODE" == "go" ]]; then
  echo "ALL DONE ($count runs)"
  $PY experiment_reports/conditional_shuffle_ablation_20260703/summarize.py || true
else
  echo "# dry-run: $count commands above (5 corruption x 3 datasets x 3 seeds). Add 'go' to execute." >&2
fi
