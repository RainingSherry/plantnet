#!/usr/bin/env bash
# Gene-program bottleneck ablation driver (Macosko, seeds 42/43/44).
# See README.md for the falsifiable design. Runs sequentially on one GPU.
set -euo pipefail

cd "$(dirname "$0")/../.."
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python
RUN=experiment_reports/program_bottleneck_ablation_20260703/run_ablation.py
OUT=experiment_reports/program_bottleneck_ablation_20260703/runs
DATA=methods/DeepLearning/scMAEs/benchmark_data/Macosko.h5ad
GPU="${1:-1}"

common="--data_path $DATA --dataset_name Macosko --n_clusters 12 --label_key resolved_label \
  --epochs 80 --variance_weight 0.02 --force_gate 1.0 --gpu $GPU --program_dim 32"

run () {  # name  program_mode  split_mode  type_dim  prog_extra  program_weight  seed
  local name=$1 pmode=$2 split=$3 tdim=$4 pextra=$5 pw=$6 seed=$7
  local dir="$OUT/${name}_seed${seed}"
  if [[ -f "$dir/summary.json" ]]; then echo "skip $dir (done)"; return; fi
  echo "=== $name seed$seed (program=$pmode split=$split type_dim=$tdim pw=$pw) ==="
  $PY $RUN $common --save_dir "$dir" \
    --program_mode "$pmode" --split_mode "$split" --type_dim "$tdim" \
    --prog_extra "$pextra" --program_weight "$pw" --seed "$seed"
}

for s in 42 43 44; do
  # A0: DEC+floor baseline (no program head)
  run a0_baseline        none     none  128 32 0.0  "$s"
  # A1: program as full-latent aux regularizer (two weights bracket the null)
  run a1_prog_w02         nmf      none  128 32 0.02 "$s"
  run a1_prog_w10         nmf      none  128 32 0.10 "$s"
  # A1-shuffle: row-permuted target, same loss form (weight matched to a1_prog_w10)
  run a1_shuffle_w10      shuffled none  128 32 0.10 "$s"
  # A2-fixed: split latent, clustering shrinks to 96 dims (mechanism reverse-test)
  run a2_fixed_w10        nmf      fixed  96 32 0.10 "$s"
  # A2-extra: split but z_type keeps full 128; program uses 32 EXTRA dims
  run a2_extra_w10        nmf      extra 128 32 0.10 "$s"
done

echo "ALL DONE"
$PY experiment_reports/program_bottleneck_ablation_20260703/summarize.py || true
