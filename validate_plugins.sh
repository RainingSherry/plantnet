#!/bin/bash
# Validation script: compare original PlantSPADE-LGCL vs plugin-enhanced version
# Datasets: SRP182008 and SRP171040
# Plugins: FCR + PolaLinearAttention + BiMamba2 + CTRGC

set -e

ROOT="/home/luolie/biopipeline/dimension-reduction/plantnet"
DATA_DIR="/data/luolie/biopipeline/dimension-reduction/plantnet/data"
RESULT_DIR="$ROOT/results/PlantSPADE_LGCL_protocol_plugins"

mkdir -p "$RESULT_DIR"

SEED=1
EPOCHS=80
LATENT_DIM=32

run_baseline() {
    local dataset=$1
    local name=$2
    local save_dir="$RESULT_DIR/${name}_baseline_seed${SEED}"
    echo "============================================"
    echo "Running BASELINE on $dataset (seed=$SEED)"
    echo "============================================"
    python "$ROOT/methods/DeepLearning/PlantSPADE_LGCL/run_plantspade.py" \
        --data_path "$DATA_DIR/${dataset}.h5ad" \
        --save_dir "$save_dir" \
        --dataset_name "$name" \
        --seed "$SEED" \
        --gpu 0 \
        --latent_dim "$LATENT_DIM" \
        --epochs "$EPOCHS" \
        --n_top_genes 2000 \
        --layers 2 \
        --temperature 0.2 \
        --contrastive_weight 0.05 \
        --module_weight 0.001 \
        --num_modules 16 \
        --module_top_k 30 \
        --negative_sampler random_zero \
        --n_clusters 0 \
        --leiden_fixed_resolution 1.0 \
        --eval_neighbors 15 \
        --use_support_attention false \
        --use_trainable_attention_refiner false \
        --use_gated_fusion false
}

run_with_plugins() {
    local dataset=$1
    local name=$2
    local save_dir="$RESULT_DIR/${name}_plugins_seed${SEED}"
    echo "============================================"
    echo "Running PLUGIN-ENHANCED on $dataset (seed=$SEED)"
    echo "Plugins: FCR + PolaLinearAttention + BiMamba2 + CTRGC"
    echo "============================================"
    python "$ROOT/methods/DeepLearning/PlantSPADE_LGCL/run_plantspade.py" \
        --data_path "$DATA_DIR/${dataset}.h5ad" \
        --save_dir "$save_dir" \
        --dataset_name "${name}_plugins" \
        --seed "$SEED" \
        --gpu 0 \
        --latent_dim "$LATENT_DIM" \
        --epochs "$EPOCHS" \
        --n_top_genes 2000 \
        --layers 2 \
        --temperature 0.2 \
        --contrastive_weight 0.05 \
        --module_weight 0.001 \
        --num_modules 16 \
        --module_top_k 30 \
        --negative_sampler random_zero \
        --n_clusters 0 \
        --leiden_fixed_resolution 1.0 \
        --eval_neighbors 15 \
        --use_support_attention false \
        --use_trainable_attention_refiner false \
        --use_gated_fusion false \
        --use_fcr false \
        --fcr_weight 0.05 \
        --use_pola_attention false \
        --pola_num_heads 8 \
        --pola_alpha 4.0 \
        --pola_attention_weight 0.05 \
        --use_mamba false \
        --mamba_d_state 64 \
        --use_ctr_gc true \
        --ctr_gc_rel_reduction 8
}

# SRP182008 (K=15, 13,514 cells)
run_baseline "SRP182008" "SRP182008"
run_with_plugins "SRP182008" "SRP182008"

# SRP171040 (K=12, 33,956 cells)
run_baseline "SRP171040" "SRP171040"
run_with_plugins "SRP171040" "SRP171040"

echo "============================================"
echo "ALL RUNS COMPLETED"
echo "Results saved to: $RESULT_DIR"
echo "============================================"
