#!/bin/bash
set -e

ROOT="/home/luolie/biopipeline/dimension-reduction/plantnet"
DATA_DIR="/data/luolie/biopipeline/dimension-reduction/plantnet/data"
GPU="${1:-0}"
SEED="${2:-1}"
LATENT_DIM=128
EPOCHS=80

# --- BiSSM1D-only on SRP182008 ---
echo "============================================"
echo "Running BiSSM1D-only on SRP182008 (seed=$SEED)"
echo "Plugins: BiSSM1D only"
echo "============================================"
CTRGC_SAVE="$ROOT/results/PlantSPADE_LGCL_protocol_plugins/SRP182008_mamba_alpha05_only_seed${SEED}"
mkdir -p "$CTRGC_SAVE"
export CUDA_VISIBLE_DEVICES="$GPU"
export OPENBLAS_NUM_THREADS=16
export OMP_NUM_THREADS=16
export MKL_NUM_THREADS=16
python "$ROOT/methods/DeepLearning/PlantSPADE_LGCL/run_plantspade.py" \
    --data_path "$DATA_DIR/SRP182008.h5ad" \
    --save_dir "$CTRGC_SAVE" \
    --dataset_name "SRP182008_ctr_gc_only" \
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
    --pola_attention_weight 0.05 \
    --use_mamba true \
    --mamba_d_state 64 \
    --ssm_alpha 0.05 \
    --use_ctr_gc false \
    --ctr_gc_rel_reduction 4

echo "============================================"
echo "BiSSM1D-only run complete!"
echo "============================================"
