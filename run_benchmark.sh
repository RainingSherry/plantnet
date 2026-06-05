#!/bin/bash
# Launch 8 experiments: scMAE and NeighborMix_scMAE on 4 datasets
# GPUs 0 and 7 are FORBIDDEN. Use GPUs 1-6.

ROOT=/home/luolie/biopipeline/dimension-reduction/plantnet
DATA_DIR=$ROOT/data

declare -A CLUSTERS
CLUSTERS[SRP171040]=12
CLUSTERS[SRP182008]=15
CLUSTERS[SRP235541]=18
CLUSTERS[hrvatin_geo_maintype_counts]=8

declare -A DATASETS
DATASETS[SRP171040]=$DATA_DIR/SRP171040.h5ad
DATASETS[SRP182008]=$DATA_DIR/SRP182008.h5ad
DATASETS[SRP235541]=$DATA_DIR/SRP235541.h5ad
DATASETS[hrvatin_geo_maintype_counts]=$DATA_DIR/hrvatin_geo_maintype_counts.h5ad

RUN_SCMAE() {
    local gpu=$1
    local dataset=$2
    local data_path=${DATASETS[$dataset]}
    local nclust=${CLUSTERS[$dataset]}
    local save_dir=$ROOT/results/scMAE/$dataset
    echo "[scMAE GPU $gpu] Starting $dataset (n_clusters=$nclust)"
    CUDA_VISIBLE_DEVICES=$gpu python $ROOT/methods/DeepLearning/scMAE/run.py \
        --data_path "$data_path" \
        --save_dir "$save_dir" \
        --n_clusters $nclust \
        --epochs 80 \
        --batch_size 256 \
        --lr 0.001 \
        --mask_prob 0.4 \
        --hidden_size 128 \
        --seed 42 \
        --gpu 0 \
        2>&1 | tee "$save_dir/stdout.log"
    local exit=$?
    echo "[scMAE GPU $gpu] Done: $dataset (exit=$exit)"
    return $exit
}

RUN_NEIGHBORMIX() {
    local gpu=$1
    local dataset=$2
    local data_path=${DATASETS[$dataset]}
    local nclust=${CLUSTERS[$dataset]}
    local save_dir=$ROOT/results/NeighborMix_scMAE/$dataset
    echo "[NeighborMix GPU $gpu] Starting $dataset (n_clusters=$nclust)"
    CUDA_VISIBLE_DEVICES=$gpu python $ROOT/methods/DeepLearning/NeighborMix_scMAE/run.py \
        --data_path "$data_path" \
        --save_dir "$save_dir" \
        --n_clusters $nclust \
        --epochs 80 \
        --batch_size 256 \
        --lr 0.001 \
        --mask_ratio 0.4 \
        --alpha 0.9 \
        --neighbor_k 5 \
        --mix_neighbors 4 \
        --mix_weight 0.5 \
        --consistency_weight 0.02 \
        --seed 42 \
        --gpu 1 \
        2>&1 | tee "$save_dir/stdout.log"
    local exit=$?
    echo "[NeighborMix GPU $gpu] Done: $dataset (exit=$exit)"
    return $exit
}

echo "=== WAVE 1: scMAE on GPUs 1-4 ==="
# GPU1: SRP171040, GPU2: SRP182008, GPU3: SRP235541, GPU4: hrvatin
RUN_SCMAE 1 SRP171040 &
RUN_SCMAE 2 SRP182008 &
RUN_SCMAE 3 SRP235541 &
RUN_SCMAE 4 hrvatin_geo_maintype_counts &
wait
echo "=== WAVE 1 COMPLETE ==="

echo "=== WAVE 2: NeighborMix_scMAE on GPUs 1-4 ==="
RUN_NEIGHBORMIX 1 SRP171040 &
RUN_NEIGHBORMIX 2 SRP182008 &
RUN_NEIGHBORMIX 3 SRP235541 &
RUN_NEIGHBORMIX 4 hrvatin_geo_maintype_counts &
wait
echo "=== WAVE 2 COMPLETE ==="

echo "ALL EXPERIMENTS COMPLETED"
