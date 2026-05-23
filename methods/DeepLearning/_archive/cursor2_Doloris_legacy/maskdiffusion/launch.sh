#!/bin/bash
# Launcher for cursor2_Doloris maskdiffusion experiments
# Usage: bash launch.sh

BASE="/home/luolie/biopipeline/dimension-reduction/plantnet"
PYTHON="python"

launch() {
    local name="$1"
    local data="$2"
    local save="$3"
    local gpu="$4"

    echo "=== Launching $name on GPU $gpu ==="
    setsid $PYTHON -u "$BASE/methods/DeepLearning/cursor2_Doloris/maskdiffusion/run.py" \
        --data_path "$data" \
        --save_dir "$save" \
        --n_top_genes 2000 \
        --latent_dim 32 \
        --epochs 150 \
        --warmup_epochs 30 \
        --batch_size 256 \
        --lr 1e-3 \
        --weight_decay 1e-4 \
        --mask_loss_weight 0.2 \
        --recon_loss_weight 0.8 \
        --diffusion_loss_weight 0.1 \
        --cluster_loss_weight 0.0 \
        --diffusion_steps 100 \
        --hidden_dim 256 \
        --diffusion_hidden_dim 256 \
        --dropout 0.1 \
        --gpu "$gpu" \
        --seed 42 \
        --eval_interval 10 \
        --cluster_methods kmeans,leiden \
        > "$save/train.log" 2>&1 &

    echo "  PID: $! (log: $save/train.log)"
}

mkdir -p "$BASE/results/scspade_cursor2/Mouse_Pancreas_1"
mkdir -p "$BASE/results/scspade_cursor2/SRP182008"

launch "Mouse_Pancreas_1" \
    "$BASE/data/Mouse_Pancreas_1.h5ad" \
    "$BASE/results/scspade_cursor2/Mouse_Pancreas_1" \
    2

launch "SRP182008" \
    "$BASE/data/SRP182008.h5ad" \
    "$BASE/results/scspade_cursor2/SRP182008" \
    3

echo "=== All jobs launched. Monitor with: tail -f results/scspade_cursor2/{dataset}/train.log ==="
