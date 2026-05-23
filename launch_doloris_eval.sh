#!/bin/bash
###############################################################################
# Unified evaluation launcher for Doloris model family
# Evaluates 6 models × 2 datasets = 12 runs on GPUs 1-6 (6 in first wave, 6 in second)
###############################################################################
set -u

PROJECT_DIR="/home/luolie/biopipeline/dimension-reduction/plantnet"
RESULTS_DIR="$PROJECT_DIR/results"
DATA_DIR="$PROJECT_DIR/data"
PYTHON="/data/luolie/conda/envs/scclubench-main/bin/python"

AVAILABLE_GPUS=(1 2 3 4 5 6)
N_GPUS=${#AVAILABLE_GPUS[@]}  # 6

# --------------------------------------------------------------------------
# Build command: Doloris/GraphDiffusion on SRP182008 (GPU 1)
# --------------------------------------------------------------------------
build_doloris_graphdiff() {
    local dataset="$1" gpu="$2"
    local data_path="$DATA_DIR/${dataset}.h5ad"
    local save_dir="$RESULTS_DIR/Doloris_GraphDiffusion/${dataset}"
    local label_col
    [[ "$dataset" == "SRP182008" ]] && label_col="Celltype" || label_col="cell_type"
    mkdir -p "$save_dir"
    "$PYTHON" "$PROJECT_DIR/methods/DeepLearning/Doloris/GraphDiffusion/run.py" \
        --data_path "$data_path" \
        --save_dir "$save_dir" \
        --label_col "$label_col" \
        --n_top_genes 2000 \
        --gene_feature_dim 64 \
        --gene_hidden_dim 64 \
        --gene_output_dim 64 \
        --cell_dim 64 \
        --epochs 80 \
        --lr 0.001 \
        --weight_decay 0.0001 \
        --corr_threshold 0.35 \
        --coexpr_top_k 20 \
        --support_weight_mode log1p_count \
        --seed 42 \
        --gpu "$gpu" 2>&1
}

# --------------------------------------------------------------------------
# Build command: cursor_Doloris/GraphDiffusion on SRP182008 (GPU 2)
# --------------------------------------------------------------------------
build_cursor_graphdiff() {
    local dataset="$1" gpu="$2"
    local data_path="$DATA_DIR/${dataset}.h5ad"
    local save_dir="$RESULTS_DIR/cursor_Doloris_GraphDiffusion/${dataset}"
    mkdir -p "$save_dir"
    cd "$PROJECT_DIR/methods/DeepLearning/cursor_Doloris/GraphDiffusion"
    "$PYTHON" train/train_plantdiffcluster.py \
        --data_path "$data_path" \
        --save_dir "$save_dir" \
        --n_hvg 1500 \
        --epochs 80 \
        --batch_size 64 \
        --lr 0.001 \
        --weight_decay 0.0001 \
        --gene_dim 64 \
        --hidden_dim 256 \
        --embed_dim 128 \
        --time_embed_dim 128 \
        --n_layers 2 \
        --pooling_strategy attention \
        --pooling_topk 50 \
        --graph_type coexpression \
        --support_strategy log1p \
        --dropout_rate 0.0 \
        --use_diffusion \
        --use_mask_predictor \
        --num_timesteps 500 \
        --ddim_steps 20 \
        --beta_schedule cosine \
        --cluster_strategy gmm \
        --cell_type_num "$([ "$dataset" = "SRP182008" ] && echo 15 || echo 13)" \
        --device cuda \
        --seed 42 2>&1
}

# --------------------------------------------------------------------------
# Build command: Doloris/DiffusionBridge (GPU 3)
# --------------------------------------------------------------------------
build_doloris_bridge() {
    local dataset="$1" gpu="$2"
    local data_path="$DATA_DIR/${dataset}.h5ad"
    local save_dir="$RESULTS_DIR/Doloris_DiffusionBridge/${dataset}"
    mkdir -p "$save_dir"
    cd "$PROJECT_DIR/methods/DeepLearning/Doloris/DiffusionBridge/train"
    "$PYTHON" run_srp182008_pipeline.py \
        --data-path "$data_path" \
        --output-dir "$save_dir" \
        --latent-dim 64 \
        --hidden-dim 256 \
        --diffusion-steps 50 \
        --batch-size 256 \
        --source-epochs 20 \
        --target-epochs 20 \
        --bridge-epochs 30 \
        --lr 0.001 \
        --recon-weight 1.0 \
        --prior-weight 0.001 \
        --zero-weight 0.25 \
        --teacher-weight 1.0 \
        --cluster-weight 1.0 \
        --gaussian-weight 0.001 \
        --support-weight 0.5 \
        --entropy-weight 0.001 2>&1
}

# --------------------------------------------------------------------------
# Build command: cursor_Doloris/DiffusionBridge (GPU 4)
# --------------------------------------------------------------------------
build_cursor_bridge() {
    local dataset="$1" gpu="$2"
    local data_path="$DATA_DIR/${dataset}.h5ad"
    local save_dir="$RESULTS_DIR/cursor_Doloris_DiffusionBridge/${dataset}"
    mkdir -p "$save_dir"
    cd "$PROJECT_DIR/methods/DeepLearning/cursor_Doloris/DiffusionBridge"
    "$PYTHON" run_pipeline.py \
        --data-path "$data_path" \
        --output-dir "$save_dir" \
        --n-hvg 2000 \
        --latent-dim 64 \
        --hidden-dim 256 \
        --time-embed-dim 128 \
        --diffusion-steps 50 \
        --dropout 0.0 \
        --batch-size 256 \
        --source-epochs 20 \
        --target-epochs 20 \
        --bridge-epochs 30 \
        --lr 0.001 \
        --teacher-mode pca_graph \
        --warmup-epochs 5 \
        --teacher-weight 1.0 \
        --cluster-weight 1.0 \
        --gaussian-weight 0.001 \
        --support-weight 0.5 \
        --entropy-weight 0.001 \
        --support-topk 256 \
        --support-blend 0.2 2>&1
}

# --------------------------------------------------------------------------
# Build command: Doloris/maskdiffusion (GPU 5)
# --------------------------------------------------------------------------
build_doloris_mask() {
    local dataset="$1" gpu="$2"
    local data_path="$DATA_DIR/${dataset}.h5ad"
    local save_dir="$RESULTS_DIR/Doloris_maskdiffusion/${dataset}"
    mkdir -p "$save_dir"
    "$PYTHON" "$PROJECT_DIR/methods/DeepLearning/Doloris/maskdiffusion/run.py" \
        --data_path "$data_path" \
        --save_dir "$save_dir" \
        --n_top_genes 2000 \
        --latent_dim 16 \
        --hidden_dim 256 \
        --diffusion_hidden_dim 128 \
        --diffusion_steps 100 \
        --mask_epochs 20 \
        --embedding_epochs 50 \
        --batch_size 256 \
        --lr 0.001 \
        --mask_threshold 0.5 \
        --cluster_method leiden \
        --seed 42 \
        --gpu "$gpu" 2>&1
}

# --------------------------------------------------------------------------
# Build command: cursor_Doloris/maskdiffusion (GPU 6)
# --------------------------------------------------------------------------
build_cursor_mask() {
    local dataset="$1" gpu="$2"
    local data_path="$DATA_DIR/${dataset}.h5ad"
    local save_dir="$RESULTS_DIR/cursor_Doloris_maskdiffusion/${dataset}"
    mkdir -p "$save_dir"
    cd "$PROJECT_DIR/methods/DeepLearning/cursor_Doloris/maskdiffusion"
    "$PYTHON" run.py \
        --data_path "$data_path" \
        --save_dir "$save_dir" \
        --n_clusters 0 \
        --n_top_genes 1000 \
        --latent_dim 64 \
        --diffusion_hidden_dims 1024,512 \
        --diffusion_steps 100 \
        --epochs 100 \
        --batch_size 256 \
        --lr 0.001 \
        --warmup_epochs 10 \
        --mask_loss_weight 0.2 \
        --diffusion_loss_weight 0.2 \
        --recon_loss_weight 0.8 \
        --cluster_loss_weight 1.0 \
        --gpu "$gpu" \
        --seed 42 2>&1
}

# --------------------------------------------------------------------------
# Run a single job
# --------------------------------------------------------------------------
run_single() {
    local run_id="$1" model_fn="$2" dataset="$3" gpu="$4"
    local save_dir="$RESULTS_DIR/${model_fn}/${dataset}"
    local log_file="$save_dir/run.log"

    echo "=========================================="
    echo "  Run ID: $run_id"
    echo "  Model:  $model_fn"
    echo "  Dataset: $dataset"
    echo "  GPU:    $gpu"
    echo "  Save:   $save_dir"
    echo "=========================================="

    local start_time=$(date +%s)
    mkdir -p "$save_dir"
    echo "[START] $(date)" > "$log_file"
    echo "[GPU] $gpu" >> "$log_file"
    echo "[DATASET] $dataset" >> "$log_file"

    # Do NOT set CUDA_VISIBLE_DEVICES here -- the --gpu N argument inside each
    # script directly selects which GPU to use (ordinal N).
    # Setting CUDA_VISIBLE_DEVICES would remap GPU N to ordinal 0, causing
    # "invalid device ordinal" errors for scripts that use --gpu N.

    case "$model_fn" in
        Doloris_GraphDiffusion)       build_doloris_graphdiff "$dataset" "$gpu" ;;
        cursor_Doloris_GraphDiffusion) build_cursor_graphdiff "$dataset" "$gpu" ;;
        Doloris_DiffusionBridge)       build_doloris_bridge "$dataset" "$gpu" ;;
        cursor_Doloris_DiffusionBridge) build_cursor_bridge "$dataset" "$gpu" ;;
        Doloris_maskdiffusion)         build_doloris_mask "$dataset" "$gpu" ;;
        cursor_Doloris_maskdiffusion)  build_cursor_mask "$dataset" "$gpu" ;;
    esac >> "$log_file" 2>&1

    local ec=$?
    local end_time=$(date +%s)
    local elapsed=$((end_time - start_time))

    echo "[END] $(date)" >> "$log_file"
    echo "[EXIT] $ec (${elapsed}s)" >> "$log_file"

    if [ $ec -ne 0 ]; then
        echo "  WARNING: Run $run_id FAILED (exit $ec, ${elapsed}s)"
        echo "{\"status\":\"error\",\"exit_code\":$ec,\"elapsed_s\":$elapsed}" > "$save_dir/status.json"
    else
        echo "  SUCCESS: Run $run_id (${elapsed}s)"
        echo "{\"status\":\"success\",\"elapsed_s\":$elapsed}" > "$save_dir/status.json"
    fi
    return $ec
}

# --------------------------------------------------------------------------
# Main: launch all 12 runs
# --------------------------------------------------------------------------
export CUDA_DEVICE_ORDER="PCI_BUS_ID"
export PYTHONPATH="$PROJECT_DIR:$PYTHONPATH"
export TMPDIR="${TMPDIR:-/data/tmp}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

echo "=========================================="
echo "  Doloris Model Family Evaluation"
echo "  6 models × 2 datasets = 12 runs"
echo "  GPUs: ${AVAILABLE_GPUS[*]}"
echo "  Python: $PYTHON"
echo "=========================================="
echo ""

# (run_id, model_fn, dataset, gpu)
declare -a TASKS=(
    "1|Doloris_GraphDiffusion|SRP182008|1"
    "2|cursor_Doloris_GraphDiffusion|SRP182008|2"
    "3|Doloris_DiffusionBridge|SRP182008|3"
    "4|cursor_Doloris_DiffusionBridge|SRP182008|4"
    "5|Doloris_maskdiffusion|SRP182008|5"
    "6|cursor_Doloris_maskdiffusion|SRP182008|6"
    "7|Doloris_GraphDiffusion|Mouse_Pancreas_1|1"
    "8|cursor_Doloris_GraphDiffusion|Mouse_Pancreas_1|2"
    "9|Doloris_DiffusionBridge|Mouse_Pancreas_1|3"
    "10|cursor_Doloris_DiffusionBridge|Mouse_Pancreas_1|4"
    "11|Doloris_maskdiffusion|Mouse_Pancreas_1|5"
    "12|cursor_Doloris_maskdiffusion|Mouse_Pancreas_1|6"
)

PIDS=()
for task in "${TASKS[@]}"; do
    IFS='|' read -r run_id model_fn dataset gpu <<< "$task"
    run_single "$run_id" "$model_fn" "$dataset" "$gpu" &
    PIDS+=($!)
    echo "  Launched run $run_id ($model_fn / $dataset) on GPU $gpu (PID ${PIDS[-1]})"
    sleep 3
done

echo ""
echo "All jobs launched. Waiting for completion..."
echo "  Total jobs: ${#PIDS[@]}"
echo ""

FAILURES=0
for i in "${!PIDS[@]}"; do
    pid="${PIDS[$i]}"
    wait $pid
    ec=$?
    if [ $ec -ne 0 ]; then
        echo "  Run $((i+1)) (PID $pid) FAILED ($ec)"
        FAILURES=$((FAILURES+1))
    else
        echo "  Run $((i+1)) (PID $pid) SUCCESS"
    fi
done

echo ""
echo "=========================================="
echo "  All jobs finished"
echo "  Failures: $FAILURES / ${#PIDS[@]}"
echo "=========================================="

exit $FAILURES
