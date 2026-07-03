# NeighborMix-MAE

NeighborMix-MAE is a lightweight masked-expression modeling line for plant single-cell representation learning.

The training corruption is:

```text
x_mix = alpha * x_i + (1 - alpha) * weighted_neighbors(x_i)
mask(x_mix) -> reconstruct x_i
```

Enhanced cells are used only during training. Evaluation extracts embeddings from the original cells.

## Phase 1

Quick validation on the three main Arabidopsis root datasets:

```bash
python methods/DeepLearning/NeighborMix_MAE/scripts/run_experiments.py \
  --phase phase1 \
  --gpus 1,2,3,4,5,6 \
  --jobs 1
```

## Phase 2

Minimal one-factor ablations:

```bash
python methods/DeepLearning/NeighborMix_MAE/scripts/run_experiments.py \
  --phase phase2 \
  --gpus 1,2,3,4,5,6 \
  --jobs 1
```

Physical GPUs `0` and `7` are rejected by the launcher and runner.

