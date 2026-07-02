# GatedNeighborMix-scMAE

Reliability-gated fusion of the three mechanisms proven effective across the
scMAE improvement search (see `scMAEs/参考文献/03_scMAE改良失败模式与结构迁移问题总结.md`
and the全benchmark结果 NeighborMix analysis).

## Core idea

One per-cell reliability field `r_i ∈ [0,1]` modulates every smoothing
mechanism at once, so they retreat together on rare/boundary cells and hand
those cells back to the self-anchored scMAE backbone.

- **scMAE backbone** (always on, all cells): swap-noise + BCE mask discriminator
  + weighted reconstruction of the real cell. Self-anchored safety net.
- **NeighborMix branch** (gated): `alpha_i = 1 − (1−alpha_min)·r_i`. Core cells
  (r≈1) mix strongly toward neighbors; rare/boundary cells (r≈0) → alpha→1 =
  mixing OFF = pure scMAE. Pseudo target is the REAL cell (anchor-recovery).
- **DEC cluster centers** (gated): confidence-gated KL further multiplied by
  `r_i`, so hard cluster pull only acts on reliable core cells.

`r_i = neighbor_agreement · local_density · membership_confidence` — three
orthogonal signals (see `reliability.py`).

## Files

- `model.py` — encoder + mask predictor + decoder + trainable DEC centers + soft membership
- `reliability.py` — the per-cell reliability field
- `loss.py` — gated scMAE + gated NeighborMix pseudo + gated DEC KL
- `run.py` — training/eval loop, PCA-KNN graph, smoke/screen stages

## Usage (smoke on Melanoma_5K)

```
python methods/DeepLearning/GatedNeighborMix_scMAE/run.py \
  --data_path methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad \
  --save_dir methods/DeepLearning/GatedNeighborMix_scMAE/runs/smoke/Melanoma_5K/seed42 \
  --dataset_name Melanoma_5K --label_key resolved_label --n_clusters 9 --smoke --no_cuda
```

## Core testable hypothesis

Gated fusion should KEEP NeighborMix's Melanoma gain (+0.042 ARI vs scMAE)
AND remove its Macosko loss (−0.073 ARI), because mixing switches off exactly
where the neighbor graph is unreliable.
