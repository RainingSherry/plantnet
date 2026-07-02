# rank02_maskfeat_gene_module_semantic_full

Independent-full scMAE candidate inspired by MaskFeat.

## Theory Basis

MaskFeat trains a masked model to predict feature descriptors of masked regions
rather than raw pixels alone. For scRNA-seq, this candidate maps that idea to
gene-module features computed from log-normalized expression: module mean,
zero-rate, and rank-profile.

## scMAE Gap

This model addresses the semantic-target gap. It keeps scMAE mask prediction and
masked expression reconstruction, then adds a semantic head trained only on
modules touched by the mask.

## NeighborMix Relation

Independent and complementary. No cell mixing is used; `mixed_cell_fraction` is
always `0.0`.

## Example

```bash
CUDA_VISIBLE_DEVICES=1 python methods/DeepLearning/scMAEs/rank02_maskfeat_gene_module_semantic_full/run.py \
  --data_path methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad \
  --save_dir methods/DeepLearning/scMAEs/runs/scMAEs_independent_smoke/rank02_maskfeat_gene_module_semantic_full/Melanoma_5K/seed42 \
  --dataset_name Melanoma_5K \
  --label_key resolved_label \
  --n_clusters 9 \
  --epochs 3 \
  --seed 42 \
  --gpu 1 \
  --no_save_h5ad \
  --smoke
```

