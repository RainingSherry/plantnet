# rank01_scdiva_absorbing_diffusion_full

Independent-full scMAE candidate inspired by ScDiVa, the highest-ranked item in
`参考文献/02_整理索引.csv` and the first method described in
`00_scMAE改良方法整理总报告.md`.

## Theory Basis

ScDiVa argues that single-cell dropout is better matched by masked discrete
diffusion than by Gaussian noise or autoregressive gene ordering. This candidate
keeps the original scMAE training spine and adds the parts that are practical for
local benchmark screening:

- absorbing mask corruption with a sampled diffusion time `t`;
- masked expression reconstruction on log-normalized expression;
- absorbing mask prediction;
- gene-specific expression quantile token prediction on masked positions.

This is not a full ScDiVa reproduction: it does not use the original large
pre-training corpus, 12-layer Transformer, RoPE, or foundation-model transfer
protocol. The indexed GitHub URL returned 404 / anonymous access unavailable, so
the implementation is reconstructed from the local paper PDF and the project
summary report.

## scMAE Gap

This model addresses the `mask/target` gap. Original scMAE replaces selected
entries with values from another cell and reconstructs expression. Here, selected
entries transition to an absorbing state, and the model must recover both
continuous expression magnitude and a discrete gene-specific token target.

## NeighborMix Relation

Independent and complementary. This candidate does not mix cells and does not use
NeighborMix. Diagnostic `mixed_cell_fraction` is always `0.0`.

## Example

```bash
CUDA_VISIBLE_DEVICES=1 python methods/DeepLearning/scMAEs/rank01_scdiva_absorbing_diffusion_full/run.py \
  --data_path methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad \
  --save_dir results/scMAEs_independent_smoke/rank01_scdiva_absorbing_diffusion_full/Melanoma_5K/seed42 \
  --dataset_name Melanoma_5K \
  --label_key resolved_label \
  --input_mode auto \
  --n_top_genes 1000 \
  --target_sum 10000 \
  --scale_input true \
  --n_clusters 9 \
  --epochs 3 \
  --batch_size 256 \
  --lr 0.001 \
  --weight_decay 0.00001 \
  --seed 42 \
  --gpu 1 \
  --no_save_h5ad \
  --smoke
```

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `training_history.json`,
`diagnostics.json`, `summary.json`, and `args.json`. If evaluation is enabled it
also writes `eval_fixed.csv`, `eval_metrics.json`, and `metrics.json`.

Smoke and screen rows are written only to:

- `methods/DeepLearning/scMAEs/新模型独立快筛单次结果.csv`
- `methods/DeepLearning/scMAEs/新模型独立快筛汇总结果.csv`

