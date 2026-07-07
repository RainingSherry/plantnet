# rank24_flow_matching_latent_denoise_full

Independent-full scMAE candidate adapted from **Flow Matching for Generative Modeling**.

## Theory basis

Read sources:

- `参考文献/00_scMAE改良方法整理总报告.md`, rank 24 section.
- `参考文献/02_整理索引.csv`, rank 24 row.
- `参考文献/01_PDF论文_按推荐程度排序/024_高_Flow_Matching_for_Generative_Modeling.pdf`.
- GitHub: `https://github.com/atong01/conditional-flow-matching`, README and `torchcfm/conditional_flow_matching.py`.

Flow Matching regresses a neural vector field to a conditional probability path. This implementation uses the target conditional OT-style path only in scMAE latent space:

`z_t = t z_clean + (1 - (1 - sigma_min)t) noise`

and predicts the conditional vector field from `z_t` and `t`.

## scMAE gap addressed

This candidate addresses the `robust loss / target` gap while preserving:

- mask prediction;
- masked expression reconstruction.

The flow head is a lightweight auxiliary regularizer. Generated latent samples are not evaluated as real cells.

`scaled_expr` is only used as encoder input when `--scale_input true`; reconstruction targets use `log_expr`. No count likelihood is used.

## NeighborMix relation

NeighborMix is not used. The relation is independent and potentially complementary. Diagnostics set `mixed_cell_fraction=0.0`.

## Smoke

```bash
CUDA_VISIBLE_DEVICES=1 python /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/rank24_flow_matching_latent_denoise_full/run.py --data_path /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad --save_dir /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/runs/scMAEs_independent_smoke/rank24_flow_matching_latent_denoise_full/Melanoma_5K/seed42 --dataset_name Melanoma_5K --label_key resolved_label --input_mode auto --n_top_genes 1000 --target_sum 10000 --scale_input true --n_clusters 9 --epochs 3 --batch_size 256 --lr 0.001 --weight_decay 0.00001 --seed 42 --gpu 1 --no_save_h5ad --smoke
```
