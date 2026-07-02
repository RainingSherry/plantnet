# rank22_ijepa_context_target_full

Independent-full scMAE candidate adapted from **I-JEPA: Self-Supervised Learning from Images with Joint-Embedding Predictive Architecture**.

## Theory basis

Read sources:

- `参考文献/00_scMAE改良方法整理总报告.md`, section 22.
- `参考文献/02_整理索引.csv`, rank 22 row.
- `参考文献/01_PDF论文_按推荐程度排序/022_高_I-JEPA_Self-Supervised_Learning_from_Images_with_Joint-Embedding_Predictive_Architecture.pdf`.
- GitHub: `https://github.com/facebookresearch/ijepa`, README, `src/masks/multiblock.py`, and `src/models/vision_transformer.py`.

I-JEPA predicts target block representations from context block representations. The teacher/target encoder is updated by EMA, and target blocks are taken from the output of the target encoder rather than from corrupted input. This scRNA adaptation uses coexpression-ordered gene-module patches instead of image patches.

## scMAE gap addressed

This candidate addresses the `teacher` and `semantic target` gaps while preserving scMAE:

- mask prediction;
- masked expression reconstruction;
- EMA target encoder patch representations;
- context-to-target JEPA predictor.

`scaled_expr` is only used as encoder input when `--scale_input true`; `log_expr` is used for masked expression reconstruction. No NB/ZINB count likelihood is used.

## NeighborMix relation

NeighborMix is not used. The relation is independent and potentially complementary. Diagnostics set `mixed_cell_fraction=0.0`.

## Smoke

```bash
CUDA_VISIBLE_DEVICES=1 python /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/rank22_ijepa_context_target_full/run.py --data_path /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad --save_dir /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/runs/scMAEs_independent_smoke/rank22_ijepa_context_target_full/Melanoma_5K/seed42 --dataset_name Melanoma_5K --label_key resolved_label --input_mode auto --n_top_genes 1000 --target_sum 10000 --scale_input true --n_clusters 9 --epochs 3 --batch_size 256 --lr 0.001 --weight_decay 0.00001 --seed 42 --gpu 1 --no_save_h5ad --smoke
```
