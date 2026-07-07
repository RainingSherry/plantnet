# rank21_audiomae_module_patch_full

Independent-full scMAE candidate adapted from **AudioMAE: Masked Autoencoders that Listen**.

## Theory basis

Read sources:

- `参考文献/00_scMAE改良方法整理总报告.md`, section 21.
- `参考文献/02_整理索引.csv`, rank 21 row.
- `参考文献/01_PDF论文_按推荐程度排序/021_高_AudioMAE_Masked_Autoencoders_that_Listen.pdf`.
- GitHub: `https://github.com/facebookresearch/AudioMAE`, README and `models_mae.py`.

AudioMAE treats an audio spectrogram as locally correlated patches, masks a high proportion of patches, encodes only visible patches, restores the full patch order with learned mask tokens, and decodes masked patches. The report warns that genes have no natural order, so this implementation first builds a coexpression anchor order from `log_expr`, then applies module-patch masking on the ordered HVGs.

## scMAE gap addressed

This model targets the `mask` and `semantic target` gaps in vanilla scMAE. It keeps:

- mask prediction;
- masked expression reconstruction.

It adds:

- coexpression-ordered gene module patches;
- AudioMAE-style high-ratio patch masking with visible-only encoder;
- shifted local decoder blocks;
- module spectral targets: patch mean, patch standard deviation, and local adjacent-difference energy from `log_expr`.

`scaled_expr` is only used for the encoder input when `--scale_input true`; `log_expr` is used for reconstruction and module targets.

## NeighborMix relation

NeighborMix is not used in this candidate. The relation is independent and potentially complementary. Diagnostics therefore set `mixed_cell_fraction=0.0`; no forced cell mixing is performed.

## Commands

Smoke:

```bash
CUDA_VISIBLE_DEVICES=1 python /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/rank21_audiomae_module_patch_full/run.py --data_path /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad --save_dir /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/runs/scMAEs_independent_smoke/rank21_audiomae_module_patch_full/Melanoma_5K/seed42 --dataset_name Melanoma_5K --label_key resolved_label --input_mode auto --n_top_genes 1000 --target_sum 10000 --scale_input true --n_clusters 9 --epochs 3 --batch_size 256 --lr 0.001 --weight_decay 0.00001 --seed 42 --gpu 1 --no_save_h5ad --smoke
```

Screen commands use the same script on Melanoma_5K, Quake_10x_Spleen, and Macosko with seed 42 and epochs 80. Screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.
