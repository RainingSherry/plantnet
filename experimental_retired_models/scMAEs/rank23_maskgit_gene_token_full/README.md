# rank23_maskgit_gene_token_full

Independent-full scMAE candidate adapted from **MaskGIT: Masked Generative Image Transformer**.

## Theory basis

Read sources:

- `参考文献/00_scMAE改良方法整理总报告.md`, rank 23 section.
- `参考文献/02_整理索引.csv`, rank 23 row.
- `参考文献/01_PDF论文_按推荐程度排序/023_高_MaskGIT_Masked_Generative_Image_Transformer.pdf`.
- GitHub: `https://github.com/google-research/maskgit`, README, repository tree, `maskgit/libml/parallel_decode.py`, and `maskgit/nets/maskgit_transformer.py`.

MaskGIT trains a bidirectional Transformer to predict masked discrete tokens and uses a cosine mask schedule plus confidence-based iterative refinement. This scRNA adaptation does not use image VQ tokens. It derives gene-specific rank tokens from each gene's `log_expr` distribution.

## scMAE gap addressed

This candidate addresses the `mask` and `semantic target` gaps while retaining scMAE:

- masked expression reconstruction;
- mask prediction;
- gene-specific token prediction on masked genes;
- replaced-expression detection;
- cosine curriculum mask ratio sampled per batch.

`scaled_expr` is only used as encoder input when `--scale_input true`; token and reconstruction targets use `log_expr`.

## NeighborMix relation

NeighborMix is not used. The relation is independent and potentially complementary. Diagnostics set `mixed_cell_fraction=0.0`.

## Smoke

```bash
CUDA_VISIBLE_DEVICES=1 python /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/rank23_maskgit_gene_token_full/run.py --data_path /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/benchmark_data/Melanoma_5K.h5ad --save_dir /home/luolie/biopipeline/dimension-reduction/plantnet/methods/DeepLearning/scMAEs/runs/scMAEs_independent_smoke/rank23_maskgit_gene_token_full/Melanoma_5K/seed42 --dataset_name Melanoma_5K --label_key resolved_label --input_mode auto --n_top_genes 1000 --target_sum 10000 --scale_input true --n_clusters 9 --epochs 3 --batch_size 256 --lr 0.001 --weight_decay 0.00001 --seed 42 --gpu 1 --no_save_h5ad --smoke
```
