# rank52_anatomy_prompt_replaced_token_full

Independent-full scMAE candidate based on rank52, "Anatomically-guided Masked Autoencoder with Domain-Adaptive Prompting for multimodal cerebral aneurysm detection and segmentation".

## Theory basis

The reference report recommends mapping this paper's anatomy-guided MAE and domain-adaptive prompting ideas to gene expression through tokenization, replaced-expression detection, self-guided/curriculum masking, and gene-specific rank tokens. The index does not provide a GitHub URL, so this implementation is reconstructed from the PDF and report.

## scMAE integration

This candidate fills the mask/target gap in scMAE:

- Gene anatomy is represented by unsupervised modules from log-expression mean, variance, and zero fraction.
- `AnatomyPromptScMAE` injects module summaries and learned module prompts into the encoder embedding.
- The original scMAE mask prediction and masked expression reconstruction are retained.
- Two auxiliary targets are added: corrupted-expression detection and gene-specific rank-token classification.

## Expression semantics

`scaled_expr` is used only as encoder input when `--scale_input true`. Reconstruction targets, rank-token quantiles, gene modules, and difficulty statistics are derived from unscaled log-expression. No scaled expression is treated as count data.

## NeighborMix

NeighborMix is not used. The relationship is independent and complementary, with `mixed_cell_fraction=0.0`.

## Outputs

Each run writes `embedding_final.npy`, `labels.npy`, `training_history.json`, `diagnostics.json`, `summary.json`, `args.json`, and optional fixed-k metrics. Smoke/screen rows are written only to `新模型独立快筛单次结果.csv` and `新模型独立快筛汇总结果.csv`.
