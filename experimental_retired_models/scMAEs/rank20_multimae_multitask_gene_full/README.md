# rank20_multimae_multitask_gene_full

Independent-full scMAE candidate inspired by MultiMAE.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 20 section.
- `02_整理索引.csv`, rank 20 row.
- `020_高_MultiMAE_Multi-modal_Multi-task_Masked_Autoencoders.pdf`.
- GitHub check: `https://github.com/EPFL-VILAB/MultiMAE` returned HTTP 200.

Adaptation:
- Retains scMAE mask prediction and masked expression reconstruction.
- Uses one shared encoder and three lightweight task decoders.
- Predicts masked log-expression, gene-specific log-expression quantile tokens, and masked module-level log-expression summaries.
- Uses scaled expression only for encoder input; token/module targets are computed from log-expression.

NeighborMix relation:
- Independent and potentially complementary.
- This candidate does not mix cells; `mixed_cell_fraction=0.0`.

Smoke/screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.

