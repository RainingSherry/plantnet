# rank10_celler_longtail_protection_full

Independent-full scMAE candidate based on Celler's long-tail ideas.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 10 Celler section.
- `02_整理索引.csv`, Celler row.
- `010_高_Celler_A_Genomic_Language_Model_for_Long-Tailed_Single-Cell_Annotation.pdf`.
- GitHub checks: `AI4science-ym/HiCeller` returned HTTP 200; the report URL `ckqqqq/PyCeller` returned HTTP 404.

Theoretical basis:
- Celler treats gene expression values as language-model-like tokens and uses masked non-zero expression prediction.
- Celler introduces Gaussian Inflation loss to increase attention to long-tail categories.
- Celler introduces Hard Data Mining for confusing high-logit alternatives.

Adaptation to scMAE:
- Retains scMAE mask prediction and masked expression reconstruction.
- Adds gene expression token prediction from log-expression quantile bins.
- Adds unsupervised Celler-style prototype regularization using model confidence rather than true labels.
- Uses density and entropy to estimate rare/boundary risk, boosting reconstruction for rare-risk cells and reducing prototype pressure on boundary cells.

NeighborMix relation:
- Independent and potentially complementary.
- This candidate does not mix cells. `mixed_cell_fraction` is always `0.0`.

This is a screen candidate only. Smoke/screen outputs must not be appended to `全benchmark结果.csv`.

