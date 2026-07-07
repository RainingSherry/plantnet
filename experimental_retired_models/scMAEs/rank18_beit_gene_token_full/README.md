# rank18_beit_gene_token_full

Independent-full scMAE candidate inspired by BEiT.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 18 section.
- `02_整理索引.csv`, rank 18 row.
- `018_高_BEiT_BERT_Pre-Training_of_Image_Transformers.pdf`.
- GitHub check: `https://github.com/microsoft/unilm/tree/master/beit` returned HTTP 200.

Theoretical basis:
- BEiT masks input patches and predicts discrete visual tokens rather than raw pixels.
- For scRNA-seq, this candidate maps the target side to gene-specific log-expression quantile tokens.

Adaptation:
- Retains scMAE mask prediction and masked expression reconstruction.
- Adds masked gene-token prediction from log-expression quantile bins.
- Adds replaced-expression detection for masked positions that are swapped from another cell.
- Uses scaled expression only as encoder input; token targets are computed from unscaled log-expression.

NeighborMix relation:
- Independent and potentially complementary.
- This candidate does not mix cells; `mixed_cell_fraction=0.0`.

Smoke/screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.

