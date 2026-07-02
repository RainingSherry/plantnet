# rank19_data2vec_latent_teacher_full

Independent-full scMAE candidate based on data2vec.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 19 section.
- `02_整理索引.csv`, rank 19 row.
- `019_高_data2vec_A_General_Framework_for_Self-supervised_Learning_in_Speech,_Vision_and_Language.pdf`.
- GitHub check: `https://github.com/facebookresearch/fairseq/tree/main/examples/data2vec` returned HTTP 200.

Adaptation:
- Retains scMAE mask prediction and masked expression reconstruction.
- Uses a student network on a strongly masked expression view.
- Uses an EMA teacher on the clean/weak view to provide normalized continuous latent targets.
- Adds a small variance regularizer and reports collapse diagnostics.

NeighborMix relation:
- Independent and potentially complementary.
- This candidate does not mix cells; `mixed_cell_fraction=0.0`.

Smoke/screen results are candidate evidence only and must not be appended to `全benchmark结果.csv`.

