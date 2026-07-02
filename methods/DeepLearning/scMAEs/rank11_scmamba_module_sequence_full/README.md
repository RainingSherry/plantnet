# rank11_scmamba_module_sequence_full

Independent-full scMAE candidate inspired by scMamba.

Read sources:
- `00_scMAE改良方法整理总报告.md`, rank 11 scMamba section.
- `02_整理索引.csv`, scMamba row.
- `011_高_scMamba_Scalable_Foundation_Model_for_Single-Cell_Multi-Omics_Integration.pdf`.
- GitHub check: `https://github.com/23AIBox/scMamba` returned HTTP 200.

Theoretical basis:
- scMamba uses patch-based cell tokenization, sequence/state-space style encoders, and cell-level representation learning.
- The report warns that genes do not have a natural order in this benchmark, so random HVG order is not acceptable.

Adaptation to scMAE:
- Genes are ordered by unsupervised co-expression anchor assignment from log-expression.
- Reordered genes are grouped into fixed patches that behave as module tokens.
- A gated depthwise sequence mixer approximates the state-space/mamba dependency path without adding a new dependency.
- The model keeps mask prediction and masked expression reconstruction, and adds a module-level reconstruction target.

NeighborMix relation:
- Independent and potentially complementary.
- This candidate does not mix cells; `mixed_cell_fraction=0.0`.

Smoke/screen outputs are candidate evidence only and must not be appended to `全benchmark结果.csv`.

