# CAAM-scMAE implementation status

This package implements CAAM-scMAE according to `BDD/CAAM_scMAE_BDD_plan/`.

Formal benchmark should register only `caam_scmae` with `--variant full`. The internal variants `control`, `axial`, and `advmask` are for ablation only.

Implemented BDD-critical contracts:

- Training batches do not return labels.
- Matched donor corruption is label-free and forbids self-donor.
- Replacement values are gene-wise and come from the same gene in donor cells.
- Random and adversarial masks use fixed-budget logic subject to eligibility.
- Budget deficit rate is recorded and checked with fail-fast.
- The adversarial generator only outputs mask scores/masks, not replacement values.
- Generator training uses the straight-through mask path and checks real generator gradient norm.
- Student parameters are frozen during generator updates.
- Axial context indices are label-free, fixed, saved, and reused.
- Context attention applies self-exclusion.
- Embedding extraction writes rows back by original cell index.
- `metrics.json` includes both `kmeans_known_k` and `leiden_fixed`.

Remaining caveats before formal benchmark:

- Do not mark `caam_scmae` as `VERIFIED/PASS` until smoke and contract tests pass.
- Formal integration belongs to Phase 11 and must follow the benchmark addendum.
- Research-level choices such as dataset split, main metric, and whether to continue after weak ablations require researcher approval.
