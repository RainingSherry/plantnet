# Protocol-analysis manuscript workspace

This directory holds the current evidence-driven CAAM/scMAE manuscript route.

Current route:

```text
protocol_analysis / diagnostic paper
```

Current guardrail:

```text
Do not present development triage, attention smoke, or known-K ARI as validation or sealed-test evidence.
```

BDD route and validation gate:

```text
BDD/CAAM_scMAE_BDD_plan/CAAM_PUBLICATION_DECISION.md
BDD/CAAM_scMAE_BDD_plan/PROTOCOL_ANALYSIS_VALIDATION_PLAN.md
paper/论文工作台/Reviewer_Risk_and_Rebuttal_Matrix_2026-06-26.md
paper/论文工作台/Protocol_Analysis_Literature_Matrix_2026-06-26.md
paper/protocol_analysis_manuscript/REPRODUCIBILITY_CHECKLIST.md
```

Validation is not approved yet. The validation plan freezes what may be tested later and forbids using validation data to choose corruption, mask policy, architecture, loss, clustering resolution, or claims.

The reviewer risk matrix lists the likely objections and the evidence or experiments needed before submission.

The protocol-analysis literature matrix maps each relevant paper to the manuscript claim it supports or weakens, the experiment design it suggests, and the safe gap sentence to use in writing.

`main.tex` now uses this matrix to frame Related Work around protocol sensitivity, leakage risks, sparse masked objectives, foundation/context motivation, and learned-mask diagnostics. It should not be edited back into an AdvMask/Axial positive-method narrative unless new validated evidence supports that route.

The Results prose now summarizes the generated Phase 13, Phase 14, attention-smoke, fixed-Leiden, post-hoc label-diagnostic, resource-smoke, feature-space-smoke, mask-ratio-smoke, and mechanism-decision evidence. It intentionally frames AdvMask and the current Axial encoder as downgraded mechanisms rather than positive contributions.

`main.tex` also includes a validation-freeze section. This is a protocol boundary, not execution approval: validation and sealed test have not been run, and validation must not be used to tune corruption type, mask policy, architecture, loss, clustering resolution, or claims.

The reproducibility checklist records generated artifact counts, source artifact roots, regeneration commands, metric labels, label-leakage controls, resource-summary and mechanism-decision boundaries, and remaining submission-readiness gaps. It is a manuscript support file, not validation evidence.

Regenerate local tables and figures from existing artifacts:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_protocol_figures.py
```

Regenerate post-hoc label diagnostics from the same existing artifacts:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_posthoc_label_diagnostics.py
```

Regenerate trainable-parameter and resource-manifest summaries from the same existing artifacts:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_resource_summary.py
```

Regenerate the development mechanism-decision matrix from the same existing artifacts:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_mechanism_decision_tables.py
```

Regenerate the feature-space development smoke summary from existing `/tmp/caam_feature_space_smoke/dev_20260626_gpu` artifacts:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_feature_space_smoke_summary.py
```

Regenerate the mask-ratio development smoke summary from existing `/tmp/caam_mask_ratio_smoke/dev_20260626_gpu` artifacts:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_mask_ratio_smoke_summary.py
```

Run or resummarize the instrumented development resource smoke:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
OMP_NUM_THREADS=1 \
MKL_NUM_THREADS=1 \
OPENBLAS_NUM_THREADS=1 \
NUMBA_NUM_THREADS=1 \
python paper/protocol_analysis_manuscript/scripts/run_instrumented_resource_smoke.py --gpu 1
```

Use `--summarize-only` to regenerate the manuscript table from existing `/tmp/caam_resource_smoke/dev_20260626` artifacts without rerunning training.

Run the manuscript claim-boundary audit:

```bash
python paper/protocol_analysis_manuscript/scripts/audit_claim_boundaries.py
```

Generated outputs:

```text
generated/data/
generated/tables/
generated/figures/
generated/artifact_manifest.json
generated/posthoc_label_diagnostics/
generated/resource_summary/
generated/instrumented_resource_smoke/
generated/feature_space_smoke/
generated/mask_ratio_smoke/
generated/mechanism_decisions/
generated/claim_audit/
```

The generated figure/table set includes both known-K K-means summaries and fixed-Leiden non-oracle summaries. Known-K metrics remain development diagnostics and must not be presented as fully unsupervised clustering evidence.

The post-hoc label diagnostics compute dominant-cluster recall, cluster purity, label F1, rare-label summaries, worst-label recovery, and label-to-cluster flow tables from existing `labels.npy` and clustering prediction artifacts. Labels are used only after training for diagnosis and must not guide corruption, masks, context construction, model selection, early stopping, validation, or manuscript claims.

The generated flow heatmaps are row-normalized label-to-cluster diagnostics for seed 42. Cluster IDs are arbitrary within each run, so the plots diagnose fragmentation and concentration only; they are not cell-type maps and should not be used as marker-gene or biological identity claims.

The resource summary reads existing `run_manifest.json`, `runtime.json`, and `metrics.json` files. It supports trainable-parameter and device-metadata claims for the full development artifact set. The instrumented resource smoke adds wall-clock time and GPU total-memory delta for three development datasets under seed 42 and three epochs; it is not a submission-scale runtime benchmark.

The feature-space smoke summary reads the paired Quake development runs comparing HVG 2000 with dense full-gene input under the same control protocol. It supports only the current dense-MLP route decision to keep HVG 2000 as the default; it does not validate HVG 2000 or reject sparse/gene-token full-gene models.

The mask-ratio smoke summary reads nine Quake development runs: mask ratios 0.2, 0.4, and 0.6 crossed with seeds 42, 2024, and 3407. It exports a table plus a sensitivity plot. It supports only a protocol-sensitivity statement; it must not be used to tune the frozen validation mask ratio.

The mechanism-decision matrix converts the current development evidence into an auditable route decision: carry scMAE-style shuffle forward, downgrade AdvMask, do not use the existing Axial encoder as a rescue path, and treat parameter-matched MLP comparisons as diagnostics. It is not validation evidence.

The claim-boundary audit checks that the manuscript does not promote development-only negative mechanism evidence into positive method, validation, synergy, or unknown-\(K\) claims. It also checks generated artifact counts, resource-smoke scope text, feature-space smoke boundaries, and mask-ratio smoke boundaries against current manifests.

The generator reads:

```text
results/CAAM_scMAE_correction/corruption_triad/formal
results/CAAM_scMAE_correction/advmask_triage/formal
/tmp/caam_attention_context_smoke/dev_20260626
/tmp/caam_resource_smoke/dev_20260626
/tmp/caam_feature_space_smoke/dev_20260626_gpu
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu
```

The `/tmp` attention, resource-smoke, feature-space-smoke, and mask-ratio-smoke paths are intentionally treated as development-only evidence. If any disappears, rerun the smoke or replace it with a reviewed, explicitly scoped artifact path before relying on the generated figures or tables.
