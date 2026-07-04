# Reproducibility and Artifact Checklist

Date: 2026-06-26

Status: manuscript support checklist. This is not validation evidence.

## Evidence Scope

Current evidence type:

```text
development-only protocol-analysis evidence
```

Not included:

```text
validation datasets
sealed test datasets
formal benchmark registration of internal variants
publication-level method validation
```

Current generated artifact count:

```text
65 files under paper/protocol_analysis_manuscript/generated
```

## Source Artifact Roots

The manuscript generators read existing artifacts from:

```text
results/CAAM_scMAE_correction/corruption_triad/formal
results/CAAM_scMAE_correction/advmask_triage/formal
/tmp/caam_attention_context_smoke/dev_20260626
/tmp/caam_resource_smoke/dev_20260626
/tmp/caam_feature_space_smoke/dev_20260626_gpu
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu
```

Run counts recorded by generated manifests:

```text
Phase 13 corruption triad: 27 runs
Phase 14 AdvMask triage: 18 runs
Attention/context smoke: 9 runs
Instrumented resource smoke: 12 runs
Feature-space smoke: 2 runs
Mask-ratio smoke: 9 runs
```

The `/tmp` attention, resource-smoke, feature-space-smoke, and mask-ratio-smoke paths are development-only and ephemeral. Before submission, they must be regenerated into reviewed artifact paths or excluded from any reproducibility package that claims stable artifacts.

## Regeneration Commands

Protocol figures and summary tables:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_protocol_figures.py
```

Post-hoc label diagnostics:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_posthoc_label_diagnostics.py
```

Resource summaries:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_resource_summary.py
```

Mechanism-decision matrix:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_mechanism_decision_tables.py
```

Feature-space smoke summary:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_feature_space_smoke_summary.py
```

Mask-ratio smoke summary:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
python paper/protocol_analysis_manuscript/scripts/build_mask_ratio_smoke_summary.py
```

Instrumented resource smoke:

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

Resummarize existing instrumented resource-smoke artifacts without rerunning training:

```bash
PATH=/data/luolie/conda/envs/scclubench-main/bin:$PATH \
MPLCONFIGDIR=/tmp/matplotlib \
NUMBA_CACHE_DIR=/tmp/numba-cache \
OMP_NUM_THREADS=1 \
MKL_NUM_THREADS=1 \
OPENBLAS_NUM_THREADS=1 \
NUMBA_NUM_THREADS=1 \
python paper/protocol_analysis_manuscript/scripts/run_instrumented_resource_smoke.py --summarize-only --gpu 1
```

Claim-boundary audit:

```bash
python paper/protocol_analysis_manuscript/scripts/audit_claim_boundaries.py
```

Static checks used during manuscript drafting:

```bash
python -m py_compile paper/protocol_analysis_manuscript/scripts/build_protocol_figures.py
python -m py_compile paper/protocol_analysis_manuscript/scripts/build_posthoc_label_diagnostics.py
python -m py_compile paper/protocol_analysis_manuscript/scripts/build_resource_summary.py
python -m py_compile paper/protocol_analysis_manuscript/scripts/build_mechanism_decision_tables.py
python -m py_compile paper/protocol_analysis_manuscript/scripts/build_feature_space_smoke_summary.py
python -m py_compile paper/protocol_analysis_manuscript/scripts/build_mask_ratio_smoke_summary.py
python -m py_compile paper/protocol_analysis_manuscript/scripts/run_instrumented_resource_smoke.py
python -m py_compile paper/protocol_analysis_manuscript/scripts/audit_claim_boundaries.py
```

## Generated Output Families

Protocol-analysis summaries:

```text
generated/data/
generated/tables/
generated/figures/
generated/artifact_manifest.json
```

Post-hoc label diagnostics:

```text
generated/posthoc_label_diagnostics/data/
generated/posthoc_label_diagnostics/tables/
generated/posthoc_label_diagnostics/figures/
generated/posthoc_label_diagnostics/artifact_manifest.json
```

Resource summaries:

```text
generated/resource_summary/data/resource_runs.csv
generated/resource_summary/data/resource_summary.csv
generated/resource_summary/data/resource_summary_by_dataset.csv
generated/resource_summary/tables/resource_summary.tex
generated/resource_summary/artifact_manifest.json
```

Instrumented resource smoke:

```text
generated/instrumented_resource_smoke/instrumented_resource_smoke.csv
generated/instrumented_resource_smoke/instrumented_resource_smoke.tex
generated/instrumented_resource_smoke/instrumented_resource_smoke_summary.tex
generated/instrumented_resource_smoke/artifact_manifest.json
```

Feature-space smoke:

```text
generated/feature_space_smoke/feature_space_smoke.csv
generated/feature_space_smoke/feature_space_smoke.json
generated/feature_space_smoke/feature_space_smoke.md
generated/feature_space_smoke/feature_space_smoke.tex
```

Mask-ratio smoke:

```text
generated/mask_ratio_smoke/mask_ratio_smoke.csv
generated/mask_ratio_smoke/mask_ratio_smoke.json
generated/mask_ratio_smoke/mask_ratio_smoke.md
generated/mask_ratio_smoke/mask_ratio_smoke.tex
generated/mask_ratio_smoke/mask_ratio_smoke_sensitivity.png
```

Mechanism-decision matrix:

```text
generated/mechanism_decisions/mechanism_decisions.csv
generated/mechanism_decisions/mechanism_decisions.json
generated/mechanism_decisions/mechanism_decisions.md
generated/mechanism_decisions/mechanism_decisions.tex
```

Claim-boundary audit:

```text
generated/claim_audit/claim_boundary_audit.md
generated/claim_audit/claim_boundary_audit.json
```

Main metric labels:

```text
kmeans_known_k = known-K development diagnostic
leiden_fixed = fixed-resolution non-oracle diagnostic
```

Known-K metrics must not be described as fully unsupervised unknown-K evidence.

## Leakage Controls

Labels are not allowed in:

```text
training
corruption
mask selection
donor selection
gene module construction
context selection
model selection
early stopping
clustering resolution selection
validation tuning
```

Labels are allowed only post hoc for:

```text
numeric label recovery diagnostics
dominant-cluster recall
cluster purity
label F1
rare-label summaries
label-to-cluster flow tables
row-normalized heatmaps
```

Current label diagnostics do not support marker-gene, named cell-type, or biological discovery claims.

## Resource Summary Boundary

The resource summary reads existing development artifacts only:

```text
run_manifest.json
runtime.json
metrics.json
```

Supported claims from the full resource summary:

```text
trainable parameter counts
generator parameter counts
device metadata
known-K and fixed-Leiden ARI paired with those resource fields
```

Supported claims from the development instrumented smoke:

```text
wall-clock seconds for the three development datasets, seed 42, 3 epochs
GPU total-memory delta from pre-run baseline to sampled peak
CPU RSS process-tree peak from psutil sampling
```

Not supported by current resource artifacts:

```text
submission-scale wall-clock runtime
submission-scale GPU memory
submission-scale CPU memory
throughput
energy or cost
```

Those missing submission-scale fields must not be inferred from logs or hardware metadata. They require future explicitly instrumented runs under the frozen protocol.

## Feature-Space Smoke Boundary

The feature-space smoke reads two existing development runs only:

```text
/tmp/caam_feature_space_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__hvg2000__control__seed42__epochs3
/tmp/caam_feature_space_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__full_gene_stress__control__seed42__epochs3
```

Supported claim:

```text
Under the current dense-MLP control protocol, a single Quake development smoke favors HVG 2000 over dense full-gene input on both parameter cost and clustering diagnostics.
```

Not supported:

```text
HVG 2000 is validated.
Full-gene input is always worse.
Sparse full-gene models or gene-token models will fail.
This result is publication-level validation.
```

## Mask-Ratio Smoke Boundary

The mask-ratio smoke reads nine existing development runs only:

```text
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p2__control__seed42__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p2__control__seed2024__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p2__control__seed3407__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p4__control__seed42__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p4__control__seed2024__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p4__control__seed3407__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p6__control__seed42__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p6__control__seed2024__epochs3
/tmp/caam_mask_ratio_smoke/dev_20260626_gpu/Quake_Smart-seq2_Lung__mask0p6__control__seed3407__epochs3
```

Supported claim:

```text
A single-dataset nine-run Quake development smoke shows that nominal mask ratio changes global effective perturbation and downstream metrics under the current corrected control protocol, with metric-dependent preferred ratios.
```

Not supported:

```text
mask_ratio 0.2 is generally best.
mask_ratio 0.6 is generally best.
mask_ratio 0.4 is validated.
mask_ratio should be changed before validation.
This result is publication-level validation.
```

## Validation and Sealed-Test Boundary

Validation has not been run.

Validation datasets reserved by the freeze plan:

```text
data/processed_scmae/Young.h5ad
data/processed_scmae/Baron.h5ad
data/SRP182008.h5ad
```

Validation must not be used to choose:

```text
corruption_type
n_top_genes
mask_ratio
architecture
loss weights
early stopping
clustering resolution
whether to reintroduce AdvMask or Axial
main manuscript claim
```

Sealed test is currently disallowed. It can be considered only after one approved validation pass under a frozen protocol and a separate explicit approval.

## Submission-Readiness Gaps

Before a serious submission, the manuscript still needs:

```text
1. approved one-pass validation under the frozen protocol;
2. stable reviewed artifact path for attention/context smoke or replacement evidence;
3. marker/rare-cell interpretation plan executed only post hoc;
4. broader wall-clock and memory instrumentation under the frozen protocol;
5. final claim audit against known-K and validation boundaries;
6. complete manuscript compile in an environment with LaTeX installed.
```

## Safe Claim Boundary

Safe:

```text
Development evidence suggests that effective corruption diagnostics and active learned masking can diverge from clustering representation quality under a corrected scMAE-style protocol.
```

Unsafe:

```text
CAAM-scMAE improves clustering.
AdvMask improves clustering.
The current Axial encoder improves clustering.
Axial plus AdvMask has synergy.
Development evidence validates a publication-ready method.
Known-K ARI proves fully unsupervised clustering performance.
```
