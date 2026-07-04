# Reviewer Risk and Rebuttal Matrix

Date: 2026-06-26

Status: manuscript strategy memo. This is not experimental evidence.

Current route:

```text
protocol_analysis / diagnostic paper
```

Current hard boundary:

```text
No validation or sealed test has been run.
Do not present development triage as publication-level evidence.
Do not claim AdvMask, current Axial, or Axial + AdvMask synergy as validated contributions.
```

## 1. Reviewer Risk Overview

The current manuscript is no longer a positive CAAM method paper. Its strongest route is a protocol-analysis paper about masked autoencoding for sparse scRNA-seq clustering. This is more honest, but it changes the likely reviewer objections.

The main reviewer risk is not "the model is weak"; that is already acknowledged. The main risk is whether the negative and diagnostic evidence is broad, rigorous, and useful enough to justify a paper.

## 2. Risk Matrix

| Risk | Likely reviewer criticism | Current evidence | Current weakness | Needed response or experiment | Claim boundary |
|---|---|---|---|---|---|
| R1: Novelty | "This is just debugging scMAE, not a paper." | Phase12-16 show protocol corrections, corruption triad, AdvMask failure, Axial failure. | The manuscript must explain why these failures reveal a general benchmark/protocol issue. | Frame contribution around protocol sensitivity, effective corruption diagnostics, and staged gates; compare to scMAE/scNAME/scCluBench protocol assumptions. | Do not claim a new state-of-the-art model. |
| R2: Negative result value | "AdvMask and Axial failed; why publish?" | AdvMask has nonzero gradients but no seed-stable clustering gain; current Axial loses to parameter-matched MLP. | Negative results need broader interpretation and reproducible artifacts. | Show that active mask learning and better corruption diagnostics can diverge from clustering utility; include tables/figures and artifact manifest. | Present as diagnostic evidence, not method superiority. |
| R3: Development-only evidence | "Three development datasets are too few." | Phase13 uses 3 datasets x 3 corruptions x 3 seeds; Phase14 uses 3 datasets x 2 variants x 3 seeds. | Validation is not run. | Run validation only after explicit approval under frozen protocol; report if qualitative conclusions hold. | Current manuscript is not validation-ready until validation plan is executed. |
| R4: Known-K dependence | "Known-K K-means ARI is oracle-like." | Reports also include fixed Leiden metrics; generated non-oracle tables/figures now exist. | Current primary development gate still uses known-K ARI. | Make known-K labeling explicit; discuss fixed-Leiden/non-oracle tables and whether conclusions agree. | Never call known-K fully unsupervised. |
| R5: Biological relevance | "Metrics do not prove biological utility." | Labels are post hoc; no marker/rare-cell analysis yet. | Biological interpretation missing. | Add marker-overlap, rare-cell recovery, cluster-label Sankey, per-label recall, and whether diagnostics predict rare-label degradation. | Do not make biological discovery claims yet. |
| R6: Baseline fairness | "Compared variants are not fair or parameter-matched." | Attention smoke includes parameter-matched MLP; Phase14 control/advmask share MLP. | Broader baseline set not yet integrated into this protocol-analysis route. | Report parameter counts and keep variants internal; compare final frozen protocol to scMAE/scNAME/Leiden/Louvain only after validation plan approval. | Do not compare internal variants as formal methods. |
| R7: Protocol overfitting | "You chose scmae_shuffle because it looked best on development." | Phase13 explicitly chooses scmae_shuffle from development. | This can be seen as tuning if reused for validation claims. | Treat scmae_shuffle as development-selected and freeze before validation; validation tests generalization, not selection. | Do not use validation to change corruption_type. |
| R8: Effective corruption diagnostic | "Effective corruption rate is obvious and not enough." | nonzero-aware donor improves effective corruption but not stable clustering. | Need clearer diagnostic interpretation. | Show scatter/table: effective corruption vs ARI; argue diagnostic is necessary but insufficient. | Do not claim effective corruption predicts quality alone. |
| R9: AdvMask implementation doubt | "Maybe AdvMask failed because implementation is wrong." | Warmup bug was found and fixed; generator_grad_norm positive on all 3 datasets; mask entropy/gini non-collapsed. | Still only 3-epoch triage. | Present the warmup correction transparently; include gradient/mask diagnostics; avoid claiming all learned masking fails. | Claim only current constrained AdvMask did not pass gate. |
| R10: Axial implementation doubt | "Maybe Axial failed because hyperparameters were bad." | Current Axial underperforms standard and parameter-matched MLP in seed-42 smoke. | Single-seed smoke only. | Treat current Axial as downgraded, not universally refuted; require new mechanism design and fresh gate if revisited. | Do not claim all attention fails. |
| R11: Reproducibility | "Can others reproduce the pipeline?" | Scripts and artifact manifests exist; generated tables/figures are reproducible from artifacts. | Paper assets currently include `/tmp` attention smoke path. | Move or regenerate attention smoke into a stable reviewed artifact path before submission. | Do not rely on ephemeral `/tmp` paths in final reproducibility package. |
| R12: Journal fit | "This is too engineering-focused for top-tier venue." | Protocol-analysis route can be valuable if tied to benchmark pitfalls and single-cell sparsity. | Need stronger literature positioning and validation breadth. | Strengthen related work around scCluBench, scMAE, scNAME, masked modeling, benchmark leakage/fairness, and negative result value. | Target high-quality bioinformatics/computational biology first, not top-tier claim yet. |

## 3. Rebuttal-ready Answers

### Q1. Why is this publishable if the proposed mechanisms failed?

Answer:

```text
The paper is not framed as a failed method paper. It studies why plausible masked-autoencoder mechanisms can look active while failing to improve clustering representations. The contribution is a corrected protocol and diagnostic analysis showing that feature-space choices, corruption semantics, and active mask learning must be evaluated separately before method claims are made.
```

Evidence needed:

```text
Phase12 protocol correction.
Phase13 corruption triad.
Phase14 AdvMask gradient-positive but effect-size-negative result.
Attention smoke showing parameter-matched MLP beats current Axial.
Generated figures/tables from protocol_analysis_manuscript/generated.
```

### Q2. Why not just tune AdvMask or Axial more?

Answer:

```text
Tuning after seeing development failures would make the validation protocol untrustworthy. The project now freezes the protocol-analysis route. New attention or masking designs would need a separate BDD branch, fresh gate, and no reuse of validation/sealed datasets for selection.
```

Evidence needed:

```text
PROTOCOL_ANALYSIS_VALIDATION_PLAN.md
CAAM_PUBLICATION_DECISION.md
risk_and_stop_criteria.md
```

### Q3. Are conclusions dependent on known-K K-means?

Answer:

```text
Known-K K-means is explicitly labeled as a development metric. The manuscript must report fixed-resolution Leiden metrics as non-oracle support and must not present known-K metrics as fully unsupervised clustering.
```

Needed addition:

```text
Fixed-Leiden generated tables/figures now exist, but the manuscript still needs prose explaining where they support or complicate the known-K trends.
```

### Q4. Do the diagnostics matter biologically?

Answer:

```text
Not established yet. Biological interpretation is a required next step. Labels and markers will be used only post hoc to test whether protocol choices affect rare-cell recovery, marker overlap, and cluster-label alignment.
```

Needed addition:

```text
Development-only biological interpretation script and validation-ready frozen biological analysis plan.
```

### Q5. How do you prevent validation leakage?

Answer:

```text
Validation is not approved or run yet. The validation plan freezes feature space, corruption choices, mask policy, architecture, metrics, and claim boundaries before any validation data are used. Validation failure must be recorded rather than tuned away.
```

Evidence needed:

```text
PROTOCOL_ANALYSIS_VALIDATION_PLAN.md
experiment_split_decision.md
```

## 4. Minimum Experiments Before Submission

Required before a serious submission:

```text
1. Validation under the frozen protocol, after explicit approval.
2. Fixed-Leiden/non-oracle metric discussion alongside known-K metrics.
3. Biological interpretation analysis: marker overlap, rare-cell recovery, Sankey/contingency.
4. Runtime and parameter-count table.
5. Stable artifact paths for attention smoke or a reviewed replacement artifact.
6. Literature-backed related work matrix for masked corruption, scRNA-seq clustering benchmarks, and negative diagnostic analyses.
```

Not required unless a new method route is reopened:

```text
1. New AdvMask tuning.
2. New Axial tuning.
3. Full CAAM 2x2 synergy factorial.
4. Formal benchmark registration of internal variants.
```

## 5. Claims To Use Carefully

Safe wording:

```text
Development evidence suggests that active mask learning and higher effective corruption can fail to improve clustering representations under a corrected scMAE-style protocol.
```

```text
The current Axial implementation did not pass a parameter-matched development smoke, so attention/context modeling remains a future design question rather than a supported contribution.
```

```text
Known-K metrics are used as development diagnostics; fixed Leiden and biological analyses are required before broader unsupervised-clustering claims.
```

Unsafe wording:

```text
CAAM-scMAE improves clustering.
AdvMask improves clustering.
Axial improves clustering.
Attention is useless for scRNA-seq.
Effective corruption predicts clustering quality.
Development evidence validates the method.
Known-K K-means proves unsupervised performance.
```

## 6. Immediate Next Writing Tasks

1. Add a `Reviewer Risk` subsection to the Discussion or Supplement.
2. Add prose interpretation for fixed-Leiden summaries in Results and Discussion.
3. Create a biological interpretation development plan/script that uses labels only post hoc.
4. Strengthen related work around protocol sensitivity and benchmark fairness.
5. Move final paper wording away from "CAAM-scMAE" unless referring to the historical project name.
