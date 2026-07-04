# Context-aware Adversarial Axial Masked Autoencoding for Single-cell RNA-seq Clustering

Working manuscript skeleton v0

Updated: 2026-06-26

## Manuscript Status

This file is the main entry point for the CAAM-scMAE manuscript. It links the current draft sections and records which claims are already supported by literature and which claims still require experiments.

Current section drafts:

```text
Introduction: Introduction_English_draft_v0.md
Related Work: Related_Work_English_draft_v0.md
Methods: Methods_English_draft_v0.md
Experiments: Experiments_English_draft_v0.md
Protocol: Experiments_Protocol_v0.md
Phase 14 AdvMask decision: Phase14_AdvMask_Triage_Decision_2026-06-26.md
Attention/context smoke decision: Attention_Context_Smoke_2026-06-26.md
Phase 16 BDD publication decision: ../../BDD/CAAM_scMAE_BDD_plan/CAAM_PUBLICATION_DECISION.md
Protocol-analysis validation freeze plan: ../../BDD/CAAM_scMAE_BDD_plan/PROTOCOL_ANALYSIS_VALIDATION_PLAN.md
Reviewer risk/rebuttal matrix: Reviewer_Risk_and_Rebuttal_Matrix_2026-06-26.md
Dataset audit: Dataset_Audit_SRP182008.md
Stage 1 plan: Stage1_scMAE_compatible_experiment_plan.md
Stage 1 diagnostics: Stage1_diagnostics_result_2026-06-26.md
Stage 1 training smoke: Stage1_training_smoke_result_2026-06-26.md
References: references_core.bib
```

## Current Research Status

Phase 14 development triage does not support AdvMask as a main contribution. AdvMask produced real nonzero generator gradients and did not show mask concentration or embedding collapse, but its mean ARI improvement over the random-mask control was only `0.002403`, below the predefined seed-variation threshold. The current recommendation is therefore:

```text
drop_or_downgrade_advmask
```

The manuscript should not currently claim adversarial masking improves clustering, and it should not claim Axial + AdvMask synergy. The strongest current route is a protocol-analysis / diagnostic paper about masked autoencoding protocol choices in sparse scRNA-seq clustering, unless a separately justified context/attention module later passes parameter-matched development tests.

The seed-42 attention/context smoke further weakens the current Axial route. Under random masking and `scmae_shuffle`, the existing Axial encoder underperformed both the standard MLP and a parameter-matched MLP on all three development datasets:

```text
mean ARI: control = 0.561970
mean ARI: axial = 0.191264
mean ARI: mlp_parammatched = 0.584531
```

This does not prove that all attention mechanisms are unsuitable for scRNA-seq. It does mean the current CAAM Axial implementation should also be downgraded, and any future TabPFN-like/context model should be treated as a new mechanism design rather than a continuation of the current full CAAM story.

The Phase 16 BDD decision is now:

```text
chosen_route = protocol_analysis
phase16_gate_result = pass_for_protocol_analysis
phase16_gate_result = fail_for_positive_method_paper
```

This means the current manuscript should be developed as a protocol-analysis / diagnostic paper. A positive method paper would require a new mechanism that passes a fresh development gate and later validation.

Validation is planned but not approved. The validation freeze plan forbids using validation datasets to choose corruption type, mask policy, architecture, loss, resolution, or final claims.

Reviewer risk memo is now available:

```text
Reviewer_Risk_and_Rebuttal_Matrix_2026-06-26.md
```

Its current high-risk points are novelty, negative-result value, development-only evidence, known-K dependence, biological interpretation, baseline fairness, and validation leakage.

## Title Candidates

Original method-paper title, now not supported by current evidence:

```text
Context-aware Adversarial Axial Masked Autoencoding for Single-cell RNA-seq Clustering
```

Current conservative working title:

```text
Protocol-aware Masked Autoencoding for Single-cell RNA-seq Clustering
```

Alternative analysis-paper title:

```text
When Does Masked Corruption Help Single-cell Clustering?
```

Historical fallback, now also unsupported by the current Axial smoke:

```text
Context-aware Axial Masked Autoencoding for Single-cell RNA-seq Clustering
```

Fallback if only protocol diagnostics remain supported:

```text
Masked Corruption Protocols Matter in Single-cell RNA-seq Clustering
```

## Abstract Draft

Single-cell RNA sequencing (scRNA-seq) enables transcriptomic profiling at single-cell resolution, but clustering scRNA-seq data remains challenging due to high dimensionality, sparsity, noise, and biological heterogeneity. Masked autoencoding has recently emerged as a promising self-supervised strategy for clustering-oriented representation learning. However, the empirical behavior of masked autoencoders can depend strongly on protocol choices that are often treated as implementation details, including feature-space selection, corruption semantics, effective corruption under sparsity, and mask-selection policy.

We therefore reframe the CAAM-scMAE project as a staged protocol analysis rather than a validated positive method. Under the corrected HVG-based scMAE-style protocol, scMAE-style gene-wise shuffle remains the strongest corruption choice among the tested alternatives. A constrained adversarial mask selector trains with nonzero generator gradients but fails to improve clustering beyond seed-level variation. A follow-up attention/context smoke shows that the existing Axial encoder underperforms both the standard MLP and a parameter-matched MLP under random masking. These negative results are treated as evidence for pruning unsupported mechanisms, not as publication-level validation.

We design a scCluBench-style evaluation protocol covering classical clustering, deep embedding, graph-based, contrastive, and masked autoencoder baselines. In addition to ACC, NMI, and ARI, the protocol evaluates rare cell recovery, marker-overlap annotation, embedding distinguishability, mask diagnostics, runtime, and memory usage. This framework is intended to test whether context-aware and informative masked reconstruction can produce robust, biologically meaningful single-cell clustering embeddings.

Current evidence constraints:

```text
[Completed development triage] AdvMask does not pass the Phase 14 effect-size gate.
[Completed development smoke] Existing Axial encoder does not pass a seed-42 parameter-matched smoke.
[Forbidden current claim] AdvMask improves clustering.
[Forbidden current claim] Current Axial improves clustering.
[Forbidden current claim] Axial + AdvMask synergy.
[Pending] Whether a new context/attention design, not the current Axial implementation, is worth a new BDD branch.
[To be inserted] Rare cell and marker-overlap analysis.
[To be inserted] Runtime and memory analysis.
```

## Keywords

```text
single-cell RNA-seq
cell clustering
masked autoencoder
self-supervised learning
attention
adversarial masking
representation learning
```

## 1. Introduction

Draft source:

```text
Introduction_English_draft_v0.md
```

Current logic:

```text
1. scRNA-seq clustering is central to cell-type discovery.
2. scRNA-seq data are high-dimensional, sparse, noisy, and heterogeneous.
3. scMAE shows masked reconstruction is useful for clustering.
4. Existing scMAE-style models do not explicitly model row/column context and use random masks.
5. TabPFN motivates two-dimensional table-aware attention.
6. CAAM-scMAE proposes bi-axial context encoding and constrained informative masking.
```

Claims supported by literature:

```text
scMAE uses masked reconstruction and mask prediction.
TabPFN uses two-dimensional table-aware attention.
scCluBench supports the need for standardized multi-dataset clustering evaluation.
Single-cell foundation models are useful but not necessarily task-optimized for clustering.
```

Claims requiring experiments:

```text
CAAM-scMAE improves clustering metrics.
Bi-axial context improves embedding distinguishability.
Full model remains scalable.
```

Claims contradicted or unsupported by current development evidence:

```text
Constrained adversarial mask selection improves clustering under the corrected Phase 14 protocol.
AdvMask should be retained as a main method module.
```

## 2. Related Work

Draft source:

```text
Related_Work_English_draft_v0.md
```

Subsections:

```text
2.1 Clustering and benchmarking in scRNA-seq analysis
2.2 Classical and graph-based clustering methods
2.3 Deep generative and autoencoder-based representation learning
2.4 Masked and contrastive learning for single-cell clustering
2.5 Foundation models and task-specific clustering objectives
2.6 Positioning of CAAM-scMAE
```

Positioning statement:

```text
CAAM-scMAE is positioned between masked autoencoding, graph-based clustering, and single-cell foundation modeling.
It avoids a fixed graph, avoids supervised labels and pseudo-labels during pretraining, extends scMAE with bi-axial context, and keeps the objective task-specific for clustering.
```

## 3. Methods

Draft source:

```text
Methods_English_draft_v0.md
```

Subsections:

```text
3.1 Problem formulation
3.2 scMAE-compatible masked reconstruction
3.3 Bi-axial context encoder
3.4 Constrained adversarial mask selector
3.5 Training schedule
3.6 Evaluation outputs
```

Core equations:

```text
X'_{ij}=X_{\pi_j(i),j}
\widetilde{X}=(1-M)\odot X+M\odot X'
L = L_rec + gamma L_mask
min_theta max_phi L(theta, phi), phi in legal mask family Phi
```

Important terminology:

```text
Use "constrained adversarial mask selector", not generic "GAN".
Use "TabPFN-inspired bi-axial context", not "using TabPFN".
Use "scMAE-compatible baseline", not "scMAE reproduction" unless implementation exactly matches original.
```

## 4. Experiments

Draft source:

```text
Experiments_English_draft_v0.md
```

Detailed protocol source:

```text
Experiments_Protocol_v0.md
```

Local development dataset audit:

```text
Dataset_Audit_SRP182008.md
```

Stage 1 baseline plan:

```text
Stage1_scMAE_compatible_experiment_plan.md
```

Stage 1 diagnostics result:

```text
Stage1_diagnostics_result_2026-06-26.md
```

Stage 1 training smoke result:

```text
Stage1_training_smoke_result_2026-06-26.md
```

Development dataset:

```text
SRP182008.h5ad
13,514 cells × 53,678 genes
sparsity: 97.60%
source: Zhang et al. 2019, Molecular Plant
primary evaluation labels: Celltype and Seurat_clusters
```

Initial diagnostics-only observation:

```text
On a 2,048-cell × 2,000-HVG subset with random 30% masking, 48.3% of masked entries remained zero-to-zero after gene-wise shuffling, and only 47.9% of masked entries changed value.
```

Planned stages:

```text
Stage 1: scMAE-compatible reproduction
Stage 2: encoder ablation
Stage 3: mask strategy ablation
Stage 4: formal benchmark
```

Required metrics:

```text
ACC
NMI
ARI
macro-F1 after best mapping
marker-overlap annotation
rare cell recovery
silhouette / ASW
embedding collapse score
mask diagnostics
runtime and memory
```

Required artifacts:

```text
config.json
metrics.json
embedding_final.npy
cluster_labels.npy
training_history.json
runtime.json
param_count.json
mask_diagnostics.json
```

## 5. Results

Status:

```text
Not yet written. Requires experiments.
```

Expected result subsections:

```text
5.1 Overall clustering performance
5.2 Ablation of encoder architecture
5.3 Ablation of mask strategy
5.4 Rare cell type recovery and biological annotation
5.5 Embedding distinguishability and collapse analysis
5.6 Mask diagnostics
5.7 Runtime and memory
```

Do not write positive result claims before evidence.

Current evidence available:

```text
Mask diagnostics support the motivation that nominal random mask ratio can overestimate effective corruption under high sparsity.
A two-epoch scMAE-compatible MLP smoke test confirms that the cached HVG training pipeline, reconstruction loss, mask prediction loss, embedding export, and metric export run end to end.
No clustering improvement evidence yet.
```

## 6. Discussion

Status:

```text
Not yet written. Depends on results.
```

Planned discussion points:

```text
Why context-aware masked reconstruction helps or does not help.
When learned mask selection is beneficial.
Relationship to graph-based clustering.
Relationship to foundation models.
Failure modes under sparsity or large gene dimension.
Biological interpretation and rare cell discovery.
```

## 7. Limitations

Draft points:

```text
1. CAAM-scMAE may require careful control of mask selector stability.
2. Bi-axial context may increase memory compared with MLP encoders.
3. Context cell selection can affect results and must remain label-free.
4. The model is designed for clustering, not universal single-cell representation learning.
5. Full benchmark conclusions require diverse datasets and multiple seeds.
```

## 8. Reproducibility Checklist

Every reported experiment should include:

```text
dataset name
preprocessing details
number of cells and genes
sparsity
number of labels used only for evaluation
random seed
model variant
parameter count
runtime
peak memory
all metrics
saved embeddings
saved clustering labels
mask diagnostics
```

## 9. Reference File

Use:

```text
references_core.bib
```

Current BibTeX coverage:

```text
scMAE
TabPFN
Geneformer
scGPT
scFoundation
CellFM
scCluBench
SC3 / Louvain / Leiden / Seurat
DEC / DESC / scDeepCluster / scDCC / scziDesk
scNAME / AttentionAE-sc / scGNN / scDSC / scCDCG
DCA / scVI / scvi-tools
GAN / WGAN-GP / MAE / BERT
single-cell foundation model survey and evaluation
Arabidopsis root scRNA-seq source dataset
scPlantDB plant atlas database
```
