# Phase 14 AdvMask Triage Decision

Date: 2026-06-26

Authoring status: research decision memo for manuscript planning. This is development evidence, not a publication claim.

## 1. Evidence Source

Phase 14 tested whether a constrained adversarial mask selector is stronger than random masking under the corrected scMAE-style protocol.

Protocol:

```text
datasets = Quake_Smart-seq2_Lung, Mouse_Pancreas_1, Limb_Muscle
corruption_type = scmae_shuffle
variants = control, advmask
seeds = 42, 2024, 3407
epochs = 3
input_mode = log1p
n_top_genes = 2000
scale_input = false
validation/sealed test = not used
Axial/full model = not run
```

Important engineering correction:

```text
The first Phase 14 attempt was invalid because epochs=3 but student_warmup_epochs=20.
That made AdvMask identical to random masking.
The runner was fixed to set student_warmup_epochs=1 for 3-epoch triage.
The final report below uses the corrected run.
```

Report:

```text
methods/DeepLearning/CAAM_scMAE/benchmark/PHASE14_ADVMASK_TRIAGE_REPORT.md
```

## 2. Result Summary

Primary metric:

```text
kmeans_known_k.ari
```

Aggregate ARI:

| dataset | control ARI | advmask ARI | advmask - control |
|---|---:|---:|---:|
| Limb_Muscle | 0.872681 | 0.864421 | -0.008260 |
| Mouse_Pancreas_1 | 0.330536 | 0.345744 | 0.015208 |
| Quake_Smart-seq2_Lung | 0.490818 | 0.491079 | 0.000261 |

Gate diagnostics:

```text
positive_ari_dataset_corruptions = 2 / 3
mean_ari_delta = 0.002403
mean_seed_std_reference = 0.013035
0.5 * mean_seed_std_reference = 0.006518
effect_size_gate_pass = false
generator_grad_norm_positive = 3 / 3
mask_concentration_flags = []
embedding_collapse_flags = []
gate_result = fail
recommendation = drop_or_downgrade_advmask
```

## 3. Interpretation

The AdvMask implementation is not trivially broken after the warmup fix: generator gradients are nonzero on all three development datasets, and the observed mask diagnostics do not show top-gene collapse or embedding collapse in this 3-epoch triage.

However, the clustering improvement is too small to survive the predefined seed-variation criterion. Two datasets have positive ARI deltas, but one is essentially zero and the average delta is below half of the seed standard-deviation reference. Therefore, Phase 14 does not support retaining AdvMask as a main method component.

This is an informative negative result: learned adversarial mask selection can be made to train, but under the corrected HVG + scMAE-style shuffle protocol it does not provide a stable clustering benefit over random masking.

## 4. Answers After This Stage

### 4.1 What Does The Result Support?

The result supports the following conservative claims:

```text
1. The corrected protocol can distinguish implementation activity from scientific usefulness.
2. AdvMask can produce real generator gradients under the fixed 3-epoch triage schedule.
3. A learned mask selector does not automatically improve clustering representation quality.
4. Random masking remains a strong and safer baseline under scMAE-style gene-wise shuffle.
```

### 4.2 What Does The Result Not Support?

The result does not support:

```text
1. AdvMask as a primary contribution.
2. A claim that learned masks are more informative than random masks for clustering.
3. AdvMask + Axial synergy.
4. Running the full model to rescue the AdvMask hypothesis.
5. A positive CAAM-scMAE method-paper narrative centered on adversarial masking.
```

### 4.3 Which Modules Should Be Kept?

Keep for now:

```text
1. corrected HVG protocol
2. scMAE-style gene-wise shuffle corruption
3. corruption diagnostics: zero_to_zero_rate, effective_corruption_rate, mean_abs_delta
4. random-mask MLP control as the main development baseline
5. AdvMask code as supplementary/negative-diagnostic infrastructure only
```

### 4.4 Which Modules Should Be Deleted Or Downgraded?

Downgrade:

```text
AdvMask should be downgraded from main method component to supplementary negative result.
```

Do not delete immediately:

```text
The code remains useful to document a tested negative hypothesis and to prevent later reintroducing the same idea without evidence.
```

### 4.5 Most Likely Paper Route Now

Current best route:

```text
protocol_analysis
```

Working route sentence:

```text
Rather than proposing a larger masked autoencoder by default, this study analyzes how feature-space choice, corruption semantics, and mask-selection policy affect clustering-oriented masked autoencoding for sparse scRNA-seq data.
```

The method-paper route is not currently supported. It would require new evidence from an independently justified context/attention module, not from AdvMask.

### 4.6 If Written As A Method Paper

Current evidence is insufficient. A method-paper contribution would need:

```text
1. a module that beats the corrected random-mask scMAE-style control beyond seed variation;
2. parameter-matched controls;
3. validation datasets after development protocol freeze;
4. biological interpretation that does not degrade;
5. no dependence on known-K or oracle clustering.
```

### 4.7 If Written As An Analysis Paper

Supported early findings:

```text
1. Full-gene/default protocol choices can distort masked AE conclusions.
2. Effective corruption diagnostics matter but do not guarantee clustering improvement.
3. Matched donor is not automatically better than scMAE-style shuffle.
4. Nonzero-aware donor improves corruption diagnostics but not stable clustering.
5. AdvMask trains but does not improve clustering beyond seed variation.
```

Potential analysis-paper contribution:

```text
Masked autoencoder performance in scRNA-seq clustering is strongly shaped by protocol details, and more complex corruption or mask-selection mechanisms can fail even when their diagnostics appear active.
```

### 4.8 Missing Experiments Before Submission

Still missing:

```text
1. independent validation datasets after freezing the corrected protocol;
2. biological interpretation analysis under the selected route;
3. runtime and memory comparison against selected baselines;
4. marker-overlap or rare-cell diagnostics;
5. final comparison against scMAE/scNAME/Leiden/Louvain and selected deep baselines;
6. if attention is revisited, parameter-matched MLP controls and no AdvMask dependency.
```

### 4.9 Claims That Cannot Be Written

Forbidden claims after Phase 14:

```text
1. "AdvMask improves clustering."
2. "Adversarial masking is the main contribution."
3. "CAAM demonstrates synergy between Axial and AdvMask."
4. "The full CAAM model is validated."
5. "The 3-epoch development trend proves publication-level improvement."
```

### 4.10 Next Smallest Informative Step

Recommended next step:

```text
Write a Phase 16 preliminary publication-decision document that records protocol_analysis as the current working route, while marking final route selection as pending validation approval.
```

Do not run next without explicit decision:

```text
1. validation datasets;
2. sealed test;
3. Axial/full model;
4. any attempt to retune AdvMask based on the Phase 14 ARI table.
```

## 5. Manuscript Reframing

The manuscript should be reframed away from:

```text
Context-aware Adversarial Axial Masked Autoencoding
```

and toward:

```text
Protocol-aware Masked Autoencoding for Single-cell RNA-seq Clustering
```

or:

```text
When Does Masked Corruption Help Single-cell Clustering?
```

Draft gap sentence:

```text
Although masked autoencoding has become a promising objective for single-cell clustering, it remains unclear how feature-space selection, corruption semantics, and learned masking policies affect representation quality under sparse scRNA-seq protocols.
```

Draft negative-result sentence:

```text
Our development triage shows that a constrained adversarial mask selector can produce nonzero generator gradients without collapse, yet still fail to improve clustering beyond seed-level variation, emphasizing that active mask learning is not sufficient evidence of useful representation learning.
```

## 6. Operational Decision

```text
Do not start Phase 15 under the current BDD.
Do not run Axial/full as a rescue experiment for AdvMask.
Keep PR #12 as a record of the Phase 14 negative gate.
Prepare a Phase 16 decision draft before any new mechanism experiment.
```
