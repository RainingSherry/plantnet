# Claim Boundary Audit

Status: generated manuscript-support audit. This is not validation evidence.

## Evidence Snapshot

- Generated artifact files: 65
- Phase 13 runs: 27
- Phase 14 runs: 18
- Attention/context smoke runs: 9
- Instrumented resource-smoke runs: 12
- Feature-space smoke runs: 2
- Mask-ratio smoke runs: 9
- Resource-smoke datasets: Quake_Smart-seq2_Lung, Mouse_Pancreas_1, Limb_Muscle
- Resource-smoke seed/epochs: 42 / 3
- GPU memory mode: total_gpu_memory_delta

## Instrumented Resource Summary

| condition | known-K ARI mean | fixed-Leiden ARI mean | wall sec mean | GPU delta MiB mean |
|---|---:|---:|---:|---:|
| control | 0.563205 | 0.420215 | 20.45 | 1449.3 |
| advmask | 0.564355 | 0.416609 | 19.12 | 2094.7 |
| axial | 0.198101 | 0.101700 | 35.33 | 2868.7 |
| mlp_parammatched | 0.576788 | 0.404586 | 17.44 | 1450.7 |

## Generated Mechanism Decisions

| study | mechanism | decision | claim boundary |
|---|---|---|---|
| phase13_corruption_triad | matched_donor | diagnostic_only | development evidence only; not validation |
| phase13_corruption_triad | nonzero_aware_donor | diagnostic_only | development evidence only; not validation |
| phase13_corruption_triad | scmae_shuffle | carry_forward | development evidence only; not validation |
| phase14_advmask_triage | advmask | drop_or_downgrade | active generator is not evidence of useful clustering representation |
| attention_context_smoke | current_axial_encoder | do_not_use_as_rescue_path | not evidence that all attention mechanisms fail |
| instrumented_resource_smoke | advmask | cost_without_stable_gain | resource smoke only; not submission-scale runtime benchmark |
| instrumented_resource_smoke | axial | cost_without_stable_gain | resource smoke only; not submission-scale runtime benchmark |
| instrumented_resource_smoke | control | diagnostic_control | resource smoke only; not submission-scale runtime benchmark |
| instrumented_resource_smoke | mlp_parammatched | diagnostic_control | resource smoke only; not submission-scale runtime benchmark |

## Feature-Space Smoke Boundary

| role | genes | student params | known-K ARI | fixed-Leiden ARI |
|---|---:|---:|---:|---:|
| hvg2000 | 2000 | 9080688 | 0.485571 | 0.463516 |
| full_gene_stress | 23341 | 1101675865 | 0.443727 | 0.389245 |
- Full-gene minus HVG known-K ARI: -0.041844.
- Full-gene minus HVG fixed-Leiden ARI: -0.074271.
- Full-gene/HVG student parameter ratio: 121.32x.
- This is a single-dataset dense-MLP development smoke, not validation and not a universal claim against full-gene modeling.

## Mask-Ratio Smoke Boundary

| mask ratio | runs | masked effective mean | global effective mean | known-K ARI mean | known-K ARI std | fixed-Leiden ARI mean | fixed-Leiden ARI std |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 3 | 0.191674 | 0.038335 | 0.462024 | 0.048981 | 0.408229 | 0.048769 |
| 0.4 | 3 | 0.191684 | 0.076674 | 0.490818 | 0.006156 | 0.430972 | 0.061284 |
| 0.6 | 3 | 0.191593 | 0.114956 | 0.470444 | 0.016125 | 0.480897 | 0.027309 |
- Best mean known-K ARI in this smoke: mask ratio 0.4.
- Best mean fixed-Leiden ARI in this smoke: mask ratio 0.6.
- Known-K ARI mean range: 0.028794; fixed-Leiden ARI mean range: 0.072668.
- Maximum known-K seed standard deviation: 0.048981; maximum fixed-Leiden seed standard deviation: 0.061284.
- This is a single-dataset mask-ratio sensitivity diagnostic, not a tuning decision for validation.

## Supported Claims

- Development evidence supports a protocol-analysis route rather than a positive CAAM method route.
- HVG 2000 remains the current dense-MLP feature-space default under development evidence.
- Nominal mask ratio can change global effective perturbation without producing monotonic downstream clustering behavior.
- Mask-ratio figures are development diagnostics and must be read with the table-level claim boundary.
- scMAE-style shuffle remains the safest corruption choice to carry forward in the frozen development protocol.
- Effective corruption diagnostics and downstream clustering quality can diverge.
- AdvMask is active but unsupported as a main clustering-improvement mechanism under current development evidence.
- The current Axial encoder is unsupported as a rescue path and should not be treated as evidence against all attention mechanisms.
- Parameter-matched MLP controls are necessary before interpreting architecture effects.

## Unsupported Claims

- CAAM-scMAE improves clustering.
- AdvMask improves clustering.
- The current Axial encoder improves clustering.
- Axial plus AdvMask has synergy.
- Development evidence is publication-level validation.
- Known-K ARI proves fully unsupervised unknown-K clustering quality.
- HVG 2000 is validated by the feature-space smoke.
- Full-gene or gene-token feature-space models are generally inferior.
- A different mask ratio should be selected from the single-dataset smoke.
- Mask ratio 0.6 is generally best because it has the highest fixed-Leiden mean in this smoke.
- The frozen validation mask ratio should be changed based on this smoke.

## Audit Findings

No blocking claim-boundary findings.

## Safe Next-Step Boundary

Validation is still not run. Any validation pass must use the frozen protocol and must not tune corruption, mask policy, architecture, loss, clustering resolution, or manuscript claims.
