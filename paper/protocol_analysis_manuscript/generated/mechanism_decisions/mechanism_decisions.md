# Mechanism Decision Matrix

Status: generated development-evidence decision matrix. This is not validation evidence.

## Decisions

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

## Key Paired Deltas

- Phase 14 AdvMask mean known-K ARI delta: 0.002403.
- Phase 14 AdvMask mean fixed-Leiden ARI delta: -0.005289.
- Attention smoke Axial mean known-K ARI delta vs control: -0.370706.
- Resource smoke AdvMask mean known-K ARI delta vs control: 0.001150.
- Resource smoke Axial mean known-K ARI delta vs control: -0.365104.

## Safe Route

Continue as protocol-analysis / diagnostic paper unless a new mechanism later passes a fresh development gate. Do not claim AdvMask, current Axial, full CAAM synergy, publication validation, or unknown-K clustering superiority from current evidence.
