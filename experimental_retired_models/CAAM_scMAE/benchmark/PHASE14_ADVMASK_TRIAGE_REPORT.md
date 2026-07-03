# Phase 14 AdvMask Triage Report

## 1. Run summary

- Complete run count: 18
- Primary metric: `kmeans_known_k.ari`
- Scope: MLP encoder only; control random mask vs AdvMask selector.

## 2. Aggregate results

| dataset | corruption | variant | ARI mean | ARI std | NMI mean | ACC mean | F1 mean | generator grad | mask entropy | mask gini | effective corruption |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Limb_Muscle | scmae_shuffle | advmask | 0.864421 | 0.010631 | 0.854959 | 0.867315 | 0.755041 | 1.000000 | 0.994606 | 0.037225 | 0.135452 |
| Limb_Muscle | scmae_shuffle | control | 0.872681 | 0.006596 | 0.856510 | 0.871920 | 0.756648 | 0.000000 | nan | 0.083480 | 0.131123 |
| Mouse_Pancreas_1 | scmae_shuffle | advmask | 0.345744 | 0.030113 | 0.578930 | 0.467480 | 0.426124 | 1.000000 | 0.998070 | 0.034962 | 0.123250 |
| Mouse_Pancreas_1 | scmae_shuffle | control | 0.330536 | 0.015155 | 0.574736 | 0.466773 | 0.429627 | 0.000000 | nan | 0.072464 | 0.120097 |
| Quake_Smart-seq2_Lung | scmae_shuffle | advmask | 0.491079 | 0.009561 | 0.713218 | 0.563842 | 0.424215 | 1.000000 | 0.998342 | 0.034919 | 0.194552 |
| Quake_Smart-seq2_Lung | scmae_shuffle | control | 0.490818 | 0.006156 | 0.711321 | 0.565036 | 0.423786 | 0.000000 | nan | 0.057461 | 0.191684 |

## 3. AdvMask minus control

`{"Limb_Muscle|scmae_shuffle": {"advmask_minus_control.kmeans_known_k.acc": -0.004604758250191798, "advmask_minus_control.kmeans_known_k.ari": -0.008260467233481839, "advmask_minus_control.kmeans_known_k.f1_macro": -0.0016069480811660153, "advmask_minus_control.kmeans_known_k.nmi": -0.0015508846361692585, "advmask_minus_control.leiden_fixed.acc": 0.03300076745970837, "advmask_minus_control.leiden_fixed.ari": 0.011188086097303818, "advmask_minus_control.leiden_fixed.f1_macro": 0.030715088092634457, "advmask_minus_control.leiden_fixed.nmi": 0.003909431474454772}, "Mouse_Pancreas_1|scmae_shuffle": {"advmask_minus_control.kmeans_known_k.acc": 0.0007069635913750094, "advmask_minus_control.kmeans_known_k.ari": 0.015207930208560594, "advmask_minus_control.kmeans_known_k.f1_macro": -0.0035028522447445587, "advmask_minus_control.kmeans_known_k.nmi": 0.004194050946878969, "advmask_minus_control.leiden_fixed.acc": -0.007069635913750427, "advmask_minus_control.leiden_fixed.ari": -0.022553600736936796, "advmask_minus_control.leiden_fixed.f1_macro": 0.007809605476586601, "advmask_minus_control.leiden_fixed.nmi": -0.007963836165327476}, "Quake_Smart-seq2_Lung|scmae_shuffle": {"advmask_minus_control.kmeans_known_k.acc": -0.0011933174224343368, "advmask_minus_control.kmeans_known_k.ari": 0.0002612480895806679, "advmask_minus_control.kmeans_known_k.f1_macro": 0.00042955058307497795, "advmask_minus_control.kmeans_known_k.nmi": 0.0018972765323518326, "advmask_minus_control.leiden_fixed.acc": -0.0015910898965791898, "advmask_minus_control.leiden_fixed.ari": -0.004502930398125882, "advmask_minus_control.leiden_fixed.f1_macro": 0.0030224634230088254, "advmask_minus_control.leiden_fixed.nmi": -0.003103519792817311}}`

## 4. Phase gate

- gate_result: `fail`
- recommendation: `drop_or_downgrade_advmask`
- positive_ari_dataset_corruptions: `["Mouse_Pancreas_1|scmae_shuffle", "Quake_Smart-seq2_Lung|scmae_shuffle"]`
- mean_ari_delta: `0.002403`
- mean_seed_std_reference: `0.013035`
- effect_size_gate_pass: `False`
- generator_grad_norm_positive: `["Limb_Muscle|scmae_shuffle", "Mouse_Pancreas_1|scmae_shuffle", "Quake_Smart-seq2_Lung|scmae_shuffle"]`
- mask_concentration_flags: `[]`
- embedding_collapse_flags: `[]`

## 5. Remaining risks

- Phase 14 must not run Axial or full model.
- Three-epoch trends must not be written as final publication claims.
- Validation and sealed test datasets are not used for this mechanism decision.
