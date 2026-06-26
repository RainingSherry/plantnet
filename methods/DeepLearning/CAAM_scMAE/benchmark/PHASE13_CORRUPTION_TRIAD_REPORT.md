# Phase 13 Corruption Triad Report

## 1. Smoke validation results

- Smoke run count: 27
- Smoke runs are implementation checks only and are not used as scientific evidence.

## 2. Formal Phase 13 results

- Formal run count: 27
- Primary metric: `kmeans_known_k.ari`

| dataset | corruption | ARI mean | ARI std | NMI mean | ACC mean | F1 mean | zero_to_zero mean | effective mean | mean_abs_delta |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Limb_Muscle | matched_donor | 0.910976 | 0.123817 | 0.928150 | 0.919758 | 0.903425 | 0.876842 | 0.123154 | 0.070328 |
| Limb_Muscle | nonzero_aware_donor | 0.967668 | 0.002784 | 0.935411 | 0.980217 | 0.968942 | 0.537621 | 0.462379 | 0.271318 |
| Limb_Muscle | scmae_shuffle | 0.777303 | 0.002639 | 0.846102 | 0.803274 | 0.740506 | 0.868886 | 0.131095 | 0.083068 |
| Mouse_Pancreas_1 | matched_donor | 0.481717 | 0.050108 | 0.699877 | 0.565394 | 0.450210 | 0.886609 | 0.113389 | 0.065304 |
| Mouse_Pancreas_1 | nonzero_aware_donor | 0.428367 | 0.030744 | 0.599331 | 0.471191 | 0.339248 | 0.554989 | 0.445011 | 0.271980 |
| Mouse_Pancreas_1 | scmae_shuffle | 0.655529 | 0.111788 | 0.716704 | 0.650230 | 0.447556 | 0.879959 | 0.120005 | 0.074899 |
| Quake_Smart-seq2_Lung | matched_donor | 0.412963 | 0.057211 | 0.694227 | 0.478918 | 0.428614 | 0.813957 | 0.186042 | 0.079121 |
| Quake_Smart-seq2_Lung | nonzero_aware_donor | 0.532312 | 0.032470 | 0.672495 | 0.590692 | 0.347014 | 0.361645 | 0.638355 | 0.254778 |
| Quake_Smart-seq2_Lung | scmae_shuffle | 0.543149 | 0.027800 | 0.741984 | 0.609785 | 0.482120 | 0.808278 | 0.191657 | 0.085424 |

## 3. Corruption recommendation for Phase 14

- Recommendation: `scmae_shuffle`
- Reason: matched donor is weaker than scMAE shuffle on at least 2/3 datasets
- Differences: `{"nonzero_aware_minus_scmae_shuffle": {"Limb_Muscle": 0.19036462182477099, "Mouse_Pancreas_1": -0.22716206489600238, "Quake_Smart-seq2_Lung": -0.010836930060872385}, "scmae_shuffle_minus_matched_donor": {"Limb_Muscle": -0.13367293697799842, "Mouse_Pancreas_1": 0.1738118517510418, "Quake_Smart-seq2_Lung": 0.13018537527210433}}`
- Nonzero-aware assessment: `{"assessment": "improves diagnostics but not consistently clustering", "beats_scmae_on_primary_metric": ["Limb_Muscle"], "improves_effective_corruption_rate": ["Limb_Muscle", "Mouse_Pancreas_1", "Quake_Smart-seq2_Lung"]}`

## 4. Remaining risks

- Phase 13 is limited to MLP encoder plus random mask; it does not validate AdvMask, Axial, or full CAAM.
- A corruption that improves mask diagnostics without improving clustering must not be claimed as a main method contribution.
- Results must stay development-only; validation and sealed test are not used in Phase 13.
