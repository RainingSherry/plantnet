# NeighborMix x std-floor ablation summary

| run | mix | varw | ARI | NMI | eff_dim | std_min | std_med | dims_std>1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dec | none | 0.0 | 0.3423 | 0.5972 | 60.4 | 0.275 | 0.709 | 34 |
| dec_floor | none | 0.02 | 0.6960 | 0.6459 | 114.2 | 1.033 | 1.060 | 128 |
| neighbormix_dec | neighbor | 0.0 | 0.3302 | 0.5982 | 35.2 | 0.119 | 0.557 | 28 |
| neighbormix_dec_floor | neighbor | 0.02 | 0.4734 | 0.5755 | 105.6 | 1.009 | 1.052 | 128 |

Interpretation guide:

- `none + varw=0`: DEC-only control.
- `none + varw=0.02`: std-floor intervention.
- `neighbor + varw=0`: NeighborMix auxiliary denoising without std-floor.
- `neighbor + varw=0.02`: test whether NeighborMix and std-floor are complementary.
