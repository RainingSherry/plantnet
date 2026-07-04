# NeighborMix x std-floor ablation summary

## Arm-level aggregate

| mix | varw | n | ARI mean+/-sd | NMI mean+/-sd | eff_dim mean+/-sd | std_min mean+/-sd | dims_std>1 mean+/-sd |
|---|---:|---:|---:|---:|---:|---:|---:|
| neighbor | 0.0 | 3 | 0.3278 +/- 0.0057 | 0.5976 +/- 0.0111 | 40.9 +/- 5.6 | 0.080 +/- 0.049 | 31.0 +/- 4.4 |
| neighbor | 0.02 | 3 | 0.4322 +/- 0.0919 | 0.5768 +/- 0.0024 | 113.1 +/- 6.5 | 1.015 +/- 0.007 | 128.0 +/- 0.0 |
| none | 0.0 | 3 | 0.3298 +/- 0.0108 | 0.5949 +/- 0.0022 | 65.6 +/- 7.5 | 0.281 +/- 0.010 | 35.3 +/- 1.5 |
| none | 0.02 | 3 | 0.7018 +/- 0.0051 | 0.6546 +/- 0.0079 | 105.7 +/- 12.4 | 1.023 +/- 0.014 | 128.0 +/- 0.0 |

## Per-run details

| run | mix | varw | ARI | NMI | eff_dim | std_min | std_med | dims_std>1 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| dec | none | 0.0 | 0.3423 | 0.5972 | 60.4 | 0.275 | 0.709 | 34 |
| dec_floor | none | 0.02 | 0.6960 | 0.6459 | 114.2 | 1.033 | 1.060 | 128 |
| dec_floor_seed43 | none | 0.02 | 0.7057 | 0.6615 | 111.4 | 1.007 | 1.042 | 128 |
| dec_floor_seed44 | none | 0.02 | 0.7037 | 0.6563 | 91.5 | 1.029 | 1.069 | 128 |
| dec_seed43 | none | 0.0 | 0.3237 | 0.5929 | 74.3 | 0.293 | 0.800 | 35 |
| dec_seed44 | none | 0.0 | 0.3234 | 0.5947 | 62.2 | 0.276 | 0.731 | 37 |
| neighbormix_dec | neighbor | 0.0 | 0.3302 | 0.5982 | 35.2 | 0.119 | 0.557 | 28 |
| neighbormix_dec_floor | neighbor | 0.02 | 0.4734 | 0.5755 | 105.6 | 1.009 | 1.052 | 128 |
| neighbormix_dec_floor_seed43 | neighbor | 0.02 | 0.4963 | 0.5796 | 117.2 | 1.014 | 1.061 | 128 |
| neighbormix_dec_floor_seed44 | neighbor | 0.02 | 0.3269 | 0.5753 | 116.5 | 1.023 | 1.057 | 128 |
| neighbormix_dec_seed43 | neighbor | 0.0 | 0.3213 | 0.5862 | 41.1 | 0.025 | 0.592 | 29 |
| neighbormix_dec_seed44 | neighbor | 0.0 | 0.3319 | 0.6084 | 46.3 | 0.097 | 0.607 | 36 |

Interpretation guide:

- `none + varw=0`: DEC-only control.
- `none + varw=0.02`: std-floor intervention.
- `neighbor + varw=0`: NeighborMix auxiliary denoising without std-floor.
- `neighbor + varw=0.02`: test whether NeighborMix and std-floor are complementary.
