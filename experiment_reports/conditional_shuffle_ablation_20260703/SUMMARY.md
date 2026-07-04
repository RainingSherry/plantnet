# Conditional / nuisance-matched shuffle ablation summary

## Arm-level aggregate (per dataset)

| dataset | corruption | n | ARI mean+/-sd | NMI mean+/-sd | eff_dim | aligned_eff_dim | eff_change | donor_pool_mean | n_bins |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Macosko | zero | 3 | 0.7018 +/- 0.0051 | 0.6546 +/- 0.0079 | 105.7 +/- 12.4 | 1.66 +/- 0.03 | 0.400 +/- 0.000 | 44808 | 1 |
| Macosko | swap_global | 3 | 0.2986 +/- 0.0340 | 0.5552 +/- 0.0123 | 126.3 +/- 0.7 | 2.71 +/- 0.01 | 0.024 +/- 0.000 | 44808 | 1 |
| Macosko | swap_lib | 3 | 0.8672 +/- 0.0028 | 0.7467 +/- 0.0133 | 126.6 +/- 1.0 | 2.40 +/- 0.40 | 0.023 +/- 0.000 | 4481 | 10 |
| Macosko | swap_ndet | 3 | 0.8640 +/- 0.0014 | 0.7257 +/- 0.0009 | 122.3 +/- 5.8 | 2.03 +/- 0.19 | 0.023 +/- 0.000 | 4481 | 10 |
| Macosko | swap_zerolib | 3 | 0.7586 +/- 0.1804 | 0.6992 +/- 0.0571 | 127.0 +/- 0.9 | 2.24 +/- 0.58 | 0.023 +/- 0.000 | 1948 | 23 |
| Melanoma_5K | zero | 3 | 0.6483 +/- 0.0030 | 0.7165 +/- 0.0052 | 126.1 +/- 0.8 | 1.42 +/- 0.08 | 0.400 +/- 0.000 | 4513 | 1 |
| Melanoma_5K | swap_global | 3 | 0.6553 +/- 0.0059 | 0.7221 +/- 0.0004 | 126.2 +/- 0.3 | 1.40 +/- 0.05 | 0.041 +/- 0.000 | 4513 | 1 |
| Melanoma_5K | swap_lib | 3 | 0.6022 +/- 0.0933 | 0.7055 +/- 0.0283 | 126.8 +/- 0.4 | 1.32 +/- 0.04 | 0.040 +/- 0.000 | 451 | 10 |
| Melanoma_5K | swap_ndet | 3 | 0.6579 +/- 0.0068 | 0.7196 +/- 0.0035 | 127.0 +/- 0.1 | 1.32 +/- 0.04 | 0.039 +/- 0.000 | 451 | 10 |
| Melanoma_5K | swap_zerolib | 3 | 0.6597 +/- 0.0038 | 0.7162 +/- 0.0007 | 126.8 +/- 0.4 | 1.36 +/- 0.02 | 0.039 +/- 0.000 | 181 | 25 |
| Quake_10x_Spleen | zero | 3 | 0.9201 +/- 0.0009 | 0.8307 +/- 0.0089 | 127.6 +/- 0.1 | 1.10 +/- 0.01 | 0.400 +/- 0.000 | 9552 | 1 |
| Quake_10x_Spleen | swap_global | 3 | 0.9226 +/- 0.0006 | 0.8418 +/- 0.0032 | 127.8 +/- 0.0 | 1.16 +/- 0.01 | 0.034 +/- 0.000 | 9552 | 1 |
| Quake_10x_Spleen | swap_lib | 3 | 0.9149 +/- 0.0085 | 0.8168 +/- 0.0243 | 127.7 +/- 0.1 | 1.16 +/- 0.01 | 0.033 +/- 0.000 | 955 | 10 |
| Quake_10x_Spleen | swap_ndet | 3 | 0.9216 +/- 0.0021 | 0.8383 +/- 0.0076 | 127.7 +/- 0.1 | 1.18 +/- 0.01 | 0.033 +/- 0.000 | 955 | 10 |
| Quake_10x_Spleen | swap_zerolib | 3 | 0.9167 +/- 0.0023 | 0.8108 +/- 0.0042 | 127.4 +/- 0.3 | 1.11 +/- 0.01 | 0.032 +/- 0.000 | 382 | 25 |

## Delta vs zero (ARI)

| dataset | corruption | ARI mean | delta | verdict (delta>=+0.02) |
|---|---|---:|---:|---|
| Macosko | zero | 0.7018 | +0.0000 | tie |
| Macosko | swap_global | 0.2986 | -0.4031 | down |
| Macosko | swap_lib | 0.8672 | +0.1654 | up |
| Macosko | swap_ndet | 0.8640 | +0.1622 | up |
| Macosko | swap_zerolib | 0.7586 | +0.0568 | up |
| Melanoma_5K | zero | 0.6483 | +0.0000 | tie |
| Melanoma_5K | swap_global | 0.6553 | +0.0070 | tie |
| Melanoma_5K | swap_lib | 0.6022 | -0.0461 | down |
| Melanoma_5K | swap_ndet | 0.6579 | +0.0096 | tie |
| Melanoma_5K | swap_zerolib | 0.6597 | +0.0115 | tie |
| Quake_10x_Spleen | zero | 0.9201 | +0.0000 | tie |
| Quake_10x_Spleen | swap_global | 0.9226 | +0.0025 | tie |
| Quake_10x_Spleen | swap_lib | 0.9149 | -0.0052 | tie |
| Quake_10x_Spleen | swap_ndet | 0.9216 | +0.0015 | tie |
| Quake_10x_Spleen | swap_zerolib | 0.9167 | -0.0034 | tie |

## Delta vs swap_global (ARI)

| dataset | corruption | ARI mean | delta | verdict (delta>=+0.02) |
|---|---|---:|---:|---|
| Macosko | zero | 0.7018 | +0.4031 | up |
| Macosko | swap_global | 0.2986 | +0.0000 | tie |
| Macosko | swap_lib | 0.8672 | +0.5685 | up |
| Macosko | swap_ndet | 0.8640 | +0.5653 | up |
| Macosko | swap_zerolib | 0.7586 | +0.4599 | up |
| Melanoma_5K | zero | 0.6483 | -0.0070 | tie |
| Melanoma_5K | swap_global | 0.6553 | +0.0000 | tie |
| Melanoma_5K | swap_lib | 0.6022 | -0.0532 | down |
| Melanoma_5K | swap_ndet | 0.6579 | +0.0026 | tie |
| Melanoma_5K | swap_zerolib | 0.6597 | +0.0044 | tie |
| Quake_10x_Spleen | zero | 0.9201 | -0.0025 | tie |
| Quake_10x_Spleen | swap_global | 0.9226 | +0.0000 | tie |
| Quake_10x_Spleen | swap_lib | 0.9149 | -0.0076 | tie |
| Quake_10x_Spleen | swap_ndet | 0.9216 | -0.0010 | tie |
| Quake_10x_Spleen | swap_zerolib | 0.9167 | -0.0058 | tie |

## Per-run details

| run | dataset | corruption | seed | ARI | NMI | eff_dim | aligned | eff_change | pool_mean | pool_min |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Macosko__swap_global__seed42 | Macosko | swap_global | 42 | 0.2967 | 0.5543 | 126.5 | 2.73 | 0.024 | 44808 | 44808 |
| Macosko__swap_global__seed43 | Macosko | swap_global | 43 | 0.2656 | 0.5434 | 126.8 | 2.71 | 0.024 | 44808 | 44808 |
| Macosko__swap_global__seed44 | Macosko | swap_global | 44 | 0.3336 | 0.5680 | 125.5 | 2.71 | 0.024 | 44808 | 44808 |
| Macosko__swap_lib__seed42 | Macosko | swap_lib | 42 | 0.8703 | 0.7314 | 127.8 | 1.94 | 0.023 | 4481 | 4480 |
| Macosko__swap_lib__seed43 | Macosko | swap_lib | 43 | 0.8658 | 0.7533 | 125.9 | 2.58 | 0.023 | 4481 | 4480 |
| Macosko__swap_lib__seed44 | Macosko | swap_lib | 44 | 0.8654 | 0.7554 | 126.0 | 2.69 | 0.023 | 4481 | 4480 |
| Macosko__swap_ndet__seed42 | Macosko | swap_ndet | 42 | 0.8645 | 0.7267 | 127.8 | 1.97 | 0.023 | 4481 | 4480 |
| Macosko__swap_ndet__seed43 | Macosko | swap_ndet | 43 | 0.8650 | 0.7251 | 116.2 | 2.23 | 0.023 | 4481 | 4480 |
| Macosko__swap_ndet__seed44 | Macosko | swap_ndet | 44 | 0.8624 | 0.7252 | 122.9 | 1.87 | 0.023 | 4481 | 4480 |
| Macosko__swap_zerolib__seed42 | Macosko | swap_zerolib | 42 | 0.8664 | 0.7365 | 127.7 | 1.89 | 0.023 | 1948 | 2 |
| Macosko__swap_zerolib__seed43 | Macosko | swap_zerolib | 43 | 0.5503 | 0.6335 | 126.0 | 2.91 | 0.023 | 1948 | 2 |
| Macosko__swap_zerolib__seed44 | Macosko | swap_zerolib | 44 | 0.8590 | 0.7275 | 127.2 | 1.93 | 0.023 | 1948 | 2 |
| Macosko__zero__seed42 | Macosko | zero | 42 | 0.6960 | 0.6459 | 114.2 | 1.68 | 0.400 | 44808 | 44808 |
| Macosko__zero__seed43 | Macosko | zero | 43 | 0.7057 | 0.6615 | 111.4 | 1.67 | 0.400 | 44808 | 44808 |
| Macosko__zero__seed44 | Macosko | zero | 44 | 0.7037 | 0.6563 | 91.5 | 1.63 | 0.400 | 44808 | 44808 |
| Melanoma_5K__swap_global__seed42 | Melanoma_5K | swap_global | 42 | 0.6527 | 0.7224 | 126.5 | 1.38 | 0.041 | 4513 | 4513 |
| Melanoma_5K__swap_global__seed43 | Melanoma_5K | swap_global | 43 | 0.6621 | 0.7221 | 125.9 | 1.46 | 0.041 | 4513 | 4513 |
| Melanoma_5K__swap_global__seed44 | Melanoma_5K | swap_global | 44 | 0.6512 | 0.7216 | 126.2 | 1.37 | 0.041 | 4513 | 4513 |
| Melanoma_5K__swap_lib__seed42 | Melanoma_5K | swap_lib | 42 | 0.6592 | 0.7230 | 127.1 | 1.32 | 0.040 | 451 | 451 |
| Melanoma_5K__swap_lib__seed43 | Melanoma_5K | swap_lib | 43 | 0.4945 | 0.6728 | 126.4 | 1.29 | 0.040 | 451 | 451 |
| Melanoma_5K__swap_lib__seed44 | Melanoma_5K | swap_lib | 44 | 0.6528 | 0.7206 | 127.0 | 1.36 | 0.040 | 451 | 451 |
| Melanoma_5K__swap_ndet__seed42 | Melanoma_5K | swap_ndet | 42 | 0.6656 | 0.7235 | 127.0 | 1.29 | 0.039 | 451 | 451 |
| Melanoma_5K__swap_ndet__seed43 | Melanoma_5K | swap_ndet | 43 | 0.6525 | 0.7176 | 127.0 | 1.32 | 0.039 | 451 | 451 |
| Melanoma_5K__swap_ndet__seed44 | Melanoma_5K | swap_ndet | 44 | 0.6557 | 0.7176 | 127.1 | 1.36 | 0.039 | 451 | 451 |
| Melanoma_5K__swap_zerolib__seed42 | Melanoma_5K | swap_zerolib | 42 | 0.6603 | 0.7165 | 126.5 | 1.36 | 0.039 | 181 | 2 |
| Melanoma_5K__swap_zerolib__seed43 | Melanoma_5K | swap_zerolib | 43 | 0.6632 | 0.7153 | 126.6 | 1.38 | 0.039 | 181 | 2 |
| Melanoma_5K__swap_zerolib__seed44 | Melanoma_5K | swap_zerolib | 44 | 0.6557 | 0.7167 | 127.2 | 1.33 | 0.039 | 181 | 2 |
| Melanoma_5K__zero__seed42 | Melanoma_5K | zero | 42 | 0.6449 | 0.7196 | 126.7 | 1.36 | 0.400 | 4513 | 4513 |
| Melanoma_5K__zero__seed43 | Melanoma_5K | zero | 43 | 0.6492 | 0.7194 | 126.4 | 1.39 | 0.400 | 4513 | 4513 |
| Melanoma_5K__zero__seed44 | Melanoma_5K | zero | 44 | 0.6507 | 0.7105 | 125.2 | 1.52 | 0.400 | 4513 | 4513 |
| Quake_10x_Spleen__swap_global__seed42 | Quake_10x_Spleen | swap_global | 42 | 0.9219 | 0.8382 | 127.9 | 1.15 | 0.034 | 9552 | 9552 |
| Quake_10x_Spleen__swap_global__seed43 | Quake_10x_Spleen | swap_global | 43 | 0.9230 | 0.8428 | 127.8 | 1.16 | 0.034 | 9552 | 9552 |
| Quake_10x_Spleen__swap_global__seed44 | Quake_10x_Spleen | swap_global | 44 | 0.9228 | 0.8445 | 127.8 | 1.17 | 0.034 | 9552 | 9552 |
| Quake_10x_Spleen__swap_lib__seed42 | Quake_10x_Spleen | swap_lib | 42 | 0.9211 | 0.8330 | 127.8 | 1.16 | 0.033 | 955 | 955 |
| Quake_10x_Spleen__swap_lib__seed43 | Quake_10x_Spleen | swap_lib | 43 | 0.9052 | 0.7889 | 127.7 | 1.15 | 0.033 | 955 | 955 |
| Quake_10x_Spleen__swap_lib__seed44 | Quake_10x_Spleen | swap_lib | 44 | 0.9186 | 0.8285 | 127.6 | 1.17 | 0.033 | 955 | 955 |
| Quake_10x_Spleen__swap_ndet__seed42 | Quake_10x_Spleen | swap_ndet | 42 | 0.9235 | 0.8459 | 127.8 | 1.17 | 0.033 | 955 | 955 |
| Quake_10x_Spleen__swap_ndet__seed43 | Quake_10x_Spleen | swap_ndet | 43 | 0.9221 | 0.8384 | 127.8 | 1.19 | 0.033 | 955 | 955 |
| Quake_10x_Spleen__swap_ndet__seed44 | Quake_10x_Spleen | swap_ndet | 44 | 0.9193 | 0.8307 | 127.6 | 1.16 | 0.033 | 955 | 955 |
| Quake_10x_Spleen__swap_zerolib__seed42 | Quake_10x_Spleen | swap_zerolib | 42 | 0.9143 | 0.8062 | 127.3 | 1.10 | 0.032 | 382 | 1 |
| Quake_10x_Spleen__swap_zerolib__seed43 | Quake_10x_Spleen | swap_zerolib | 43 | 0.9172 | 0.8119 | 127.2 | 1.11 | 0.032 | 382 | 1 |
| Quake_10x_Spleen__swap_zerolib__seed44 | Quake_10x_Spleen | swap_zerolib | 44 | 0.9188 | 0.8144 | 127.8 | 1.11 | 0.032 | 382 | 1 |
| Quake_10x_Spleen__zero__seed42 | Quake_10x_Spleen | zero | 42 | 0.9211 | 0.8409 | 127.7 | 1.10 | 0.400 | 9552 | 9552 |
| Quake_10x_Spleen__zero__seed43 | Quake_10x_Spleen | zero | 43 | 0.9199 | 0.8243 | 127.5 | 1.11 | 0.400 | 9552 | 9552 |
| Quake_10x_Spleen__zero__seed44 | Quake_10x_Spleen | zero | 44 | 0.9193 | 0.8268 | 127.8 | 1.11 | 0.400 | 9552 | 9552 |

## Decision rule

- **S1/S2/S3 > S0 (swap_global)**：ARI +>=0.02 且多种子稳定 -> 限制供体池确实逼模型学更细结构，是从 scMAE 内部改对了。
- **S1/S2/S3 ~= S0**：nuisance 匹配无效 -> swap 捷径不是瓶颈。
- **swap_global(S0) vs zero**：先确认 swap-noise 本身相对零填充的效应方向。
- **donor_pool_mean / pool_min 过小（个位数）**：swap 退化成近似恒等，该 arm 结果不可信 —— 需调大 n_nuisance_bins 对应的箱宽（减少箱数）。
- **eff_change 远低于 mask_prob**：供体与目标高度相同（池太窄或过稀疏），同样是退化信号。
