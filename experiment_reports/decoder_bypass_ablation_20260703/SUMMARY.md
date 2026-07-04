# Decoder-bypass 消融汇总

## 按 (数据集 × decoder_mode) 聚合

| 数据集 | decoder_mode | n | ARI mean±sd | NMI mean±sd | eff_dim | aligned_eff_dim |
|---|---|---:|---:|---:|---:|---:|
| Macosko | concat (D0=赢家对照) | 3 | 0.7018 +/- 0.0051 | 0.6546 +/- 0.0079 | 105.7 +/- 12.4 | 1.66 +/- 0.03 |
| Macosko | none | 3 | 0.4708 +/- 0.1895 | 0.6100 +/- 0.0392 | 87.0 +/- 22.8 | 1.76 +/- 0.08 |
| Macosko | lowrank | 3 | 0.6272 +/- 0.2125 | 0.6519 +/- 0.0641 | 94.6 +/- 20.0 | 1.71 +/- 0.04 |
| Melanoma_5K | concat (D0=赢家对照) | 3 | 0.6483 +/- 0.0030 | 0.7165 +/- 0.0052 | 126.1 +/- 0.8 | 1.42 +/- 0.08 |
| Melanoma_5K | none | 3 | 0.7153 +/- 0.0520 | 0.7364 +/- 0.0184 | 125.5 +/- 0.5 | 1.38 +/- 0.03 |
| Melanoma_5K | lowrank | 3 | 0.6781 +/- 0.0512 | 0.7263 +/- 0.0154 | 125.9 +/- 0.9 | 1.35 +/- 0.05 |
| Quake_10x_Spleen | concat (D0=赢家对照) | 3 | 0.9201 +/- 0.0009 | 0.8307 +/- 0.0089 | 127.6 +/- 0.1 | 1.10 +/- 0.01 |
| Quake_10x_Spleen | none | 3 | 0.9208 +/- 0.0008 | 0.8337 +/- 0.0012 | 127.6 +/- 0.1 | 1.10 +/- 0.00 |
| Quake_10x_Spleen | lowrank | 3 | 0.9203 +/- 0.0009 | 0.8317 +/- 0.0057 | 127.6 +/- 0.1 | 1.10 +/- 0.01 |

## 相对 concat (D0 原始 scMAE decoder) 的 ARI delta

| 数据集 | decoder_mode | ARI mean | delta vs D0 | 判定(|Δ|≥0.02) |
|---|---|---:|---:|---|
| Macosko | concat | 0.7018 | +0.0000 | ≈持平 |
| Macosko | none | 0.4708 | -0.2310 | ↓变差 |
| Macosko | lowrank | 0.6272 | -0.0746 | ↓变差 |
| Melanoma_5K | concat | 0.6483 | +0.0000 | ≈持平 |
| Melanoma_5K | none | 0.7153 | +0.0671 | ↑改善 |
| Melanoma_5K | lowrank | 0.6781 | +0.0298 | ↑改善 |
| Quake_10x_Spleen | concat | 0.9201 | +0.0000 | ≈持平 |
| Quake_10x_Spleen | none | 0.9208 | +0.0007 | ≈持平 |
| Quake_10x_Spleen | lowrank | 0.9203 | +0.0002 | ≈持平 |

## 每-run 明细

| run | ARI | NMI | eff_dim | aligned | std_min | base_loss |
|---|---:|---:|---:|---:|---:|---:|
| Macosko__concat__seed42 | 0.6960 | 0.6459 | 114.2 | 1.68 | 1.033 | 0.5120 |
| Macosko__concat__seed43 | 0.7057 | 0.6615 | 111.4 | 1.67 | 1.007 | 0.5116 |
| Macosko__concat__seed44 | 0.7037 | 0.6563 | 91.5 | 1.63 | 1.029 | 0.5046 |
| Macosko__lowrank__seed42 | 0.5223 | 0.6099 | 77.7 | 1.69 | 1.026 | 0.5252 |
| Macosko__lowrank__seed43 | 0.8717 | 0.7258 | 116.6 | 1.67 | 1.033 | 0.5013 |
| Macosko__lowrank__seed44 | 0.4874 | 0.6201 | 89.5 | 1.76 | 1.033 | 0.5157 |
| Macosko__none__seed42 | 0.3640 | 0.5830 | 106.5 | 1.79 | 1.012 | 0.5344 |
| Macosko__none__seed43 | 0.6896 | 0.6549 | 62.0 | 1.67 | 1.017 | 0.5121 |
| Macosko__none__seed44 | 0.3588 | 0.5921 | 92.4 | 1.82 | 1.024 | 0.5326 |
| Melanoma_5K__concat__seed42 | 0.6449 | 0.7196 | 126.7 | 1.36 | 1.060 | 0.4651 |
| Melanoma_5K__concat__seed43 | 0.6492 | 0.7194 | 126.4 | 1.39 | 1.042 | 0.4655 |
| Melanoma_5K__concat__seed44 | 0.6507 | 0.7105 | 125.2 | 1.52 | 1.054 | 0.4663 |
| Melanoma_5K__lowrank__seed42 | 0.6446 | 0.7149 | 126.6 | 1.31 | 1.035 | 0.4671 |
| Melanoma_5K__lowrank__seed43 | 0.6525 | 0.7202 | 126.2 | 1.33 | 1.062 | 0.4648 |
| Melanoma_5K__lowrank__seed44 | 0.7370 | 0.7439 | 124.8 | 1.40 | 1.040 | 0.4673 |
| Melanoma_5K__none__seed42 | 0.6554 | 0.7154 | 126.0 | 1.40 | 1.011 | 0.4697 |
| Melanoma_5K__none__seed43 | 0.7427 | 0.7441 | 125.1 | 1.35 | 1.053 | 0.4708 |
| Melanoma_5K__none__seed44 | 0.7479 | 0.7496 | 125.4 | 1.39 | 1.049 | 0.4664 |
| Quake_10x_Spleen__concat__seed42 | 0.9211 | 0.8409 | 127.7 | 1.10 | 1.033 | 0.4477 |
| Quake_10x_Spleen__concat__seed43 | 0.9199 | 0.8243 | 127.5 | 1.11 | 1.038 | 0.4490 |
| Quake_10x_Spleen__concat__seed44 | 0.9193 | 0.8268 | 127.8 | 1.11 | 1.045 | 0.4482 |
| Quake_10x_Spleen__lowrank__seed42 | 0.9210 | 0.8375 | 127.6 | 1.10 | 1.033 | 0.4482 |
| Quake_10x_Spleen__lowrank__seed43 | 0.9206 | 0.8316 | 127.4 | 1.11 | 1.040 | 0.4483 |
| Quake_10x_Spleen__lowrank__seed44 | 0.9193 | 0.8260 | 127.6 | 1.10 | 1.046 | 0.4485 |
| Quake_10x_Spleen__none__seed42 | 0.9200 | 0.8323 | 127.5 | 1.10 | 1.033 | 0.4485 |
| Quake_10x_Spleen__none__seed43 | 0.9211 | 0.8345 | 127.5 | 1.10 | 1.037 | 0.4480 |
| Quake_10x_Spleen__none__seed44 | 0.9215 | 0.8343 | 127.6 | 1.10 | 1.049 | 0.4482 |

## 判据

- **D1(none)/D2(lowrank) 相对 D0(concat) ARI +≥0.02 且多种子多数据集稳定** → 原始 decoder 的 G 维 mask 旁路是结构捷径，去掉/压缩它能得到更适合聚类的 embedding。
- **D1/D2 ≈ D0** → mask 旁路无害无益，原 scMAE decoder 结构对聚类不是瓶颈。
- **D1/D2 < D0** → mask 旁路对重构/聚类是有益的，不应移除。
- concat 应复现 DEC+std-floor 赢家（Macosko ARI≈0.70）作为 sanity check。
