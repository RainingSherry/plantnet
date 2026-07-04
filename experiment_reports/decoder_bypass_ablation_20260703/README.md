# Decoder-bypass 消融（Phase 1）

路线图 `experimental_retired_models/scmae_structural_pivot_20260703/ROADMAP.md` 的
Phase 1。这是转向「从 scMAE 内部结构动手」后的第一个、也是最便宜的实验。

## 假设

原始 scMAE 的 decoder 输入是 `[latent, mask_logits]`，即把 G 维的 mask 预测**整个**
喂给 decoder。这对重构有利（decoder 知道哪里被 corrupt）。但我们的目标是聚类，不是
重构。假设：

```text
这个 G 维 mask 旁路 = 一条结构捷径，
让 encoder 不必在 latent 里编码「哪里被扰动 / 该细胞的表达结构」，
因此 latent 对聚类可能变弱。
```

## 三个 arm（唯一变量 = decoder 如何接收 mask）

| arm | decoder | 含义 |
|---|---|---|
| **D0 `concat`** | `decoder([latent, mask_logits])` | 原始 scMAE = 复现 DEC+std-floor 赢家（对照）|
| **D1 `none`** | `decoder(latent)` | mask 只进 BCE 判别 loss，不进 decoder |
| **D2 `lowrank`** | `decoder([latent, pool(mask_logits)→r])` | mask 只以 r=16 维摘要进 decoder |

其它一切与赢家**完全一致**：encoder / mask_predictor / DEC centers / std-floor
(`variance_weight=0.02`) / `force_gate=1.0` / mask_prob=0.4 / 80 epochs。
`mask_logits` 在三种模式下都保留（BCE mask 判别器 loss 不变），只改它是否/如何进
decoder。这样能干净归因到「mask 旁路」这一个变量。

## 矩阵与运行

3 decoder_mode × 3 数据集 (Macosko k=12 / Melanoma_5K k=9 / Quake_10x_Spleen k=5)
× 3 种子 (42/43/44) = 27 runs。并行调度到 GPU 1-6（禁用 0/7）：

```bash
python experiment_reports/decoder_bypass_ablation_20260703/dispatch.py
python experiment_reports/decoder_bypass_ablation_20260703/summarize.py
```

`dispatch.py` 每 GPU 同时 1 个 run，已完成的自动跳过（断点续跑）。

## 判据（见 SUMMARY.md 自动 delta 表）

```text
D1/D2 相对 D0  ARI +≥0.02 且多种子多数据集稳定 -> mask 旁路是结构捷径，去掉更好
D1/D2 ≈ D0                                      -> 旁路无害无益，decoder 非瓶颈
D1/D2 < D0                                      -> 旁路有益，不应移除
```

sanity check：`concat` 在 Macosko 应复现赢家 ARI≈0.70。

## 防泄露

不按测试 ARI 选任何超参；label 只在最终评测用一次；三数据集固定协议、多种子，
避免单数据集调参（Macosko 种子 sd≈±0.087）。
