# Conditional-shuffle 消融 —— 结论（2026-07-03）

45/45 runs 完成。骨干 = scMAE + DEC + std-floor 赢家；唯一变量 = corruption 供体来源。

## 结果（ARI mean±sd，3 seeds；括号为相对同数据集 zero-mask 赢家的 delta）

| arm | Macosko (赢家0.702) | Melanoma (0.648) | Quake (0.920) |
|---|---|---|---|
| zero（复现赢家） | 0.702±0.004 | 0.648±0.002 | 0.920±0.001 |
| swap_global (S0) | 0.299±0.028 (**−0.40**) | 0.655±0.005 (+0.01) | 0.923±0.001 (flat) |
| **swap_lib (S1)** | **0.867±0.002 (+0.16)** | 0.602±0.076 (−0.05) | 0.915±0.007 (flat) |
| **swap_ndet (S2)** | **0.864±0.001 (+0.16)** | 0.658±0.006 (+0.01) | 0.922±0.002 (flat) |
| swap_zerolib (S3) | 0.759±0.147 (双峰) | 0.660±0.003 (+0.01) | 0.917±0.002 (flat) |

所有 swap arm 的 `eff_change`≈0.023（z-scaled 稀疏数据上 swap 只改变~2.3%位置），
zero-mask 的 eff_change=0.40。arm 之间 eff_change 相同 → 差异纯粹来自**供体来自谁**。

## 四个结论

1. **原始 scMAE 的全局 swap-noise 在细粒度 Macosko 上是灾难**（0.30 vs zero-mask 0.70）。
   这是对原 scMAE corruption 设计的一个干净负结果：细粒度聚类需要 zero-masking，
   而非全局 swap。（Melanoma/Quake 上 S0 无害，因为它们不需要精细区分。）

2. **nuisance 匹配 swap（S1/S2）在 Macosko 上给出 +0.16 的稳定巨大增益**（三种子 sd<0.003），
   且在**两种独立聚类头下都成立**（KMeans-knownK 0.870 vs 0.696；label-free Leiden
   0.894 vs 0.780）→ 是真实的 embedding 改善，不是 KMeans-known-K 的评测假象。

3. **但完全不泛化**：Melanoma 和 Quake 上所有 swap arm 都在种子噪声内持平。
   这是一个 Macosko 特有现象。

4. **显然的机制解释被证伪**：NMI(nuisance箱,标签) 在 Macosko(lib0.11/ndet0.13) 与
   Melanoma(lib0.10/ndet0.14) 几乎相同，但只有 Macosko 受益。所以"同箱供体≈同类型
   供体"不成立（否则 Melanoma 也该涨）。Macosko 为何独特 = **未知**。

## 纪律性判定

按路线图 Section 4 的多数据集纪律（禁止单数据集调参；Macosko 种子 sd±0.087），
**这不是通用赢家，不作为论文头条。稳健跨数据集主线仍是 DEC + std-floor。**

## 植物数据检验（Phase 2b，runs_plant/）—— 决定性负结果

项目真实目标是植物 scRNA。检验 matched-swap 是否迁移：

| arm | SRP182008 k=15 (zero=0.390) | CRA002977_1 k=7 (zero=0.683) |
|---|---|---|
| swap_global | 0.397 (+0.01) | 0.508 (**−0.175**) |
| swap_lib | 0.389 (−0.00) | 0.537 (**−0.146**) |
| swap_ndet | 0.384 (−0.01) | 0.544 (**−0.139**) |

**Macosko 的 +0.16 完全不迁移。** SRP182008 全持平（像 Melanoma）；CRA002977_1 上
所有 swap（含 matched）比 zero-mask 差 ~0.15。跨 5 数据集（Macosko/Melanoma/Quake/
SRP182008/CRA002977_1）总结：**只有 Macosko 从 matched-swap 受益；其余持平或变差。
zero-mask + DEC + std-floor 在所有数据集（含两个植物集）都是稳健选择。**

## 全局结论（"攻击 corruption" 分支闭合）

1. corruption 设计（zero-mask / global-swap / matched-swap）对结果有巨大、**数据集
   依赖**的影响，但**没有任何一种 corruption 在所有数据集上都胜过 zero-mask**。
2. "matched > global" 的次序在多个数据集一致（限制供体池总比全局好），但"swap 是否
   胜过 zero-mask"是数据集特异的——只有 Macosko 上 matched-swap 胜出。
3. Macosko 为何独特仍未知（nuisance≈type 已被 NMI 检验否证）。
4. **对项目（植物）无直接价值**：植物数据上 matched-swap 不帮忙甚至有害。
5. 稳健跨数据集主线不变：**scMAE + DEC + 每维 std-floor（zero-mask）**。

## 保留的干净科学产出

- 原始 scMAE 全局 swap-noise 在细粒度数据（Macosko/CRA002977_1）上劣于 zero-mask，
  可作为"corruption 设计对细粒度聚类的影响"的干净消融证据。
- Macosko 的 ~0.87 盆地是真实的（两种独立扰动可达、且 label-free Leiden 头验证），
  是一个孤立但真实的现象，机制未解——可作为 future work 或诊断论文的一个观察点。
