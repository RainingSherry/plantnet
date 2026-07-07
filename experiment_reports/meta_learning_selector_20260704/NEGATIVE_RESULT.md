# 负结果章：scRNA 聚类的算法选择 / regime routing（2026-07-04）

严格重做（用户辩驳后修正了三处方法学问题）。这份是论文"适用边界 / 为何不做算法选择"一章的
方法与结论。核心一句话：

> **portfolio 的 oracle 空间真实且显著（Common-18 强池 VBS−SBS=+0.073，CI 不跨 0），
> 但当前 benchmark 的任务数（n≈18–23）与方法覆盖度，不足以让任何可兑现的 label-free
> selector 兑现这个空间——跨 3 类模型 × 2 套特征 × LODO，最好仅捕获 +14% 且不显著。**

## 1. 相对第一版的三处修正（认账）

1. 第一版"kill-shot"只在 **PCA/VF 两方法池**（oracle 仅 +0.0225）上测 1-NN → 推不出
   "17 方法 portfolio 不工作"。**测错了池子。**
2. 第一版 coverage≥18 段出现 **VBS<SBS**（数学上不可能）→ 评估 universe 不一致的 bug。
3. 第一版用全体 z-score（泄露）、true-K 当特征、只 1-NN + 12 特征 → 只是 sanity check。

修正 = 严格矩形矩阵（VBS≥SBS 断言通过）+ 多 selector + 训练集统计标准化 + known-K/
deployment 分版 + paired bootstrap CI。

## 2. 严格矩形矩阵（评估对象）

| 池 | 方法 × 数据集 | SBS | VBS | gap | 门槛(≥0.04) |
|---|---|---|---|---|---|
| Full-23 small | 2 × 23 (VarFloor/PCA) | 0.5938 | 0.6162 | +0.0225 | ✗ |
| **Common-18 strong** | 10 × 18 | scMAE 0.7030 | 0.7760 | **+0.0729** | ✓ |

Common-18 强池方法 = DEC/Leiden/Louvain/SC3/scDCC/scMAE/NeighborMix/scMAE-w/o/
VarFloor/PCA（都覆盖同一 18 集，严格矩形）。冠军散在 7 方法：NeighborMix5/Leiden3/
PCA3/scDCC3/VarFloor2/SC3/scMAE。oracle(VBS) 的 paired bootstrap CI [+0.034,+0.121]
**不跨 0 → headroom 真实显著**。

## 3. selector 套件结果（Common-18，LODO）

判定门槛（用户设定）：捕获 ≥30–50% oracle gap 且 paired bootstrap CI 不跨 0。

| selector | 特征版 | 捕获 gap | CI 跨 0 |
|---|---|---|---|
| ORACLE(上限) | — | 100% | 不跨(显著) |
| always-SBS | — | 0% | — |
| 1-NN | deployment(去真K) | −28% | 跨 |
| ridge-perf-regression | deployment | +1% | 跨 |
| decision-tree(d3) | deployment | −143% | 显著更差 |
| 1-NN | known-K | **+14%(最好)** | **跨(不显著)** |
| ridge-perf-regression | known-K | −8% | 跨 |
| decision-tree(d3) | known-K | −129% | 显著更差 |

**无一达标。** 最好的 1-NN(known-K) 仅 +14% 且不显著；多个 selector 显著差于 always-SBS。

## 4. 结论与根因

- oracle portfolio 空间真实显著（+0.073）——"没有单一 SOTA、不同数据需不同方法"成立。
- 但**可兑现的 label-free selector 无法兑现它**。根因：**benchmark 任务数太少**——regime
  routing 需几百 task（AutoML 经验），这里 18 集 / 7 赢家 / LODO 17 训练任务 → 学不动 +
  严重过拟合。**不是想法错，是数据稀缺卡死。**

## 5. 对主线的价值（不与方案三冲突，反而喂它）

这是一个**论文级的高质量负结果**，构成 benchmark 论文的"适用边界 / regime routing 可行性"
一章：量化说明"为何该转向可靠性 + GPU 加速，而非算法选择或造新方法"。

## 6. 未尽之处（诚实标注，预期不翻盘）

1. 未补跑 23×17 完整矩阵（可把 universe 从 18 抬回 23，但 n≈20 根本瓶颈不变）。
2. 未做 leave-one-tissue/platform-out（更严）。
3. 未做不确定性 router（conformal abstention / top-k）——LODO 下 3 类已失败，abstention
   仅降覆盖、top-k 换指标，翻盘概率低。

脚本：`compute_metafeatures.py`(23集 label-free meta-feature)、`proper_selector.py`
(严格矩形 + 多 selector + 双 bootstrap CI)。数据：`metafeatures.csv`。
