# Granularity_scMAE_experiments

本文件夹汇集 2026-07（分支 `Granularity`）围绕 scMAE 聚类改良的**探索性实验网络**。原始网络（`../scMAE`、`../NeighborMix_scMAE`、`../CutAware_NeighborMix_scMAE`、`../scMAEs/rank*` 等）**不在此处，未被改动**。

## 目标
非 benchmark。为 scMAE 寻找能带来**跨数据集稳定提升**的结构机制，面向 CS 一区论文。

## 结构
```
Granularity_scMAE_experiments/
├── README.md                    # 本文件（总览 + benchmark 参数说明）
├── EXPERIMENT_LOG.md            # 详细实验流程（背景/思路/踩坑/结果）
├── AdaptiveSwitch_scMAE/        # ★ 赢家：DEC + 每维方差下限
├── AdaptiveGranularity_scMAE/   # 负结果线：自适应粒度聚类目标（含 rank29 解剖 runs/）
├── GatedNeighborMix_scMAE/      # 负结果线：可靠性门控 NeighborMix（含 rank13 忠实复现）
└── ReliableRecon_scMAE/         # 负结果线：局部可靠性精度加权重构
```
每个子网络内有独立 `README.md`（当前模型的流程与思路）。

## 一句话结论
**当前赢家 = `AdaptiveSwitch_scMAE`**，配置 `--var_mode hinge --variance_weight 0.02 --force_gate 1.0`。
机制定位应写成诊断型：在 scMAE+DEC 组合中，DEC 锐化会在细粒度数据上诱发**每维方差坍缩**；VICReg 式 std-floor 是验证该诊断的有效干预项，但 std-floor 本身不 novel。
多 seed（42/2024/3407）：Melanoma ARI 0.648±0.003 / Quake 0.920±0.001 / Macosko 0.576±0.087（最好 0.695）。首个三数据集同时达标的配置。

### 2026-07-02 纯 DEC 普适性检查
`../dec_with_floor/runs/universality_test/` 在 Macosko 上跑了纯 DEC baseline vs DEC+floor：

- 默认早停 (`tol=0.001`)：两组在 epoch 1 收敛，ARI 都是 0.239，基本只评估到 KMeans 初始化后的状态，不可作为 floor 证据。
- 强制训练 (`tol=0.0`, 300 epochs)：baseline 与 floor 逐元素完全相同，ARI 都是 0.457；DEC+floor 日志中 `Floor: 0.000000`，最终 embedding 每维 std 最小值约 2.27，32/32 维有效。
- 结论：当前纯 DEC 实现**没有出现方差坍缩，也没有触发 floor**。因此不能声称“floor 普适救活纯 DEC”；更稳妥的论文 claim 是“scMAE+DEC 的诊断性失败模式 + 多数据集 benchmark + std-floor 作为诊断干预”。

## 运行环境
```
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python   # numpy1.26 / torch2.4 / anndata0.10
# GPU 1-6（禁 0/7）
```
数据：`../scMAEs/benchmark_data/{Melanoma_5K,Quake_10x_Spleen,Macosko}.h5ad`

---

## 为什么当前实验没有用"统一的 benchmark 参数"？（重要说明）

这是**有意的、分阶段的证据管理策略**，来自 `../scMAEs/参考文献/00_scMAE改良方法整理总报告.md` 第 7 条的证据分级：

> legacy-lightweight（灵感筛选）→ independent-full（候选证据）→ formal-verified（多 seed、多数据集、可进论文）。

本轮处于 **independent-full 快筛阶段**，目的是**廉价、快速地在假设空间里做机制探索与证否**，不是产出投稿级数字。具体差异与理由：

1. **只用 3 个数据集、主要 seed 42**：快筛要的是"这个机制方向是否有戏"的信号，不是统计显著性。3 个数据集刻意选成"人格互补"（Quake 易/均衡差、Melanoma 异质、Macosko 细粒度），能最快暴露机制冲突。跑满 15 数据集 × 5 seed 再筛，成本高 20 倍且没必要——大多数假设在第一个数据集就被证否了。

2. **评估口径与"全benchmark"不同**：本 harness = `benchmark_data/*.h5ad` + `scMAE_family` 预处理 + KMeans(known k)。而 `../scMAEs/全benchmark结果.csv` 用的是另一套 3-seed 预处理。**两者数字不可直接比**（例：纯 scMAE 在本 harness 上 Quake ARI≈0.50，在全benchmark 里 0.92）。所以本轮所有候选都**只与同 harness 的参照（rank13/rank29/纯 backbone）比**，这是内部一致的、公平的相对比较。

3. **与已有网络的对比用它们各自的默认超参**（如 rank13 cluster_weight=0.35、rank29 的 6 个损失权重），因为要复现它们"已知的行为"作为锚点，而不是重新调参。忠实复现已由 `GatedNeighborMix_scMAE/runs/repro_check`（Quake 0.9118 = rank13 原文）验证。

4. **消融内部严格控制变量**：同一份代码里改单一开关（`--force_gate`、`--var_mode`、`--variance_weight`、`--reliability_lambda`），保证归因干净——这一点是**统一的**。

**结论**：不统一的是"数据集数量/seed 数量/跨 harness 比较"，这是快筛阶段的成本取舍；统一的是"同 harness 内部、单变量消融"。真正的统一 benchmark（10–15 数据集、5 seed、对标 scDCC/scDeepCluster/scNAME/DESC）属于下一阶段 **formal-verified**，见 EXPERIMENT_LOG.md 末尾投稿路线。赢家配置已做了 3-seed 的初步稳健性验证，是进入 formal 阶段的候选。
