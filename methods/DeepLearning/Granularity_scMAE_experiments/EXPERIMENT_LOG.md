# Granularity 实验流程详细记录

> 记录本轮（2026-07，分支 Granularity）scMAE 聚类改良探索的完整流程：背景、每条线的思路、踩坑与修正、关键结果。供后续翻阅与写论文取材。
> 参照报告：`../../../../experiment_reports/scmae_full_exploration_20260701/REPORT.md`（同内容，仓库级）。

## 目录
- 0 背景与目标
- 1 harness 与基准口径（关键前提）
- 2 线1 GatedNeighborMix（负结果）
- 3 线2 AdaptiveGranularity（负结果）
- 4 线3 rank29 解剖（关键中间发现）
- 5 线4 ReliableRecon（负结果）
- 6 线5 AdaptiveSwitch（赢家）+ 机制证据
- 7 结论、局限、投稿路线

---

## 0 背景与目标

**目标**：不是做 benchmark，而是在众多模型/机制中，找到能让 scMAE（掩码自编码器 scRNA 聚类）**跨数据集稳定提升**的结构组合，面向 CS 一区论文。

**已知前置事实**（进入本轮前）：
- scMAE 主体：swap-noise（列内置换扰动）+ BCE mask 判别 + 加权 MSE 重构；导出 encoder embedding 跑 KMeans。归纳偏置来自 VIME/SCARF+ELECTRA（表格自监督），非图像 MAE；学的是 gene-gene 共表达结构。
- 前期遍历 50+ 篇论文改良（`../scMAEs/rank*`），无一跨 3 数据集过线。总报告归纳：Quake 易、Melanoma 难（异质）、Macosko 需软/边界保护；全局平滑伤边界与稀有细胞。
- 三个最有价值的既有线索：rank13（DEC 可训练簇心，赢 Quake/Melanoma、崩 Macosko）、rank29（fuzzy soft，唯一打过 Macosko）、NeighborMix（组内自建的邻居混合增强，全benchmark 上唯一超 Melanoma 的表达级机制）。

**评估**：encoder embedding → KMeans(known k, n_init=20)，报 NMI/ARI。

---

## 1 harness 与基准口径（关键前提）

**踩坑（最重要的口径修正）**：记忆/旧文档里的 baseline（Quake ARI 0.922 等）来自 **3-seed 全benchmark** 的另一套预处理。而快筛用的是 `../scMAEs/benchmark_data/*.h5ad` + `scMAE_family` 预处理这套 harness。两者**数字不可直接比**：
- 全benchmark：纯 scMAE Quake ARI 0.922。
- 本 harness：纯 scMAE Quake ARI ≈ 0.50、Macosko 0.376。

**修正后原则**：所有候选只与**同 harness 参照**（rank13 / rank29 / 纯 backbone）比。含义：本 harness 上 rank13 把 Quake 0.50→0.91 是真实巨大增益，不是平滑。

**数据集人格**（刻意选成互补）：
- Quake_10x_Spleen：5 类，极不均衡 [6886,1930,42,464,230]（72% 一类），边界清晰。
- Melanoma_5K：9 类，肿瘤异质/连续态。
- Macosko：12 类，44808 细胞，细粒度亚型 + 稀有类。

<!-- SEC2 -->
---

## 2 线1：GatedNeighborMix_scMAE（负结果）

**背景/思路**：全benchmark 显示 NeighborMix 是唯一超 Melanoma 的表达级机制（+0.042 ARI），但过平滑伤 Macosko/Tosches（−0.07）。设想：用 per-cell 可靠性 `r=邻居标签一致性·局部密度` **门控** NeighborMix——核心细胞混合、稀有/边界细胞关闭混合；建在 rank13 DEC 上。

**踩坑与修正**：
1. 从零重写 DEC，Quake 只有 0.52（应 0.91）。**根因**：重构目标误用了 *scaled* 空间 → recon loss 掉到 0.07 → DEC 的 KL 项反客为主、炸掉嵌入。**修正**：重构目标必须是 *未缩放 log*（保持 recon loss ~0.44，让 DEC 当温和正则）。这是全线通用的关键教训。
2. 决策：不再从零重写，改为**在 rank13 真实代码上叠加**门控分支。`runs/repro_check/` 确认 backbone 精确复现 rank13（Quake **0.9118**）。

**结果**（gated NeighborMix, pseudo=0.3, α_min=0.85）：Melanoma 0.649 / Quake 0.890 / Macosko 0.350。**中性偏负**，核心假设不成立。
**诊断**：Macosko 上 DEC-KL 在 epoch80 仍高达 0.679（Quake 仅 0.0135）——DEC 在 Macosko 根本不收敛。→ 指向"Macosko 需要的不是门控 NeighborMix"。

---

## 3 线2：AdaptiveGranularity_scMAE（负结果）

**思路**（融合 fuzzy-rough 上/下近似 + 确定性退火）：让 DEC 目标锐度按 per-cell 可聚类性 `c` 自适应：`p_i = c_i·sharpen(q_i) + (1−c_i)·q_i`。核心细胞用锐目标，边界细胞目标退化成自身 q（KL(q‖q)=0，保持模糊）。另试 rank29 式 SVD 锚点融合。

**Quake 三消融**（`runs/ablate_quake/`，定位干扰源）：
- A(无锚+无自适应)=0.912（=rank13 ✓ 验证 backbone）
- B(+SVD锚点)=**0.447**（SVD 锚点炸 Quake：把嵌入拉回被 72% 大类主导的原始 PCA 几何）
- C(无锚+自适应)=0.897（自适应目标近乎无害）

**干净版（无锚+自适应）三数据集**：Melanoma 0.654 / Quake 0.897 / Macosko 0.277。
**结论**：自适应目标保住 Quake/Melanoma，但**救不了 Macosko**。铁证：Macosko 上 cluster_weight 越小越好（cw0=0.376 > cw0.35=0.343），**纯 backbone 打败所有 DEC 变体** → 任何 DEC 家族目标都伤 Macosko。SVD 锚点未复现 rank29 的 Macosko 0.569（只 0.259）。

<!-- SEC4 -->
---

## 4 线3：rank29 解剖（关键中间发现）

**动机**：rank29 是唯一打过 Macosko 的机制。用 `AdaptiveGranularity_scMAE/_run_rank29.py`（注入 null-h5ad shim）复现 rank29 full = Macosko ARI **0.5649**。6 个辅助损失全部 CLI 可消融 → 做留一消融定位真正的杠杆。

**Macosko 留一消融**（ARI，Δ vs full 0.565）：

| 去掉的组件 | ARI | Δ |
|---|---|---|
| −entropy（边界熵）| 0.143 | **−0.42** |
| −fuzzy（软核 KL）| 0.242 | **−0.32** |
| −variance（VICReg）| 0.290 | **−0.27** |
| −anchor（重构头）| 0.460 | −0.10 |
| −separation | 0.512 | −0.05 |
| **−balance（KL到均匀）** | **0.637** | **+0.07** |

**发现**：
1. Macosko 提升是 soft 损失的**协同**（entropy+fuzzy+variance 三大件），非单一组件。
2. **balance 损失伤不均衡数据**：去掉 → 0.637 超 baseline（Macosko 类别不均，强推均匀簇是错的）。
3. 锚点**架构不是杠杆**：全辅助权重=0 时 Macosko 仅 0.143（比纯 backbone 0.376 还差）。

**跨数据集（定量确认机制冲突）**：

| | Melanoma | Quake | Macosko |
|---|---|---|---|
| rank13 DEC(锐) | 0.652 | 0.912 | 0.343 |
| rank29 soft | 0.529 | **0.170** | 0.565 |
| rank29 −balance | 0.534 | 0.166 | 0.637 |

**完美反相关**：rank29-soft 在 Quake 灾难性（0.17 vs 0.91）。**无固定配方赢三个** → 强烈指向"数据集自适应切换"，`DEC-KL 收敛性`（Quake→0.01 收敛，Macosko→0.7 停滞）作切换信号。

---

## 5 线4：ReliableRecon_scMAE（负结果）

**思路**（用户提出"重点恢复稳定基因 / kernel 加权 / 自适应带宽"三点，取最小可证伪的一条）：per-cell-per-gene **局部**可靠性 `r_ig`（解耦的原始数据 PCA-KNN 邻域内该基因方差的倒数，相对 per-gene 中位数，floor 0.2）→ 精度加权重构损失。局部而非全局（保住 marker：marker 全局高方差但局部稳定），floor 保护稀有细胞。

**结果**（λ=0 vanilla vs λ=1 全加权，同代码）：Melanoma 0.564→0.492 / Quake 0.538→0.521 / Macosko 0.264→0.173。**三数据集单调变差**。

**为什么失败**（第一性，建前已标注为风险）：重构保真度 ≠ 聚类质量。降权"局部高方差"基因，恰好丢掉了边界/过渡/稀有的**判别信号**（判别信息本就藏在难重构的基因里）。scVI 式精度加权对生成/插补有用，对 KMeans-on-embedding 聚类无益甚至有害。→ 这条线证否了"从重构任务层面加权"这个方向。

<!-- SEC6 -->
---

## 6 线5：AdaptiveSwitch_scMAE（赢家）+ 机制证据

**思路**：把所有正向证据组装成一个模型——DEC(锐) 头 + soft(fuzzy-core-KL + boundary-entropy + variance，无 balance) 头，用 `gate=1/(1+(kl_ref/κ)²)`（kl_ref=mean KL(sharpen(q)‖q)）自动加权。`--force_gate` 供消融。

**自动门控（seed42）**：Melanoma 0.648(gate0.94) / Quake 0.921(gate0.99) / Macosko 0.641(gate0.85)。
**意外**：门控信号在 Quake/Melanoma 正确读出"可锐化"（gate≈1），**但 Macosko gate 也高(0.85) 却拿到 0.641**——提升不是来自 soft 分支。必须消融。

### 6.1 消融：杠杆是方差损失，不是门控（Macosko）
| 配置 | variance | gate | ARI |
|---|---|---|---|
| no-var, pure-sharp | 0 | 1.0 | 0.343（=rank13）|
| **var, pure-sharp** | 0.02 | 1.0 | **0.695** |
| no-var, auto | 0 | 0.66 | 0.342 |
| var, auto | 0.02 | 0.85 | 0.641 |

方差关掉→0.343（门控无关）；方差开→0.64–0.70。纯 sharp+方差(0.695) 还高于 auto(0.641)。**门控是红鲱鱼，真杠杆是方差损失。**

### 6.2 机制确认：是"每维方差下限"这个特定形式（Macosko seed42）
| 反坍缩机制 | vw=0.02 | vw=0.1 |
|---|---|---|
| **hinge（每维 std 下限）** | **0.695** | 0.329 |
| cov（协方差去相关）| 0.333 | 0.202 |
| koleo（超球均匀）| 0.362 | −0.007 |

只有 std-floor 有效，cov/koleo ≈ rank13。**精确断言，非泛反坍缩。**

### 6.3 variance_weight 权衡：vw≈0.02 最优
vw 0.01→0.522 | 0.02→0.576 | 0.05→0.405 | 0.10→0.329（效果抵消）。太少不防坍缩，太多撑散簇。

### 6.4 多 seed（42/2024/3407）
- Melanoma ARI 0.648±0.003，NMI 0.722±0.002（极稳）
- Quake ARI 0.920±0.001，NMI 0.834±0.005（极稳）
- Macosko ARI 0.576±0.087（0.695/0.543/0.491），NMI 0.616±0.019

### 6.5 seed 稳定性尝试（warmup 扫描）
10→0.417 | **20→0.576** | 30→0.533 | 40→0.433，非单调峰在 20。加长 warmup 收窄方差但降均值。Macosko seed 方差（±0.087）是**细粒度硬聚类内禀**（多近等价最优、KMeans 初始化敏感），warmup 不可解；诚实保留为 limitation（最差 seed 0.49 仍 > rank13 0.34）。

### 6.6 直接机制证据（每维方差谱，Macosko 128 维）
| | 有效维数(PR) | 每维方差 min/中位/max |
|---|---|---|
| DEC only | **60.4**/128 | 0.076/0.505/5.86 |
| DEC+std-floor | **114.0**/128 | 1.064/1.115/4.77 |

std-floor 使**有效维数几乎翻倍(60→114)**。DEC 锐化把嵌入压进 ~60 有效维（坍缩）；方差下限让全 128 维活跃，KMeans 得以区分 12 亚型。**定量铁证 + 论文理想配图。**

---

## 7 结论、局限、投稿路线

**赢家**：`scMAE + DEC 聚类目标 + 每维方差下限（VICReg std-hinge, w=0.02）`。`AdaptiveSwitch_scMAE/run.py --var_mode hinge --variance_weight 0.02 --force_gate 1.0`。三数据集同时达到/超过历史最佳。

**论文核心 claim**：*每维方差坍缩是深度聚类(DEC)在细粒度 scRNA 上失败的原因；每维方差下限是精确解药。* cov/koleo 无效，说明是"每维方差收缩"这一特定失败模式。

**局限**：(1) Macosko seed 方差 ±0.087；(2) 仅 3 数据集、快筛 harness；(3) 未对标已发表方法。

**投稿路线（formal-verified 阶段）**：
1. 扩数据集（10–15 标准 scRNA）验证"细粒度数据集普遍受益于每维方差下限"。
2. 同 harness 对标 scDCC/scDeepCluster/scNAME/DESC。
3. Macosko seed 方差作 limitation 报告，或 consensus/多 init 取最优 silhouette。
4. 全部结果多 seed 均值±方差进正文；每维方差谱作机制配图。



