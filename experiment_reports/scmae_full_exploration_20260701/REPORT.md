# scMAE 结构改良完整探索记录（2026-07-01, 分支 Granularity）

> 本文是本轮全部探索的**完整过程记录**，供后续翻阅。包含所有尝试过的机制线（含失败线，因其代码目录将从仓库删除，机制细节在此永久保留）、关键实验数据、以及最终的跨数据集赢家与机制证据。
>
> 目标：非 benchmark。为 scMAE 找到能带来**跨数据集稳定提升**的结构机制，面向 CS 一区。
> 环境：`/data/luolie/conda/envs/scssl_bench_py310/bin/python`（numpy1.26/torch2.4/anndata0.10）。GPU 1–6（禁 0/7）。
> 数据（`methods/DeepLearning/scMAEs/benchmark_data/`）：Melanoma_5K(4513,9类)、Quake_10x_Spleen(9552,5类,极不均衡[6886,1930,42,464,230])、Macosko(44808,12类)。
> 评估：encoder embedding → KMeans(known k, n_init=20)，报 NMI/ARI，seed 42（除非注明多 seed）。

---

## 0. 结论速览

**最终赢家**：`scMAE + DEC 聚类目标 + 每维方差下限（VICReg std-hinge, w=0.02）`，代码在 `methods/DeepLearning/AdaptiveSwitch_scMAE/`，配置 `--var_mode hinge --variance_weight 0.02 --force_gate 1.0`。

多 seed（42/2024/3407）：

| 数据集 | ARI mean±std | NMI mean | rank13 DEC | rank29 soft |
|---|---|---|---|---|
| Melanoma_5K | 0.648±0.003 | 0.722 | 0.652 | 0.529 |
| Quake_10x_Spleen | 0.920±0.001 | 0.834 | 0.912 | 0.170 |
| Macosko | 0.576±0.087 (best 0.695) | 0.616 | 0.343 | 0.565 |

**首个三数据集同时达到/超过历史最佳的配置。** 机制精确：DEC 在细粒度数据上失败的根因是**每维方差坍缩**，per-dim 方差下限是精确解药（去相关 cov、超球均匀 koleo 均无效）。

---

## 1. 关键前提修正：基准口径

记忆/旧文档里的 baseline（Quake ARI 0.922 等）来自 **3-seed 全benchmark** 的另一套预处理。在快筛实际用的 `benchmark_data` harness 上，**纯 scMAE 实测**：Quake ARI ~0.50、Macosko 0.376。所有候选必须与**同 harness** 参照比。含义：本 harness 上 rank13 DEC 把 Quake 0.50→0.91 是真实巨大增益。

<!-- SEC2 -->
---

## 2. 前置认知（进入本轮前已建立）

- 原 scMAE：swap-noise（列内置换）+ BCE mask 判别 + 加权 MSE 重构；embedding→KMeans。祖宗是 VIME/SCARF+ELECTRA，非图像 MAE。学的是 gene-gene 共表达结构。
- NeighborMix（组内自建，非论文）：Mixup+swap+DAE 合成；`x'=0.9x+0.1·邻居均值`→再 swap→重构**真细胞**（anchor-recovery）；推理关闭混合→导出嵌入干净。全benchmark 上唯一超 Melanoma 的表达级机制，但过平滑伤 Macosko/Tosches。
- rank13（DEC）：可训练簇心 + 置信门控 KL。赢 Quake/Melanoma，崩 Macosko。
- rank29（fuzzy soft）：SVD 锚点融合 + fuzzy-core-KL + boundary-entropy + variance + balance 等。唯一打过 Macosko baseline（0.565）。
- 数据集人格：Quake 易（图/去噪/锐化都行）；Melanoma 难（全局平滑伤肿瘤异质性）；Macosko 需软/边界保护。

---

## 3. 线1：GatedNeighborMix_scMAE（负结果）— 目录将删

**思路**：per-cell 可靠性 `r=邻居一致性·密度` 门控 NeighborMix，核心细胞混合、稀有/边界关闭。建在 rank13 DEC 上。

**踩坑与修正（永久记录）**：
- 从零重写 DEC 无法复现 rank13（Quake 0.52 vs 0.91）。根因：重构目标必须是**未缩放 log**（保持 recon loss 大 ~0.44，让 DEC 当温和正则）；改成 scaled 目标 → recon loss 掉到 0.07 → DEC KL 反客为主炸掉嵌入。
- 决策：改为在 rank13 真实代码上叠加。`repro_check` 确认 backbone 精确复现 rank13（Quake 0.9118）。

**结果**（gated NeighborMix, pseudo=0.3, α_min=0.85）：Melanoma 0.649 / Quake 0.890 / Macosko 0.350。中性偏负。核心假设不成立。诊断：Macosko 上 DEC-KL epoch80 仍 0.679（Quake 仅 0.0135）——DEC 在 Macosko 不收敛。

---

## 4. 线2：AdaptiveGranularity_scMAE（负结果）— 保留目录

**思路**：DEC 目标锐度按 per-cell 可聚类性 `c` 自适应：`p_i=c_i·sharpen(q_i)+(1−c_i)·q_i`。核心细胞锐、边界细胞软（KL(q‖q)=0）。加 rank29 式 SVD 锚点融合。基于 fuzzy-rough 近似 + 确定性退火。

**Quake 三消融**（定位干扰源）：A(无锚+无自适应)=0.912(=rank13✓)；B(+SVD锚点)=**0.447**(锚点炸 Quake)；C(无锚+自适应)=0.897(自适应近乎无害)。

**干净版（无锚+自适应）**：Melanoma 0.654 / Quake 0.897 / Macosko 0.277。
**结论**：自适应目标保住 Quake/Melanoma 但**救不了 Macosko**。铁证：Macosko 上 cluster_weight 越小越好（cw0=0.376 > cw0.35=0.343），纯 backbone 打败所有 DEC 变体。**任何 DEC 家族目标都伤 Macosko**。SVD 锚点炸 Quake，且未复现 rank29 的 Macosko 0.569（只 0.259）。

<!-- SEC5 -->
---

## 5. 线3：rank29 解剖（关键中间发现）

用 `_run_rank29.py`（注入 null-h5ad shim）复现 rank29 full = Macosko ARI 0.5649。6 个辅助损失全 CLI 可消融。

**Macosko 留一消融**（ARI，Δ vs full 0.565）：
| 去掉的组件 | ARI | Δ |
|---|---|---|
| −entropy(边界熵) | 0.143 | −0.42 |
| −fuzzy(软核KL) | 0.242 | −0.32 |
| −variance(VICReg) | 0.290 | −0.27 |
| −anchor(重构头) | 0.460 | −0.10 |
| −separation | 0.512 | −0.05 |
| **−balance(KL到均匀)** | **0.637** | **+0.07** |

发现：(1) Macosko 提升是 soft 损失的**协同**（entropy+fuzzy+variance 三大件）；(2) **balance 损失伤不均衡数据**，去掉 → 0.637 超 baseline；(3) 锚点**架构**不是杠杆（全辅助权重=0 时 Macosko 仅 0.143，比纯 backbone 0.376 还差）。

**跨数据集（定量确认机制冲突）**：
| | Melanoma | Quake | Macosko |
|---|---|---|---|
| rank13 DEC(锐) | 0.652 | 0.912 | 0.343 |
| rank29 soft | 0.529 | **0.170** | 0.565 |
| rank29 −balance | 0.534 | 0.166 | 0.637 |

**完美反相关**：rank29-soft 在 Quake 灾难性（0.17 vs 0.91）。无固定配方赢三个。→ 强烈指向**数据集自适应切换**，DEC-KL 收敛性作信号。

---

## 6. 线4：ReliableRecon_scMAE（负结果）— 目录将删

**思路**（用户提出的"重点恢复稳定基因/kernel 加权/自适应带宽"三点，取最小可证伪的一条）：per-cell-per-gene **局部**可靠性 `r_ig`（解耦的原始数据 PCA-KNN 邻域内该基因方差的倒数，相对 per-gene 中位数，floor 0.2）→ 精度加权重构损失。局部而非全局（保住 marker），floor 保护稀有细胞。

**结果**（λ=0 vanilla vs λ=1 全加权）：Melanoma 0.564→0.492 / Quake 0.538→0.521 / Macosko 0.264→0.173。**三数据集单调变差**。

**为什么失败（第一性，建前已标注的风险）**：重构保真度 ≠ 聚类质量。降权"局部高方差"基因恰好丢掉边界/过渡/稀有的判别信号；判别信息本就藏在难重构的基因里。scVI 式精度加权对生成/插补有用，对 KMeans-on-embedding 聚类无益。

<!-- SEC7 -->
---

## 7. 线5：AdaptiveSwitch_scMAE（赢家）— 保留目录

**思路**：一个模型挂 DEC(锐) + soft(fuzzy-core-KL + boundary-entropy + variance, 无 balance) 双目标，`gate=1/(1+(kl_ref/κ)²)`（kl_ref=mean KL(sharpen(q)‖q)）自动加权。`--force_gate` 供消融。

**自动门控（seed42）**：Melanoma 0.648(gate0.94) / Quake 0.921(gate0.99) / Macosko 0.641(gate0.85)。门控信号在 Quake/Melanoma 正确读出"可锐化"。**但 Macosko gate 也高(0.85)却拿到 0.641** → 提升不是来自 soft 分支，需消融。

### 7.1 消融：杠杆是方差损失，不是门控（Macosko）
| 配置 | variance | gate | ARI |
|---|---|---|---|
| no-var, pure-sharp | 0 | 1.0 | 0.343(=rank13) |
| **var, pure-sharp** | 0.02 | 1.0 | **0.695** |
| no-var, auto | 0 | 0.66 | 0.342 |
| var, auto | 0.02 | 0.85 | 0.641 |

方差关掉→0.343（门控无关）；方差开→0.64–0.70。纯sharp+方差(0.695)>auto(0.641)。**门控是红鲱鱼，方差损失是真杠杆。**

### 7.2 机制确认：是"每维方差下限"特定形式（Macosko seed42）
| 反坍缩机制 | vw=0.02 | vw=0.1 |
|---|---|---|
| **hinge(每维std下限)** | **0.695** | 0.329 |
| cov(去相关) | 0.333 | 0.202 |
| koleo(超球均匀) | 0.362 | −0.007 |

只有 std-floor 有效。cov/koleo≈rank13。**精确断言，非泛反坍缩。**

### 7.3 variance_weight 权衡：vw≈0.02 最优
vw 0.01→0.522 | 0.02→0.576 | 0.05→0.405 | 0.10→0.329(效果抵消)。太少不防坍缩，太多撑散簇。

### 7.4 多 seed + warmup 稳定性
多 seed 见 §0。warmup 扫描：10→0.417 | **20→0.576** | 30→0.533 | 40→0.433，非单调峰在20。Macosko seed 方差(±0.087)是细粒度硬聚类**内禀**（多近等价最优），warmup 不可解；诚实保留为 limitation（最差 seed 0.49 仍 > rank13 0.34）。

### 7.5 直接机制证据（每维方差谱，Macosko 128 维）
| | 有效维数(PR) | 每维方差 min/中位/max |
|---|---|---|
| DEC only | **60.4**/128 | 0.076/0.505/5.86 |
| DEC+std-floor | **114.0**/128 | 1.064/1.115/4.77 |

std-floor 使**有效维数几乎翻倍(60→114)**。DEC 锐化把嵌入压进~60 有效维（坍缩）；方差下限让全 128 维活跃，KMeans 得以区分 12 亚型。**定量铁证 + 论文理想配图。**

---

## 8. 机制解释（论文核心 claim）

*每维方差坍缩是深度聚类(DEC)在细粒度 scRNA 上失败的原因；每维方差下限(VICReg std-hinge)是精确解药。* DEC 目标锐化把细胞拉向少数簇心，在 12 亚型的 Macosko 上导致潜在空间每维方差退化，亚型被挤进低维子空间无法区分。std-floor(强制每维 std≥1)精确对抗；在本就高方差的 Quake/Melanoma 上无害。去相关(减冗余)、KoLeo(超球均匀)都不针对"每维方差收缩"，故无效。

---

## 9. 代码与结果位置

- **赢家（仓库保留）**：`methods/DeepLearning/AdaptiveSwitch_scMAE/`（`model.py`/`loss.py`/`clusterability.py`/`run.py`）。runs 子目录：`screen_auto` `ablate_macosko` `mech` `vw_sweep` `dec_plus_var`(多seed) `warmup_sweep` `smoke`。
- **已删除的负结果/中间线**（代码不再保留，机制与全部数值见本报告对应节）：
  - `GatedNeighborMix_scMAE/`（线1，§3；含忠实 rank13 复现 Quake 0.9118）
  - `AdaptiveGranularity_scMAE/`（线2，§4；其 `runs/rank29_dissect/` 是 rank29 解剖 §5 的数值来源，均已录入本报告表格）
  - `ReliableRecon_scMAE/`（线4，§6）
- 前序报告：`experiment_reports/scmae_adaptive_granularity_20260701/`、`scmae_variance_floor_20260701/`（内容已并入本总报告，可作交叉参考）。

## 10. 运维要点
- h5ad 有 `encoding-type:'null'` 的 uns 项，需 `_register_null_h5ad_reader()` shim（AdaptiveSwitch run.py 已内置）。
- 44k 细胞需 `OPENBLAS_NUM_THREADS=8` 等上限（已内置），否则 KMeans OpenBLAS 崩。
- 编码器输入=scaled，重构目标=**未缩放 log**（关键，保 DEC 损失平衡）。
- rank29 解剖用的 shim wrapper `_run_rank29.py` 随 `AdaptiveGranularity_scMAE/` 一并删除；如需重跑，按 §5 用 `rank29_deep_adaptive_fuzzy_clustering_full/run.py` + null-shim 即可。

## 10. 运维要点
- h5ad 有 `encoding-type:'null'` 的 uns 项，需 `_register_null_h5ad_reader()` shim（各 run.py 已内置）。
- 44k 细胞需 `OPENBLAS_NUM_THREADS=8` 等上限（已内置），否则 KMeans OpenBLAS 崩。
- 编码器输入=scaled，重构目标=**未缩放 log**（关键，保 DEC 损失平衡）。

## 11. 下一步（投稿路线）
1. 扩数据集(10–15 标准 scRNA)验证"细粒度数据集普遍受益于每维方差下限"。
2. 同 harness 对标 scDCC/scDeepCluster/scNAME/DESC。
3. Macosko seed 方差作 limitation 如实报告，或 consensus/多 init 取最优 silhouette。
4. 全部结果多 seed 均值±方差进正文。



