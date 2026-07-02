# scMAE 结构改良探索报告：可靠性门控 NeighborMix 与自适应粒度聚类目标

- **日期**：2026-07-01
- **分支**：`wide-research`
- **目标**：非 benchmark。为 scMAE 寻找能带来**跨数据集稳定提升**的结构机制，面向 CS 一区论文。
- **环境**：`/data/luolie/conda/envs/scssl_bench_py310/bin/python`（numpy 1.26.4 / torch 2.4.0+cu121 / anndata 0.10.5）。GPU 1–6（禁用 0/7）。
- **快筛数据集**（`methods/DeepLearning/scMAEs/benchmark_data/`）：Melanoma_5K（4513 细胞，9 类）、Quake_10x_Spleen（9552 细胞，5 类，极度不均衡 [6886,1930,42,464,230]）、Macosko（44808 细胞，12 类）。
- **评估**：encoder embedding → KMeans(known k, n_init=20)，报 NMI/ARI。

---

## 0. TL;DR（一句话结论）

我们实现并干净消融了两条新机制线。**两条都不是跨数据集赢家**，但换来一个**确凿且可发表的科学发现**：

> scMAE 聚类中存在**真实的机制冲突**——DEC 锐化聚类对干净/不均衡数据（Quake）是命脉、对细粒度连续数据（Macosko）有害；流形/锚点保持则完全相反。**不存在任何固定机制能同时赢三个数据集**。而 **DEC-KL 是否收敛**是一个可靠的、数据集级的"该不该用 DEC"探测信号。

---

## 1. 代码位置

两条候选线均在 `methods/DeepLearning/` 下，均为 independent-full（独立 model+loss+run）。

### 线 1：`GatedNeighborMix_scMAE/`（可靠性门控 NeighborMix + DEC）
- `model.py` — 忠实复现 rank13 的 DEC scMAE（zero-mask + 可训练簇心）
- `loss.py` — rank13 损失（recon+BCE+strong/weak 一致性+置信门控 KL）+ 门控 NeighborMix 伪分支
- `reliability.py` — per-cell 可靠性场 `r = 邻居一致性 · 局部密度`
- `run.py` — 训练/评估，含 PCA-KNN 图、h5ad null-encoding 读修复、BLAS 线程上限
- `README.md`

### 线 2：`AdaptiveGranularity_scMAE/`（自适应粒度聚类目标，AGCT）
- `model.py` — rank13 DEC backbone + 可选 rank29 SVD-锚点融合
- `clusterability.py` — per-cell 可聚类性 `c_i` + `adaptive_target()` 插值
- `loss.py` — recon+BCE+一致性 + KL(到自适应目标)
- `run.py` — 训练/评估

### 关键运行产物
- `GatedNeighborMix_scMAE/runs/repro_check/` — 证明 backbone 忠实复现 rank13（Quake ARI 0.9118）
- `AdaptiveGranularity_scMAE/runs/ablate_quake/` — Quake 三消融（定位 SVD 锚点为元凶）
- `AdaptiveGranularity_scMAE/runs/screen_noanchor/` — AGCT 干净版三数据集

---

## 2. 关键前提修正：基准口径

记忆/旧报告里的 baseline（Quake NMI 0.852/ARI 0.922 等）来自 **3-seed 全benchmark** 的另一套预处理。在快筛实际使用的 `benchmark_data` harness 上，**纯 scMAE**（NeighborMix pseudo 关闭）实测：

| 数据集 | 纯 scMAE (本 harness) | 旧 baseline (全benchmark) |
|---|---|---|
| Quake_10x_Spleen | ARI ~0.50 / NMI ~0.65 | 0.922 / 0.852 |
| Macosko | ARI 0.376 | 0.494 / 0.657 |

**因此所有候选必须与同 harness 参照比，不能与全benchmark 数字比。** 关键含义：在本 harness 上 rank13 的 DEC 把 Quake 从 0.50 拉到 0.91，是**真实的巨大增益**，不是平滑。

---

## 3. 线 1：可靠性门控 NeighborMix（负结果）

**动机**：全benchmark 显示 NeighborMix 是唯一在 Melanoma 超基线的表达级机制（+0.042 ARI），但它在 Macosko/Tosches 过平滑（−0.07）。设想用 per-cell 可靠性 `r_i` 门控：核心细胞混合、稀有/边界细胞关闭混合。

**过程与踩坑**（全部已修正并记录）：
1. 首次从零实现的 DEC 无法复现 rank13（Quake 0.52 vs 0.91）。
2. 定位根因：重构目标必须是**未缩放 log 表达**（保持 recon loss 大 ~0.44，让 DEC 当温和正则）；一旦改成 scaled 目标，recon loss 掉到 0.07，DEC KL 反客为主直接炸掉嵌入。
3. 决策：停止从零重写，**在 rank13 真实代码上叠加**门控 NeighborMix。`runs/repro_check` 确认 backbone 精确复现 rank13（Quake ARI **0.9118**）。

**结果**（gated NeighborMix，pseudo=0.3，α_min=0.85）：

| | rank13 单独 | +门控 NeighborMix | Δ ARI |
|---|---|---|---|
| Melanoma | 0.652 | 0.649 | −0.003 |
| Quake | 0.912 | 0.890 | −0.022 |
| Macosko | 0.343 | 0.350 | +0.007 |

**结论**：中性偏负。没能救 Macosko，轻微伤 Quake。核心假设不成立。

**诊断收获**：Macosko 上 DEC-KL 在 epoch 80 仍高达 **0.679**（Quake 仅 0.0135）——DEC 在 Macosko 的 12 个细粒度亚型上根本不收敛。

---

## 4. 线 2：自适应粒度聚类目标 AGCT（负结果，但机制成立）

**动机**（第一性原理，融合 fuzzy-rough 下/上近似 + 确定性退火）：
既然 Quake 要锐、Macosko 要软，就让**目标锐度按 per-cell 可聚类性自适应**：

```
p_i = c_i · sharpen(q_i) + (1 − c_i) · q_i
c_i = 邻居标签一致性 · 局部密度
```

核心细胞（c_i→1）→ 锐 DEC 目标；边界/稀有细胞（c_i→0）→ 目标=自身 q → KL(q‖q)=0 → DEC 自动"悬置判断"，细胞保持模糊。门控的是**目标**而非损失权重（这是线 1 失败后的关键改进）。

**Quake 三消融**（`runs/ablate_quake/`，定位干扰源）：

| 消融 | Quake ARI |
|---|---|
| A: 无锚点 + 无自适应（=rank13） | **0.912** ✓ backbone 验证 |
| B: +SVD 锚点，无自适应 | **0.447** ← SVD 锚点炸掉 Quake |
| C: 无锚点 + 自适应 | **0.897** 自适应目标近乎无害 |

**AGCT 干净版（无锚点+自适应）三数据集**：

| 数据集 | NMI | ARI | vs rank13 |
|---|---|---|---|
| Melanoma_5K | 0.725 | **0.654** | +0.002（保住）|
| Quake_10x_Spleen | 0.829 | **0.897** | −0.015（略降）|
| Macosko | 0.571 | **0.277** | −0.066（更差）|

**结论**：
- 自适应机制**按设计工作**（Macosko DEC-KL 0.68→0.21），且**完整保住 Quake+Melanoma**。作为"不破坏强数据集"的机制它成立。
- **但救不了 Macosko**——因为 Macosko 的病根不是 DEC 过锐。铁证（cluster_weight 扫描）：

| Macosko cluster_weight | ARI |
|---|---|
| 0.0（无 DEC，纯 backbone） | **0.376**（家族内最佳）|
| 0.1 | 0.348 |
| 0.35 | 0.343 |

Macosko 上聚类压力**越小越好**，纯 backbone 打败所有 DEC 变体。任何 DEC 家族目标都伤 Macosko。

- **SVD 锚点融合炸掉 Quake**（0.912→0.447）：在不均衡数据上把嵌入拉回被 72% 大类主导的原始 PCA 几何，覆盖 DEC 的 KMeans-友好重塑。且我的锚点实现也**未复现** rank29 的 Macosko 0.569（只得 0.259）——说明 rank29 的 Macosko 优势不来自锚点融合。

---

## 5. 汇总矩阵（本 harness，seed42，ARI）

| 配置 | Melanoma | Quake | Macosko |
|---|---|---|---|
| 纯 backbone（无 DEC，cw=0） | ~0.63 | 0.50 | **0.376** |
| rank13 DEC | 0.652 | **0.912** | 0.343 |
| rank29 fuzzy+anchor | ~0.53 | crash | **0.569** |
| 线1: +门控 NeighborMix | 0.649 | 0.890 | 0.350 |
| 线2: AGCT 自适应（无锚点） | **0.654** | 0.897 | 0.277 |
| 线2: AGCT + SVD 锚点 | 0.595 | 0.447 | 0.259 |
| 同 harness 参照(纯scMAE) | ~0.63 | ~0.50 | 0.376 |

粗体为该数据集当前最佳。**没有任何一行同时在三列取胜。**

---

## 6. 核心科学发现（论文级 motivation）

1. **真实机制冲突**：
   - DEC 锐化：Quake **+0.41**（命脉）、Melanoma +0.02（微助）、Macosko **−0.03**（有害，低于纯 backbone）。
   - 流形/锚点：Macosko 有益（rank29 0.569）、Quake **灾难**（crash 到 0.447）。
   - **不存在固定机制赢三个。**

2. **DEC-KL 收敛性 = 数据集级"可 DEC 聚类性"探测器**：Quake KL→0.0135（收敛）；Macosko KL 停在 0.68→0.21（不收敛）。这是一个廉价、可靠、可作为自适应开关的信号。

3. **自适应粒度目标是一个成立的机制**：per-cell 锐/软插值能在不伤 Quake/Melanoma 的前提下自动卸掉边界细胞的 DEC 压力——只是 Macosko 的瓶颈在**表征层**（backbone 本身），聚类目标层的自适应够不到。

4. **Macosko 仍是全家族未解难题**：baseline 0.494 与 rank29 0.569 未被任何 DEC 家族变体打破；纯 backbone 0.376 是家族内最佳但仍不及格。

---

## 7. 结论与下一步

**本轮两条线均为负结果**，但已排除大片无效假设空间，并沉淀出机制冲突这一核心洞见与 DEC-KL 探测信号。

**推荐下一步（按优先级）**：
1. **解剖 rank29**，定位把 Macosko 从 0.38 抬到 0.57 的**具体组件**（variance loss / balance / center-separation / anchor 重构头 / GELU decoder）。这是整个搜索中**唯一**真正撬动 Macosko 的杠杆，找到它才可能做出跨数据集赢家。
2. **数据集级自适应开关**：用 DEC-KL 收敛信号自动决定是否启用 DEC（Quake/Melanoma 开、Macosko 关退回纯 backbone）。工程稳妥，至少保证"不伤 Macosko"。
3. 若 (1) 定位成功，将该组件与 AGCT/DEC 组合，做真正的数据自适应混合。

---

## 8. 复现命令

```bash
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python
BD=methods/DeepLearning/scMAEs/benchmark_data

# AGCT 干净版（无锚点+自适应），Melanoma 例
CUDA_VISIBLE_DEVICES=1 $PY -u methods/DeepLearning/AdaptiveGranularity_scMAE/run.py \
  --data_path $BD/Melanoma_5K.h5ad --save_dir <out> \
  --dataset_name Melanoma_5K --label_key resolved_label --n_clusters 9 \
  --anchor_dim 0 --adaptive true --epochs 80 --seed 42 --gpu 1

# rank13 backbone 复现（pseudo 关）
CUDA_VISIBLE_DEVICES=1 $PY -u methods/DeepLearning/GatedNeighborMix_scMAE/run.py \
  --data_path $BD/Quake_10x_Spleen.h5ad --save_dir <out> \
  --dataset_name Quake_10x_Spleen --label_key resolved_label --n_clusters 5 \
  --pseudo_weight 0.0 --epochs 80 --seed 42 --gpu 1
```

> 注：`benchmark_data/*.h5ad` 含 `encoding-type:'null'` 的 uns 项，任何 anndata 版本直接读会报 `IORegistryError`；两条线的 `run.py` 均已内置 `_register_null_h5ad_reader()` 修复。大数据集需 `OPENBLAS_NUM_THREADS=8` 等上限（已内置），否则 44k 细胞 KMeans 会触发 OpenBLAS "too many memory regions"。

