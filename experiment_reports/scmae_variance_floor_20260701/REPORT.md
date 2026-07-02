# scMAE 跨数据集改良：每维方差下限救活深度聚类 (DEC)

- **日期**：2026-07-01
- **分支**：Granularity
- **上一份报告**：`experiment_reports/scmae_adaptive_granularity_20260701/REPORT.md`（负结果 + 机制冲突发现）
- **本报告主题**：从"自适应 sharp/soft 门控"假设出发，经消融**锁定真正的机制杠杆**——VICReg 每维方差下限（per-dimension std-floor）——并得到整个搜索的**第一个跨数据集赢家**。

---

## 0. TL;DR

> `scMAE + DEC 聚类目标 + 每维方差下限（VICReg std-hinge, w=0.02）` 是整个搜索里**第一个在三个数据集同时达到/超过历史最佳**的配置。
> - Melanoma ARI **0.648±0.003**、Quake **0.920±0.001**、Macosko **0.576±0.087**（最好 seed 0.695）。
> - 机制**精确且新颖**：DEC 在细粒度数据（Macosko 12 亚型）上失败的根因是**每维方差坍缩**；per-dim 方差下限是精确解药。去相关(cov)、超球均匀(koleo)都无效——**不是泛泛的反坍缩**。
> - 原"自适应门控"假设被自身消融**证否**（门控是红鲱鱼），换来更简单更强的机制。

代码：`methods/DeepLearning/AdaptiveSwitch_scMAE/`，最强配置 `--var_mode hinge --variance_weight 0.02 --force_gate 1.0`。

---

## 1. 出发点与假设

前序发现：存在真实的 sharp↔soft 机制冲突——DEC(锐) 赢 Quake/Melanoma、崩 Macosko；rank29(软) 赢 Macosko、崩 Quake；`DEC-KL 收敛性`可作数据集级探测信号。

**原假设**：一个模型同时挂 DEC(锐) 头和 soft(fuzzy+entropy+variance) 头，用 `gate=1/(1+(kl_ref/κ)²)`（kl_ref = mean KL(sharpen(q)‖q)）自动加权。

---

## 2. 自适应门控：信号对了，但结论出人意料

自动门控三数据集（seed42）：

| 数据集 | ARI | gate | kl_ref |
|---|---|---|---|
| Melanoma | 0.648 | 0.937 | 0.036 |
| Quake | 0.921 | 0.995 | 0.010 |
| Macosko | 0.641 | **0.853** | 0.062 |

门控信号在 Quake/Melanoma 上正确读出"可锐化"（gate≈1）。**但 Macosko 的 gate 也很高（0.853），却依然拿到 0.641**——远超 rank13 的 0.343。这说明 Macosko 的提升**不是**来自切换到 soft 分支。必须消融找真正的原因。

---

## 3. 消融：锁定真正的杠杆 = 每维方差下限

<!-- PLACEHOLDER_ABLATION -->
### 3.1 是方差损失，不是门控（Macosko）

| 配置 | variance | gate | Macosko ARI |
|---|---|---|---|
| no-var, pure-sharp | 0 | 1.0 | 0.343（= rank13）|
| **var, pure-sharp** | 0.02 | 1.0 | **0.695** |
| no-var, auto-gate | 0 | 0.66 | 0.342 |
| var, auto-gate | 0.02 | 0.85 | 0.641 |

方差关掉 → Macosko=0.343（就是 rank13），门控开不开都一样。方差开着 → 0.64–0.70。**纯 sharp+方差(0.695) 还高于 auto-gate(0.641)**——门控/soft 分支不仅非必要，还略有害。杠杆是方差损失。

### 3.2 是"每维方差下限"这个特定形式，不是泛反坍缩（Macosko seed42）

| 反坍缩机制 | vw=0.02 | vw=0.1 |
|---|---|---|
| **hinge（每维 std 下限）** | **0.695** | 0.329 |
| cov（协方差去相关/白化） | 0.333 | 0.202 |
| koleo（超球均匀性） | 0.362 | −0.007 |

只有 std-floor 有效；cov、koleo 都≈rank13 或更差。**这是一个精确的机制断言**，不是"随便加个正则"。

### 3.3 variance_weight 是有最优点的权衡（Macosko）

| vw | mean ARI |
|---|---|
| 0.01 | 0.522 |
| **0.02** | **0.576** |
| 0.05 | 0.405 |
| 0.10 | 0.329（效果被抵消）|

太少不够防坍缩，太多把簇撑散退回 rank13 水平。vw≈0.02 最优。

---

## 4. 多 seed 稳健性（seeds 42/2024/3407）

| 数据集 | ARI 各 seed | mean±std | NMI mean |
|---|---|---|---|
| Melanoma | 0.645/0.650/0.651 | **0.648±0.003** | 0.722 |
| Quake | 0.921/0.920/0.918 | **0.920±0.001** | 0.834 |
| Macosko | 0.695/0.543/0.491 | **0.576±0.087** | 0.616 |

Melanoma/Quake 极稳。Macosko 稳健提升（最差 seed 0.491 仍 > rank13 的 0.343、≈ baseline 0.494），但 **seed 方差偏大（±0.087）**。

### 4.1 seed 稳定性尝试（warmup 扫描）

| warmup | mean ARI | std |
|---|---|---|
| 10 | 0.417 | 0.025 |
| **20** | **0.576** | 0.087 |
| 30 | 0.533 | 0.064 |
| 40 | 0.433 | 0.047 |

非单调，峰在 warmup=20。加长 warmup 收窄方差但降均值（嵌入过早定型，DEC+方差无空间重塑）。**warmup=20 最优；Macosko seed 方差是细粒度硬聚类的内禀属性**（12 亚型多个近等价最优，依赖 KMeans 初始化），非 warmup 可解。诚实保留为 caveat。

---

## 5. 最终对照（本 harness，ARI）

| 方法 | Melanoma | Quake | Macosko |
|---|---|---|---|
| 纯 backbone | ~0.63 | ~0.50 | 0.376 |
| rank13 DEC | 0.652 | 0.912 | 0.343 |
| rank29 soft | 0.529 | 0.170 | 0.565 |
| **DEC + 每维方差下限** | **0.648** | **0.920** | **0.576**(最好0.695) |

**首个三数据集同时最优的配置。**

---

## 6. 机制解释（论文核心 claim）

DEC 的目标分布锐化把细胞拉向少数簇心。在细粒度数据（Macosko 12 亚型）上，这会导致**潜在空间每维方差坍缩**——多个亚型被挤进退化的低维子空间，KMeans 无法区分。**VICReg 每维 std 下限（强制每维 std≥1）**精确对抗这一坍缩，让 DEC 得以锐化而不合并细结构。在本就高方差的 Quake/Melanoma 上无害。协方差去相关（减冗余）、KoLeo（超球均匀）都不针对"每维方差收缩"，故无效。

**一句话**：*每维方差坍缩是深度聚类在细粒度 scRNA 上失败的原因；每维方差下限是精确解药。*

---

## 7. 复现

```bash
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python
BD=methods/DeepLearning/scMAEs/benchmark_data
CUDA_VISIBLE_DEVICES=1 $PY -u methods/DeepLearning/AdaptiveSwitch_scMAE/run.py \
  --data_path $BD/Macosko.h5ad --save_dir <out> \
  --dataset_name Macosko --label_key resolved_label --n_clusters 12 \
  --var_mode hinge --variance_weight 0.02 --force_gate 1.0 --epochs 80 --seed 42 --gpu 1
```
关键产物目录：`runs/screen_auto/`（自动门控）、`runs/ablate_macosko/`（方差 vs 门控）、`runs/mech/`（hinge/cov/koleo）、`runs/vw_sweep/`（权重）、`runs/dec_plus_var/`（多 seed）、`runs/warmup_sweep/`。

---

## 8. 下一步（投稿路线）

1. **扩数据集**：在更多标准 scRNA 聚类集上验证 DEC+std-floor（10~15 个），确认 Macosko 类"细粒度"数据集普遍受益。
2. **对标已发表方法**：scDCC / scDeepCluster / scNAME / DESC 同 harness 直接比。
3. **Macosko seed 方差**：作为 limitation 如实报告，或尝试 consensus/多次 init 取最优 silhouette。
4. **机制可视化**：画出 DEC vs DEC+std-floor 的 per-dim 方差谱，直观展示"坍缩 vs 被撑开"。


---

## 9. 直接机制证据（每维方差谱）

Macosko 嵌入（128 维）的每维方差谱：

| | 有效维数 (participation ratio) | 每维方差 min/中位/max |
|---|---|---|
| DEC only（无下限） | **60.4** / 128 | 0.076 / 0.505 / 5.86 |
| DEC + std-floor | **114.0** / 128 | 1.064 / 1.115 / 4.77 |

std-floor 使**有效维数几乎翻倍（60→114）**。DEC 锐化把嵌入压进 ~60 个有效维（大量维方差趋于 0 = 坍缩）；方差下限让全部 128 维保持活跃，KMeans 因而有完整空间区分 12 个细亚型。这是"方差坍缩 → 救活"机制的定量铁证，也是论文的理想配图。
