# GatedNeighborMix_scMAE（负结果线）

## 一句话
per-cell 可靠性门控 NeighborMix（核心细胞混合、稀有/边界关闭）+ DEC 聚类——保住 Melanoma，但**未救 Macosko**（0.350）。含**忠实 rank13 复现**（Quake 0.9118）。

## 当前模型的流程

```
输入 x(scaled) ──> encoder(g→256→128→128) ──> latent z
                                              │
          ┌───────────────────────────────────┼────────────────────────┐
          ▼                                   ▼                        ▼
  mask_predictor(z)                  student_q(z, 簇心)         decoder([z, mask])
          │                          = Student-t 软分配              │
          ▼                                   │                        ▼
     mask_logits                              │                    重构 x'
                                              │
        per-cell 可靠性 r_i ←─────────────────┤
        (邻居标签一致性 × 局部密度 × 成员置信)
                │                             │
                ├──> gate NeighborMix:        ├──> gate DEC:
                │    alpha_i = 1-(1-αₘᵢₙ)·r_i │    KL 系数 × r_i
                │    r≈1核心→强混(α→αₘᵢₙ)     │    r≈1核心→强拉
                │    r≈0边界→纯scMAE(α→1)      │    r≈0边界→KL≈0
                ▼                             ▼
   NeighborMix pseudo 损失            置信门控 KL(p‖q)
   x_pseudo = alpha·x + (1-alpha)·邻居均值
   重构目标 = 真细胞 x (anchor-recovery)
```

训练损失：
```
L = scmae(recon + mask)                          # 主体，永远开（安全网）
  + pseudo_weight · L_pseudo(gated NeighborMix)  # r_i 门控
  + cluster_scale · cluster_weight · r_i·KL(p‖q) # r_i 门控 DEC
```

- 可靠性 r_i（`reliability.py`）：邻居标签一致性（伪标签 consistency）× 局部密度（KNN 距离倒数）× 成员置信（max(q)）。
- NeighborMix 门控：核心细胞 r≈1 → alpha→0.85（强混）；边界/稀有 r≈0 → alpha→1.0（混合关闭，退回纯 scMAE）。
- DEC 门控：KL 损失系数 × r_i，边界细胞 KL≈0（不被强拉向簇心）。

## 思路与为什么失败

**设计意图**：全benchmark 显示 NeighborMix 是唯一超 Melanoma 的表达级机制（+0.042 ARI），但过平滑伤 Macosko/Tosches（−0.07）。假设：用可靠性**门控**混合——核心细胞混、边界细胞不混 → 保住 Melanoma 增益 + 去掉 Macosko 损失。

**踩坑与修正**（关键教训，全线通用）：
1. 从零重写 DEC，Quake 只有 0.52（应 0.91）。**根因**：重构目标误用 *scaled* 空间 → recon loss 掉到 0.07 → DEC KL 反客为主炸掉嵌入。**修正**：重构目标必须是 *未缩放 log*（保持 recon loss ~0.44，让 DEC 当温和正则）。
2. 决策：改为**在 rank13 真实代码上叠加**门控分支。`runs/repro_check/` 确认 backbone 精确复现 rank13（Quake **0.9118** = 原文 0.912）。

**结果**（gated NeighborMix, pseudo=0.3, α_min=0.85，seed42 ARI）：
- Melanoma: **0.649** (略低于 rank13 的 0.652，未达到 NeighborMix 全开的 +0.042)
- Quake: 0.890 (低于 rank13 的 0.912)
- Macosko: **0.350** (略高于 rank13 的 0.343，但远低于 baseline 0.494)

**中性偏负**，核心假设不成立。

**诊断**：Macosko 上 DEC-KL 在 epoch80 仍高达 0.679（Quake 仅 0.0135）——DEC 在 Macosko **根本不收敛**，门控只是"减轻伤害"但救不了根本问题（后来在 AdaptiveSwitch 里证实是"每维方差坍缩"）。

**结论**：此线证明"门控 NeighborMix"不是杠杆。但它的副产物——**忠实 rank13 复现**——成为后续所有 DEC 变体的验证锚点。

## 关键结果（seed42, ARI）
| 数据集 | GatedNeighborMix | rank13 DEC(参照) |
|---|---|---|
| Melanoma | 0.649 | 0.652 |
| Quake | 0.890 | 0.912 |
| Macosko | 0.350 | 0.343 |

**repro_check**（Quake，验证 backbone）：本实现 0.9118 vs rank13 原文 0.912 ✓

## 文件
- `model.py`（encoder + mask + decoder + DEC 簇心 + soft membership）
- `reliability.py`（per-cell r_i：邻居一致性 × 密度 × 置信）
- `loss.py`（gated scMAE + gated NeighborMix pseudo + gated DEC KL）
- `run.py`（训练/评估循环，PCA-KNN 图，smoke/screen 阶段）

关键参数：`--pseudo_weight`（NeighborMix 权重）、`--neighbor_alpha_min`（核心细胞混合强度）、`--cluster_weight`（DEC 权重）。
