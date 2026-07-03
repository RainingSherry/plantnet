# AdaptiveGranularity_scMAE（负结果线）

## 一句话
自适应聚类目标锐度（per-cell 可聚类性门控）+ 可选 SVD 锚点融合——保住 Quake/Melanoma，但**救不了 Macosko**（0.277）。

## 当前模型的流程

```
输入 x(scaled) ──> encoder(g→256→128→128) ──> latent z
                                              │
                     (可选) SVD锚点融合 ←─────┤
                          z' = 0.7·z + 0.3·SVD(raw_PCA 投影回 128)
                                              │
                                              ▼
                                      student_q(z', 簇心)
                                   = Student-t 软分配 q
                                              │
                      计算 per-cell 可聚类性 c_i ←──────┐
                      (邻居标签一致性 × 局部密度)        │
                                              │         │
                                              ▼         │
                        自适应目标 p_i = c_i·sharpen(q_i) + (1-c_i)·q_i
                        （核心细胞锐、边界细胞保持模糊 q）
                                              │
                                              ▼
                                    置信门控 KL(p‖q)
```

训练损失：
```
L = scmae(recon + mask)                          # 主体
  + cluster_scale · cluster_weight · KL(p‖q)     # warmup 后 ramp
```

- 锚点融合（`--anchor_weight`）：把 latent 拉向原始 PCA 几何（保留线性判别信息）。
- 可聚类性 c_i：fuzzy-rough 式上/下近似 + 确定性退火（epoch 0→50 从 0.1→0.9 退火）。
- 核心细胞 c≈1 → 纯锐目标 sharpen(q)；边界细胞 c≈0 → 目标退化成 q（KL(q‖q)=0，无约束）。

## 思路与为什么失败
**设计意图**：DEC 全局锐化伤边界/稀有 → 让**边界细胞不被强推**（自适应目标 p=q 时 KL=0）。

**Quake 三消融**（`runs/ablate_quake/`）定位干扰源：
- A(无锚+无自适应)=**0.912**（= rank13 ✓）
- B(+SVD锚点)=**0.447**（锚点炸 Quake：72% 大类主导的 PCA 几何毁掉 DEC 学到的细微边界）
- C(无锚+自适应)=0.897（自适应目标近乎无害）

**跨三数据集（无锚+自适应）**：Melanoma 0.654 / Quake 0.897 / Macosko **0.277**。

**失败根因**（铁证）：Macosko 上 `cluster_weight` 越小越好：
- cw=0 (纯 backbone): ARI **0.376**
- cw=0.35 (rank13 default): ARI 0.343
- cw+自适应: ARI 0.277

→ **任何 DEC 家族目标都伤 Macosko**。自适应只能"减轻伤害"（把边界细胞退回 q），但救不了核心细胞被锐化后的**表征坍缩**（后来在 AdaptiveSwitch 里证实是"每维方差坍缩"）。

**rank29 解剖副产物**：`runs/rank29_dissect/` 完整记录 rank29 留一消融（entropy/fuzzy/variance/balance 等），见 ../EXPERIMENT_LOG.md §4。

## 关键结果（seed42, ARI）
| 数据集 | AdaptiveGranularity(无锚) | rank13 DEC | 纯 backbone |
|---|---|---|---|
| Melanoma | 0.654 | 0.652 | ~0.63 |
| Quake | 0.897 | 0.912 | ~0.50 |
| Macosko | 0.277 | 0.343 | 0.376 |

**结论**：此线证明"调聚类目标的锐度"救不了 Macosko，问题在**表征层**（后来被 variance-floor 解决）。

## 文件
- `model.py`（backbone + 可选锚点融合 + 簇心）
- `loss.py`（KL 聚类损失）
- `clusterability.py`（per-cell c_i 估计：邻居一致性 × 密度 + 退火）
- `_run_rank29.py`（wrapper，注入 null-h5ad shim 复现 rank29；全 6 个损失 CLI 可消融）
- `run.py`（主程序）

关键参数：`--cluster_weight`、`--anchor_weight`（0=无锚点）、`--adaptive_target`（bool，开/关自适应）。
