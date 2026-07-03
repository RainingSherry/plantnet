# ReliableRecon_scMAE（负结果线）

## 一句话
纯 backbone（无 DEC）+ per-cell-per-gene **局部可靠性**精度加权重构——三数据集**单调变差**，证否了"从重构任务加权"这个方向。

## 当前模型的流程

```
输入 x(scaled) ──> encoder(g→256→128→128) ──> latent z
                                              │
              ┌───────────────────────────────┴─────────────────┐
              ▼                                                 ▼
      mask_predictor(z)                               decoder([z, mask_logits])
              │                                                 │
              ▼                                                 ▼
         mask_logits                                       重构 x'
              │                                                 │
              └──────────────────> 精度加权重构损失 <───────────┘
                                   w_ig = w_pos · ((1-λ) + λ·r_ig)
                                   L_recon = Σ w_ig · SmoothL1(x'_ig, target_ig)

r_ig = per-cell-per-gene 局部可靠性（一次性预计算，固定）:
  1. 在 **解耦的原始数据 PCA-KNN 图**（非 live embedding）上
  2. 对每个细胞 i，看基因 g 在其 KNN 邻居中的方差
  3. 相对 per-gene 中位数归一化 → r_ig ∈ [floor, 1]
  4. 局部稳定的基因 r 高（可信），局部波动大的 r 低（技术噪声/dropout）
```

训练损失：
```
L = (1-mask_weight) · L_recon(λ-weighted)  +  mask_weight · BCE(mask_logits, mask)
```

- `λ=0` 退化成 vanilla scMAE；`λ=1` 全精度加权。
- **局部 vs 全局**（设计关键）：marker 基因全局高方差但在其所属细胞类型的 KNN 内**局部稳定**，故 r 高（被保护）。全局方差会误杀 marker。
- floor=0.2 保护稀有细胞（即使所有基因都略不稳，也不至于被完全降权）。

## 思路与为什么失败

**设计意图**：Macosko baseline 低 → 猜测是 backbone 重构质量不够 → 用局部可靠性加权，让模型**重点恢复稳定基因、容忍噪声基因**，提升整体表征。

**结果**（λ=0 vs λ=1，同代码，seed42 ARI）：
| 数据集 | λ=0 vanilla | λ=1 全加权 | Δ |
|---|---|---|---|
| Melanoma | 0.564 | 0.492 | −0.072 |
| Quake | 0.538 | 0.521 | −0.017 |
| Macosko | 0.264 | 0.173 | −0.091 |

**三数据集单调变差**。

**失败根因**（第一性，建前已标注为风险）：
1. **重构保真度 ≠ 聚类质量**。scMAE 的 encoder 学的是 gene-gene 共表达结构（哪些基因一起变化），不是"重构得最精确"。
2. 降权"局部高方差"基因，恰好丢掉了**边界/过渡/稀有细胞的判别信号**——这些信号藏在难重构的基因里（细胞类型之间的 transition 本就是"不稳定"）。
3. scVI 式精度加权对**生成任务**(插补/去噪)有用，因为目标就是"像真数据"；对 **KMeans-on-embedding 聚类**无益甚至有害，因为判别边界来自"异常/不稳定"。

**结论**：此线证否了"从重构任务层面救 Macosko"。问题在聚类目标（后来被 variance-floor 解决）。

## 关键结果（seed42, ARI）
| 数据集 | λ=0 vanilla | λ=1 全加权 |
|---|---|---|
| Melanoma | 0.564 | 0.492 |
| Quake | 0.538 | 0.521 |
| Macosko | 0.264 | 0.173 |

## 文件
- `model.py`（纯 scMAE backbone，无簇心）
- `loss.py`（λ-blended 精度加权重构 + mask BCE）
- `reliability_recon.py`（`compute_local_reliability`：解耦 PCA-KNN → per-gene 局部方差 → r_ig field）
- `run.py`（主程序，`--reliability_lambda` 控制 λ）

复现：`--reliability_lambda 0` = vanilla，`--reliability_lambda 1` = 全加权。
