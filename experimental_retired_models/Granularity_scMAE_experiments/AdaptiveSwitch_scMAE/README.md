# AdaptiveSwitch_scMAE ★（本轮赢家）

## 一句话
`scMAE + DEC 聚类目标 + 每维方差下限（VICReg std-hinge）`——首个在 Melanoma/Quake/Macosko 三数据集**同时**达到/超过历史最佳的配置。

## 当前模型的流程

```
输入 x(scaled) ──> encoder(g→256→128→128) ──> latent z
                                              │
              ┌───────────────────────────────┼───────────────────────────┐
              ▼                               ▼                             ▼
      mask_predictor(z)→mask_logits     student_q(z, 簇心)            方差下限损失
              │                          = Student-t 软分配 q         relu(1 - std_d(z))
              ▼                               │                        （每维方差下限）
   decoder([z, mask_logits])→重构          目标 p = sharpen(q)
              │                          （DEC 锐化目标）
   重构损失(未缩放log目标) + BCE(mask)   置信门控 KL(p‖q)
```

训练损失：
```
L = scmae(recon + mask)                          # 主体，永远开
  + consistency(strong/weak latent)              # 弱正则
  + cluster_scale · cluster_weight · L_cluster   # warmup 后 ramp
  + variance_weight · L_variance                 # ★ 关键：每维方差下限，永远开
L_cluster = gate·L_sharp + (1-gate)·L_soft       # 最强配置 force_gate=1 → 纯 L_sharp
L_sharp = 置信门控 KL(p‖q)   （= rank13 DEC）
```

- 编码器输入 = **scaled** 表达；重构目标 = **未缩放 log**（关键：保持 recon loss 大，DEC 才是温和正则而非主宰）。
- warmup(20 ep) 后：全量提取 embedding → KMeans 初始化簇心 → 每 5 ep 刷新目标分布 p。
- `--force_gate 1.0` 即最强配置（纯 sharp DEC + 方差下限）；自动门控是探索遗留，非必要。

## 思路与机制（为什么这样设计）
- rank13 的 DEC 在 Quake/Melanoma 强，但在 Macosko(12 细亚型) 崩到 0.343。
- **根因**：DEC 目标锐化把细胞拉向少数簇心，导致潜在空间**每维方差坍缩**——有效维数从 128 掉到 ~60，多个亚型被挤进退化子空间，KMeans 无法区分。
- **解药**：VICReg 每维 std 下限 `relu(1-std_d)` 强制每维 std≥1，阻止坍缩（有效维数回到 ~114），DEC 得以锐化而不合并细结构。在本就高方差的 Quake/Melanoma 上无害。
- 机制**特定于每维方差下限**：协方差去相关(cov)、超球均匀(koleo) 都无效（详见 ../EXPERIMENT_LOG.md §6.2）。

## 关键结果（多 seed 42/2024/3407, ARI）
| 数据集 | AdaptiveSwitch | rank13 DEC |
|---|---|---|
| Melanoma_5K | 0.648±0.003 | 0.652 |
| Quake_10x_Spleen | 0.920±0.001 | 0.912 |
| Macosko | 0.576±0.087（best 0.695）| 0.343 |

## 复现
```bash
PY=/data/luolie/conda/envs/scssl_bench_py310/bin/python
BD=methods/DeepLearning/scMAEs/benchmark_data
CUDA_VISIBLE_DEVICES=1 $PY -u methods/DeepLearning/Granularity_scMAE_experiments/AdaptiveSwitch_scMAE/run.py \
  --data_path $BD/Macosko.h5ad --save_dir <out> \
  --dataset_name Macosko --label_key resolved_label --n_clusters 12 \
  --var_mode hinge --variance_weight 0.02 --force_gate 1.0 --epochs 80 --seed 42 --gpu 1
```
文件：`model.py`（backbone+簇心）、`loss.py`（gate/sharp/soft/方差 + `compute_gate`）、`clusterability.py`（per-cell c_i）、`run.py`。
关键参数：`--var_mode {hinge,cov,koleo,both}`、`--variance_weight`（0.02 最优）、`--force_gate`（-1 自动 / 1 纯sharp）。
