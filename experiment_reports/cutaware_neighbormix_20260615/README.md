# CutAware NeighborMix scMAE Experiment 2026-06-15

本目录归档本次把 scGNN / scDSC / AttentionAE-sc / scCDCG 思路接入 NeighborMix-scMAE 的代码、轻量输出和结论。

## 代码位置

- 实际方法代码副本：`code/CutAware_NeighborMix_scMAE/`
- 工作区主代码位置：`methods/DeepLearning/CutAware_NeighborMix_scMAE/`
- 完整本地输出位置：`results/experimental/cutaware_neighbormix_20260615/`

## 已实现的网络改动

- scGNN 思路：先用 scMAE/PCA 表征构图，支持 KNN、可靠性边权、剪枝和图刷新。
- scDSC 思路：在 scMAE 编码器后加入 cluster head，并提供聚类/图结构训练目标。
- AttentionAE-sc 思路：保留 attention fusion probe，用来诊断表达 assignment 和图 assignment 是否一致。
- scCDCG/DCGC 思路：新增 `canm_cut_reweighted_mix`，先用无监督图分割识别候选跨簇边，再在 NeighborMix 混合前显式下调这些边权。

## 代表性结果

5 个代表数据集、seed 42、80 epoch，未使用 GPU 0/7。

| 方法 | n_runs | mean ARI | min ARI | max ARI |
| --- | ---: | ---: | ---: | ---: |
| canm_cut_reweighted_mix | 5 | 0.632087 | 0.373287 | 0.939177 |
| canm_diagnostic_only | 5 | 0.605034 | 0.410052 | 0.744993 |
| canm_cut_ot | 5 | 0.256310 | 0.000000 | 0.733063 |
| rg_none | 8 | 0.631446 | 0.366215 | 0.906919 |
| rg_reliability | 8 | 0.626196 | 0.344666 | 0.942741 |

逐数据集结果见 `tables/canm_cut_reweighted_mix_per_dataset.csv`，总表见 `tables/cutaware_vs_rg_detail.csv` 和 `tables/cutaware_vs_rg_summary.csv`。

## 结论

直接把 NCut/OT loss 加到 encoder 上会不稳定：`canm_cut_ot`、`canm_mix_plus_cut` 和 warm 版本在多个数据集上出现 ARI=0 或 NaN/collapse。

把 scCDCG 的“cut-informed”思想改成 NeighborMix 边权机制后，训练稳定性明显改善：`canm_cut_reweighted_mix` 在 5 个代表数据集上平均 ARI 为 0.632087，明显高于 direct cut/OT 路线，并且高于 RG phase1 `rg_reliability` 的代表集均值。

但它还没有全面超过 expression-only/scMAE-style baseline：相对 `canm_diagnostic_only`，主要提升来自 Wang，Macosko/SRP182008/Tosches/Melanoma_5K 是接近或小幅下降。因此当前判断是：cut-informed 边权是比直接 cut loss 更合理的移植方式，但还需要学习式 edge reliability 或训练中图刷新，不能只依赖初始 PCA-KMeans 切边。

## 大文件说明

本归档不包含 `*.pt`, `*.npy`, `*.h5`, `*.h5ad`。这些完整输出在本机 `results/experimental/cutaware_neighbormix_20260615/`，GitHub 只提交轻量摘要和代码。
