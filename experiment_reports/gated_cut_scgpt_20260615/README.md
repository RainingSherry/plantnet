# Gated CutAware NeighborMix + scGPT Run 2026-06-15

本目录归档本轮进一步优化、近期文献源码检索、门控迁移和 scGPT 运行结果。

## 新增代码

- `methods/DeepLearning/CutAware_NeighborMix_scMAE/model.py`
  - 新增 `edge_gate`、`edge_gate_scores(...)`。
- `methods/DeepLearning/CutAware_NeighborMix_scMAE/mixing.py`
  - 新增 `make_gated_neighbor_mixed_batch(...)`。
- `methods/DeepLearning/CutAware_NeighborMix_scMAE/run.py`
  - 新增 `canm_gated_cut_mix` 和 `canm_gated_cut_warm`。
- `methods/Foundation/scGPT/run_plantnet.py`
  - 新增 PlantNet-safe scGPT runner。
- `methods/Foundation/scGPT/run_plantnet_suite.py`
  - 新增 5 数据集批量运行脚本，禁止 GPU 0/7。

## 运行范围

代表性数据集：

- `SRP182008`
- `Melanoma_5K`
- `Macosko`
- `Tosches`
- `Wang`

全部使用 seed 42、80 epoch。未使用 GPU 0/7。Macosko 在并发预处理阶段多次触发本地 OpenBLAS `rc=-11`，均用单进程、线程数 1 补跑完成。

## 结果

| 方法 | mean ARI | 结论 |
| --- | ---: | --- |
| `canm_cut_reweighted_mix` | 0.632087 | 当前最稳；cut-informed 边权直接进 NeighborMix 比 direct cut loss 更有效。 |
| `canm_gated_cut_mix` | 0.572877 | 门控会放大早期噪声；Tosches/Macosko 明显下降。 |
| `canm_gated_cut_warm` | 0.571910 | warm 和去 cluster gate loss 后仍未改善。 |
| `scGPT_runner_pca_fallback` | 0.550407 | 当前环境不是实际 scGPT transformer；是安全 PCA fallback。 |

逐数据集表：

- `tables/gated_warm_scgpt_comparison.csv`
- `tables/cutaware_vs_rg_summary.csv`
- `tables/scgpt_plantnet_summary.csv`

## scGPT 状态

官方 scGPT 需要 `scgpt` 包和 checkpoint。当前环境：

- `scgpt` 包未安装；
- `methods/Foundation/scGPT/` 有 `vocab.json` 和 `args.json`，但没有 `best_model.pt`；
- 本机只找到 `/home/luolie/biopipeline/PRESCRIBE/scLLM_weights/scGPT/embedding.pkl`，不是可加载 checkpoint；
- PyPI dry-run 会安装旧 `torchtext`、`scvi-tools<1.0` 和 CUDA wheel，风险较高，未污染当前 base 环境。

因此本轮 scGPT 全部执行了安全 runner，并在 `summary.json` 里记录：

```text
used_scgpt_transformer = false
fallback_reason = scgpt package unavailable: No module named 'scgpt'
```

## 当前判断

门控不是当前最优方向。对这些数据，最有效的改法仍是明确 cut-informed edge prior，也就是 `canm_cut_reweighted_mix`。如果继续做 gate，需要让 gate 预测“跨簇边概率”而不是让伪重构 loss 自己学，否则它会重新变成一个平滑器。
