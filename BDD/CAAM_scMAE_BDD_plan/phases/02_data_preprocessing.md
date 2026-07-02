# Phase 2: 数据与 preprocessing

## 目标

实现 CAAM 的输入语义，确保 corruption 发生在 log1p、非负、未 gene-wise scale 的表达上，并修正主 benchmark 的 feature space（特征空间）协议：**主协议使用 label-free HVG，不再默认使用 full-gene input**。

## 关键修订

早期 BDD 要求 benchmark_mode 下：

```text
input_mode = log1p
n_top_genes = 0
scale_input = false
```

该设定现在仅保留为 **full-gene stress test（全基因压力测试）**。主 benchmark / quick ablation / development ablation 的默认 feature space 应为：

```text
input_mode = log1p
n_top_genes = 2000   # 可通过 CLI 显式改为 3000 或 0
scale_input = false
```

`n_top_genes=0` 只能表示：

```text
full-gene stress test
```

不得作为主表协议、默认 development ablation 协议或论文主结果协议。

## 操作步骤

1. 实现或修正 `data/preprocessing.py`。
2. 实现或修正 `data/dataset.py`，训练 Dataset 只返回：

```python
{
    "index": int,
    "x": Tensor[G],
    "batch_code": Optional[int],
    "library_size": float,
    "zero_ratio": float,
}
```

3. benchmark_mode 下必须使用 label-free feature selection（无标签特征选择）：

```text
input_mode = log1p
n_top_genes = 2000 by default
scale_input = false
canonical_input = true
```

4. 如果上游 benchmark 已经提供 HVG matrix（高变基因矩阵），则允许：

```text
n_top_genes = 0
feature_space = external_hvg
```

但必须在 `resolved_config.yaml` 和 artifact 中明确记录：

```text
feature_space_source = external_hvg
```

5. 如果 CAAM 自己执行 HVG selection（高变基因选择），必须满足：

```text
label-free
deterministic under seed
no class-balanced subsample
no label/cell_type/n_clusters access
save selected_gene_indices.npy
save selected_gene_names.txt if available
```

6. standalone raw count 下执行：

```text
normalize_total -> log1p -> HVG selection
```

7. benchmark_mode 下不得重复：

```text
normalize_total
log1p
scale
class-balanced subsample
```

8. evaluation labels 单独保存为 evaluation bundle，不进入 DataLoader。
9. 若 `scale_input=true`，日志必须警告 zero value 语义改变；主协议默认不得 scale。
10. 所有运行必须记录：

```text
input_mode
n_top_genes
actual_n_genes_after_selection
feature_space_source
selected_gene_indices_path
scale_input
```

## 验收条件

```text
training Dataset 不含 label/cell_type/true_cluster/n_clusters
主 benchmark 默认 n_top_genes=2000，而不是 0
n_top_genes=0 只用于 full-gene stress test 或 external_hvg 输入
standalone raw count 输出非负 log1p 表达
scale_input 默认 false
selected_gene_indices.npy 可复现
相同 seed + 相同 data_path + 相同 n_top_genes 得到相同 feature space
```

## 必须检查的上一阶段/历史问题

Codex 在执行本阶段或后续 correction 阶段前必须检查：

```text
1. registry.py 是否仍在 benchmark_mode 下强制 n_top_genes=0；若是，必须修正。
2. method_manifest.yaml 是否通过 extra_args 强制 --n_top_genes 0；若是，必须移除或改为 2000。
3. validate_formal_smoke.py 是否仍要求 preprocessing.n_top_genes == 0；若是，必须改为兼容主协议 2000 和 full-gene stress test。
4. 任何已有 smoke PASS 若依赖旧 n_top_genes=0 协议，协议改变后必须重新跑 smoke 或在报告中标记 old-protocol smoke。
```
