# Phase 2: 数据与 preprocessing

## 目标

实现 CAAM 的输入语义，确保 corruption 发生在 log1p、非负、未 gene-wise scale 的表达上。

## 操作步骤

1. 实现 `data/preprocessing.py`。
2. 实现 `data/dataset.py`，训练 Dataset 只返回：

```python
{
    "index": int,
    "x": Tensor[G],
    "batch_code": Optional[int],
    "library_size": float,
    "zero_ratio": float,
}
```

3. benchmark_mode 下复用 canonical input：

```text
input_mode = log1p
n_top_genes = 0
scale_input = false
```

不得重新：

```text
normalize_total
log1p
HVG selection
scale
class-balanced subsample
```

4. standalone raw count 下执行：

```text
normalize_total -> log1p -> HVG selection
```

5. evaluation labels 单独保存为 evaluation bundle，不进入 DataLoader。
6. 若 `scale_input=true`，日志必须警告 zero value 语义改变。

## 验收条件

```text
training Dataset 不含 label/cell_type/true_cluster/n_clusters
benchmark_mode 不重复 preprocessing
standalone raw count 输出非负 log1p 表达
scale_input 默认 false
```

