# Phase 3: DonorCandidateProvider 与 corruption

## 目标

实现 label-free matched gene-wise donor replacement 与 eligibility。

## 操作步骤

1. 实现 `data/donor_candidates.py`。
2. 基于以下字段构造 donor pool：

```text
batch_code
library_size
zero_ratio
```

3. donor 优先级：

```text
1. same batch + same library-size bin + same zero-ratio bin
2. same batch + nearest bins
3. global + same bins
4. global random
```

4. 保证每个 donor：

```text
r != i
```

5. 对每个位置独立采样 donor：

```text
V_ij = X_rij,j
```

不得让一个 target cell 的所有 masked genes 共享同一 donor row。

6. 实现 eligibility：

```text
E_ij = not isclose(V_ij, X_ij)
```

7. MaskSelector 只能在 `E=1` 的位置选择 mask。
8. 记录 budget deficit；超阈值时 fail-fast。
9. 实现 corruption 输出契约：

```python
{
    "x_tilde": ...,
    "replacement": ...,
    "selected_mask": ...,
    "effective_mask": ...,
    "eligibility": ...,
    "donor_indices": ...,
    "budget_deficit": ...
}
```

## 验收条件

```text
replacement 来自同一 gene
donor r != i
不同 gene 独立 donor
不使用 label
ineligible 不被 mask
deficit 被记录
```

