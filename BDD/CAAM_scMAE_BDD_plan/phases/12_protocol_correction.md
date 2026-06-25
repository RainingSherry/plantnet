# Phase 12: Correction 1 — Protocol correction

## 目标

修正 CAAM 的主实验协议，使其更接近 scMAE/scNAME/scCluBench 背景下的公平比较：

```text
1. 主协议使用 HVG feature space，而不是 full-gene input。
2. zero-to-zero corruption 不再导致默认 fail-fast。
3. strict effective budget 变成显式开关。
4. formal smoke 与 validator 适配新协议。
```

## previous_phase_check

执行本阶段前，Codex 必须检查并在报告中写明：

```text
1. caam_scmae 当前是否仍在 method_manifest.yaml 中保持 authenticity=PENDING。
2. default_in_formal 是否仍为 false。
3. 是否存在 control/axial/advmask 的 formal manifest entry。
4. registry.py 是否仍在 benchmark_mode 下强制 n_top_genes=0。
5. validate_formal_smoke.py 是否仍要求 n_top_genes==0。
6. risk_and_stop_criteria.md 是否已经将 budget deficit 默认降级为 diagnostic。
7. 当前工作区是否有非 CAAM 未提交改动；不得提交它们。
```

若 1-3 不满足，不得继续。

## 修改要求

### 1. Benchmark feature space

`benchmark_mode` 默认：

```text
input_mode = log1p
n_top_genes = 2000
scale_input = false
```

CLI 显式传入时允许：

```text
--n_top_genes 3000
--n_top_genes 0
```

其中：

```text
n_top_genes=0 => full-gene stress test 或 external_hvg input
```

必须记录：

```text
feature_space_source
actual_n_genes_after_selection
selected_gene_indices_path
selected_gene_names_path if available
```

### 2. Strict effective budget

新增或修正配置：

```text
strict_effective_budget: false by default
```

当 `strict_effective_budget=false` 时：

```text
budget_deficit_rate 不触发 fail-fast
zero_to_zero_rate 不触发 fail-fast
effective_corruption_rate 不触发 fail-fast
```

当 `strict_effective_budget=true` 时，允许旧 BDD 的 deficit gate。

### 3. Diagnostics

每个 run 必须输出：

```text
corruption_stats.json
```

至少包含：

```text
corruption_type
mask_ratio
n_top_genes
actual_n_genes
zero_to_zero_rate
effective_corruption_rate
budget_deficit_rate
mean_abs_delta
mean_abs_delta_masked
strict_effective_budget
```

### 4. Validator

`validate_formal_smoke.py` 不得硬编码要求 `n_top_genes == 0`。它应检查：

```text
n_top_genes in {2000, 3000, 0}
if n_top_genes == 0: feature_space_source in {full_gene_stress, external_hvg}
scale_input == false
benchmark_mode == true
variant == full
```

### 5. Manifest smoke status

如果 formal behavior 改变，则必须：

```text
1. 先把当前 smoke PASS 视为 old-protocol smoke；
2. 重新跑 GPU 1-seed smoke；
3. 重新跑 GPU 3-seed smoke；
4. validate_formal_smoke.py 通过后，才允许保持 PASS。
```

不要改 authenticity，不要改 default_in_formal。

## 测试要求

新增或更新：

```text
test_protocol_correction.py
```

最少测试：

```text
1. benchmark_mode 默认 n_top_genes=2000。
2. CLI --n_top_genes 0 不被覆盖。
3. strict_effective_budget 默认 false。
4. budget_deficit_rate 存在但不导致 fail-fast。
5. selected_gene_indices 在同 seed 下可复现。
6. validate_formal_smoke.py 接受 n_top_genes=2000。
```

## 实验要求

本阶段只跑 smoke，不做研究结论：

```text
1. 1-seed GPU formal smoke
2. 3-seed GPU formal smoke
3. validate_formal_smoke.py
```

不得触碰 validation/sealed test。

## 输出报告

新增：

```text
methods/DeepLearning/CAAM_scMAE/benchmark/PHASE12_PROTOCOL_CORRECTION_REPORT.md
```

报告必须包含：

```text
old_protocol_behavior
new_protocol_behavior
smoke commands
smoke status
whether smoke PASS remains valid
files changed
tests run
```

## 提交策略

只提交：

```text
CAAM code
CAAM tests
CAAM BDD/report
method_manifest smoke 字段，仅当重新 smoke gate 满足
```

不得提交：

```text
results/
data/smoke/*.h5ad
非 CAAM 遗留改动
```
