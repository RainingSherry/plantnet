# Phase 11: 正式 benchmark 接入

## 目标

将最终 CAAM-scMAE 作为一个正式方法接入仓库正式 benchmark。只接入 Model C。

## 操作步骤

1. 修改 `methods/method_manifest.yaml`，只注册：

```text
caam_scmae
```

2. `caam_scmae` 调用：

```text
methods/DeepLearning/CAAM_scMAE/run.py --variant full
```

3. 修改或复用 `envs/runtime_registry.yaml` 中的 Python runtime。
4. 确认 `scripts/run_formal_benchmark.py` 能按正式体系调用 CAAM runner。
5. 不修改：

```text
methods/DeepLearning/PlantSPADE_LGCL/scripts/run_single.py
```

6. Benchmark 模式传入：

```text
--benchmark_mode true
--skip_eval true
--input_mode log1p
--n_top_genes 0
--scale_input false
```

7. CAAM runner 输出：

```text
embedding_final.npy
artifact_manifest.json
resolved_config.yaml
runtime.json
```

8. formal benchmark 使用现有评测体系读取 embedding 并输出正式结果。

## 验收条件

```text
正式主方法列表只出现 caam_scmae
正式输出目录只出现 caam_scmae
caam_scmae 对应 variant full
artifact_manifest status == complete
config_hash 匹配
required_files 全部存在
```

## 不允许

```text
不注册 caam_scmae_control
不注册 caam_scmae_axial
不注册 caam_scmae_advmask
不把 internal ablation 结果混入正式主表
不使用 label-selected resolution
不将 known-K 结果称为 unknown-K
```

