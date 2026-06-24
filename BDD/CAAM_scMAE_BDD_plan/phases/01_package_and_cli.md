# Phase 1: 独立包与 CLI

## 目标

创建 CAAM-scMAE 独立包骨架，不接入正式 benchmark，不实现复杂模型。

目标路径：

```text
methods/DeepLearning/CAAM_scMAE/
```

## 操作步骤

1. 创建包目录与 `__init__.py`。
2. 创建 `run.py`，支持基础 CLI 参数。
3. 创建 `registry.py`，统一解析 `variant`。
4. 创建 `configs/`，至少包含：

```text
control.yaml
axial.yaml
advmask.yaml
full.yaml
benchmark_main.yaml
benchmark_methods.yaml
benchmark_datasets.yaml
```

5. 创建 BDD 指定子目录：

```text
data/
corruption/
mask_generator/
models/
losses/
trainers/
diagnostics/
evaluation/
benchmark/
tests/
docs/
```

6. 实现配置合并：

```text
默认值 < YAML 配置 < CLI 显式参数
```

7. 每次运行保存：

```text
resolved_config.yaml
runtime.json
environment.txt
git_commit.txt
run_manifest.json
```

## 验收条件

```text
python methods/DeepLearning/CAAM_scMAE/run.py --help 可执行
--variant 只接受 control/axial/advmask/full
resolved_config.yaml 可保存
不创建仓库顶层重复 configs/data/tests/scripts
```

## 不允许

```text
不接入 PlantSPADE-LGCL runner
不修改 formal benchmark
不实现标签路径
```

