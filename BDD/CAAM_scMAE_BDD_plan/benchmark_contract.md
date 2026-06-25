# CAAM-scMAE benchmark contract

## 正式接入方式

正式 benchmark 入口：

```text
scripts/run_formal_benchmark.py
methods/method_manifest.yaml
envs/runtime_registry.yaml
```

禁止把 CAAM 主入口接入：

```text
methods/DeepLearning/PlantSPADE_LGCL/scripts/run_single.py
```

该 runner 属于 PlantSPADE-LGCL 协议，不作为 CAAM-scMAE 主 benchmark 入口。

## 正式方法注册

正式主表只注册：

```text
caam_scmae
```

调用方式：

```text
methods/DeepLearning/CAAM_scMAE/run.py --variant full
```

不得注册到正式主表：

```text
caam_scmae_control
caam_scmae_axial
caam_scmae_advmask
caam_scmae_scmae_shuffle
caam_scmae_matched_donor
caam_scmae_nonzero_aware
```

这些模型或 corruption settings 只用于 CAAM 内部消融、correction pipeline 或 supplementary analysis。

## 主 feature-space 协议

主 benchmark / quick ablation / development ablation 默认使用：

```text
input_mode = log1p
n_top_genes = 2000
scale_input = false
```

`n_top_genes=0` 只允许用于：

```text
1. full-gene stress test
2. 上游已提供 external_hvg matrix，并且 artifact 明确记录 feature_space_source=external_hvg
```

不得把 full-gene stress test 结果放入主表，不得用 full-gene stress test 决定模型机制是否成立。

## 内部消融

内部 ablation runner：

```text
methods/DeepLearning/CAAM_scMAE/benchmark/run_ablation.py
```

建议输出目录：

```text
results/CAAM_scMAE_ablation/
```

历史 2x2 factorial 内部消融包含：

```text
Model 0: control
Model A: axial
Model B: advmask
Model C: full
parameter-matched MLP
```

Correction Pipeline 中，必须先执行 corruption triad 和 AdvMask triage，再决定是否恢复 2x2 factorial。

正式 benchmark 输出目录只应出现：

```text
caam_scmae
```

## CLI 契约

最低支持：

```bash
python methods/DeepLearning/CAAM_scMAE/run.py \
  --data_path PATH \
  --save_dir PATH \
  --n_clusters K \
  --seed SEED \
  --gpu GPU \
  --variant control|axial|advmask|full
```

必须支持：

```text
--config
--data_path
--save_dir
--dataset_name
--method_name
--variant
--n_clusters
--seed
--gpu
--no_cuda
--benchmark_mode
--input_mode
--n_top_genes
--target_sum
--scale_input
--skip_eval
--resume
--overwrite
--epochs
--corruption_type scmae_shuffle|matched_donor|nonzero_aware_donor
--strict_effective_budget true|false
```

配置优先级：

```text
默认值 < YAML 配置 < CLI 显式参数
```

最终配置保存为：

```text
resolved_config.yaml
```

## Corruption contract

必须支持三种 corruption：

```text
A. scmae_shuffle: gene-wise shuffle，接近 scMAE/scNAME 的经验分布替换。
B. matched_donor: 当前匹配 batch/library_size/zero_ratio 的 donor 替换。
C. nonzero_aware_donor: 在 donor 值中优先选择能产生非零或数值变化的位置。
```

共同约束：

```text
replacement value 必须来自同一 gene
不得由 generator 生成 replacement value
不得使用 label/cell_type/n_clusters 选择 donor
必须记录 zero_to_zero_rate、effective_corruption_rate、mean_abs_delta、budget_deficit_rate
```

`budget_deficit_rate` 默认是 diagnostic metric，不得默认 fail-fast。只有显式设置：

```text
--strict_effective_budget true
```

才允许因为 effective budget deficit 超阈值而退出。

## GPU 语义

Standalone 模式：

```text
--gpu 表示物理 GPU
允许 1,2,3,4,5,6
禁止 0,7
```

Benchmark 模式：

```text
父进程设置 CUDA_VISIBLE_DEVICES=<physical_gpu>
子进程内部固定使用 cuda:0
--gpu 只作为逻辑 GPU 或元信息记录
不得把 --gpu 误认为物理 GPU
```

必须写入 `runtime.json`：

```json
{
  "physical_gpu": 3,
  "cuda_visible_devices": "3",
  "logical_device": "cuda:0"
}
```

## n_clusters 限制

`n_clusters` 只能用于评测元信息，不得影响：

```text
训练
掩码生成
context 选择
gene module 生成
donor selection
loss
early stopping
```

必须测试：

```text
相同 seed、相同数据、不同 n_clusters
首个 batch mask 相同
首个 forward 输出相同
首个 student loss 相同
```

## Unknown-K 主协议

主表只保留：

```text
leiden_fixed
```

不得把以下内容放入主表：

```text
silhouette-selected resolution
label-selected resolution sweep
oracle sweep
```

如做 sweep，只能作为 supplementary/oracle，不得用于主结果或调参。

## Artifact manifest

每次完整运行必须保存：

```text
artifact_manifest.json
```

只有同时满足以下条件，正式 benchmark 才允许 skip training：

```text
status == complete
config_hash 匹配
required_files 全部存在
```

只存在 `embedding_final.npy` 不足以认定运行完整。

## 协议变更后的 smoke 规则

如果任何修改改变以下行为：

```text
n_top_genes 默认值
corruption_type 默认值
strict_effective_budget 默认值
artifact schema
formal manifest extra_args
```

则旧 smoke 只能标记为 old-protocol smoke。必须重新运行：

```text
1-seed GPU formal smoke
3-seed GPU formal smoke
validate_formal_smoke.py
```

之后才允许维持或恢复 `smoke: PASS`。
