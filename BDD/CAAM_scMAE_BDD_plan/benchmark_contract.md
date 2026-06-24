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
```

这些模型只用于 CAAM 内部消融。

## 内部消融

内部 ablation runner：

```text
methods/DeepLearning/CAAM_scMAE/benchmark/run_ablation.py
```

建议输出目录：

```text
results/CAAM_scMAE_ablation/
```

内部消融必须包含：

```text
Model 0: control
Model A: axial
Model B: advmask
Model C: full
parameter-matched MLP
```

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
```

配置优先级：

```text
默认值 < YAML 配置 < CLI 显式参数
```

最终配置保存为：

```text
resolved_config.yaml
```

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

