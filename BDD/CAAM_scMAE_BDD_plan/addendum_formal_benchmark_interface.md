# CAAM-scMAE formal benchmark interface addendum

本补充文件用于修正 CAAM-scMAE 与当前仓库正式 benchmark 的接口关系。执行优先级高于 `benchmark_contract.md` 中的泛化表述。

## 1. Formal benchmark 只注册最终方法

正式 benchmark 主表只注册一个方法：

```text
caam_scmae = Model C = Axial encoder + Adversarial mask selector
```

不得在正式 `methods/method_manifest.yaml` 主方法列表中注册：

```text
caam_scmae_control
caam_scmae_axial
caam_scmae_advmask
```

这些 variant 只用于 CAAM 内部消融。

## 2. CAAM run.py 必须兼容当前 formal runner 的真实传参

当前正式 benchmark 由：

```text
scripts/run_formal_benchmark.py
methods/method_manifest.yaml
```

调度。`run_formal_benchmark.py` 调用单个方法时主要传入：

```text
--data_path
--save_dir
--n_clusters
--seed
--epochs
--gpu / --no_cuda
```

因此 `methods/DeepLearning/CAAM_scMAE/run.py` 必须兼容这些参数。

以下参数必须支持，但不能设为 formal benchmark 下的必填项：

```text
--dataset_name
--method_name
--variant
--benchmark_mode
```

formal benchmark 下：

```text
--variant full
--benchmark_mode true
--method_name caam_scmae
```

应通过 `methods/method_manifest.yaml` 的 `extra_args` 传入。若未传入，`run.py` 默认 variant 可为 `full`，dataset_name 可从 `data_path` 或 `save_dir` 推断。

## 3. method_manifest.yaml 注册格式

正式注册建议为：

```yaml
- key: caam_scmae
  name: CAAM-scMAE
  full_name: Constrained Adversarial Axial Masked Autoencoder
  path: methods/DeepLearning/CAAM_scMAE/run.py
  category: DeepLearning
  authenticity: PENDING
  smoke: UNKNOWN
  gpu_policy: PASS
  default_in_formal: false
  framework: PyTorch
  required_artifacts:
    - metrics.json
    - embedding_final.npy
    - labels.npy
    - args.json
    - artifact_manifest.json
  extra_args:
    - --variant
    - full
    - --benchmark_mode
    - "true"
    - --method_name
    - caam_scmae
```

开发期不得直接伪装成：

```text
authenticity: VERIFIED
smoke: PASS
```

只有 smoke test 和 artifact contract 通过后，才允许改为 VERIFIED/PASS。

## 4. runtime registry 不可假定已经存在

`envs/runtime_registry.yaml` 可能不存在。实现时必须支持三种情况：

```text
1. method_manifest.yaml 中 runtime.python 显式指定 Python 路径
2. envs/runtime_registry.yaml 存在并能解析 runtime env
3. 二者均不存在时，早期开发使用当前 sys.executable
```

不得假设 `envs/runtime_registry.yaml` 已经存在。

## 5. GPU 语义必须兼容两种运行模式

如果 `CUDA_VISIBLE_DEVICES` 已经由父进程设置，例如：

```bash
CUDA_VISIBLE_DEVICES=3
```

则 CAAM 子进程内部应使用：

```text
cuda:0
```

并记录：

```json
{
  "physical_gpu": 3,
  "cuda_visible_devices": "3",
  "logical_device": "cuda:0"
}
```

如果 `CUDA_VISIBLE_DEVICES` 不存在，则：

```text
--gpu 表示物理 GPU
必须为 1,2,3,4,5,6
禁止 0,7
```

不得把 benchmark 隔离模式下的逻辑 GPU 与 standalone 模式下的物理 GPU 混淆。

## 6. metrics.json 必须能被当前 formal runner 解析

`metrics.json` 至少必须包含：

```json
{
  "kmeans_known_k": {
    "acc": 0.0,
    "nmi": 0.0,
    "ari": 0.0,
    "f1_macro": 0.0
  }
}
```

可以额外包含：

```json
{
  "leiden_fixed": {...},
  "diagnostics": {...}
}
```

但不能省略 `kmeans_known_k`，否则当前 `run_formal_benchmark.py` 的 summary 可能无法正确解析主指标。

## 7. artifact_manifest.json 是 CAAM 自己的完整性门槛

每次完整运行必须保存：

```text
artifact_manifest.json
```

其中至少包含：

```json
{
  "status": "complete",
  "config_hash": "...",
  "required_files": [
    "metrics.json",
    "embedding_final.npy",
    "labels.npy",
    "args.json",
    "artifact_manifest.json"
  ],
  "variant": "full",
  "seed": 42,
  "data_path": "...",
  "embedding_shape": [0, 0]
}
```

当前 formal runner 只会检查 `method_manifest.yaml.required_artifacts` 是否存在；不会自动检查 `config_hash`。因此：

```text
config_hash 一致性必须由 CAAM run.py 的 resume/skip 逻辑自己检查。
```

## 8. 失败时必须返回非零 exit code

如果 CAAM 发生以下情况：

```text
loss NaN/Inf
label leakage test failure
generator gradient failure
context self-exclusion failure
artifact incomplete
embedding shape invalid
fail-fast condition triggered
```

`run.py` 必须：

```text
1. 写入内部 failure 信息
2. 不得留下完整 required_artifacts 组合
3. 返回非零 exit code
```

否则 `scripts/run_formal_benchmark.py` 可能误判为 success。

## 9. 正式 benchmark 与内部 ablation 必须分离

正式 benchmark：

```text
methods/method_manifest.yaml 只注册 caam_scmae -> --variant full
```

内部消融：

```text
methods/DeepLearning/CAAM_scMAE/benchmark/run_ablation.py
```

内部消融输出建议：

```text
results/CAAM_scMAE_ablation/
```

不得将 control/axial/advmask 作为正式主表方法写入 canonical benchmark 目录。

## 10. 执行要求

Codex 在进入 Phase 11 formal benchmark 接入前，必须先完成本 addendum 的接口约束。Phase 1-10 的内部研发可以先进行，但 `run.py` 的 CLI 和 artifact contract 必须从一开始兼容本文件。
