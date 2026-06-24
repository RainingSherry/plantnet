# CAAM-scMAE BDD implementation plan

本目录整理 CAAM-scMAE 最新 BDD 与补充说明，作为后续实现时的阶段化执行参考。

## 适用范围

目标包：

```text
methods/DeepLearning/CAAM_scMAE/
```

对外正式方法名：

```text
caam_scmae
```

正式 benchmark 只注册最终模型：

```text
caam_scmae = Model C = Axial encoder + Adversarial mask selector
```

内部仍必须实现：

```text
Model 0: Controlled-scMAE
Model A: Axial-scMAE
Model B: AdvMask-scMAE
Model C: CAAM-scMAE
```

Model 0/A/B/C 用于 CAAM 内部消融，不进入正式 benchmark 主方法列表。

## 文档优先级

执行时按以下优先级理解需求：

```text
补充说明 > CAAM-scMAE 完整 BDD v1.0 > 早期 BDD > 早期模型设想说明
```

补充说明已经修正 benchmark 接入方式：

```text
不修改 PlantSPADE-LGCL runner 作为 CAAM 主入口
正式 benchmark 通过 scripts/run_formal_benchmark.py
正式方法注册在 methods/method_manifest.yaml
runtime 注册在 envs/runtime_registry.yaml
```

## 目录内容

```text
execution_order.md
  11 个开发阶段的总顺序与验收边界

benchmark_contract.md
  正式 benchmark、内部 ablation、GPU、unknown-K、artifact 契约

test_matrix.md
  BDD 要求的测试矩阵与阶段归属

risk_and_stop_criteria.md
  高风险点、人工核查项、fail-fast 与止损条件

phases/
  每个阶段的具体操作步骤、产物和验收条件
```

## 总原则

1. 不新增 BDD 未声明的 loss、corruption、value generator 或训练标签路径。
2. 不把 Model 0/A/B/C 混成不可拆模型。
3. 不让 generator 生成 replacement value。
4. 不让真实标签进入训练、donor、context、gene module、early stopping 或 resolution selection。
5. 不把 known-K 结果称作 unknown-K 或完全无监督。
6. 不跳过 shape、label leakage、donor、budget、gradient、context self-exclusion、reproducibility 测试。
7. 遇到 BDD 未覆盖且会改变研究机制的问题，写 TODO 并停止，不自行发明机制继续。

