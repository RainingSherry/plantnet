# CAAM Correction Pipeline BDD Index

## 为什么需要本修订

CAAM-v1 的工程 smoke 已通过，但研究判断不能继续沿用早期假设。已有结果显示：

```text
1. full-gene benchmark_mode 使 14k-23k genes 输入成为默认，造成参数量和 mask budget 问题；这不是 scMAE/scNAME 的典型主协议。
2. strict effective budget 把 zero-to-zero corruption 当作硬失败，而 scMAE/scNAME 的 shuffle corruption 并不要求每个 mask 位置都发生数值变化。
3. matched donor 目前只是研究假设，不是已验证优势。
4. subsample_2k 早期 5-model ablation 为 no_positive_interaction。
5. 三个 development datasets 在旧协议下均 blocked_by_budget_deficit，但该结论主要说明旧协议过严，不能直接判死模型。
```

因此，本目录新增 Phase 12-16，用于一次性整理并执行较优先的问题，而不是继续碎片化修补。

## 总体策略

```text
先修协议，再判断机制。
先判断 corruption，再判断 AdvMask。
先判断 AdvMask，再让 Axial 重新进入。
先用 development datasets 冻结方案，再碰 validation。
sealed test 最后只运行一次。
```

## 分阶段执行顺序

### Phase 12: Protocol correction

目标：修正主输入协议和 budget 语义。

必须完成：

```text
1. benchmark_mode 默认 n_top_genes=2000，而不是 0。
2. n_top_genes=0 只作为 full-gene stress test。
3. strict_effective_budget 默认 false。
4. budget_deficit_rate / zero_to_zero_rate / effective_corruption_rate 作为 diagnostics。
5. formal smoke validator 兼容新协议。
```

### Phase 13: Corruption triad

目标：比较三种 corruption，不再默认坚持 matched donor。

必须比较：

```text
A. scmae_shuffle: scMAE-style gene-wise shuffle
B. matched_donor: matched donor shuffle
C. nonzero_aware_donor: nonzero/change-aware donor shuffle
```

第一轮只跑：

```text
MLP encoder + random mask
```

不要跑 Axial，不要跑 full。

### Phase 14: AdvMask triage

目标：判断 AdvMask 是否有稳定价值。

只比较：

```text
control vs advmask
```

使用 Phase 13 中最好的 corruption 或前两名 corruption。

### Phase 15: Axial re-entry

目标：只有在 AdvMask 通过后，才允许 Axial 重回主实验。

必须证明：

```text
1. Axial 优于 parameter-matched MLP，或至少不是参数量导致。
2. full 模型优于 axial 与 advmask。
3. Delta_AB > 0，且 paired CI 支持。
```

### Phase 16: Publication decision

目标：决定论文路线。

可能路线：

```text
A. 正向方法论文：新 corruption / AdvMask / scalable design 稳定提高聚类与生物解释。
B. 协议分析论文：主要贡献是发现 masked corruption 在 scRNA-seq 中受 HVG/zero-to-zero/feature space 强烈影响。
C. 停止 CAAM：若 corrected protocol 仍无稳定增益。
```

## 每个阶段都必须检查上一阶段

Codex 执行每份 BDD 前，必须先写一个 `previous_phase_check` 小节，检查上一阶段是否满足：

```text
1. 所有要求文件是否存在。
2. resolved_config / artifact schema 是否包含新增字段。
3. 是否误改 formal manifest、authenticity、default_in_formal。
4. 是否提交 results/ 或 data/smoke/*.h5ad。
5. 是否使用 sealed test 进行调参或机制选择。
6. 是否根据 ARI/NMI 自动改模型结构。
```

若上一阶段未满足，不得继续下一阶段。

## 不允许的行为

```text
1. 不得把 corruption variants 注册为 formal methods。
2. 不得把 control/axial/advmask 注册进主表。
3. 不得继续用 full-gene input 作为主表协议。
4. 不得默认 strict effective budget。
5. 不得在 AdvMask 未通过前跑 full 模型并解读 synergy。
6. 不得用 validation/test 选择 corruption_type、mask_ratio、HVG 数量或是否保留 Axial。
7. 不得用 known-K 结果包装为 unknown-K。
```
