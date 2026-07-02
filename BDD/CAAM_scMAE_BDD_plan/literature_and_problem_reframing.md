# Literature and Problem Reframing for CAAM Correction

## 目的

本文件不是综述论文，而是用于约束 Codex 和研究执行：后续修改必须围绕领域中真实存在的问题，而不是继续维护一个工程上漂亮但科学上站不住的 CAAM-v1。

## 关键文献事实

### scMAE

scMAE 使用 processed gene expression matrix，而不是在正文中明确声明 full-gene raw input。其核心 corruption 是：

```text
1. 对每个 gene column 在 cells 之间 random shuffle。
2. 生成 Bernoulli mask。
3. 用 shuffled value 替换 masked entries。
4. 使用 mask predictor 预测哪些位置被扰动。
5. 使用 weighted MSE 同时重建 masked 和 unmasked entries。
```

因此，scMAE-style corruption 不要求每个 masked position 都发生数值变化。zero-to-zero 是可能存在的 label noise，而不是 hard failure。

### scNAME

scNAME 使用 mask estimation、ZINB loss、soft K-means loss 和 neighborhood contrastive loss。它也使用类似 empirical distribution / feature replacement 的思想，并强调高维 scRNA-seq 中无信息基因会加剧 curse of dimensionality。

### scCluBench

scCluBench 的启发不是“所有方法必须 full-gene 输入”，而是：

```text
1. 数据集、预处理、评估协议需要标准化。
2. 方法需要在多数据集、多指标、可解释下游任务中比较。
3. 高稀疏、高维和大规模数据是 scRNA 聚类真实挑战。
4. 表征坍缩、过平滑和 clustering objective 脱节是重要问题。
```

因此 CAAM 主协议应使用统一且合理的 HVG feature space；full-gene input 只能作为 scalability stress test。

### sciLaMA / paired gene-cell designs

sciLaMA 说明外部 gene embeddings 和 paired cell-gene decoder 可以避免单纯 full-gene MLP 的一些表达瓶颈。对 CAAM 的启发是：若未来重写 decoder，应优先考虑 scalable decoder 或 low-rank cell-gene decoder，而不是让参数量随 gene number 近似 O(G^2) 爆炸。

## CAAM-v1 问题重构

早期 CAAM-v1 叙事：

```text
Axial encoder + adversarial mask selector 产生协同提升。
```

当前应冻结该叙事，改为先回答：

```text
在 scRNA-seq masked autoencoding 中，corruption 机制、HVG feature space 和 zero-to-zero noise 如何影响 clustering representation？
```

## 三个 correction hypotheses

### H1: 主输入空间假设

```text
HVG feature space 比 full-gene input 更适合作为主 benchmark 协议。
```

验证：HVG=2000/3000 与 full-gene stress test 分开报告。

### H2: Corruption 假设

```text
matched donor 不一定优于 scMAE-style gene-wise shuffle；nonzero-aware donor 可能提高 effective corruption，但也可能制造 artificial shortcuts。
```

验证：A/B/C corruption triad，只改变 corruption_type。

### H3: AdvMask 假设

```text
AdvMask 只有在最佳 corruption 下稳定优于 random mask，才值得保留。
```

验证：control vs advmask，不跑 Axial。

## 论文路线候选

### 路线 A：方法论文

要求：

```text
1. corrected protocol 后，新 corruption 或 AdvMask 稳定优于 scMAE-style baseline。
2. 多个 development/validation datasets 均有提升。
3. biological interpretation 不下降。
4. runtime 和参数量可接受。
```

### 路线 B：协议分析论文

如果没有稳定新方法，但能系统证明：

```text
1. full-gene vs HVG 对 masked AE 影响很大；
2. zero-to-zero rate 与 mask learning/reconstruction/clustering 有关系；
3. matched donor 和 scMAE-style shuffle 在不同 sparsity 下有边界；
```

则可以写成 masked corruption protocol analysis。

### 路线 C：停止 CAAM

如果 correction 后仍无稳定增益，停止 CAAM，不再添加模块。

## 禁止的论文表述

在完成 Phase 15 之前，不得写：

```text
CAAM-scMAE outperforms scMAE
Axial and AdvMask are synergistic
matched donor is more biologically appropriate and empirically better
nonzero-aware corruption solves dropout
```

只能写：

```text
We evaluate whether ...
We test the effect of ...
We observe that ...
```
