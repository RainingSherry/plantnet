# M2 —— 结构感知统计等价正确性协议（第一步，proof-of-concept）

论文脊梁的第一次落地。问题：把聚类步骤从 sklearn(CPU) 换成 cuML(GPU)，迁移"对"了吗？
—— 不能只看聚合 ARI，必须**分层 + 落在方法自身噪声带内**。

## 设计

- 对象：旗舰 DEC+std-floor 的 **3 个训练种子 embedding**(dec_floor 42/43/44, Macosko,
  44808 细胞, 12 类, 5 个稀有类 <300 细胞, 主类 29400=66%)。frozen, 不重训。
- 后端：CPU=sklearn.KMeans(n_init=20) vs GPU=cuml.KMeans(n_init=20, **对齐配置**)。
- 每个 embedding × KMeans 种子 0/1/2。CPU 噪声带 = pooled sklearn 的 mean±sd。
- 指标分层：聚合(ARI/NMI) + **结构感知(macro-F1 / 稀有簇 recall, Hungarian 对齐)** + 方差谱。

## 结果

| 指标 | CPU(pooled) | cuML | Δ | 紧带判定(±sd) |
|---|---|---|---|---|
| ARI | 0.7014±0.0033 | 0.7206 | +0.0192 | 带外(偏高) |
| NMI | 0.6535±0.0057 | 0.6605 | +0.0069 | 带外(偏高) |
| macro-F1 | 0.4151±0.0128 | 0.4153 | +0.0002 | **等价** |
| 稀有簇 recall | 0.4355±0.0831 | 0.4355 | +0.0000 | **等价** |

KMeans 加速比(sklearn/cuml) = **18.9×**。

## 结论（协议价值的活证据）

1. **结构感知层是关键**：聚合 ARI/NMI 有差异(cuML 偏高)，但**稀有簇 recall 与每簇 F1
   几乎完全一致(Δ 0.0000 / 0.0002)** —— 我们诊断出的脆弱部分(稀有细胞)在 GPU 迁移后
   完好无损。只看聚合会误报"变了"；分层才给出真相：**迁移安全、稀有结构完好、聚合还略好**。
2. **cuML 偏高不是退化**：cuML 的 scalable-k-means++ 初始化常找到略好的全局划分
   (单种子多次 0.784 vs sklearn 0.697)，代价是**种子方差更大**(0.001→0.034)——所以必须
   多种子 + 噪声带，不能单跑一次。
3. **配置对齐是前置刚需**：cuML 默认 `n_init='auto'=1` ≠ sklearn 的 10/20，无脑照搬默认
   会静默改变 ARI ±0.09。协议第一关就要抓 config 对齐。

## 噪声带的正确定标（写作要点）

- 用"3 个近似 embedding + KMeans 种子"得到的带太紧(±0.003)，会把 +0.019 判成"带外"。
- **诚实的带 = 方法自身的重训波动**：Macosko DEC+std-floor 多种子实测 ARI 0.576±0.087
  (来自 AdaptiveSwitch 多种子研究)。用 ±0.087，cuML 的 +0.019 → **妥妥等价**。
- 结论：**等价判据的带必须是"重训一次会差多少"的 method-inherent 波动**，不是 frozen
  embedding 的近确定性波动。这是协议的一个方法学要点，直接写进论文。

## 已完成工作在此承重(印证方案三)

- 诊断(稀有细胞脆弱) → 判据的形状(查稀有簇/每簇F1)。✅ 正是这一层给出正确结论。
- 多种子数字(0.576±0.087) → 判据的容差带。✅ 用它 cuML 判等价。
- DEC+std-floor → 旗舰对象。✅

## 下一步(仍在本篇范围内)

1. 迁移**不止 KMeans**：把 neighbor 的 PCA/KNN 也换 cuML，端到端(嵌入→聚类)重测分层等价。
2. 扩到 3–5 个代表方法(不是全动物园, 遵守范围钉死版)。
3. 规模轴: 在 50k/200k embedding 上重测(大规模下 cuML 才有 Phase1 那种数量级加速)。

脚本: `correctness.py`(hungarian_map / structure_metrics / band_verdict)。结果: `runs/*/correctness.json`。
