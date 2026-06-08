# NeighborMix-scMAE Dataset Analysis

## 核心问题

> **NeighborMix-scMAE 的收益是否由数据集的局部邻域可靠性决定？**

## 1. 七数据集 NeighborMix 效果汇总

| Dataset | delta_ari | effect_group | 说明 |
|---|---|---|---|
| **Quake_Smart-seq2_Lung** | +0.0404 | **positive** | NeighborMix 提升 ARI |
| **hrvatin_geo_maintype_counts** | +0.0372 | **positive** | NeighborMix 明显提升 |
| **SRP182008** | +0.0544 | **positive** | NeighborMix 明显提升 |
| **SRP235541** | +0.0018 | neutral | 基本持平 |
| **Pollen** | +0.0015 | neutral | 基本持平 |
| **Wang** | -0.0012 | neutral | 基本持平（边界情况） |
| **SRP171040** | -0.0333 | **negative** | NeighborMix 反而下降 |

## 2. 数据集性质对比

| Dataset | cells | genes | clusters | zero% | knn_purity_k10 | cross_type_k10 | silhouette_pca | imbalance | effect |
|---|---|---|---|---|---|---|---|---|---|
| Wang | 9,519 | 14,561 | 2 | 85.3% | 0.9942 | 0.0058 | 0.4365 | 3.55x | neutral |
| hrvatin | 48,266 | 25,187 | 8 | 94.2% | 0.9951 | 0.0049 | 0.1144 | 26.6x | **positive** |
| Pollen | 301 | 21,721 | 11 | 64.1% | 0.8658 | 0.1342 | 0.2264 | 7.71x | neutral |
| Quake | 1,676 | 23,341 | 11 | 89.1% | 0.8013 | 0.1987 | -0.0685 | 27.7x | **positive** |
| SRP235541 | 27,798 | 53,678 | 18 | 93.7% | 0.8901 | 0.1099 | 0.0861 | 9.34x | neutral |
| SRP171040 | 33,956 | 53,678 | 12 | 95.2% | 0.9184 | 0.0816 | 0.0800 | 13.2x | **negative** |
| SRP182008 | 13,514 | 53,678 | 15 | 97.6% | 0.7397 | 0.2603 | -0.1199 | 9.66x | **positive** |

## 3. NeighborMix 效果与数据集性质的相关性

Spearman 相关系数（对 delta_ari）：

| Property | r | p-value | 解释 |
|---|---|---|---|
| **silhouette_pca** | **-0.607** | 0.148 | silhouette 越低，NeighborMix 越有效 |
| knn_purity_k10 | -0.571 | 0.180 | 纯度越高，反而效果越弱 |
| cross_type_edge_ratio_k10 | +0.571 | 0.180 | 跨类边越多，NeighborMix 越有效 |
| class_imbalance_ratio | +0.464 | 0.294 | 类别越不平衡，NeighborMix 越有效 |
| zero_fraction | +0.321 | 0.482 | 弱正相关 |
| n_cells | -0.071 | 0.879 | 无相关 |
| n_genes | +0.222 | 0.632 | 无相关 |

> **注：相关性均未达统计显著（p > 0.05），因为只有 7 个数据点，功效不足。**

## 4. 关键发现

### 4.1 反直觉的负相关性：silhouette vs delta_ari

**发现**：silhouette_pca 与 delta_ari 呈负相关（r = -0.607）。即：

- **全局聚类分离度差（低 silhouette）的数据集，NeighborMix 反而更有效**
- **全局分离度好（高 silhouette）的数据集，NeighborMix 增益小**

这与"邻域可靠性决定 NeighborMix 有效性"的假设 **相反**。可能解释：

- 高 silhouette 数据集（如 Wang，sil=0.44）本身线性可分，不需要邻域混合
- 低 silhouette 数据集（如 SRP182008，sil=-0.12）需要借助邻域信息来弥补全局结构的不足

### 4.2 KNN 纯度的悖论

**发现**：KNN 纯度与 delta_ari 呈弱负相关（r = -0.57）。最高纯度的数据集（Wang, 0.99；hrvatin, 0.99）中，Wang 为 neutral，hrvatin 为 positive。

**可能的解释**：高纯度邻域意味着可靠的局部信号，但 Wang 是二分类（极简任务），NeighborMix 无法进一步提升；hrvatin 的高纯度来自其可靠的生物学结构，NeighborMix 成功利用了这一优势。

### 4.3 SRP171040：NeighborMix 为什么会失效？

SRP171040（delta_ari = -0.033）是唯一明确 negative 的数据集。

关键性质：
- knn_purity_k10 = 0.9184（高）
- cross_type_edge_ratio_k10 = 0.0816（低）
- silhouette_pca = 0.08（接近零，类边界模糊）
- zero_fraction = 0.9518（极高稀疏度）

**解释**：SRP171040 的邻域结构**看起来可靠**（高纯度），但 silhouette 接近零表明**全局类别分离度差**。NeighborMix 在这类数据上可能错误地混合了边界细胞，引入噪声而非增益。

## 5. NeighborMix 适用边界判定

| Dataset | Effect | 判定依据 | 邻域结构 | 全局分离 |
|---|---|---|---|---|
| hrvatin | **positive** | delta_ari=+0.037 | 极高纯度（0.995） | 中等（sil=0.11） |
| SRP182008 | **positive** | delta_ari=+0.054 | 低纯度（0.74） | 差（sil=-0.12） |
| Quake | **positive** | delta_ari=+0.040 | 中等纯度（0.80） | 差（sil=-0.07） |
| SRP235541 | neutral | delta_ari=+0.002 | 较高纯度（0.89） | 差（sil=0.09） |
| Pollen | neutral | delta_ari=+0.002 | 高纯度（0.87） | 中等（sil=0.23） |
| Wang | neutral | delta_ari=-0.001 | 极高纯度（0.994） | 好（sil=0.44） |
| SRP171040 | **negative** | delta_ari=-0.033 | 高纯度（0.92） | 差（sil=0.08） |

## 6. 下一步模型改进方向

基于以上分析，NeighborMix 的改进应聚焦于：

### 6.1 可靠性感知混合（Reliability-Aware NeighborMix）

当 KNN 纯度高但 silhouette 低时（即"局部一致但全局模糊"），NeighborMix 可能产生误混合。

**改进**：引入**邻域可靠性评分**：

```python
reliability = silhouette_pca * knn_purity_k10
mix_weight = sigmoid(alpha * reliability + beta)
```

在可靠性低的边界区域降低混合权重。

### 6.2 边界感知混合（Boundary-Aware NeighborMix）

SRP171040 的问题是边界细胞被错误混合。

**改进**：使用 HDBSCAN 或密度峰值检测识别边界细胞，仅对核心区域细胞应用 NeighborMix。

### 6.3 自适应混合权重（Adaptive Mixing Weight）

当前 NeighborMix 使用固定权重。根据数据集性质自适应调节：

- 高稀疏度（zero_fraction > 0.95）→ 降低混合权重
- 低 silhouette → 根据 silhouette 动态调节
- 类别极度不平衡（imbalance > 20x）→ 考虑分层采样后混合

### 6.4 跨类型邻居检测与过滤

发现 cross_type_edge_ratio 与 delta_ari 呈正相关（r=0.57）。这意味着：

- 跨类型边多的数据集 NeighborMix 反而更好（看似矛盾）
- 但这可能是因为 SRP182008（跨类型边 0.26）提供了强信号

**改进**：在构建 KNN 图时，优先选择**类内最近邻**作为混合伙伴，避免跨类型邻居的干扰。

## 7. 结论

**NeighborMix 的有效性并非由单一因素决定**，而是数据集局部邻域可靠性与全局聚类结构的交互结果：

1. **高纯度 + 低 silhouette → NeighborMix 有效**（hrvatin、Quake）：局部信息弥补全局不足
2. **低纯度 + 差 silhouette → NeighborMix 有效**（SRP182008）：需要邻域信息才能分类
3. **高纯度 + 高 silhouette → Neutral**（Wang）：已是强基准，增益空间小
4. **高纯度 + 低 silhouette → Negative**（SRP171040）：局部高置信但全局不可分，混合引入噪声

**核心假设部分成立**：邻域可靠性确实影响 NeighborMix 效果，但需要结合全局分离度（silhouette）才能完整解释。单纯的 knn_purity 或 cross_type_edge_ratio 都不能可靠预测 NeighborMix 的有效性。
