# scMAE 改良失败模式与结构迁移问题总结

本文总结当前 independent-full 候选模型快筛中暴露出的主要问题，重点分析为什么大量论文机制迁移到 scMAE 后仍未形成稳定的跨数据集提升。

## 1. 当前快筛现象

截至目前，已实现并筛选了大量 independent-full 候选模型，但仍未出现满足有效门槛的模型。

有效门槛定义为：同一个模型在 Melanoma_5K、Quake_10x_Spleen、Macosko 三个数据集中，至少两个数据集的 ARI 或 NMI 不低于当前 scMAE baseline，并且没有 collapse。

目前观察到的主要模式是：

- Quake_10x_Spleen 最容易提升，多个 graph、diffusion、GRN、causal mask 类方法都能超过 baseline。
- Melanoma_5K 最难突破，当前最佳结果仍略低于 baseline。
- Macosko 只有少数方法在 ARI 上过线，且不能同步提升其他数据集。
- 目前没有任何候选模型能在两个数据集同时过线。

这说明问题不是单纯的训练失败，而是机制迁移后缺乏跨数据集稳定性。

## 2. 评价门槛本身较硬

当前 baseline 已经较强：

- Melanoma_5K：NMI 0.735414，ARI 0.668029
- Quake_10x_Spleen：NMI 0.851730，ARI 0.922275
- Macosko：NMI 0.657465，ARI 0.494268

在 Melanoma_5K 上，最接近的候选距离 baseline 只差约 0.008 NMI 和 0.016 ARI，但仍未过线。

在 Quake_10x_Spleen 上，多数有效机制只带来小幅提升，说明该数据集边界较清晰，轻量去噪或图一致性即可获益。

在 Macosko 上，个别方法 ARI 能过线，但 NMI 仍不足，说明模型可能改善了部分大类结构，却没有稳定改善整体类别信息。

## 3. 数据集结构差异导致机制不稳定

### 3.1 Quake_10x_Spleen

Quake_10x_Spleen 的细胞类型边界相对清楚，局部邻域结构较可靠。因此以下机制容易带来提升：

- graph consistency
- diffusion denoising
- soft graph edge confidence
- pseudo-condition denoising
- GRN/causal dependency consistency

这些机制本质上都在增强相似细胞之间的一致性。

### 3.2 Melanoma_5K

Melanoma_5K 可能包含更强的肿瘤异质性、连续状态、边界细胞和非离散谱系结构。

在这种数据中，全局平滑类机制容易产生副作用：

- 把边界细胞拉向大类；
- 弱化恶性细胞内部异质性；
- 混淆激活状态和细胞类型；
- 强化表达相似但 label 不同的邻域。

因此很多候选能接近 baseline，但无法超过 baseline。

### 3.3 Macosko

Macosko 细胞数大、类别复杂，可能依赖细粒度 subtype、局部 marker 和稀有类信号。

HVG 1000 加全局自监督目标可能更容易学习主轴变化，却不一定保留小类判别信息。全局去噪或全局 graph regularization 也可能压平稀有 subtype。

## 4. 论文机制与 scMAE 主体耦合过浅

许多候选只是将论文方法改造成辅助模块：

- diffusion denoiser；
- edge predictor；
- causal predictor；
- token prediction head；
- graph consistency head；
- teacher/student consistency head。

这些模块训练时 loss 可以下降，但最终 KMeans 使用的仍是 scMAE encoder embedding。如果辅助任务没有强力反向塑造 encoder 几何，最终聚类指标不会明显改变。

典型表现是：某些候选在训练日志中 auxiliary loss 明显下降，但最终 embedding 指标与普通 scMAE 相近。

这说明新机制被“接上了”，但没有真正成为表征形成过程的一部分。

## 5. 迁移时丢失原论文关键上下文

许多论文机制依赖特定前提，迁移到 scMAE 后被简化成 pseudo 信号。

例如：

- diffusion 方法原本依赖完整生成过程、采样策略和强 denoising prior；
- GRN/GAN 方法原本依赖真实 TF 列表、真实或高质量 GRN；
- graph 方法原本依赖可靠邻接图；
- fuzzy/rough 方法原本依赖明确的下近似、上近似和边界集合；
- NLP/CV mask 方法原本依赖稳定 token 或 patch 语义。

迁移后这些前提经常变成：

- pseudo graph；
- pseudo TF；
- pseudo condition；
- coexpression prior；
- SVD/KMeans pseudo labels；
- gene token proxy。

这些 pseudo 结构本身可能有噪声，甚至与真实 label 结构不一致。

## 6. scRNA-seq 不等同于图像 patch 或自然语言 token

许多图像和 NLP 方法迁移时，默认基因维度可以像 patch 或 token 一样处理。但基因表达矩阵有特殊性：

- 基因顺序没有自然空间语义；
- 基因之间没有固定局部连续关系；
- 表达值是连续、稀疏、零膨胀且批次敏感的；
- 同一基因在不同细胞类型中的作用不同；
- HVG 选择会改变输入词表；
- gene-gene 关系具有条件依赖性，非固定静态关系。

因此，直接套用 patch mask、token prediction、row/column attention 等结构，可能形式上合理，但归纳偏置并不完全匹配 scRNA-seq。

## 7. 全局平滑损害 boundary 和 rare-cell

当前很多候选机制本质上都在做平滑：

- graph neighbor consistency；
- diffusion denoising；
- teacher/student consistency；
- GRN edge consistency；
- causal regulator-target consistency；
- soft graph smoothing。

这些机制在结构清楚的数据集上有效，但在复杂数据中可能损害：

- 边界细胞；
- 稀有细胞；
- 连续谱系；
- 肿瘤异质性；
- 激活态和亚型之间的细小差异。

虽然当前 diagnostics 中有 `collapse_warning=false`，但这只说明 embedding 没有整体塌缩，不代表边界和稀有类得到保护。

需要区分：

- 不 collapse；
- 小类不被吞并；
- 边界细胞不过度平滑；
- 稀有 subtype 不被均值化。

目前诊断更多是事后记录，还没有充分进入训练控制。

## 8. pseudo 结构与 label 结构不一定一致

许多候选使用表达相似性或相关性构建 pseudo structure：

- KNN graph；
- coexpression GRN；
- pseudo TF-target mask；
- pseudo condition；
- soft edge confidence；
- SVD/KMeans pseudo label。

但 clustering label 的真实结构可能由多个因素决定：

- 细胞谱系；
- cell type；
- cell state；
- cell cycle；
- tumor/normal 差异；
- rare marker；
- batch effect；
- 技术噪声；
- 平台差异。

表达相关性强不一定意味着对 label 聚类有利。某些 pseudo structure 甚至会强化错误邻域或错误边。

## 9. scMAE 主目标过强，辅助目标影响不足

scMAE 必须保留：

- mask prediction；
- masked expression reconstruction。

这保证了原型模型主体不被破坏，但也带来优化问题：主 reconstruction loss 经常主导训练。

辅助 loss 权重如果太小：

- loss 下降但 embedding 几何几乎不变；
- 新机制成为旁路任务；
- 最终 KMeans 指标不明显变化。

辅助 loss 权重如果太大：

- reconstruction 可能变差；
- embedding variance 可能下降；
- graph/denoise 可能过平滑；
- Melanoma/Macosko 边界可能被破坏。

因此，当前结构缺少稳定的多目标平衡机制。

## 10. 缺少跨数据集自适应控制

同一个固定配置很难同时适配三类数据集。

理想情况下，模型应根据数据诊断动态调整：

- mask ratio；
- graph strength；
- denoise strength；
- edge confidence threshold；
- teacher momentum；
- rare-cell veto；
- boundary penalty；
- auxiliary loss weight。

目前大多数候选使用固定超参，导致：

- Quake 上平滑有益；
- Melanoma 上平滑过强；
- Macosko 上 rare subtype 被削弱。

## 11. 诊断没有形成训练闭环

当前已经输出大量 diagnostics：

- `edge_survival`
- `neighbor_purity_proxy`
- `mixed_cell_fraction`
- `boundary_entropy`
- `rare_risk_fraction`
- `embedding_variance`
- `cluster_mass_min/max`
- `collapse_warning`

但这些多是事后记录，没有充分反向控制训练过程。

理想闭环应为：

1. 训练中监控 boundary、rare-cell、embedding variance、cluster mass。
2. 根据诊断动态调整 mask、mix、graph edge、teacher target 和 loss weight。
3. 当 boundary 或 rare-cell 风险升高时，自动降低平滑强度或关闭混合。
4. 当 embedding variance 降低时，增加 variance/contrastive 约束。

目前多数模型缺少这个闭环。

## 12. 最有希望的已有线索

尽管没有模型达到 2/3 有效门槛，但已有一些方向显示局部优势：

### 12.1 Melanoma_5K

最接近 baseline 的是：

- `rank13_masked_sc_cluster_target_full`
- `rank07_dinobloom_self_distill_full`
- `rank56_scinfovae_mi_zinb_full`
- `rank59_soft_graph_clustering_full`

这说明 Melanoma 可能更需要 cluster target、teacher/self-distillation、稳健 latent regularization，而不是强全局平滑。

### 12.2 Quake_10x_Spleen

过线较多的是：

- `rank60_qdiffusion_latent_denoising_full`
- `rank61_scdiffusion_pseudocond_ddpm_full`
- `rank64_groundgan_causal_mask_full`
- `rank59_soft_graph_clustering_full`
- `rank63_planet_grn_attention_full`

说明 Quake 对 graph/denoise/GRN/causal consistency 友好。

### 12.3 Macosko

最突出的是：

- `rank29_deep_adaptive_fuzzy_clustering_full`

这说明 Macosko 可能更需要 fuzzy/adaptive clustering 和 boundary-aware 策略，而不是单纯重建或全局图平滑。

## 13. 后续设计建议

后续不宜继续单纯“逐篇论文模块迁移”，因为这很可能继续产生单数据集亮点而非跨数据集稳定提升。

更建议围绕现有有效信号做结构重组：

1. 以 `rank13` 的 cluster target 机制作为 Melanoma 主线。
2. 引入 `rank29` 的 fuzzy/adaptive boundary 机制保护 Macosko subtype。
3. 轻量吸收 `rank60/rank64` 在 Quake 上有效的 denoise/causal consistency。
4. 增加 boundary/rare-cell veto，避免全局平滑吞并小类。
5. 将 diagnostics 从事后记录变成训练时 controller。
6. 使用动态 loss weight，而不是固定辅助权重。

一个更合理的下一代候选应具备：

- scMAE 主体；
- mask prediction + masked expression reconstruction；
- cluster-aware semantic target；
- fuzzy boundary/rare-cell protection；
- adaptive auxiliary strength；
- optional graph/denoise consistency；
- diagnostics-driven controller。

## 14. 总结

当前失败的根本原因不是代码实现或运行流程单点问题，而是结构迁移层面的系统性问题：

- 论文机制与 scMAE encoder 耦合不够深；
- pseudo structure 与真实 label structure 不稳定一致；
- 全局平滑损害边界和稀有类；
- scMAE reconstruction 主目标压制辅助目标；
- 缺少跨数据集自适应；
- diagnostics 尚未进入训练闭环。

下一阶段应从“继续迁移新论文”转向“基于已有失败模式重构一个 boundary-aware、rare-cell-aware、adaptive-controller 的 scMAE 变体”。
