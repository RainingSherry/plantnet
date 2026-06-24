# CAAM-scMAE natural language model overview

本文档用自然语言说明 CAAM-scMAE 的任务、目标、模型框架、理论支撑和实现重点。它不替代 BDD 中的工程约束，而是帮助开发者理解为什么这些约束必须存在。

## 1. 任务

CAAM-scMAE 面向单细胞 RNA-seq 聚类中的自监督表示学习。输入是细胞 × 基因表达矩阵，训练阶段不使用真实细胞类型标签，最终输出每个细胞的 embedding，再交给统一 benchmark 的聚类评测流程。

模型不是直接预测标签，也不是端到端训练分类器。它的目标是通过 masked autoencoding 学到更适合聚类的细胞表示。

## 2. 目标

CAAM-scMAE 解决两个问题。

第一，普通 scMAE 的 MLP encoder 把一个细胞的一整行基因表达直接压成向量，没有显式建模表达矩阵中的 gene axis 与 cell axis。

第二，普通随机 mask 不一定有信息量，模型可能只学会识别扰动痕迹，而不是学习基因上下文关系。

因此 CAAM-scMAE 的目标是：用 Axial encoder 增强二维上下文建模能力，用受约束的 adversarial mask selector 选择更难但合法的 mask 位置，用 matched gene-wise donor corruption 避免单点数值捷径。

## 3. 四模型关系

内部必须实现四个模型：

- Model 0: Controlled-scMAE = MLP encoder + random mask
- Model A: Axial-scMAE = Axial encoder + random mask
- Model B: AdvMask-scMAE = MLP encoder + adversarial mask
- Model C: CAAM-scMAE = Axial encoder + adversarial mask

Model 0/A/B/C 构成严格的 2×2 因子设计。正式 benchmark 主表只注册最终方法 caam_scmae，也就是 Model C。Model 0/A/B 只用于 CAAM 内部消融，不能进入正式主方法列表。

判断 Model C 是否真正有效，不能只看 C 是否优于 baseline，而要看交互项：

Delta_AB = Model C - Model A - Model B + Model 0

只有当 Model C 同时优于 Model A 和 Model B，并且交互项为正，才可以声称 Axial encoder 与 AdvMask selector 存在协同作用。

## 4. 总体框架

训练流程是：原始表达 X 先经过 donor provider 生成 replacement value V 和 eligibility，再由 mask selector 选择 mask M，然后构造扰动输入 X_tilde。模型用 encoder 得到细胞表示 Z，再用 mask head 预测 mask，用 decoder 重建原始表达。

推理阶段只使用 clean X 经过 encoder 得到 Z。最终用于聚类的是 Z，不是 decoder 输出，也不是 mask head 输出。

## 5. Axial encoder 的理论

单细胞表达矩阵天然是二维结构：行是细胞，列是基因。MLP encoder 主要利用同一细胞内其他基因的信息，而 Axial encoder 进一步显式建模两个方向。

Gene-axis attention 在同一个细胞内部的 gene modules 之间进行信息交互，用来学习基因程序之间的关系。

Cell-axis attention 让当前细胞访问一组固定的 context cells，用来学习同一 gene module 在细胞群体中的上下文关系。

context_indices 必须在训练前固定并保存。每个 epoch 可以刷新 context tokens，但不能重新选择 context cells。如果 query cell 出现在 context set 中，必须 self-exclusion，防止 query 读取自己的 clean context。

## 6. AdvMask selector 的理论

AdvMask selector 的作用是选择更有训练价值、更难恢复的 mask 位置。它只能输出 mask logits 和 mask，不能生成 replacement value。

如果 generator 可以自由生成 replacement value，它可能制造极端值、无效扰动、marker gene 攻击或 batch artifact，使模型学习人造噪声，而不是学习上下文结构。

因此 CAAM-scMAE 使用 constrained adversarial mask selection，而不是自由 GAN。

## 7. Matched gene-wise donor corruption

扰动值必须来自同一 gene 的 donor cell：V_ij = X_rj。

这样做的目的是尽量保持每个 gene 自身的边际分布，让单点数值不容易暴露 mask。模型若要识别和重建被扰动位置，就必须利用该细胞其他基因的上下文。

donor 必须满足 r != i，且每个 cell-gene 位置独立采样 donor，不能整行 donor。donor pool 可按 batch、library size、zero ratio 等无标签技术变量匹配，但不能使用真实标签。

## 8. Loss 逻辑

Student 的目标是降低重建损失和 mask prediction 损失。重建损失必须分 masked 和 visible 两部分，并分别按位置数归一化。Mask prediction 使用 BCEWithLogitsLoss，mask head 输出 logits。

Generator 的目标相反：选择让 student 更难的位置。generator loss 中 reconstruction loss 和 mask loss 带负号，同时加入 coverage、distortion、entropy 等约束，防止 mask 集中到少数基因、扰动过强或 mask pattern 坍缩。

## 9. Generator gradient 是硬风险

AdvMask 是否真的成立，取决于 generator 是否能从 loss 收到真实梯度。不能只检查 mask.requires_grad，必须检查 generator 参数的 grad norm 是否大于 0。

Generator step 中 student 参数要 freeze，但 student forward 不能放入 torch.no_grad。原因是 loss 仍然需要通过 student 计算图回传到 X_tilde，再回到 mask_st 和 generator。

Generator step 必须使用连续扰动公式：

X_tilde = X * (1 - mask_st) + V.detach() * mask_st

不得退回 hard torch.where(mask_hard.bool(), V, X)，否则梯度会被截断。

## 10. Context cache 与复现性

Axial 模型最容易出现两个问题：不可复现和 self-copy 泄漏。

必须保证 context_indices 训练前固定，context refresh 只刷新 token，不重新选择 cell；刷新 context cache 时使用 eval mode 和 no_grad；context key/value detach；DataLoader 使用固定 torch.Generator，第一版 num_workers=0；embedding extraction 按原始 index 写回。

同 seed 两次短训练必须检查 context_indices、gene_module_ids、first_batch_indices、first_mask_hard、first_donor_indices、first_loss 和 context_cache_checksum 是否一致或在容差内一致。

## 11. Benchmark 与 ablation 的边界

正式 benchmark 只回答：CAAM-scMAE 作为最终方法是否优于已有方法。因此正式 benchmark 只注册 caam_scmae -> --variant full。

内部 ablation 回答：Axial 是否有用，AdvMask 是否有用，二者是否协同。因此内部 ablation 必须跑 control、axial、advmask、full 和 parameter-matched MLP，但这些不进入正式主表。

## 12. 最容易偏离的地方

实现中不得：让 generator 生成 replacement value；把 donor replacement 写成整行 donor；让 label 进入 training、donor、context、gene module 或 early stopping；把 Model 0/A/B/C 混成一个模型；generator step 使用 hard torch.where；generator step 用 no_grad 包住 student forward；让 query attend 自己的 clean context；把 known-K 结果写成 fully unsupervised；把 control、axial、advmask 注册进正式 benchmark。

## 13. 一句话总结

CAAM-scMAE 的核心是：用受约束的 adversarial mask selection 构造更难但合法的自监督任务，用 axial context encoder 显式建模 gene-axis 和 cell-axis 上下文，从而学习更难被 shortcut 解决、更适合单细胞聚类的 cell embedding。
