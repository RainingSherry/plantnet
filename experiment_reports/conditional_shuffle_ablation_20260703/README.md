# 条件化 / nuisance 匹配 swap-noise 消融（Phase 2）

外部消融脚手架（位于 `experiment_reports/`，不在 `methods/`）。骨干 = 已验证的跨数据集赢家
**scMAE + DEC + 每维 std-floor**（`AdaptiveSwitch_scMAE`，`variance_weight=0.02`，`force_gate=1.0`）。
路线图见 `experimental_retired_models/scmae_structural_pivot_20260703/ROADMAP.md` 的 Section 5 Phase 2。

## 假设（这是本实验的科学核心）

```text
原始 scMAE 用 swap-noise（交换噪声）做 corruption：对每个基因列，从「全体细胞」里随机抽
一个该基因的值来替换目标位置。问题是——「这个值是否被替换过」常常能靠技术轴判断出来：
文库大小（library size）、检测基因数（n_detected）、零率（zero-rate）。这是一条捷径。

若把供体池限制在「同一 nuisance 分箱」内，就抹掉了这条捷径，模型无法再靠技术轴识别被换的
位置，只能去学更细的 gene-gene 条件结构。若这样能提升聚类，说明我们第一次从 scMAE 内部
（corruption 任务本身）把它改对了，而不是又在外面挂一层 X 派生的伪结构。
```

## swap-noise vs zero-masking（关键区别）

赢家骨干当前用的是 **zero-masking**：`AdaptiveSwitchScMAE.random_mask` 把被选中的位置直接填
`0.0`（`masked_fill(mask.bool(), 0.0)`）。原始 scMAE 论文用的是 **swap-noise**：被选位置换成
其它细胞的同基因值。Phase 2 的假设针对的正是 swap-noise 的**供体来源（donor pool）**，所以本
runner 先实现真正的 per-gene independent swap-noise（区别于 `scMAE_family.apply_scmae_noise`
那种「整行换同一个供体细胞」的近似），再让它的供体池可条件化。

唯一被改动的变量是 corruption 的**替换值来源**；`mask_prob`、DEC、std-floor、`force_gate=1`、
`variance_weight=0.02` 等一切保持与赢家一致，保证干净归因。返回的 mask 指示矩阵语义完全不变
（1=该位置被选中 corrupt），因此 `loss.py` 的 BCE mask-discriminator 和加权重构无需任何改动。

## 四个 arm + 一个复现对照（`--corruption`）

| arm | 值 | 供体池 | 含义 |
|---|---|---|---|
| 复现赢家 | `zero` | —— | 零填充 mask，= 现有赢家 `random_mask`，作为「复现赢家」对照 |
| **S0** | `swap_global` | 全体细胞 | 复现原始 scMAE 全局 swap-noise，**swap 基线** |
| **S1** | `swap_lib` | 同 library-size 分箱 | 去掉「靠文库大小判断被换位置」的捷径 |
| **S2** | `swap_ndet` | 同 n_detected 分箱 | 去掉「靠检测基因数判断」的捷径 |
| **S3** | `swap_zerolib` | 同 (zero-rate × library) 联合分箱 | 去掉「靠零率×文库联合判断」的捷径 |

## nuisance 变量与分箱（务必从 UNSCALED 原始数据算）

- **nuisance 必须从未 scale 的原始 counts 算，绝不能从编码器输入（scale 后）算** —— 否则
  library-size 这条技术轴已被总计数归一化抹掉，整个 Phase 2 假设就不成立。本 runner 重新读
  h5ad 的原始 counts（`layers[counts]` / `raw.X` / raw-looking `X`）计算：
  - `library_size` = 每细胞原始总 counts（全基因）
  - `n_detected` = 每细胞原始 counts>0 的基因数（全基因）
  - `zero_rate` = 每细胞在 HVG 空间的零率（= 编码器实际看到的稀疏度）
- 分箱用**等频分箱（按秩切分）**：单轴 `n_nuisance_bins=10`（S1/S2），联合轴每轴 `n_joint_bins=5`
  （S3，共 25 箱）。按秩切分对大量并列值稳健，保证箱大小均衡。
- **corruption 施加在编码器输入空间**：箱标签用未缩放量算，但实际交换的值取自 encoder 输入矩阵里
  同基因列、同箱内另一细胞的值（与 `model.random_mask` 施加空间一致）。

## donor-pool 监控（必须看的失败模式）

箱太窄 → 供体池太小 → swap 退化成近似恒等变换。每个 arm 都在 `summary.json` /
`neighbor_profile.json` 记录各箱平均/最小/最大供体池大小与 cell-weighted 平均。**若 `pool_min`
掉到个位数，该 arm 结果不可信**，应减少箱数（加宽箱）。

### 稀疏数据的重要现象：`eff_change << mask_prob`

runner 同时记录两个诊断量（对所有 arm 同一定义，可比）：
- `designed_mask` = 设计选中率，应 ≈ `mask_prob=0.4`（各 arm 相同，证明 corruption 预算一致）；
- `eff_change` = 实际数值改动率。

在稀疏 scRNA 上（Macosko HVG 零率中位数 ≈ **0.976**），swap-noise 的实际改动率会远低于
`mask_prob`：大量被选位置是「0 换 0」的空操作。经验公式 `eff_change ≈ mask_prob·(1−zero_rate²)`。
Macosko 上 ≈ `0.4·(1−0.976²) ≈ 0.019`，与实测 ≈ 0.023 吻合。**这不是 bug，而是 swap-noise 在
稀疏数据上的真实性质** —— 有效扰动集中在原本被检测到（非零）的基因上。`zero` 填零 arm 因为把
非零位置也改成了 0，实际改动率反而更贴近 `mask_prob`。这两种 corruption 在稀疏数据上力度天然不
同，是 arm 之间的固有差异，不是实现偏差；跨 swap arm（S0/S1/S2/S3 同机制）的相对比较仍然有效。

## 运行

```bash
# dry-run：打印完整 45 条 run 命令（5 corruption × 3 数据集 × 3 种子），GPU 可参数化
bash experiment_reports/conditional_shuffle_ablation_20260703/run_all.sh 3
# 真跑（顺序执行，跳过已完成）：
bash experiment_reports/conditional_shuffle_ablation_20260703/run_all.sh 3 go
# 或交给外部调度器逐条分发到不同 GPU（自行覆盖 --gpu）：
bash experiment_reports/conditional_shuffle_ablation_20260703/run_all.sh | parallel -j6
# 聚合：
python experiment_reports/conditional_shuffle_ablation_20260703/summarize.py
```

数据集：Macosko(k=12)、Melanoma_5K(k=9)、Quake_10x_Spleen(k=5)；种子 42/43/44；epochs=80。
GPU 只用 1-6，禁止 GPU 0/7（runner 的 `get_device` 已有保护）。

## 判据（见 SUMMARY.md 自动计算的 delta）

```text
S1/S2/S3 > S0 (swap_global)  ARI +>=0.02 且多种子稳定  -> 限制供体池真的有效，从 scMAE 内部改对了
S1/S2/S3 ~= S0                                          -> nuisance 匹配无效，swap 捷径不是瓶颈
swap_global(S0) vs zero                                 -> 先确认 swap-noise 本身相对零填充的效应方向
donor_pool_min 个位数 / eff_change 远低于同数据集其它 swap arm -> swap 退化，该 arm 不可信
```

判据锚定在 **相对 S0** 的 delta（这才是 Phase 2 的真正问题：限制供体池有没有用），
同时报告相对 `zero`（复现赢家）的 delta 作为总体参照。

## 防泄露声明（见 ROADMAP.md Section 4）

- 不用测试 ARI/NMI 选任何超参（箱数、mask_prob、loss 权重、epoch 全部固定沿用赢家默认）。
- label 只在最终评测用一次（KMeans known-K）。
- 多数据集固定协议（Macosko + Melanoma + Quake），不做单数据集调参。
- 报告 ARI 时并列 NMI；固定 K 下比较。
