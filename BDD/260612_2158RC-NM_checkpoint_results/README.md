# RC-NM v4.1 checkpoint 分析与数据索引

本目录保存本次执行 `BDD/260612_1748RC-NM.md` 后的结论和关键数据表。

执行结果来源：

`results/rc_nm_checkpoint_v4_1`

正式 run 数：

- Stage 0: 2 runs
- Stage 1A: 16 runs
- Stage 1B: 50 runs
- Smoke run: 1 run，不纳入正式判定

冻结数据集：

- Wang
- worm_neuron_cell
- Macosko
- Tosches
- Pollen

## 最终判断

最终判定为：

`B_continue_with_narrowed_claim`

关键数值：

- `analytic_RC` vs `fixed_control` 平均 delta ARI: `+0.001561`
- `analytic_RC` vs `random_delta_matched` 平均 delta ARI: `+0.007101`
- `analytic_RC` 在 2/5 个数据集均值上优于 `random_delta_matched`
- `analytic_RC` 在 5/10 个 dataset-seed pair 上优于 `fixed_control`
- random perturbation matching 只有 4/10 个 Stage 1B random-control run 通过 BDD 阈值

严格解释：

RC-NM 目前只能作为 feasibility checkpoint。现有证据不支持声称成熟方法、SOTA、顶会/顶刊就绪，或 backbone-agnostic 已经被实验验证。

允许的收窄主张：

> Reliability-controlled local shrinkage is plausible as a narrow perturbation-controlled regularizer in the scMAE NeighborMix setting, but the reliability mechanism remains weak and needs stronger random matching and cleaner mechanism separation before paper-level claims.

中文：

> 可靠性控制的局部收缩在 scMAE NeighborMix 设置中有一定可行性，但可靠性机制证据仍弱；在进入论文级主张前，需要更强的 random matching 和更干净的机制拆分。

## Stage 1A pseudo objective 结论

`analytic_RC_rec_only` 在两个 sentinel 数据集上都优于 `analytic_RC_full_pseudo`：

- Wang: rec-only 均值 `0.969439`，full-pseudo 均值 `0.967587`
- Macosko: rec-only 均值 `0.300174`，full-pseudo 均值 `0.255821`

因此 Stage 1B 主矩阵使用 `rec_only`。

## 关键问题

1. `analytic_RC` 相比 `fixed_control` 的平均提升只有 `+0.001561` ARI，效果过小。
2. `analytic_RC` 相比 `random_delta_matched` 的优势不稳定，只在 2/5 个数据集均值上占优。
3. `random_delta_matched` 的匹配质量不足，只有 4/10 个 random-control run 通过阈值。
4. Pollen seed 2024 上 `analytic_RC` 相比 `fixed_control` 明显下降，delta ARI 为 `-0.067358`。
5. Main-5 主矩阵没有同管线 `none` baseline，因此 negative transfer vs noMix 和 worst-case drop vs noMix 不能被强解释。

## 文件说明

- `FINAL_REVIEW_SUMMARY.md`: 简短最终审稿式结论。
- `rc_nm_checkpoint_all_runs.csv`: 正式 68 个 run 的完整汇总表。
- `stage1A_pseudo_objective_sentinel.csv`: Stage 1A pseudo objective sentinel 汇总。
- `stage1B_dataset_method_summary.csv`: Stage 1B 数据集 x 方法均值表。
- `stage1B_pairwise_vs_fixed_control.csv`: 每个 dataset-seed-method 相对 `fixed_control` 的差值。
- `stage1B_analytic_RC_decision_by_dataset.csv`: `analytic_RC` 在每个数据集上的判定表。
- `stage1B_ABC_decision.csv`: A/B/C 最终判定表。
- `random_matching_diagnostics.csv`: random-control 匹配质量诊断。
- `analytic_RC_diagnostics.csv`: `analytic_RC` 的 gate、edge、TV 诊断。

## 执行注意

本次运行中发现并修复了一个关键执行错误：初始 runner 未将 `--seed 2024` 传给训练脚本，导致 seed2024 目录实际使用 seed42。修复后已重跑 Stage 1A 和 Stage 1B 的 seed2024 任务，最终表格中的 seed 分布已经正确。
