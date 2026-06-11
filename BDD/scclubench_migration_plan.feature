# language: zh-CN

# 基本信息
# 项目: plantnet dimension-reduction benchmark
# 本地路径: /home/luolie/biopipeline/dimension-reduction/plantnet
# 源代码目录: OtherMode/scCluBench-main
# GitHub 源路径: https://github.com/RainingSherry/plantnet/tree/refactor/minimal-method-isolation-and-neighbormix-ablations/OtherMode/scCluBench-main
# 目标目录: methods/
# 创建日期: 2026-06-11
# 作者: Codex
# 状态: Draft
#
# 当前观察:
# - scCluBench 源模型目录包含 11 个主要模型:
#   DeepLearning: dec, desc, scDCC, scDeepCluster, scMAE, scNAME, scziDesk
#   GNN: AttentionAE-sc, scCDCG, scDSC, scGNN
# - methods/method_manifest.yaml 当前有 20 个方法条目，其中 VERIFIED=15,
#   ENV-GATED=1, PENDING_AUDITED=1, PLACEHOLDER=3。
# - results/formal/benchmark_summary_all_current.csv 当前是成功运行的 canonical 表:
#   18 个数据集、12 个方法、176 个 dataset-method 行、461 个成功 run。
# - 当前未全覆盖的典型情况:
#   attentionae_sc 8/18 datasets; scdsc 2/18 datasets;
#   scname/sczidesk 11/18 datasets; desc/sccdcg/scgnn/scdeepcluster 无 canonical 全量行。
#
# 初步原因分类:
# - 迁移完成不等于全量运行完成: SOURCE_MIGRATED、SMOKE_PASS、FULL_COVERAGE 是不同状态。
# - 运行环境冲突: TensorFlow、旧 PyTorch、Scanpy、rpy2/R 等依赖不能放在一个环境中。
# - 部分模型只被脚本安排在 scMAE 11 数据集上运行，不等价于 18 数据集 full benchmark。
# - scGNN 属于慢速/不确定模型，小数据 smoke 也可能 600s+，需要单独性能门槛。
# - scDeepCluster 仍为 ENV-GATED，不能因为代码已复制就进入正式结果。
# - 部分历史/临时运行脚本绕过统一 runner、使用旧 GPU 策略或硬编码环境，导致覆盖和日志不一致。

功能: 将 OtherMode/scCluBench-main 可靠迁移到 methods 并推进全数据集覆盖
  为了让迁移后的模型可审计、可复现、可在统一 benchmark 中比较
  作为 benchmark 维护者
  我希望每个 scCluBench 模型都经过源代码迁移、运行契约、环境隔离、smoke、子集运行、全量运行的分阶段验收

  背景:
    假如 当前项目根目录为 "/home/luolie/biopipeline/dimension-reduction/plantnet"
    并且 迁移源目录为 "OtherMode/scCluBench-main"
    并且 迁移目标目录为 "methods"
    并且 统一运行入口为 "scripts/run_formal_benchmark.py"
    并且 当前 canonical 结果表为 "results/formal/benchmark_summary_all_current.csv"

  @P0 @inventory
  场景: 建立源模型到目标模型的一一迁移清单
    假如 源目录中存在 scCluBench 模型目录
    当 维护者生成迁移清单时
    那么 每个源模型都应有唯一 target_path、method key、category 和 source_path
    并且 "AttentionAE-sc" 这类非法 Python 包名应记录规范化目标名 "AttentionAE_sc"
    并且 迁移清单应区分 proposed/ablation、本项目新增模型、scCluBench 原始模型和 placeholder
    并且 未迁移、迁移中、已迁移但未验证的模型不能被标为 FULL_COVERAGE

  @P0 @state-model
  场景: 迁移状态不得与运行覆盖状态混用
    假如 一个模型的源代码已经复制到 methods
    当 维护者更新 manifest 或文档时
    那么 该模型最多只能标记为 SOURCE_MIGRATED
    并且 只有通过 smoke run 后才能标记为 SMOKE_PASS 或 VERIFIED
    并且 只有 18 个数据集按要求 seed 全部成功后才能标记为 FULL_COVERAGE
    并且 scMAE11 子集成功不能替代 18 数据集 full benchmark 成功

  @P0 @no-othermode-runtime
  场景: 迁移后的方法运行时不得依赖 OtherMode
    假如 一个方法被纳入 methods
    当 执行该方法的 run.py
    那么 Python import、数据加载、模型定义和训练逻辑都应来自 methods 内部或标准依赖
    并且 运行命令、PYTHONPATH 和配置文件不得指向 "OtherMode/scCluBench-main"
    并且 authenticity.json 应记录 source_path 和 target_path 用于追踪，而不是运行时依赖

  @P0 @runtime
  场景: 每个模型必须声明并验证运行环境
    假如 一个模型需要 TensorFlow、旧版 PyTorch、R/rpy2 或特殊 CUDA 组合
    当 维护者将它加入 method_manifest.yaml
    那么 manifest 应声明 runtime.env_name
    并且 envs/runtime_registry.example.yaml 应包含对应逻辑环境
    并且 本地 runtime_registry.yaml 缺少该环境时 runner 应 fail fast，不应静默 fallback 到 sys.executable
    并且 environment.json 和 status.json 应记录实际 python_executable 与 runtime_env

  @P0 @authenticity
  场景: manifest、core card 和迁移文档状态保持一致
    假如 method_manifest.yaml 中某方法为 VERIFIED
    当 维护者查看 docs/model_core_cards 和 docs/migration_status.md
    那么 三处的 authenticity、smoke、runtime、known_deviations 应一致
    并且 不应出现 manifest 已 VERIFIED 但 core card 仍写 PENDING_AUDITED 的冲突
    并且 状态冲突时该模型不得进入默认正式运行列表

  @P0 @data-contract
  场景: 所有迁移模型使用统一数据输入契约
    假如 输入数据为 h5ad 或由 h5 转换得到的 canonical h5ad
    当 任何迁移模型运行时
    那么 它应通过统一数据预处理或显式记录自己的预处理策略
    并且 必须解析 resolved_label 或明确 label_key
    并且 h5 转换缓存必须校验源文件路径、大小、mtime 和转换参数
    并且 不允许用测试标签进行训练决策或 checkpoint selection

  @P0 @output-contract
  场景: 所有可汇总运行必须满足统一输出契约
    假如 一个方法运行结束且 return_code 为 0
    当 runner 校验输出目录时
    那么 至少应存在 metrics.json
    并且 若 manifest 声明 required_artifacts，则必须全部存在
    并且 metrics.json 应包含 acc、nmi、ari、f1_macro、fmi、v_measure、homogeneity、completeness
    并且 status.json 应明确 success、failed、timeout、incomplete 或 skipped
    并且 incomplete 不应被汇入 mean/std 主表

  @P0 @gpu-policy
  场景: 迁移运行脚本不得绕过统一 GPU 策略
    假如 批量运行脚本需要分配 GPU
    当 它调用任何模型
    那么 它应通过 scripts/run_formal_benchmark.py 调度
    并且 不得使用被禁止 GPU
    并且 不得在临时脚本中硬编码与 runner 冲突的 GPU 策略
    并且 CPU-only 运行必须显式传入 --no_cuda 并记录在 status.json

  @P1 @smoke
  场景: 每个迁移模型先通过最小 smoke run
    假如 一个模型完成 SOURCE_MIGRATED
    当 维护者尝试升级它的运行状态
    那么 应先在小数据集上运行 epochs=1、pretrain_epochs=1、seed=42 的 smoke
    并且 smoke 应验证命令行参数、运行环境、输出契约和无 OtherMode 依赖
    并且 smoke 失败时应记录失败类型为 env_error、data_error、model_error、output_contract_error 或 timeout

  @P1 @subset
  场景: scMAE 11 数据集子集运行只能证明子集覆盖
    假如 一个模型只在 data/scMAE 或 data/processed_scmae 的 11 个数据集上运行
    当 汇总覆盖率时
    那么 coverage_scope 应标记为 scMAE11
    并且 不能在 README、manifest 或结果表中称为 all_datasets
    并且 进入 18 数据集 full benchmark 前应创建补跑任务覆盖缺失的 7 个非 scMAE11 数据集

  @P1 @full-coverage
  场景: 18 数据集 full benchmark 的完成条件
    假如 某模型目标是全数据集比较
    当 查看 canonical 结果表
    那么 该模型应在全部 18 个数据集上都有行
    并且 每个 dataset-method 行的 n_success 应等于计划 seed 数
    并且 status_summary 中不应存在 failed、timeout、incomplete 或 unknown
    并且 全覆盖完成后才允许标记 full_coverage: true

  @P1 @coverage-matrix
  场景: 自动生成覆盖矩阵用于发现未跑满模型
    假如 results/formal/benchmark_summary_all_current.csv 存在
    当 维护者运行覆盖检查脚本
    那么 输出应包含每个方法的 dataset_count、success_run_count、planned_run_count、coverage_scope
    并且 应列出缺失数据集和缺失 seed
    并且 attentionae_sc、scdsc、scname、sczidesk 这类 partial coverage 模型应被明确标记为 incomplete_full_coverage
    并且 desc、sccdcg、scgnn、scdeepcluster 这类无 canonical 行模型应被明确标记为 not_in_current_canonical

  @P1 @failure-triage
  场景: 失败运行必须保留可诊断信息
    假如 一个迁移模型在任意数据集或 seed 上失败
    当 runner 写入 status.json 和 run.log
    那么 error 应包含失败类型和关键错误摘要
    并且 run.log 应保留完整命令、stdout、stderr
    并且 汇总脚本不应删除失败目录，除非另有 cleanup report 记录
    并且 下一次补跑应能根据 status.json 跳过已成功 seed 并只重跑失败 seed

  @P1 @scheduler
  场景: 批量补跑使用统一可恢复调度
    假如 需要补跑多个模型、数据集和 seed
    当 维护者启动批量任务
    那么 任务计划应来自 manifest 和覆盖矩阵
    并且 不应手写只覆盖部分数据集的临时脚本作为最终方案
    并且 每个任务应经过 dry-run 生成命令清单
    并且 调度器应支持 resume、skip-success 和 fail-continue

  @P1 @slow-model
  场景: 慢速模型使用单独性能门槛
    假如 模型被标记为 SLOW_INCONCLUSIVE
    当 维护者决定是否继续全量运行
    那么 应先记录小数据 smoke 的 elapsed_seconds、内存和 GPU/CPU 占用
    并且 scGNN 这类 194 cells 仍可能 600s+ 的模型应进入 slow_queue
    并且 slow_queue 模型不得阻塞其它已验证模型的 full benchmark

  @P2 @docs
  场景: 文档明确解释为何已迁移模型没有全数据集结果
    假如 用户查看 README、docs/migration_status.md 或结果目录 README
    当 某模型已迁移但没有 18 数据集结果
    那么 文档应说明具体原因类别
    并且 原因应从 env_gated、pending_smoke、slow_inconclusive、subset_only、failed_needs_triage、deliberately_excluded 中选择
    并且 文档应给出下一步动作和负责人可执行命令

  @P2 @canonical-results
  场景: 当前结果入口与迁移覆盖报告分离
    假如 用户需要比较模型效果
    当 使用结果表
    那么 应使用 benchmark_summary_all_current.csv 作为当前成功运行结果入口
    并且 应使用单独 coverage report 判断哪些迁移模型未跑满
    并且 历史表、临时表和只含子集的表不得混作 full benchmark 结论

  @acceptance
  场景: 本轮迁移 BDD 的最小验收
    假如 本 BDD 被执行
    当 维护者完成第一轮实现
    那么 应生成 scCluBench 源模型到 methods 的迁移矩阵
    并且 应生成 current coverage matrix
    并且 应列出每个未全覆盖模型的阻塞原因和下一步补跑计划
    并且 不应把 SOURCE_MIGRATED 模型自动纳入 FULL_COVERAGE 或默认正式运行
