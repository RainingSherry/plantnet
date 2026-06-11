# language: zh-CN

# 基本信息
# 项目: plantnet dimension-reduction benchmark
# 路径: /home/luolie/biopipeline/dimension-reduction/plantnet
# 创建日期: 2026-06-11
# 作者: Codex
# 状态: Executed
# 目的: 在不立即修改现有模型代码的前提下，明确后续修复 benchmark 运行、迁移和结果汇总风险的 BDD 验收标准。
# 范围: scripts/run_formal_benchmark.py, methods/method_manifest.yaml, methods/Traditional/sc3/run.py,
#       methods/evaluation.py, methods/DeepLearning/scMAE_family.py,
#       methods/DeepLearning/NeighborMix_scMAE, results/formal 汇总脚本与文档。
# 非目标: 本文件不修改模型逻辑、不重跑模型、不删除已有结果。

功能: benchmark 迁移与运行可靠性修复
  为了避免后续模型迁移、批量运行和结果汇总时再次产生错误或混乱
  作为 benchmark 维护者
  我希望关键运行入口、模型清单、评估指标和结果表格都有一致、可验证的行为

  背景:
    假如 当前项目根目录为 "/home/luolie/biopipeline/dimension-reduction/plantnet"
    并且 当前最全且清理后的结果表为 "results/formal/benchmark_summary_all_current.csv"
    并且 用户已经明确不再需要 "scanpy_standard" 模型

  @P0 @metrics @sc3
  场景: SC3 不再覆盖统一评估接口生成的正确指标
    假如 SC3 运行结束并产生聚类标签和 embedding
    当 保存 benchmark 输出时
    那么 SC3 应只使用统一且经过验证的评估函数计算 ACC 和 F1-macro
    并且 Hungarian 标签映射方向应为 "predicted cluster -> true label"
    并且 ACC 应比较编码后的真实标签与对齐后的预测标签
    并且 F1-macro 应使用对齐后的预测标签
    并且 "metrics.json" 不应被一段独立的错误指标计算逻辑覆盖

  @P0 @method-selection @scanpy-standard
  场景: 默认正式运行不会再次包含 scanpy_standard
    假如 用户未在命令行显式传入 "--methods"
    当 正式 benchmark runner 解析默认方法列表时
    那么 "scanpy_standard" 不应出现在将要运行的方法集合中
    并且 "methods/method_manifest.yaml" 中 scanpy_standard 不应是 default_in_formal true
    并且 README 与 runner help 中不应再把 scanpy_standard 作为默认正式方法展示

  @P0 @method-selection @manifest
  场景: 默认方法选择逻辑与文档声明一致
    假如 README 或 runner help 声明默认正式方法集合
    当 runner 在未指定 "--methods" 的情况下选择方法
    那么 实际选择出的集合应与声明的集合完全一致
    并且 新增 VERIFIED 方法不会仅因 default_in_formal 字段误设而进入默认正式运行
    并且 DESC、scCDCG、AttentionAE_sc 等扩展方法的默认运行资格应被显式确认

  @P1 @results @canonical-table
  场景: 当前结果只有一个推荐入口表
    假如 用户需要比较所有已运行模型在所有数据集上的表现
    当 查看 results/formal 目录中的汇总文件时
    那么 "benchmark_summary_all_current.csv" 应被标记为当前最全、最干净的推荐入口
    并且 旧的 "benchmark_summary_final.csv"、"benchmark_summary_combined.csv" 或每数据集旧 summary 不应被文档推荐为主入口
    并且 如保留旧表，应明确标注为历史表或派生表

  @P1 @results @script
  场景: 过期的 reformat_benchmark 脚本不会读取已删除中间文件
    假如 "results/formal/benchmark_summary_merged.csv" 不存在
    当 执行 "results/formal/reformat_benchmark.py"
    那么 脚本不应因硬编码缺失文件而失败
    并且 它应读取当前 canonical 表，或明确退出并提示该脚本已废弃
    并且 输出文件名不应与旧的、不完整的 summary 语义冲突

  @P1 @data-conversion @cache
  场景: h5 到 h5ad 的自动转换缓存可追踪源数据
    假如 输入数据为 ".h5" 文件
    并且 data/processed 中已存在同名 h5ad 和 meta.json
    当 runner 判断是否复用缓存时
    那么 meta.json 应记录源文件绝对路径、文件大小、mtime 或内容 hash、转换参数和 label_key
    并且 当源文件或转换参数变化时应重新转换
    并且 不应仅凭 dataset_name 相同就复用旧缓存

  @P1 @output-verification
  场景: 输出完整性校验按方法声明 required artifacts
    假如 一个方法在 manifest 中声明了所需输出文件
    当 runner 校验一次运行是否成功时
    那么 应使用该方法的 required artifacts 列表进行校验
    并且 对于只需要参与汇总的运行，至少应要求 "metrics.json" 和成功状态存在
    并且 不同方法不应被强制套用同一组输出文件导致误判 incomplete

  @P2 @results @seed-parsing
  场景: collect_results 能可靠识别 run_id 中的 seed
    假如 运行目录名形如 "method__seed42__20260611_120000"
    当 collect_results 解析 seed
    那么 应识别出 seed 为 42
    并且 即使缺少 status.json，也不应把该运行错误归为 unknown seed
    并且 seed 解析应使用稳定的正则或结构化 run_id 规则

  @P2 @neighbormix @docs
  场景: NeighborMix-scMAE 文档只描述实际生效参数
    假如 README 或方法说明列出默认超参数
    当 用户根据文档迁移或重跑 NeighborMix-scMAE
    那么 文档应说明当前有效参数是 use_pseudo、pseudo_weight、alpha、neighbor_k、mix_neighbors 和 mask_ratio
    并且 mix_weight、consistency_weight、target_mode 应被标记为兼容保留且不参与当前训练逻辑
    并且 文档应说明 pseudo branch 的重建目标是原始真实细胞

  @P2 @neighbormix @maintenance
  场景: NeighborMix-scMAE 不保留会误导维护者的重复 helper
    假如 NeighborMix-scMAE 主流程实际调用 scMAE_family 中的数据加载、设备选择、评估和保存函数
    当 开发者维护 NeighborMix-scMAE
    那么 本地未被 main 使用的重复 helper 应被删除、重命名为 legacy，或添加明确注释
    并且 后续修复应优先落在实际被 main 调用的 scMAE_family 函数中

  @P2 @scmae @mask
  场景: scMAE 系列记录实际有效 mask rate
    假如 mask_ratio 设置为 0.4
    当 apply_scmae_noise 在稀疏 scRNA 数据上执行替换扰动
    那么 训练日志或 history 应记录实际有效 mask rate
    并且 文档应说明相同值替换不会计入有效 mask
    并且 若期望严格按采样 mask 训练，应明确改为使用 should_swap 作为 mask

  @P2 @runtime @environment
  场景: environment.json 收集不会显著拖慢或干扰每次运行
    假如 runner 在每个方法运行前写 environment.json
    当 某些环境没有 torch、tensorflow 或 scanpy
    那么 版本采集失败不应影响方法运行
    并且 不相关框架的 import 不应造成明显启动延迟
    并且 每个版本采集命令应有短超时和清晰的失败记录

  @acceptance
  场景: 修复完成后的最小验收
    假如 上述 P0 和 P1 场景均通过
    当 运行一次 dry-run 和一次小数据 smoke run
    那么 默认方法列表不包含 scanpy_standard
    并且 SC3 输出的 ACC/F1 与统一 evaluation 结果一致
    并且 汇总脚本生成的 canonical 表不含 scanpy_standard
    并且 canonical 表能清楚区分成功、失败、跳过和不完整运行

# 执行记录
# 2026-06-11:
# - P0/P1/P2 场景已落地到代码和文档。
# - dry-run 验收通过: 默认方法为 9 个，不包含 scanpy_standard。
# - SC3 smoke run 通过: /tmp/plantnet_bdd_smoke2/bdd_sc3_smoke2。
# - 当前 canonical flat 表已生成: results/formal/benchmark_summary_all_current_flat.csv。
