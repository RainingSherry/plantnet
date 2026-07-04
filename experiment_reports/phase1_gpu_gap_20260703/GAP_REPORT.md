# Phase 1 — CPU↔GPU 缺口实测报告（2026-07-03，进行中）

目的：量化深度单细胞方法（旗舰 = scMAE+DEC+std-floor / AdaptiveSwitch）流水线里
「时间花在哪、哪些是 CPU-bound、GPU 化能省多少、随细胞数如何变化」，作为 GPU-迁移
benchmark（总纲领强度A）的事实基础。

计时纪律：顺序跑、线程绑 48 核，避免并行争抢污染 wall-time。每档 30 epochs。

## 分阶段耗时（秒）

| run | device | n_cells | load | neighbor(PCA+KNN) | kmeans_init | train/epoch | eval_km | GPU利用率 |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Melanoma | cuda | 4513 | 3.1 | 0.7 | 0.5 | 0.271 | 0.4 | 5.9% |
| Melanoma | cpu  | 4513 | 3.1 | 0.7 | 0.4 | 0.287 | 0.3 | — |
| Quake | cuda | 9552 | 3.6 | 2.7 | 0.5 | 0.589 | 0.4 | 6.0% |
| Quake | cpu  | 9552 | 2.7 | 2.7 | 0.5 | 0.640 | 0.3 | — |
| Macosko | cuda | 44808 | 7.1 | 46.2 | 4.7 | 2.832 | 3.9 | 5.4% |

## 三个核心发现（都指向"迁移是非平凡问题"）

### 发现1：经典步骤是 CPU-bound 且**超线性爆炸→整条 CPU 管线大规模不可行**
本 harness 的 neighbor 阶段(PCA + **cosine** KNN, k=16) 随细胞数(HVG-1000, CPU sklearn,
48线程)：**2.9s(10k) → 56.3s(50k) → 161.4s(86k)**；整条 CPU 管线(load+PCA+KNN+2×KMeans
n_init=20+训练)在 **200k 细胞 >900s 超时**。经典步骤占比 81%→92%→95%。
**GPU 迁移把这几步拉回秒级**(独立单步实测, `../../文献调研_CS一区选题_20260703/cuML加速实测.md`)：
PCA 7–15×、KMeans 20–33×、KNN 7–28×，且**规模越大越划算**(PCA/KMeans 加速比随细胞数上升)。
例：KMeans@200k CPU 64.3s → cuML 2.1s(30.7×)；PCA@200k 15.4×。
kmeans_init/eval_km 同为 CPU sklearn, 同样随规模增长。

### 发现2（反直觉，关键）：小模型的 GPU 训练**并不比 CPU 快**
train/epoch 的 CPU vs GPU 几乎相同：Melanoma 0.287 vs 0.271，Quake 0.640 vs 0.589。
原因：模型是 128 维小 MLP，GPU 计算量微不足道；瓶颈在 **CPU 数据管线**(random_mask /
DataLoader / numpy)。GPU 利用率仅 **5-6%** 直接印证 GPU 在挨饿。
=> **朴素的 `.to(cuda)` 对这类深度单细胞方法几乎零收益。** 这颠覆了"把模型丢上 GPU 就快"
的想当然，是一个诚实且更有价值的发现。

### 发现3：真正的加速需要"重构整条管线"，不是"移动模型"
时间三大块：(a) 经典步骤(PCA/KNN/KMeans, CPU, 超线性) (b) 数据管线(masking/loader, CPU,
饿死GPU) (c) 模型前反向(GPU, 但量小)。真加速 = cuML 化经典步骤 + GPU 常驻数据/融合 kernel
喂饱 GPU。**这正是一个"迁移 agent"要自己诊断+决策的非平凡问题**（移哪、怎么移、正确性守住），
而不是机械 `.to(cuda)`。→ 强力支撑 benchmark+agent 选题的意义。

## 对选题/论文的含义

- benchmark **必须分阶段计时**：因为"模型"不是时间所在；naive 全流程 speedup 会掩盖真相。
- 正确性判据(总纲领M2)仍是脊梁：cuML 版 KNN/KMeans 与 sklearn 结果不 bit-match(近似/随机)，
  必须用统计等价判据验证"迁移后聚类指标不退化"。
- "naive GPU 端口零收益 → 必须智能迁移" 本身就是论文的一个卖点级 motivating finding。

## 已补齐（两个后台 agent）

1. ✅ **cuML vs sklearn 加速比**：`rapids_bench` 环境(RAPIDS 25.04/cuML 25.04, mamba一次建成)
   实测 —— PCA 7–15×、KMeans 20–33×、KNN 7–28×，规模越大越划算。见 `cuML加速实测.md`。
2. ✅ **规模阶梯数据**：`scale_bench/`(真实 10k/50k/86k/200k Arabidopsis 拼接 + 合成 500k/1M)。
   CPU 超线性曲线已测到 200k(整条管线超时)；500k/1M CPU 已停(浪费, 结论在200k已成立)。

## 计时口径说明（诚实校对）

- 本 harness 的 neighbor 用 **cosine KNN + PCA**，且整条管线含 **2×KMeans(n_init=20)** + 训练；
  cuML 报告测的是 **孤立单步**(euclidean KNN, KMeans n_init=10)。故"200k 超时"指**整条 CPU
  管线**不可行，不是 KNN 单步——两组数字互补，别混用。
- 精确 KNN 两侧都是 O(n²)，大规模加速比收敛到 ~7×；换近似近邻(cuML cagra/NN-descent)可进一步
  拉开——这是 Phase 3/benchmark 里值得单列的一条。

## 一句话结论（Phase 1）

深度单细胞方法的 wall-time 不在"模型"上：**朴素 `.to(cuda)` ≈ 零收益**(小模型 GPU 训练不比
CPU 快、利用率 5%)；真正的时间在 **CPU 经典步骤(PCA/KNN/KMeans)**，它们超线性、大规模不可行，
而 cuML 迁移能拿到 7–33× 且规模越大越值。**"缝"已量化坐实**——迁移是非平凡的分阶段诊断问题，
正是 benchmark + agent 选题的立足点。下一步 M2：为"迁移后聚类指标不退化"定义统计等价正确性协议。

产物：`measure.py`、`run_ladder.sh`、`run_scale.sh`、`summary.csv`、`bench_cuml.py`、本报告。
