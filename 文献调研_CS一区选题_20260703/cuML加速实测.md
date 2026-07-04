# cuML (GPU) vs scikit-learn (CPU) 经典单细胞步骤加速实测

日期：2026-07-03
硬件：8× NVIDIA H100 80GB HBM3，驱动 580.159.03 / CUDA 13.0；实测固定单卡 `CUDA_VISIBLE_DEVICES=1`。
CPU 侧线程绑定 48（`OMP/OPENBLAS/MKL/NUMEXPR_NUM_THREADS=48`）。

## 1. 结论速览

在 HVG-1000 的真实单细胞数据上，cuML 单卡 H100 对三个经典步骤（PCA / KMeans / KNN 图）相对 48 线程 CPU 的 sklearn，均取得**明显加速**：

- PCA：约 7–15×
- KMeans（n_init=10）：约 20–33×（加速最大的一步）
- KNN 图（k=15）：约 7–28×

所有 12 个 (规模 × 步骤) 单元的 CPU 侧都在 600s 内跑完，无超时、无不可行项。

## 2. 主结果表：规模 × 步骤 × (CPU秒 / GPU秒 / 加速比)

矩阵：X 转 dense float32；PCA `n_components=50`，KMeans `n_clusters=15, n_init=10`，KNN `n_neighbors=15, metric=euclidean`（fit + 对全体点 kneighbors 建图）。cuML 均先预热 1 次再计时；数据→显存的搬运不计入 GPU 计时。

| 规模 (cells×genes) | 步骤 | CPU sklearn (s) | GPU cuML (s) | 加速比 |
|---|---|---:|---:|---:|
| 10k (10000×1000) | PCA | 0.154 | 0.018 | 8.4× |
| 10k | KMeans | 4.129 | 0.200 | 20.7× |
| 10k | KNN | 0.527 | 0.019 | 27.6× |
| 50k (50000×1000) | PCA | 0.226 | 0.030 | 7.5× |
| 50k | KMeans | 19.301 | 0.758 | 25.4× |
| 50k | KNN | 3.586 | 0.342 | 10.5× |
| 86k (86215×1000) | PCA | 0.295 | 0.025 | 11.7× |
| 86k | KMeans | 27.440 | 0.844 | 32.5× |
| 86k | KNN | 8.124 | 1.027 | 7.9× |
| 200k (200000×1000, 跨物种) | PCA | 0.559 | 0.036 | 15.4× |
| 200k | KMeans | 64.298 | 2.096 | 30.7× |
| 200k | KNN | 38.976 | 5.451 | 7.2× |

原始数据同时写入 `experiment_reports/phase1_gpu_gap_20260703/bench_cuml_results.csv`。

## 3. 加速比随细胞数如何变化

- **KMeans 是最稳、最大的加速点（约 20–33×）**：CPU 时间随细胞数近似线性上涨（10k→200k：4.1s→64.3s），GPU 侧涨得很慢（0.2s→2.1s），所以规模越大 GPU 越划算。这是 EM 式迭代 + n_init=10 重复的典型 GPU 友好负载。
- **PCA 的加速比随规模走高（8.4×→15.4×）**：本身两侧都很快（GPU 全程 <0.04s），绝对收益小，但大规模下 CPU 的 SVD 相对更吃力，比值上升。
- **KNN 图的加速比随规模反而回落（27.6×→7.2×）**：小数据时 GPU 近乎"免费"（0.02s）导致比值虚高；规模上来后 GPU 的暴力近邻是 O(n²) 量级，5.5s@200k，CPU（brute + 48 线程）也在同数量级，故比值收敛到 ~7×。若换成近似近邻（cuML 的 NN descent / cagra）大规模下应能进一步拉开，但本次为公平对齐用的是精确 brute。
- 总体趋势：**细胞数越大，GPU 迁移的净收益越明显**（PCA、KMeans 都随规模上扬），唯独精确 KNN 因两侧同为 O(n²) 而比值趋稳。这条"规模越大越值"的曲线正是 GPU 迁移 benchmark 想要的核心叙事。

## 4. 环境搭建（记录 + 踩坑）

**采用方案 A（mamba/conda 官方钉定组合），一次成功，未回退到方案 B。**

建环境命令：
```bash
mamba create -y -n rapids_bench -c rapidsai -c conda-forge -c nvidia \
    rapids=25.04 python=3.10 'cuda-version=12.8'
/data/luolie/conda/envs/rapids_bench/bin/pip install "rapids-singlecell==0.13.1"
```

冒烟验证（GPU 上 import + 跑一次 PCA）通过：
```bash
CUDA_VISIBLE_DEVICES=1 /data/luolie/conda/envs/rapids_bench/bin/python \
  -c "import cuml,cupy,cugraph,rapids_singlecell; print('OK', cuml.__version__)"
```

**最终各库版本（rapids_bench）：**

| 库 | 版本 |
|---|---|
| python | 3.10.20 |
| cuml | 25.04.00 |
| cudf | 25.04.00 |
| cugraph | 25.04.01 |
| cupy | 14.1.1 |
| rapids-singlecell | 0.13.1 |
| scikit-learn | 1.7.2 |
| numpy | 2.0.2 |
| scipy | 1.15.2 |
| anndata | 0.11.4 |
| scanpy | 1.11.5 |
| CUDA runtime（cupy 报告） | 12.9 |

**踩坑 / 注意事项：**

1. **只新建独立环境 `rapids_bench`，全程绝对路径，未 `conda activate`，未改动任何现有环境。** 环境体积约 12G，装在 `/data/luolie/conda/envs/rapids_bench`。
2. **mamba solve 较慢**：`rapids=25.04` 的依赖解析 + 下载花了约 8–10 分钟（几百 MB~GB），期间 env 目录一直为空属正常，耐心等即可，最终 `Transaction finished` / EXIT_CODE=0。
3. **cu12 wheel 在 CUDA 13.0 驱动上前向兼容可用**：环境是 CUDA 12.9 运行时，跑在 580.159.03（CUDA 13.0）驱动上无问题，cuML/cupy 均正常初始化并出结果。
4. **驱动 CUDA 版本 ≠ 运行时 CUDA 版本**：无系统 nvcc 不影响，RAPIDS wheel 自带运行时。
5. **numpy 2.0.2**：pip 装 rapids-singlecell 时带入 numpy 2.x，与本环境的 cuml 25.04 / sklearn 1.7.2 均兼容，冒烟与全量 benchmark 都无 `_ARRAY_API not found` 之类的 numpy2/torch 冲突（本环境无 torch，故不涉及）。
6. **数据读取**：HVG-1000 的 h5ad 里 `uns` 只有 `hvg`，本批数据未触发 `encoding-type: null` 读崩溃；脚本仍防御性注册了 null-shim（`_register_null_h5ad_reader`），对含 null uns 的数据也安全。
7. **文件名**：200k 文件实际名为 `real_200k_hvg1000_xspecies.h5ad`（任务描述里 `real_200k_xspecies` 的命名顺序与实际略有出入），脚本已按实际路径配置。

## 5. 计时方法学（可复现）

- 每个 (数据 × 步骤 × 后端) 单元在**独立子进程**中运行：GPU 显存每单元干净、互不干扰；CPU 侧有硬性 600s 墙钟超时（超时记 `CPU>600s(不可行)`，本次无触发）。
- CPU 线程绑到 48；sklearn 子进程设 `CUDA_VISIBLE_DEVICES=""` 确保不误用 GPU。
- cuML 均先预热 1 次（吸收 CUDA context / JIT / allocator 初始化），再计时第二次；`cp.cuda.Stream.null.synchronize()` 保证计时含完整 kernel 执行。
- 数据 host→device 搬运（`cp.asarray`）在计时窗口之外，反映稳态计算加速比。

计时脚本：`experiment_reports/phase1_gpu_gap_20260703/bench_cuml.py`
（`--worker` 单元模式 / `--run-all` 全量编排；结果 CSV 同目录）。
