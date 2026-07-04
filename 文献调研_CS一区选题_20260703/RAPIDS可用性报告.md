# 本机 RAPIDS / GPU 经典单细胞管线可用性报告

调研日期：2026-07-03
调研目标：确认本机能否为「经典单细胞管线（PCA / KNN 图 / KMeans / Leiden / UMAP）」提供 GPU 加速（RAPIDS / cuML / cuPy / cuGraph / rapids-singlecell），为 CPU vs GPU 耗时 benchmark 项目做可行性判断。
调研原则：只调研 + 最小验证，不激活环境、不做大规模安装或改动。

---

## 一、结论速览（TL;DR）

1. **本机没有任何 conda 环境安装了 RAPIDS 系列库**（cuml / cupy / cudf / cugraph / rapids_singlecell 全部缺失，13 个环境无一例外）。
2. **硬件与驱动完全支持 RAPIDS**：8× H100 80GB（算力 9.0），驱动 580.159.03 / CUDA 13.0。已实测现有 torch 的 cu121 wheel 在该驱动上 `torch.cuda.is_available() == True`，证明 **CUDA 12 的 wheel 可在 CUDA 13 驱动上前向兼容运行** → RAPIDS 的 `*-cu12` wheel 可跑。
3. **没有现成可直接用的环境**，但可行路径明确：**新建一个专用环境**，用 conda/mamba 或 pip 安装 `rapids-singlecell + RAPIDS 25.04 (cu12)`。不建议改动主 benchmark 环境 `scssl_bench_py310`。
4. 安装体量偏大（完整 RAPIDS 栈约 **4–6 GB** 磁盘；仅核心大件 `libcuml` 单个 wheel 就 424.6 MB，`cupy` 132.7 MB）。可控但非轻量，需预留下载时间与磁盘。

---

## 二、CUDA / 驱动信息

| 项目 | 值 |
|---|---|
| GPU | 8× NVIDIA H100 80GB HBM3 |
| 算力 (compute capability) | 9.0（实测 `torch.cuda.get_device_capability()` = (9, 0)） |
| 驱动版本 | 580.159.03 |
| 驱动支持的最高 CUDA | 13.0（`nvidia-smi` 报告） |
| 系统 CUDA toolkit（nvcc） | **无**（`nvcc not found in PATH`） |
| torch 运行时 CUDA（cu121 wheel 实测） | 12.1，`cuda.is_available()=True`，可用 |

要点：
- **系统无 nvcc 不影响 RAPIDS**。RAPIDS 的 pip wheel（`*-cu12`）自带 CUDA 运行时库（`nvidia-*-cu12` 依赖），不需要系统装 CUDA toolkit。
- **CUDA 版本匹配**：驱动是 CUDA 13，但 RAPIDS 目前主线发行是 `cu12`（CUDA 12.x）wheel。因 NVIDIA 驱动向前兼容，cu12 wheel 在 CUDA 13 驱动上可正常运行（已由现有 torch cu121 实测佐证）。**选 `cu12` 版本即可**。
- **GPU 使用限制**：按任务要求，0 号和 7 号禁用。当前 1/2/3 号最空闲（各约 5.5GB 被 llama-server 占用，利用率 0%），4/5/6 号在跑任务（利用率 54–85%）。跑 benchmark 时用 `CUDA_VISIBLE_DEVICES=1`（或 2/3）指定单卡。

---

## 三、各 conda 环境 RAPIDS 库矩阵

探测方式：直接用各环境 python 绝对路径 `import`，未激活任何环境。
根目录：`/data/luolie/conda/envs/`

图例：`OK`=已安装，`—`=缺失。

| 环境 | Python | cuml | cupy | cudf | cugraph | rapids_singlecell | sklearn | scanpy | torch | numpy |
|---|---|---|---|---|---|---|---|---|---|---|
| cellagent | 3.10.13 | — | — | — | — | — | 1.7.2 | 1.11.5 | 2.1.0+cu121 | 2.2.6 ⚠ |
| MixDiffusion | 3.8.20 | — | — | — | — | — | 1.3.2 | — | 2.4.1+cu121 | 1.24.4 |
| PhytoCluster | 3.10.20 | — | — | — | — | — | 1.7.2 | 1.10.4 | 2.5.1+cu121 | 2.2.6 |
| plantnet | 3.8.20 | — | — | — | — | — | — | — | — | — |
| scclubench-foundation | 3.10.20 | — | — | — | — | — | 1.7.2 | 1.11.5 | 2.1.0+cu121 | 2.2.6 ⚠ |
| scclubench-main | 3.9.25 | — | — | — | — | — | 1.6.1 | 1.10.0 | 2.1.2+cu118 | 1.24.3 |
| scclubench-sccdcg | 3.8.20 | — | — | — | — | — | — | — | 1.12.0+cu113 | 1.24.4 |
| scdeepcluster_legacy | — | 环境损坏：无 python 可执行文件 | | | | | | | | |
| scilama | 3.11.15 | — | — | — | — | — | 1.6.1 | 1.11.3 | 2.5.1+cu121 | 2.4.4 |
| scPlantLLM_Py_Env | 3.10.20 | — | — | — | — | — | 1.4.2 | 1.10.1 | 2.3.0+cu121 | 1.26.4 |
| scPlantLLM_R_Env | (R 环境) | 未探测（R 环境，与本任务无关） | | | | | | | | |
| scssl_bench_py310 | 3.10.20 | — | — | — | — | — | 1.7.2 | 1.9.6 | 2.4.0+cu121 | 1.26.4 |
| test_new_env | 3.10.20 | — | — | — | — | — | — | — | — | — |

⚠ 注：`cellagent` 与 `scclubench-foundation` 存在 numpy 2.x 与 torch(cu121) 的 ABI 不匹配警告（`Failed to initialize NumPy: _ARRAY_API not found`），这两个环境本身已有隐患，不适合再叠加 RAPIDS。

结论：**RAPIDS 五件套在所有环境均缺失**。当前所有「经典步骤」确实都是 CPU 的 sklearn/scanpy（scanpy 内部的 Leiden/UMAP/KNN 也是 CPU）。

---

## 四、安装可行性分析（走 pip / conda 新建路径）

### 4.1 版本匹配（关键）

`rapids-singlecell 0.13.1`（PyPI 最新）的依赖钉死如下（从 wheel METADATA 读取）：

- `Requires-Python: >=3.10, <3.14`
- extra `[rapids12]` → `cuml-cu12==25.04.*`、`cugraph-cu12==25.04.*`、`cudf-cu12==25.04.*`、`cupy-cuda12x`

即 **rapids-singlecell 0.13.1 ↔ RAPIDS 25.04 (cu12)** 是官方钉定组合。PyPI 上 `cuml-cu12` 已发布到 26.2.0，但要配 rapids-singlecell 就用 **25.04**。

已实测确认 cp310（Python 3.10）wheel 全部存在且可下载：
`cuml_cu12-25.4.0`、`cugraph_cu12-25.4.1`、`cudf_cu12-25.4.0`、`cupy_cuda12x-14.1.1`。

### 4.2 体量（实测 wheel 大小）

| 包 | 顶层 wheel 大小 | 说明 |
|---|---|---|
| cupy-cuda12x | 132.7 MB | GPU 数组基础库 |
| libcuml-cu12 | **424.6 MB** | cuML 的 C++ 大件（最大单包） |
| cuml-cu12 | 9.9 MB | 薄 wheel，实际算法在 libcuml |
| cugraph-cu12 | 3.2 MB | + 依赖 libcugraph（另计，数百 MB 级） |
| cudf-cu12 | 1.8 MB | + 依赖 libcudf（另计，数百 MB 级） |
| cuvs-cu12 | 0.73 MB | 近邻检索（KNN 用） |

顶层薄 wheel 之外还会拉入 `libcudf-cu12`、`libcugraph-cu12`、`librmm-cu12`、`pylibraft-cu12`、`libcuvs-cu12`、`dask/distributed`、`nvidia-*-cu12`（CUDA 运行时）等一批传递依赖。**完整栈落盘约 4–6 GB**，属可控但非轻量，务必预留磁盘与下载时间。

### 4.3 pip dry-run 的坑（已踩，供后续注意）

- 本机 pip 为 26.0.1，其 `pip install --dry-run --report` **并非纯解析**：为读取传递依赖的 metadata，它会**实际下载 wheel**（观测到它开始拉 424 MB 的 `libcuml`）。因此不能靠 dry-run「零下载」看全量依赖，已主动中止以省流量。
- NVIDIA 源（`https://pypi.nvidia.com`）解析很慢，pip 回溯（backtracking）会让不加钉死的解析长时间挂起（多次 7 分钟超时未完成）。**安装时务必钉死精确版本**以避免回溯。
- 本机默认 index 是清华镜像（`pypi.tuna.tsinghua.edu.cn`），RAPIDS 包在镜像里没有，需显式加 `--extra-index-url https://pypi.nvidia.com`。

---

## 五、cuML vs sklearn 计时对比

**未执行**。原因：本机当前没有任何环境装有 cuML（属于报告要求里的「若没有任何环境有 RAPIDS」分支），且任务明确要求「不要真装除非无冲突且体积可控」。因此本节留待环境搭好后补测。

建议的最小验证脚本（装好环境后运行，随机造矩阵，不读大数据）：

```python
# CUDA_VISIBLE_DEVICES=1 python bench_kmeans.py
import time, numpy as np
X = np.random.RandomState(0).rand(50000, 50).astype(np.float32)

from sklearn.cluster import KMeans as skKMeans
t = time.time()
skKMeans(n_clusters=12, n_init=10, random_state=0).fit(X)
print("sklearn KMeans:", round(time.time()-t, 2), "s")

import cupy as cp
from cuml.cluster import KMeans as cuKMeans
Xg = cp.asarray(X)
cuKMeans(n_clusters=12, n_init=10, random_state=0).fit(Xg); cp.cuda.Stream.null.synchronize()  # 预热
t = time.time()
cuKMeans(n_clusters=12, n_init=10, random_state=0).fit(Xg); cp.cuda.Stream.null.synchronize()
print("cuml KMeans:", round(time.time()-t, 2), "s")
# 同理可加 cuml.PCA vs sklearn.PCA、cuml.neighbors.NearestNeighbors vs sklearn
```

---

## 六、如何得到一个能跑 GPU 经典管线的环境（明确建议）

**推荐方案 A（首选）：新建专用 conda 环境，用 conda/mamba 装 RAPIDS。**
本机有 `conda 26.1.1` 和 `mamba`。conda/mamba 解析比 pip 在 NVIDIA 源上更干净、更少回溯，且能一并管理 CUDA 运行时。

```bash
# 不激活、不动现有环境；GPU 用 1/2/3 号
mamba create -n rapids_bench -c rapidsai -c conda-forge -c nvidia \
    rapids=25.04 python=3.10 'cuda-version=12.8'
/data/luolie/conda/envs/rapids_bench/bin/pip install rapids-singlecell
```

**推荐方案 B（备选）：纯 pip 装进一个全新 venv/conda 环境（Python 3.10 或 3.11）。**
务必钉死版本、显式加 NVIDIA 源，避免回溯：

```bash
mamba create -n rapids_bench python=3.10 -c conda-forge
/data/luolie/conda/envs/rapids_bench/bin/pip install \
    --extra-index-url https://pypi.nvidia.com \
    "rapids-singlecell[rapids12]==0.13.1" \
    "cuml-cu12==25.4.0" "cugraph-cu12==25.4.1" "cudf-cu12==25.4.0" "cupy-cuda12x"
```

**不推荐**：直接往主 benchmark 环境 `scssl_bench_py310` 里装。理由：(1) 会顺带把 scanpy 从 1.9.6 升到 ≥1.10（rapids-singlecell 要求），可能扰动现有 CPU 基线的可复现性；(2) 4–6 GB 依赖 + numpy/cupy 版本约束有污染现有环境的风险。benchmark 要「CPU vs GPU 同题对比」，更应保持 CPU 基线环境不变、GPU 另起干净环境。

**环境选择要点**：
- Python 必须 3.10–3.13（rapids-singlecell 0.13.1 的 `Requires-Python`）。现有可复用的 3.10/3.11 候选：`PhytoCluster`(3.10, scanpy1.10.4)、`scPlantLLM_Py_Env`(3.10)、`scilama`(3.11)——但同样建议新建而非改动它们。
- 选 `cu12` 版本（不是 cu11），匹配驱动的前向兼容。
- 装完做 `import cuml, cupy, cugraph, rapids_singlecell` 冒烟测试 + 上面的 KMeans 计时。

---

## 七、中文简报结论

- **现状**：本机 13 个 conda 环境无一安装 RAPIDS（cuml/cupy/cudf/cugraph/rapids_singlecell 全缺）。现有经典单细胞步骤（PCA/KNN/KMeans/Leiden/UMAP）确实都跑在 CPU 的 sklearn/scanpy 上，深度方法才在 torch GPU 上。
- **可行性**：硬件（8×H100，算力9.0）、驱动（580.159.03 / CUDA 13.0）完全支持 RAPIDS；已实测 cu12 wheel 在该驱动上可用（前向兼容）。系统无 nvcc 不影响（RAPIDS wheel 自带 CUDA 运行时）。
- **不可直接用现成环境**，但**可以装**：官方钉定组合是 `rapids-singlecell 0.13.1 + RAPIDS 25.04 (cu12) + Python 3.10/3.11`，cp310 wheel 已确认可下载。
- **建议**：新建专用环境 `rapids_bench`（方案 A conda/mamba 首选，方案 B pip 备选），不要改动主 benchmark 环境 `scssl_bench_py310`。体量约 4–6 GB（libcuml 单件 424MB、cupy 132MB），需预留磁盘与下载时间；pip 安装务必钉死版本并加 `--extra-index-url https://pypi.nvidia.com` 以避开慢速回溯。
- **cuML vs sklearn 计时**：本次未测（无现成 cuML 环境，且按要求未实装）；环境搭好后用本报告第五节脚本补测即可。
- **跑测提醒**：GPU 仅用 1–6 号（禁 0/7），当前 1/2/3 最空闲，用 `CUDA_VISIBLE_DEVICES` 指定单卡。
