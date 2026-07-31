# TSI-Denoising

> 一维线性密集台阵高频面波模态分离与三台干涉去噪工具包  
> High Frequency surface-wave mode separation, and three-station interferometry denoising for dense 1-D arrays

[中文说明](#中文说明) · [English](#english)

---

## 中文说明

### 1. 项目简介

TSI-Denoising 是一个面向科研和教学的 Python 程序包，用于处理一维线性密集台阵中的背景噪声互相关函数（ambient-noise cross-correlation，ANC）。程序将 ObsPy `Stream` 封装为经过验证的 `Wavefield` 处理对象，并提供从数据读取、预处理、频散分析、模态分离到三台干涉（three-station interferometry，TSI）迭代去噪的一套工作流。

本仓库同时包含为**第十二届地震学算法和程序培训班**（2026 年 8 月 10–12 日）准备的两个教学案例：

- **MARS DAS**：海底 DAS 单模态 Scholte 波的 TSI 去噪，
- **RR Array**：跨断层密集台阵多模态瑞利波分离和分模态去噪。

除非具体接口另有说明，程序统一使用以下单位：

| 物理量 | 单位 |
|---|---|
| 距离 | km |
| 相速度 | km/s |
| 相关时间 | s |
| 频率 | Hz |

> [!IMPORTANT]
> 原始教学数据不纳入代码仓库。运行 `python retrieve_datasets.py` 可下载并安全解压两个案例的输入数据；`processed/` 中的缓存和结果由教程生成或复用。

### 2. 仓库结构

```text
TSI-Denoising/
├── src/tsi_denoising/              # Python 程序包
│   ├── io/                         # SAC 目录读取
│   ├── mode_separation/            # 极化与相位匹配分离
│   ├── denoising/                  # 三台干涉与诊断绘图
│   ├── wavefield.py                # Wavefield 数据模型
│   ├── preprocessing.py            # 共同预处理
│   └── masw.py                     # MASW 计算与绘图
├── tutorial/
│   ├── MARS_DAS/                   # Monterey 湾海底 DAS 单模态 Scholte 波案例
│   └── RR_Array/                   # 跨 San Jacinto 断裂带密集台阵多模态瑞利波案例
├── tests/                          # 单元测试
├── pyproject.toml                  # 包元信息和依赖
├── environment.yml                 # 推荐的 Conda 环境
├── requirements.txt                # 运行依赖
└── README.md
```

### 3. 安装与数据准备

#### 6.1 获取程序

本项目通过 GitHub 分发。安装前请准备 [Git](https://git-scm.com/) 和 [Miniforge](https://github.com/conda-forge/miniforge)（或 Miniconda）。

~~~bash
git clone https://github.com/YuanYusung/TSI-Denoising.git
cd TSI-Denoising
~~~

#### 6.2 创建 Conda 环境

项目根目录的 `environment.yml` 定义了推荐环境。首次安装时执行：

~~~bash
conda env create -f environment.yml
conda activate tsi-denoising
python -m pip install -e .
~~~

环境使用 Python 3.10，并包含 NumPy、SciPy、Matplotlib、ObsPy、JupyterLab 和 ipykernel。`pip install -e .` 会以可编辑模式安装程序包；修改仓库中的 Python 源码后，无须重新安装。

#### 6.3 使用 Jupyter 与验证安装

若要在 Notebook 中明确选择该环境，请注册内核：

~~~bash
python -m ipykernel install --user --name tsi-denoising --display-name "Python (TSI-Denoising)"
~~~

打开教程 Notebook 后选择 `Python (TSI-Denoising)`。使用下列命令验证安装：

~~~bash
python -c "import tsi_denoising; print('TSI-Denoising imported successfully')"
~~~

#### 6.4 获取教学数据

原始 SAC 教学数据不随代码仓库发布。在项目根目录运行：

~~~bash
python retrieve_datasets.py
~~~

脚本从公开数据源下载归档，并只解压缺失的 `tutorial/RR_Array/input/` 和 `tutorial/MARS_DAS/input/`；已有目录不会被覆盖。两个教程会在各自的 `processed/` 目录保存可再生成的 NPZ 缓存和结果。

### 4. 核心能力

- 递归读取一个目录中的 SAC 互相关波形，并构造经过验证的 `Wavefield`；
- 统一台站名来源、台站对方向、采样信息和距离元数据；
- 对互相关函数进行正负时间对称化、速度窗 taper 和零相位带通滤波；
- 使用相移法计算频率–相速度 MASW 能量图；
- 根据 ZZ、ZR、RZ、RR 四分量极化关系分离逆进和顺进瑞利波；
- 使用参考频散曲线执行相位匹配模态分离；
- 对单个台站对生成 TSI 几何、输入波场和干涉结果诊断图；
- 对整个波场进行串行或多进程迭代 TSI 去噪；
- 使用 NPZ 文件保存和恢复波场、MASW 与去噪结果。

### 5. 处理流程

```text
SAC 台站对互相关函数
          |
          v
读取、规范化和验证
          |
          v
共同预处理与 MASW 频散诊断
          |
          +----------------------+
          |                      |
          v                      v
目标频带内单一模态主导数据          四分量、多模态数据
例如 MARS DAS                    例如 RR Array
      |                             |
      |                          极化分离
      |                             |
      |                        相位匹配模态分离
      |                             |
      +-------------+---------------+
                    |
           分模态三台干涉去噪与 QC

```


### 6. 适用范围与限制

本程序主要面向台站编号能够代表空间顺序的一维或近似一维阵列。TSI 假设目标面波主要沿阵列方向传播，并依赖不同台站组合之间稳定的几何和相位关系。

使用时需要注意：

- 二维台阵、弯曲测线或强横向不均匀介质可能不满足距离顺序和一维传播假设；
- 极化分离要求 ZZ、ZR、RZ、RR 四个分量具有相同台站对、距离、采样率和时间轴；
- ZR/RZ 符号、径向正方向和正负相关时间定义必须在数据制作阶段保持一致；
- MASW 能量峰值不等同于自动完成模态阶次判定，仍需结合理论频散、极化和空间连续性；
- TSI 不会自动消除模态交叉项，多模态数据应先进行可靠的模态分离；

### 7. 输入数据要求

#### 7.1 SAC 波形

每个 SAC 文件必须只包含一个台站对的一条一维互相关波形。`read_sac_directory()` 会递归搜索目录，默认匹配大小写不同的 `.sac`、`.SAC` 和 `.SAC_s` 等名称。

| 信息 | 首选字段 | 回退字段 | 要求 |
|---|---|---|---|
| 虚拟源名 | `stats.sac.kevnm` | `stats.source` | 非空且含数字编号 |
| 接收台站名 | `stats.sac.kstnm` | `stats.station` | 非空且含数字编号 |
| 台间距 | `stats.sac.dist` | 无 | 有限正数，单位 km |
| 相关时间起点 | `stats.sac.b` | 无 | 有限数，单位 s |

同一 `Wavefield` 内的 Trace 必须具有相同的样本数、采样间隔与 `sac.b`；数据必须是一维有限实数数组；有向台站对不能重复；不同台站名不能映射到同一个数字编号。

#### 7.2 台站编号与方向规范化

程序使用台站名最后一组连续数字确定空间顺序：

| 名称 | 解析编号 |
|---|---:|
| `RR010` | 10 |
| `LINE2_ST040` | 40 |
| `DAS00345` | 345 |

构造 `Wavefield` 或调用 `read_sac_directory()` 时，程序会删除自相关、统一台站对为较小编号指向较大编号、去除重复的相反方向记录，并在必要时交换头段与反转波形。默认还会检查距离是否随对端编号严格增加。

对于弯曲测线、二维阵列或编号不代表空间位置的数据，请先用 ObsPy 读取为 `Stream`，再以 `check_distance_order=False` 构造 `Wavefield`。这只跳过距离顺序检查，不会放宽其他数据一致性要求。

#### 7.3 相关时间轴

互相关波形应覆盖目标面波所需的正、负相关时间。`preprocess_stream()` 与 `Wavefield.preprocess()` 会按正负时间分支进行对称化；若原始数据只有正时间半支，应先在负时间端补零，并正确设置 `sac.b`。

### 8. 推荐工作流

根 README 仅说明通用流程和接口；完整的可运行代码、参数选择、质量控制与结果解释见教程。

#### 8.1 单分量、单一主导模态数据

适用于目标频带内由单一面波模态主导的数据，例如 MARS DAS 案例中的 Scholte 波。

~~~text
读取 SAC 互相关波形
        |
构造并验证 Wavefield
        |
对称化、速度窗和带通滤波
        |
时间–距离波场与 MASW 频散诊断
        |
确认目标频带内由单一模态主导
        |
单台站对三台干涉诊断
        |
全波场迭代去噪与质量控制
~~~

完整步骤见 [MARS DAS 教程](tutorial/MARS_DAS/README.md) 和 [MARS DAS Notebook](tutorial/MARS_DAS/run.ipynb)。

#### 8.2 四分量、多模态瑞利波数据

适用于具有 ZZ、ZR、RZ、RR 四分量互相关波形、且目标频带内存在多个瑞利波模态的数据。

~~~text
读取四分量 SAC 数据
        |
统一预处理与一致性检查
        |
波场与 MASW 频散诊断
        |
逆进–顺进极化分离
        |
相位匹配模态分离
        |
确认目标模态的空间与频散连续性
        |
各目标模态分别进行三台干涉去噪与质量控制
~~~

完整步骤见 [RR Array 教程](tutorial/RR_Array/README.md) 和 [RR Array Notebook](tutorial/RR_Array/run.ipynb)。

> [!IMPORTANT]
> 极化分离不等同于模态分离。逆进或顺进波场仍可能包含多个模态。将波场输入三台干涉前，应结合频散谱、参考频散曲线与空间连续性确认目标模态已可靠分离。

### 9. 关键接口与参数参考

推荐始终从 `tsi_denoising` 顶层导入公开接口。除非另有说明，距离单位为 km、速度单位为 km/s、时间单位为 s、频率单位为 Hz。

#### 9.1 波场构造、读取与预处理

| 接口 | 参数 | 说明 |
|---|---|---|
| `Wavefield(stream, *, component=None, copy=True, check_distance_order=True)` | `stream` | 必填的非空 ObsPy `Stream`；所有 Trace 必须满足第 7 节规范。 |
|  | `component` | 可选分量标签，如 `"ZZ"` 或 `"RR"`；只作为元数据。 |
|  | `copy` | 默认 `True`，复制输入 Trace；仅在调用方自行管理原始 Stream 时使用 `False`。 |
|  | `check_distance_order` | 默认 `True`，要求编号与距离顺序相容；弯曲或二维阵列可设为 `False`。 |
| `read_sac_directory(directory, pattern=None, *, component=None)` | `directory` | 必填 SAC 根目录；递归读取并返回验证、规范化后的 `Wavefield`。 |
|  | `pattern` | 一个通配符或通配符序列；默认自动匹配常见 SAC 扩展名。 |
|  | `component` | 可选分量名；省略时采用目录名。 |
| `preprocess_stream(wavefield, fmin=0.5, fmax=5.0, vmin=0.1, vmax=2.5, taper_fraction=0.05)` | `wavefield` | 必填输入波场；返回新的 `Wavefield`，不修改输入。 |
|  | `fmin`、`fmax` | 零相位带通滤波下限与上限；要求 `0 < fmin < fmax < Nyquist`。 |
|  | `vmin`、`vmax` | 距离相关速度窗范围；因果窗为 `distance / vmax` 至 `distance / vmin`。 |
|  | `taper_fraction` | 速度窗两端 taper 比例，必须位于 0 至 0.5。 |

`Wavefield.preprocess()` 使用同一组参数，但会原地替换对象内的数据并返回自身；需要保留原始波场时应使用 `preprocess_stream()`。

#### 9.2 MASW 频散分析

| 接口 | 参数 | 说明 |
|---|---|---|
| `MASW(wavefield, velocities=None, fmin=0.5, fmax=5.0, padding_factor=5, dist_threshold=0.2)` | `wavefield` | 必填、通常已预处理的输入波场；构造对象仅保存配置，须调用 `.compute()`。 |
|  | `velocities` | 相速度采样数组；省略时使用 0.2–2.5 km/s 的 231 点线性网格。 |
|  | `fmin`、`fmax` | MASW 频率范围；必须落在有效奈奎斯特频率内。 |
|  | `padding_factor` | FFT 零填充倍数，必须是不小于 1 的整数；增大仅加密频率采样。 |
|  | `dist_threshold` | 参与成像的最小台间距；小于该值的记录被排除。 |
| `MASW.compute(n_jobs=1)` | `n_jobs` | 频率行计算的进程数；从 1 开始验证，再按 CPU 与内存提高。 |
| `compute_masw(wavefield, velomin=0.2, velomax=2.5, fmin=0.5, fmax=5.0, padding_factor=5, *, velocities=None, dist_threshold=0.2, n_jobs=1)` | `velomin`、`velomax` | 未提供 `velocities` 时生成 231 点线性速度网格的范围。 |
|  | 其余参数 | 与 `MASW` 构造参数和 `.compute(n_jobs=...)` 含义相同；直接返回 `(velocity, frequency, amplitude)`。 |

#### 9.3 极化与相位匹配模态分离

| 接口 | 参数 | 说明 |
|---|---|---|
| `polarization_separate(*, zz, rz, zr, rr, swap_polarization=False, use_rr=True, vmin=0.1, vmax=2.5, taper_fraction=0.05)` | `zz`、`rz`、`zr`、`rr` | 必填且仅可按关键字传入的四个 `Wavefield`；它们必须有相同的规范化台站对、距离、采样率与时间轴。 |
|  | `swap_polarization` | 默认 `False`；若径向正方向或相关定义与教程相反，设为 `True` 以交换逆进与顺进标签。 |
|  | `use_rr` | 默认 `True`，在顺进组合中使用 RMS 归一化 RR 项；设为 `False` 时使用三分量组合。 |
|  | `vmin`、`vmax`、`taper_fraction` | 对结果应用的距离相关速度窗；含义与预处理接口相同。输入波场不会被修改。 |
| `phase_match_separate(wavefield, reference_curve=None, *, fmin=0.5, fmax=5.0, t_window=0.2, keep_positive=True, return_reference=False, masw_cache_path=None, vmin=0.1, vmax=2.5, taper_fraction=0.05)` | `wavefield` | 必填输入波场，通常为极化分离后的候选模态；函数不执行带通滤波。 |
|  | `reference_curve` | 可选 `(N, 2)` 数组，每行是频率与相速度；也接受 `(frequency, velocity)` 二元组。省略时加载可用 MASW 缓存或计算 MASW，并显示拾取界面。 |
|  | `fmin`、`fmax` | 应用参考相位的频带。 |
|  | `t_window` | 相位对齐后高斯时间窗的宽度，必须为正数。 |
|  | `keep_positive` | 默认 `True`，保留正时间支并衰减负时间支；`False` 时反向处理。 |
|  | `return_reference` | 默认 `False`；为 `True` 时返回 `(separated_wavefield, used_curve)`。 |
|  | `masw_cache_path` | 可选、已完成计算的 MASW NPZ；与 `reference_curve` 互斥。 |

#### 9.4 三台干涉去噪与结果绘图

| 接口 | 参数 | 说明 |
|---|---|---|
| `denoise_station_pair_demo(wavefield, station_pair, *, include_convolution=True, sqrt_spectrum=True, taper_output=False, fmin=0.5, fmax=4.0, window_padding=0.1, periods=(0.8, 0.3), time_limits=(-2.0, 8.0))` | `wavefield`、`station_pair` | 必填目标波场与两元素台站名序列；用于绘制代表性台站对诊断图。 |
|  | `include_convolution` | 默认 `True`，加入内侧第三台站卷积项；`False` 可对比仅使用外侧互相关的结果。 |
|  | `sqrt_spectrum` | 默认 `True`，对候选干涉谱应用平方根幅度归一化。 |
|  | `taper_output` | 默认 `False`；仅控制候选输出是否再次施加速度窗与带通，不影响输入已有的预处理。 |
|  | `fmin`、`fmax`、`window_padding` | 仅在 `taper_output=True` 时使用的输出频带，以及默认 0.1 s 的信号窗余量。 |
|  | `periods`、`time_limits` | 诊断图高斯窄带周期序列和显示相关时间范围。 |
| `denoise_wavefield_iteratively(wavefield, example_pair, *, threshold, first_iteration_convolution, max_iterations=6, sqrt_spectrum=True, taper_output=False, fmin=0.5, fmax=5.0, distance_threshold=0.0, signal_vmin=0.2, signal_vmax=2.0, window_padding=0.2, n_jobs=1)` | `wavefield`、`example_pair` | 必填目标波场与代表性台站对；每轮输出都会峰值归一化。 |
|  | `threshold` | 必填、非负的全波场相对 L2 变化阈值；变化不大于该值即停止。 |
|  | `first_iteration_convolution` | 必填布尔值；控制第一轮是否加入内侧卷积。后续各轮始终同时使用外侧互相关和内侧卷积。 |
|  | `max_iterations` | 最大完成迭代次数，默认 6。 |
|  | `sqrt_spectrum`、`taper_output`、`fmin`、`fmax` | 与单台站对诊断中的同名参数含义一致。 |
|  | `distance_threshold` | 参与候选叠加的最小台间距，默认 0 km。 |
|  | `signal_vmin`、`signal_vmax`、`window_padding` | TSI 信号窗范围及余量，默认 0.2 km/s、2.0 km/s 和 0.2 s。 |
|  | `n_jobs` | 台站对计算进程数，默认 1。 |
| `plot_denoised_result(wavefield, result, *, periods=(0.8, 0.3), time_limits=(-2.0, 8.0), jitter_duplicate_distances=True)` | `wavefield`、`result` | 必填原始波场与对应 `DenoisingResult`；绘制窄带前后波场、示例台站对和迭代历史。 |
|  | `periods`、`time_limits`、`jitter_duplicate_distances` | 控制窄带周期、显示时窗和相同距离记录的稳定小偏移。 |

`denoise_wavefield_iteratively()` 返回 `DenoisingResult`。其 `iterations`、`relative_changes`、`converged` 与 `stop_reason` 分别提供完成轮数、每轮变化、是否达到阈值及停止原因（`"threshold"` 或 `"max_iterations"`）。

### 10. NPZ 持久化

`Wavefield`、`MASW` 和 `DenoisingResult` 均提供版本化压缩 NPZ 的保存和加载接口，并以 `np.load(..., allow_pickle=False)` 安全读取。

| 对象 | 保存与加载 | 行为 |
|---|---|---|
| `Wavefield` | `.save(path, overwrite=False)` / `.load(path)` | 保存数据、台站对、距离、采样信息、分量与距离检查设置；加载时重新验证。 |
| `MASW` | `.save(path, overwrite=False)` / `.load(path)` | 保存配置、波场和已计算的频散结果；未计算对象也可保存。 |
| `DenoisingResult` | `.save(base_path, overwrite=False)` / `.load(base_path)` | 写入 `<name>_wavefield.npz` 与 `<name>_info.npz` 两个文件。 |

所有保存接口默认拒绝覆盖已有文件，父目录必须已存在。仅在确认需要替换结果时传入 `overwrite=True`。

### 11. 教程导航

- [RR Array：多模态瑞利波分离与三台干涉去噪](tutorial/RR_Array/README.md)
- [MARS DAS：单一主导 Scholte 波三台干涉去噪](tutorial/MARS_DAS/README.md)
- [RR Array Notebook](tutorial/RR_Array/run.ipynb)
- [MARS DAS Notebook](tutorial/MARS_DAS/run.ipynb)

根 README 介绍通用接口和推荐工作流；两个教程 README 进一步讨论数据背景、参数选择、物理解释、质量控制和预期图件。

### 12. 性能与可重复性

- `MASW.compute(n_jobs=...)` 和 `denoise_wavefield_iteratively(n_jobs=...)` 支持多进程。建议先以 `n_jobs=1` 验证，再根据 CPU 与内存提高进程数。
- 首次读取、预处理、MASW 与迭代去噪可能耗时较长。建议在 `processed/` 保存 NPZ，并在后续会话中复用。
- 四分量数据必须使用一致的预处理参数；修改频带、速度窗、极化公式、参考频散、距离阈值或 TSI 参数后，应重新生成相应下游缓存。
- 发表结果时，应记录程序版本、输入数据版本、频带、速度范围、相速度网格、距离阈值、迭代阈值与停止原因。

### 13. 常见问题

#### 缺少 source 或 receiver

错误信息会列出尝试过的字段。请设置 `sac.kevnm` 或 `stats.source` 作为 source，并设置 `sac.kstnm` 或 `stats.station` 作为 receiver；空字符串也会被视为缺失。

#### 台站名不含数字，或多个名称得到相同编号

程序依赖台站数字编号确定空间顺序。请使用如 `RR010`、`RR011` 的唯一名称，避免 `STA01` 与 `NODE01` 同时代表不同台站。

#### 距离顺序检查失败

首先检查 `sac.dist` 是否以 km 为单位，以及台站编号是否按空间位置排序。若阵列本来就不是直线，请自行读入 ObsPy `Stream` 后，使用 `check_distance_order=False` 构造 `Wavefield`。

#### 四分量极化分离提示台站对不匹配

ZZ、ZR、RZ、RR 必须包含相同的规范化台站对。检查是否有缺失 SAC 文件、重复 pair、不同的台站名回退字段或不一致的距离。

#### `fmax must be below the Nyquist frequency`

最高滤波频率必须严格小于采样率的一半。降低 `fmax` 或使用经过正确抗混叠处理的更高采样率数据。

#### MASW 绘图提示先调用 `compute()`

构造 `MASW(wavefield)` 只保存配置。请先调用 `.compute()`，再调用 `.plot()`。

#### 保存时提示文件已存在

默认行为用于保护已有结果。更换文件名，或在确认后传入 `overwrite=True`。

### 14. 引用

若您在研究中使用本程序包，请引用：

> Qiu, H., Niu, F., & Qin, L. (2021).  
> Denoising surface waves extracted from ambient noise recorded by 1-D linear array using three-station interferometry of direct waves.  
> *Journal of Geophysical Research: Solid Earth*, **126**, e2021JB021712.  
> <https://doi.org/10.1029/2021JB021712>

示例数据集的相关论文：

> Yuan, Y., Qiu, H., Chi, B., & Qin, L. (2026).  
> *Mitigating the Resolution–SNR Trade-Off in DAS Ambient Noise Imaging: Application to Monterey Bay*.  
> Manuscript under review at *Journal of Geophysical Research: Solid Earth*.

### 15. 作者

**袁宇嵩（Yusong Yuan）**  
中国地质大学（武汉）  
中国科学技术大学  
✉️ yuanyusong25@gmail.com

**裘鸿瑞（Hongrui Qiu）**  
中国地质大学（武汉）  
✉️ qiuhongrui@gmail.com  
✉️ qiuhongrui@cug.edu.cn

---

## English

### 1. Overview

TSI-Denoising is a research and teaching Python package for ambient-noise cross-correlation (ANC) waveforms recorded by dense 1-D arrays. It wraps ObsPy `Stream` objects in a validated `Wavefield` model and provides a connected workflow for data ingestion, preprocessing, dispersion diagnosis, modal separation, and iterative three-station interferometry (TSI) denoising.

The repository includes two teaching cases prepared for the **12th Seismological Algorithms and Programs Training Course** (10–12 August 2026):

- **MARS DAS**: TSI denoising of a dominant Scholte-wave mode in submarine DAS data;
- **RR Array**: multimode Rayleigh-wave separation and mode-by-mode TSI denoising in a dense fault-crossing array.

Unless an individual interface states otherwise, distances are in km, phase velocities in km/s, correlation times in s, and frequencies in Hz.

> [!IMPORTANT]
> Raw teaching data are not included with the source repository. Run `python retrieve_datasets.py` to download and safely extract the two tutorial inputs. Each tutorial generates or reuses NPZ caches and results under `processed/`.

### 2. Repository Structure

~~~text
TSI-Denoising/
├── src/tsi_denoising/              # Python package
│   ├── io/                         # SAC-directory ingestion
│   ├── mode_separation/            # Polarization and phase matching
│   ├── denoising/                  # TSI and diagnostic plotting
│   ├── wavefield.py                # Wavefield data model
│   ├── preprocessing.py            # Shared preprocessing
│   └── masw.py                     # MASW calculation and plotting
├── tutorial/
│   ├── RR_Array/                   # Four-component Rayleigh-wave case
│   └── MARS_DAS/                   # Submarine DAS Scholte-wave case
├── pyproject.toml                  # Package metadata and runtime dependencies
├── environment.yml                 # Recommended Conda environment
├── requirements.txt                # Runtime dependencies
└── README.md
~~~

### 3. Installation and Data Preparation

#### 3.1 Get the source code

Install [Git](https://git-scm.com/) and [Miniforge](https://github.com/conda-forge/miniforge) (or Miniconda), then clone the repository:

~~~bash
git clone https://github.com/YuanYusung/TSI-Denoising.git
cd TSI-Denoising
~~~

#### 3.2 Create the Conda environment

The repository's `environment.yml` defines the recommended environment:

~~~bash
conda env create -f environment.yml
conda activate tsi-denoising
python -m pip install -e .
~~~

The environment uses Python 3.10 and includes NumPy, SciPy, Matplotlib, ObsPy, JupyterLab, and ipykernel. Editable installation makes local source changes immediately available without reinstalling.

#### 3.3 Use Jupyter and verify the package

To select this environment explicitly in Jupyter, register its kernel:

~~~bash
python -m ipykernel install --user --name tsi-denoising --display-name "Python (TSI-Denoising)"
~~~

Select `Python (TSI-Denoising)` in a tutorial notebook. Verify the installation with:

~~~bash
python -c "import tsi_denoising; print('TSI-Denoising imported successfully')"
~~~

#### 3.4 Download the teaching data

Run the following command from the repository root:

~~~bash
python retrieve_datasets.py
~~~

The script downloads the public archive and extracts only missing `tutorial/RR_Array/input/` and `tutorial/MARS_DAS/input/` directories; existing directories are never overwritten. Reproducible intermediate products belong in each tutorial's `processed/` directory.

### 4. Main Capabilities

- Recursively read SAC cross-correlations and construct validated `Wavefield` objects;
- Normalize station-name sources, pair directions, sampling metadata, and distances;
- Symmetrize correlations and apply distance-dependent velocity tapers and zero-phase band-pass filters;
- Compute frequency–phase-velocity MASW images with phase-shift stacking;
- Separate retrograde and prograde Rayleigh-wave components from ZZ, ZR, RZ, and RR data;
- Extract a target mode by phase-matched filtering with a reference dispersion curve;
- Diagnose the geometry, input gathers, and output of TSI for one station pair;
- Iteratively denoise a complete wavefield serially or with multiple processes;
- Save and restore wavefields, MASW products, and denoising results as NPZ files.

### 5. Processing Workflow

~~~text
SAC station-pair cross-correlations
          |
          v
Read, normalize, and validate
          |
          v
Shared preprocessing and MASW diagnosis
          |
          +-------------------------------+
          |                               |
          v                               v
Single dominant mode in target band       Four-component, multimode data
for example, MARS DAS                     for example, RR Array
          |                               |
          |                         Polarization separation
          |                               |
          |                         Phase-matched separation
          |                               |
          +---------------+---------------+
                          |
              Mode-specific TSI denoising and QC
~~~

### 6. Scope and Limitations

The package targets 1-D or approximately 1-D arrays whose station numbering represents spatial order. TSI assumes that the target surface wave propagates mainly along the array and that its geometry and phase relations remain stable across station combinations.

- Two-dimensional arrays, curved lines, and strong lateral heterogeneity can violate the distance-order and 1-D propagation assumptions.
- Polarization separation requires ZZ, ZR, RZ, and RR wavefields with identical station pairs, distances, sampling rates, and time axes.
- ZR/RZ signs, radial-positive direction, and correlation-time conventions must be consistent when data are produced.
- A MASW amplitude peak does not automatically identify a modal order; use theoretical dispersion, polarization, and spatial continuity.
- TSI does not automatically remove modal cross terms. Separate multimode wavefields reliably before denoising.

### 7. Input Data Requirements

#### 7.1 SAC waveforms

Each SAC file must contain one one-dimensional cross-correlation for one station pair. `read_sac_directory()` searches recursively and, by default, matches common case variants of `.sac`, `.SAC`, and `.SAC_s`.

| Field | Preferred source | Fallback | Requirement |
|---|---|---|---|
| Virtual source name | `stats.sac.kevnm` | `stats.source` | Non-empty and contains a numeric identifier |
| Receiver name | `stats.sac.kstnm` | `stats.station` | Non-empty and contains a numeric identifier |
| Pair distance | `stats.sac.dist` | None | Finite and positive, in km |
| Correlation-time origin | `stats.sac.b` | None | Finite, in s |

All traces in one `Wavefield` must have identical sample counts, sampling intervals, and `sac.b`; waveform data must be finite one-dimensional real arrays; directed pairs must be unique; and distinct station names must not resolve to the same numeric identifier.

#### 7.2 Station numbering and pair normalization

The final continuous digit group in a station name defines its spatial number: `RR010 → 10`, `LINE2_ST040 → 40`, and `DAS00345 → 345`.

When a `Wavefield` is constructed or `read_sac_directory()` is called, the package removes autocorrelations, normalizes pairs to lower-numbered source toward higher-numbered receiver, removes duplicate reverse-direction records, and swaps metadata and reverses waveforms when needed. By default, it also enforces increasing distance with station number.

For curved lines, two-dimensional arrays, or non-spatial numbering, read the data into an ObsPy `Stream` and construct `Wavefield(..., check_distance_order=False)`. This disables only distance-order validation.

#### 7.3 Correlation-time axis

Correlations must cover the required positive and negative times. `preprocess_stream()` and `Wavefield.preprocess()` symmetrize time branches; zero-pad a missing negative branch and set `sac.b` correctly before processing.

### 8. Recommended Workflows

The root README provides the common workflow and API reference. Complete, runnable procedures, QC guidance, and interpretation are maintained in the tutorials.

#### 8.1 Single-component data with one dominant mode

This branch applies when one surface-wave mode dominates the target band, such as the Scholte wave in the MARS DAS case.

~~~text
Read SAC correlations
        |
Build and validate Wavefield
        |
Symmetrize, velocity-window, and band-pass filter
        |
Time–distance and MASW diagnosis
        |
Confirm one dominant target-band mode
        |
One-pair TSI diagnostic
        |
Iterative full-wavefield denoising and QC
~~~

See the [MARS DAS tutorial](tutorial/MARS_DAS/README.md) and [MARS DAS notebook](tutorial/MARS_DAS/run.ipynb).

#### 8.2 Four-component, multimode Rayleigh-wave data

This branch applies to ZZ, ZR, RZ, and RR correlations containing multiple Rayleigh-wave modes.

~~~text
Read four-component SAC data
        |
Shared preprocessing and consistency checks
        |
Wavefield and MASW diagnosis
        |
Retrograde–prograde polarization separation
        |
Phase-matched modal separation
        |
Confirm spatial and dispersive continuity of the target mode
        |
Denoise each target mode separately with TSI and QC
~~~

See the [RR Array tutorial](tutorial/RR_Array/README.md) and [RR Array notebook](tutorial/RR_Array/run.ipynb).

> [!IMPORTANT]
> Polarization separation is not modal separation. A retrograde or prograde wavefield can still contain multiple modes. Before TSI, verify reliable target-mode separation using dispersion images, a reference curve, and spatial continuity.

### 9. Key Interfaces and Parameters

Import public interfaces from `tsi_denoising`. Distances are km, velocities km/s, times s, and frequencies Hz unless noted otherwise.

#### 9.1 Wavefields, SAC ingestion, and preprocessing

| Interface | Parameter | Description |
|---|---|---|
| `Wavefield(stream, *, component=None, copy=True, check_distance_order=True)` | `stream` | Required non-empty ObsPy `Stream` satisfying Section 7. |
|  | `component` | Optional metadata label such as `"ZZ"` or `"RR"`; it does not alter processing. |
|  | `copy` | Default `True`; defensively copy traces. Use `False` only when the caller controls the input stream. |
|  | `check_distance_order` | Default `True`; enforce compatible number and distance ordering. |
| `read_sac_directory(directory, pattern=None, *, component=None)` | `directory` | Required SAC root directory; recursively returns a validated, normalized `Wavefield`. |
|  | `pattern` | One wildcard or an iterable of wildcards; common SAC suffixes are matched by default. |
|  | `component` | Optional component label; the directory name is used when omitted. |
| `preprocess_stream(wavefield, fmin=0.5, fmax=5.0, vmin=0.1, vmax=2.5, taper_fraction=0.05)` | `wavefield` | Required input; returns a new `Wavefield` and does not modify it. |
|  | `fmin`, `fmax` | Zero-phase band-pass limits; require `0 < fmin < fmax < Nyquist`. |
|  | `vmin`, `vmax` | Distance-dependent velocity-window bounds; the causal window spans `distance / vmax` to `distance / vmin`. |
|  | `taper_fraction` | Cosine-taper fraction at both window edges; must be between 0 and 0.5. |

`Wavefield.preprocess()` accepts the same parameters but replaces the object's data in place. Use `preprocess_stream()` to preserve the input wavefield.

#### 9.2 MASW

| Interface | Parameter | Description |
|---|---|---|
| `MASW(wavefield, velocities=None, fmin=0.5, fmax=5.0, padding_factor=5, dist_threshold=0.2)` | `wavefield` | Required, normally preprocessed input. Construction only stores configuration; call `.compute()` to calculate the image. |
|  | `velocities` | Phase-velocity sample array; defaults to 231 linear samples from 0.2 to 2.5 km/s. |
|  | `fmin`, `fmax` | MASW frequency range within the valid Nyquist interval. |
|  | `padding_factor` | Integer FFT zero-padding multiplier of at least 1; it refines frequency sampling, not physical resolution. |
|  | `dist_threshold` | Minimum pair distance included in imaging. |
| `MASW.compute(n_jobs=1)` | `n_jobs` | Number of worker processes for frequency rows. Start with 1 before increasing it. |
| `compute_masw(wavefield, velomin=0.2, velomax=2.5, fmin=0.5, fmax=5.0, padding_factor=5, *, velocities=None, dist_threshold=0.2, n_jobs=1)` | `velomin`, `velomax` | Bounds of the 231-sample default grid when `velocities` is omitted. |
|  | Remaining parameters | Match `MASW` and `.compute(n_jobs=...)`; returns `(velocity, frequency, amplitude)`. |

#### 9.3 Polarization and phase-matched separation

| Interface | Parameter | Description |
|---|---|---|
| `polarization_separate(*, zz, rz, zr, rr, swap_polarization=False, use_rr=True, vmin=0.1, vmax=2.5, taper_fraction=0.05)` | `zz`, `rz`, `zr`, `rr` | Required keyword-only `Wavefield` inputs with identical normalized pairs, distances, sampling rates, and time axes. |
|  | `swap_polarization` | Default `False`; set `True` to exchange retrograde and prograde labels when data conventions are reversed. |
|  | `use_rr` | Default `True`; include RMS-normalized RR in the prograde combination. `False` uses the three-component combination. |
|  | `vmin`, `vmax`, `taper_fraction` | Distance-dependent output taper parameters. Inputs are not modified. |
| `phase_match_separate(wavefield, reference_curve=None, *, fmin=0.5, fmax=5.0, t_window=0.2, keep_positive=True, return_reference=False, masw_cache_path=None, vmin=0.1, vmax=2.5, taper_fraction=0.05)` | `wavefield` | Required candidate-mode wavefield. This function does not band-pass filter its input. |
|  | `reference_curve` | Optional `(N, 2)` array of frequency and phase velocity, or a `(frequency, velocity)` pair. When omitted, a cached or newly computed MASW image is displayed for picking. |
|  | `fmin`, `fmax`, `t_window` | Phase-matching band and positive Gaussian time-window width. |
|  | `keep_positive` | Default `True`; retain the positive-time branch. Use `False` to retain the negative-time branch. |
|  | `return_reference` | Default `False`; when true, return `(separated_wavefield, used_curve)`. |
|  | `masw_cache_path` | Optional computed MASW NPZ. It is mutually exclusive with `reference_curve`. |

#### 9.4 TSI denoising and results

| Interface | Parameter | Description |
|---|---|---|
| `denoise_station_pair_demo(wavefield, station_pair, *, include_convolution=True, sqrt_spectrum=True, taper_output=False, fmin=0.5, fmax=4.0, window_padding=0.1, periods=(0.8, 0.3), time_limits=(-2.0, 8.0))` | `wavefield`, `station_pair` | Required target wavefield and two-station sequence; renders a representative diagnostic. |
|  | `include_convolution` | Default `True`; include inner-third-station convolutions. Set false to compare outer correlations only. |
|  | `sqrt_spectrum` | Default `True`; apply square-root spectral amplitude normalization. |
|  | `taper_output` | Default `False`; controls only repeated output tapering/filtering, not input preprocessing. |
|  | `fmin`, `fmax`, `window_padding` | Output band when tapering is enabled, and the default 0.1 s signal-window padding. |
|  | `periods`, `time_limits` | Diagnostic narrow-band periods and displayed correlation-time range. |
| `denoise_wavefield_iteratively(wavefield, example_pair, *, threshold, first_iteration_convolution, max_iterations=6, sqrt_spectrum=True, taper_output=False, fmin=0.5, fmax=5.0, distance_threshold=0.0, signal_vmin=0.2, signal_vmax=2.0, window_padding=0.2, n_jobs=1)` | `wavefield`, `example_pair` | Required target wavefield and representative pair; each output iteration is peak normalized. |
|  | `threshold` | Required non-negative whole-wavefield relative-L2 stopping threshold. |
|  | `first_iteration_convolution` | Required Boolean controlling inner convolutions only in the first iteration; later iterations always include them. |
|  | `max_iterations`, `n_jobs` | Maximum completed iterations (default 6) and worker-process count (default 1). |
|  | `distance_threshold`, `signal_vmin`, `signal_vmax`, `window_padding` | Candidate-stack minimum distance and signal-window bounds; defaults are 0 km, 0.2 km/s, 2.0 km/s, and 0.2 s. |
| `plot_denoised_result(wavefield, result, *, periods=(0.8, 0.3), time_limits=(-2.0, 8.0), jitter_duplicate_distances=True)` | `wavefield`, `result` | Plot narrow-band input/output wavefields and iteration history. The remaining parameters set periods, time limits, and duplicate-distance jitter. |

`denoise_wavefield_iteratively()` returns `DenoisingResult`. Its `iterations`, `relative_changes`, `converged`, and `stop_reason` report the completed count, change history, convergence state, and either `"threshold"` or `"max_iterations"`.

### 10. NPZ Persistence

`Wavefield`, `MASW`, and `DenoisingResult` use versioned compressed NPZ files read with `np.load(..., allow_pickle=False)`.

| Object | Save and load | Behavior |
|---|---|---|
| `Wavefield` | `.save(path, overwrite=False)` / `.load(path)` | Stores data, station pairs, distances, sampling metadata, component, and distance-order setting; loading validates again. |
| `MASW` | `.save(path, overwrite=False)` / `.load(path)` | Stores configuration, wavefield, and any computed dispersion result. |
| `DenoisingResult` | `.save(base_path, overwrite=False)` / `.load(base_path)` | Writes `<name>_wavefield.npz` and `<name>_info.npz`. |

All save interfaces require an existing parent directory and refuse overwriting by default. Pass `overwrite=True` only when replacement is intentional.

### 11. Tutorials

- [RR Array: multimode Rayleigh-wave separation and TSI denoising](tutorial/RR_Array/README.md)
- [MARS DAS: TSI denoising of a dominant Scholte-wave mode](tutorial/MARS_DAS/README.md)
- [RR Array notebook](tutorial/RR_Array/run.ipynb)
- [MARS DAS notebook](tutorial/MARS_DAS/run.ipynb)

The root README covers common interfaces and workflows. The two tutorial READMEs cover data context, parameter selection, physical interpretation, QC, and expected figures.

### 12. Performance and Reproducibility

- `MASW.compute(n_jobs=...)` and `denoise_wavefield_iteratively(n_jobs=...)` support multiprocessing. Start at `n_jobs=1`, then increase only as CPU and memory allow.
- First-pass ingestion, preprocessing, MASW, and iterative denoising can take time. Save reusable NPZ products under `processed/`.
- Four-component inputs must use identical preprocessing. Changing the frequency band, velocity window, polarization formula, reference curve, distance threshold, or TSI settings requires regenerating downstream caches.
- For a publication, record package and input-data versions, frequency and velocity ranges, the velocity grid, distance threshold, iteration threshold, and stop reason.

### 13. Troubleshooting

#### Missing source or receiver

The error lists the fields attempted. Set `sac.kevnm` or `stats.source` for the source and `sac.kstnm` or `stats.station` for the receiver. Empty strings are treated as missing.

#### A station name has no digits, or two names resolve to one number

Use unique numeric station names such as `RR010` and `RR011`. Avoid assigning distinct stations names such as `STA01` and `NODE01`, because both resolve to 1.

#### Distance-order validation fails

Check that `sac.dist` is in km and that numbers follow spatial order. For a non-linear array, read a `Stream` yourself and construct `Wavefield(..., check_distance_order=False)`.

#### Polarization separation reports mismatched pairs

ZZ, ZR, RZ, and RR must contain the same normalized pairs. Check for missing files, duplicate pairs, inconsistent fallback names, or distance mismatches.

#### `fmax must be below the Nyquist frequency`

The upper filter frequency must be strictly below half the sampling rate. Lower `fmax` or use correctly anti-aliased data sampled faster.

#### MASW plotting asks for `compute()`

`MASW(wavefield)` only stores configuration. Call `.compute()` before `.plot()`.

#### Saving reports that a file already exists

Saving protects existing outputs by default. Choose another name or explicitly set `overwrite=True`.

### 14. Citation

If you use this package in research, please cite:

> Qiu, H., Niu, F., & Qin, L. (2021).  
> Denoising surface waves extracted from ambient noise recorded by 1-D linear array using three-station interferometry of direct waves.  
> *Journal of Geophysical Research: Solid Earth*, **126**, e2021JB021712.  
> <https://doi.org/10.1029/2021JB021712>

For the example data set, please also cite:

> Yuan, Y., Qiu, H., Chi, B., & Qin, L. (2026).  
> *Mitigating the Resolution–SNR Trade-Off in DAS Ambient Noise Imaging: Application to Monterey Bay*.  
> Manuscript under review at *Journal of Geophysical Research: Solid Earth*.

### 15. Authors

**Yusong Yuan**  
China University of Geosciences (Wuhan)  
University of Science and Technology of China  
✉️ yuanyusong25@gmail.com

**Hongrui Qiu**  
China University of Geosciences (Wuhan)  
✉️ qiuhongrui@gmail.com  
✉️ qiuhongrui@cug.edu.cn
