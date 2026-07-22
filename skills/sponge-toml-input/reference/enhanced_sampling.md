# 增强采样参数

## 集体变量（CV）

| 参数 | 类型 | 说明 |
|------|------|------|
| `cv_in_file` | string | 集体变量定义文件路径 |

也可通过 `[CV]` section 内联定义：

| 参数 | 类型 | 说明 |
|------|------|------|
| `CV_type` | string | CV 类型 |
| `CV_period` | int | CV 更新频率 |
| `CV_minimal` / `CV_maximum` | float | CV 取值范围 |

## Metadynamics

通过 `[META]` 或 `[meta]` section 配置：

| 参数 | 类型 | 说明 |
|------|------|------|
| `sink` | string | Sink metadynamics 模式 |

## Voronoi milestoning 首达检测

在 CV 文件中配置：

```toml
[voronoi_detector]
CV = ["distance"]
milestone_file = "milestones.txt"
source_interface = "S_0_1"
```

`source_interface` 是无向发射界面。轨迹在它两侧 cell 之间的回穿只更新
当前侧并计数，不终止，也不重置首次到达时间；首次命中其他已声明邻接界面
才终止并导出 `voronoi_hit_<interface>` 坐标和速度。检测按 MD step 对已提交
状态采样，不做连续时间插值；当前只支持单 MPI 进程。milestone manifest 的
完整格式和初态要求见用户文档 `collective-variables.md`。

## SITS

SITS（Self-guided Integrated Tempering Sampling）参数：

| 参数 | 类型 | 说明 |
|------|------|------|
| `SITS_mode` | string | 运行模式 |
| `SITS_atom_numbers` | int | 参与 SITS 的原子数 |
| `SITS_k_numbers` | int | k 空间点数 |
| `SITS_T_low` / `SITS_T_high` | float | 温度范围（K） |
| `SITS_record_interval` | int | 记录间隔 |
| `SITS_update_interval` | int | 更新间隔 |
| `SITS_nk_fix` | int | 固定 k 空间点数 |
| `SITS_nk_in_file` | string | k 空间输入文件 |
| `SITS_pe_a` / `SITS_pe_b` | float | 势能参数 |
| `SITS_fb_interval` | int | 反馈间隔 |

`SITS_mode` 可选值：

| 值 | 说明 |
|----|------|
| `"observation"` | 观测阶段 |
| `"iteration"` | 迭代阶段 |
| `"production"` | 生产阶段 |
| `"empirical"` | 经验模式 |
| `"amd"` | AMD 模式 |
| `"gamd"` | GaMD 模式 |
