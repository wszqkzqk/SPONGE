# SPONGE H5 Bundle Track T3 报告

## Track

T3 单文件 trajectory H5。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T3 细化章节。本 track 的目标是 non-VDS trajectory H5 routing，同时保持 trajectory 和 observable cadence 相互独立。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 添加 trajectory H5 writer construction | 已完成 | `MD_INFORMATION::trajectory_output::Initial_H5_Trajectory` 在 `output_h5_trajectory_path` 启用且 `vds=false` 时创建 `HighFiveBackend` 和 `TrajectoryH5Writer`。 |
| 写入初始 H5MD metadata 和 SPONGE schema metadata | 已完成 | `Initial_H5_Trajectory` 调用 `TrajectoryH5Writer::Open_Single_File`，后者通过 `H5MDWriter` 初始化公共 H5MD/SPONGE layout。 |
| 定义 particle datasets | 已完成 | `Initial_H5_Trajectory` 调用 `Define_Particle_Datasets`，定义 position、box，以及显式 velocity/force 输出存在时的 velocity/force dataset。 |
| append position 和 box frames | 已完成 | `Main_Print` 的 trajectory cadence 内调用 `Append_H5_Trajectory_Frame`，写入 `/particles/all/position/value` 和 `/particles/all/box/edges/value`。 |
| append velocity frames | 部分完成 | 当 legacy velocity 输出显式启用并设置 `is_vel_traj` 时，H5 trajectory 同步写 velocity；默认不强制写 velocity。 |
| append force frames | 已完成 | 沿用现有显式 `frc` 开关；`is_frc_traj` 启用时，H5 trajectory 定义 `/particles/all/force/value` 并在 particle frame 中同步写入 force。 |
| 捕获 `Step_Print` 注册 metadata | 已完成 | `Initial_H5_Trajectory` 在所有 `Step_Print_Initial` 后运行，读取 `controller->outputs_key`，写入原始列名和 HDF5-safe 列名。 |
| append scalar observables | 已完成 | `Main_Print` 在 `Print_To_Screen_And_Mdout` 重置 `outputs_content` 前调用 `Append_H5_Observable_Frame`。 |
| 保持 mdout 和 trajectory cadence 分离 | 已完成 | observable frame 写入位于 `Check_Mdout_Step()` 分支；particle frame 写入位于 `Check_Trajectory_Step()` 分支。 |
| 修复 H5 bundle 模式下 legacy NULL sidecar 崩溃风险 | 已完成 | `trajectory_output::Initial` 和 `Append_Crd_Traj_File`/`Append_Box_Traj_File` 对 legacy file pointer 增加 NULL guard。 |
| finalize status | 已完成 | `Main_Clear` 调用 `Finalize_H5_Trajectory`，正常结束时 finalize 并 close writer。 |

## 架构说明

T3 runtime 接入后的数据流为：

```text
Main_Initial
  -> controller.Print_First_Line_To_Mdout()
  -> md_info.output.Initial_H5_Trajectory()

Main_Print mdout cadence
  -> modules Step_Print(...)
  -> md_info.output.Append_H5_Observable_Frame()
  -> controller.Print_To_Screen_And_Mdout()

Main_Print trajectory cadence
  -> md_info.Crd_Vel_dd_to_Device(...)
  -> legacy crd/vel/box append if explicitly active
  -> md_info.output.Append_H5_Trajectory_Frame()
      -> position/box
      -> optional velocity when explicit vel output is active
      -> optional force when explicit frc output is active

Main_Clear
  -> md_info.output.Finalize_H5_Trajectory()
```

H5 writer 结构为：

```text
trajectory_output
  -> HighFiveBackend
  -> TrajectoryH5Writer
  -> H5MDWriter
```

`Append_H5_Observable_Frame` 使用 `controller->outputs_content` 中已经格式化但尚未 reset 的值。HDF5 dataset name 由 mdout 原始列名 sanitize 得到，并为重复名称追加数字后缀。

box 写入使用 H5MD `/particles/all/box/edges/value` 的 3x3 edges 表达：PBC 时来自 `md_info->pbc.cell`，非 PBC 时退化为由 `sys.box_length` 组成的对角矩阵。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| VDS trajectory runtime | 属于 T6 runtime 子阶段。 |
| observable-only H5 runtime | 属于 T4 runtime 子阶段。 |
| HDF5 inspection 和 legacy 数值对照 | 未运行编译/测试/验证。 |

## 审查点

- `Initial_H5_Trajectory` 在 `Print_First_Line_To_Mdout` 之后调用，是为了捕获完整 `Step_Print_Initial` 列集合。
- 启用任意 H5 output 后，legacy sidecar 默认关闭；T3 同时修复了该模式下 legacy file pointer 为 NULL 时的 trajectory 初始化和 append 崩溃风险。
- 当前 velocity 和 force 只有在 legacy `vel`/`frc` 输出显式开启时才写入 H5 trajectory。后续如果希望 H5 trajectory 独立控制 velocity/force，需要新增明确策略键。
