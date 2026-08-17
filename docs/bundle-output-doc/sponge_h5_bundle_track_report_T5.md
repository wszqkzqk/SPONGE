# SPONGE H5 Bundle Track T5 报告

## Track

T5 Restart H5。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T5 细化章节。契约要求 `*.spgr.h5` 只包含一个可用于 launch 的 restart state，使用 H5MD-compatible 的单帧 particle state，并在 `/parameters/restart` 下存放 SPONGE-specific continuation 数据。本轮推进 runtime structural restart 接入。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 添加 restart H5 writer construction | 已完成 | `MD_INFORMATION::trajectory_output::Export_H5_Restart_File` 在 `output_h5_restart_path` 启用时创建 `HighFiveBackend` 和 `RestartH5Writer`。 |
| 写入恰好一个可 launch state | 已完成 | 每次 restart cadence 新建 writer 并覆盖 `output_h5_restart_path`，`RestartH5Writer::Write_Structural_State` 保证单 writer 内只能写一个 structural state。 |
| coordinate 和 box state runtime 写入 | 已完成 | `Export_H5_Restart_File` 调用 `Crd_Vel_Device_To_Host` 后写 `/particles/all/position/value` 和 `/particles/all/box/edges/value`。 |
| velocity state runtime 写入 | 已完成 | `Export_H5_Restart_File` 以 `include_velocity=true` 定义并写入 `/particles/all/velocity/value`。 |
| `/run` metadata | 已完成 | `RestartH5Writer::Write_Structural_State` 写入 `/run/current_step`、`/run/current_time` 和 `/run/state_type`。 |
| `/parameters/restart` base layout | 已完成 | `RestartH5Writer::Open`/`Ensure_Base_Layout` 创建 `/parameters/restart`、thermostat/barostat/bias 子组。 |
| NHC restart state API | 已完成 | `Export_H5_Restart_File` 将 NHC coordinate/velocity 打包并调用 `RestartH5Writer::Write_Nose_Hoover_Chain_State`。 |
| SITS state API/runtime route | 已完成 | `Export_H5_Restart_File` 接收 SITS module name 与 `nk` buffer，并调用 `RestartH5Writer::Write_Sits_State` 写 `/parameters/restart/bias/sits/<module>/nk`。 |
| Metad path-dependent text snapshots | 部分完成 | `Export_H5_Restart_File` 接收 metad module name 和已有 text export 路径，将存在的 `hills/history/edge/potential_export/direct_export` 写入 `/parameters/restart/bias/meta/<name>/...`。 |
| 避免 topology/protocol ownership | 已完成 | runtime structural restart 只写粒子状态和 run metadata，不写 topology/protocol 定义。 |
| legacy restart gating | 已完成 | `Main_Print` 中 restart cadence 先写 H5 restart；legacy `Export_Restart_File` 和 `nhc.Save_Restart_File` 仅在 `Should_Write_Legacy_Restart` 为 true 时调用。 |
| retention 行为 | 采用单文件覆盖 | H5 restart 只保留一帧，按契约覆盖 `output_h5_restart_path`；legacy `max_restart_export_count` 仅保留给显式 legacy sidecar。 |

## 架构说明

T5 runtime 接入后的 restart 输出数据流为：

```text
Main_Initial
  -> md_info.output.Initial_H5_Restart()

Main_Print restart cadence
  -> md_info.output.Export_H5_Restart_File(&controller)
  -> if legacy restart explicitly enabled or no H5 output:
       md_info.output.Export_Restart_File()
       nhc.Save_Restart_File()
```

H5 restart writer 生命周期是 per-export 的：

```text
Export_H5_Restart_File
  -> HighFiveBackend
  -> RestartH5Writer::Open
  -> Define_Structural_State(atom_count, include_velocity=true)
  -> Write_Structural_State(step,time,coordinates,box,velocities)
  -> optional Write_Nose_Hoover_Chain_State(...)
  -> optional Write_Sits_State(...)
  -> optional Write_Metad_State_Text(...)
  -> Finalize
  -> Close
```

该策略与“restart 仅保留一帧”的约定一致，避免把 restart 文件变成轨迹文件。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| Metad accumulated state 结构化二进制 schema | 当前仅写已有 text snapshots；grid/scatter/hill 内存结构的 canonical restart schema 仍需单独设计。 |
| H5 restart reader/resolver | 当前只实现输出；作为下一段 run 输入仍需 reader/resolver。 |
| HDF5 inspection 和 continuation 验证 | 未运行编译/测试/验证。 |

## 审查点

- `output_h5_restart_path` 是单文件输出路径；每次 restart cadence 覆盖写出最新 single-frame restart。
- 启用任意 H5 输出后，legacy restart sidecar 不再隐式写出；只有显式 `rst` key 或完全未启用 H5 输出时才写 legacy restart 和 NHC legacy restart。
- H5 restart 当前始终写 velocity，因为 SPONGE restart continuation 通常需要 velocity。
- NHC 和 SITS `nk` continuation state 已接入 H5 restart；metad 已接入 text snapshot 级别的 accumulated state，但尚未完成结构化二进制 restart schema。
