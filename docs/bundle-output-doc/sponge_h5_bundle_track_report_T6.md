# SPONGE H5 Bundle Track T6 报告

## Track

T6 VDS trajectory 输出。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T6 细化章节。契约要求提供 user-facing trajectory H5MD wrapper，shard 由 writer 派生，按 trajectory frame count 轮转，记录 manifest metadata，并且不允许用户配置 shard directory。本轮推进 runtime VDS 分支接入。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 派生内部 shard directory | 已完成 | T2 resolver 生成 `ResolvedOutputPlan::trajectory.derived_shard_root`，`VdsTrajectoryH5Writer` 消费该路径。 |
| runtime 选择 VDS writer | 已完成 | `Initial_H5_Trajectory` 在 `output_h5_trajectory_vds=true` 时创建 `HighFiveBackendFactory` 和 `VdsTrajectoryH5Writer`。 |
| 打开 wrapper H5MD | 已完成 | `VdsTrajectoryH5Writer::Open` 使用 user-facing `output_h5_trajectory_path` 打开 wrapper，并写入 chunk size metadata。 |
| 第一帧 trajectory 到来时打开 shard 0 | 已完成 | runtime `Append_H5_Trajectory_Frame` 调用 `VdsTrajectoryH5Writer::Append_Particle_Frame`，内部在首帧时 `Rotate_To_New_Shard`。 |
| 按 chunk size 轮转 shard | 已完成 | `VdsTrajectoryH5Writer` 使用 `current_shard_frame_count_ >= chunk_size_` 触发轮转；默认 chunk size 仍为 20。 |
| 写入 shard-local H5MD-like dataset | 已完成 | 每个 shard 使用 `TrajectoryH5Writer::Open_Single_File` 和 T3 相同 layout。 |
| mark shard complete after frame data | facade 层已完成 | `Complete_Current_Shard` 只 finalize 非空 shard 并记录 manifest entry。 |
| 维护 manifest | 已完成 | `VdsShardManifestEntry` 记录 index、path、frame range、step/time range 和 status。 |
| store manifest in wrapper | 已完成 | finalize 时 `Write_Manifest_To_Wrapper` 写 `/parameters/sponge/output/shard_manifest`。 |
| finalize 前 manifest validation | 已完成 | `VdsTrajectoryH5Writer::Finalize` 在写 wrapper manifest 前调用 complete manifest validation，不连续或 incomplete shard 会标记 wrapper failed。 |
| 支持独立 observable frame count | 部分完成 | VDS writer 维护独立 observable count；runtime 仅在至少一个 trajectory shard 打开后写 observable，避免首个 trajectory frame 前无 shard 的错误。 |
| finalize wrapper without payload copy | 已完成 | `Finalize_H5_Trajectory` 在 VDS 模式调用 `VdsTrajectoryH5Writer::Finalize`，写 manifest 并 finalize wrapper，不复制 shard payload。 |

## 架构说明

T6 runtime 接入后的 trajectory writer 选择为：

```text
output_h5_trajectory_path set, output_h5_trajectory_vds=false
  -> HighFiveBackend
  -> TrajectoryH5Writer

output_h5_trajectory_path set, output_h5_trajectory_vds=true
  -> HighFiveBackendFactory
  -> VdsTrajectoryH5Writer
      -> wrapper H5MDWriter
      -> shard TrajectoryH5Writer instances
```

runtime 数据流仍复用 T3 接入点：

```text
Main_Initial
  -> md_info.output.Initial_H5_Trajectory()

Main_Print trajectory cadence
  -> md_info.output.Append_H5_Trajectory_Frame()
      -> VdsTrajectoryH5Writer::Append_Particle_Frame()
      -> shard rotation by chunk size

Main_Clear
  -> md_info.output.Finalize_H5_Trajectory()
      -> complete current shard
      -> validate complete manifest
      -> write manifest to wrapper
      -> finalize wrapper
```

shard 路径仍由 writer 派生：

```text
<derived_shard_root>/segment_000000.spg.h5md
<derived_shard_root>/segment_000001.spg.h5md
```

用户不能通过 mdin 设置 shard directory。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| 真实 HDF5 VDS dataset materialization | 当前 wrapper 写 manifest 和 metadata；`WriterBackend` 尚未提供创建 HDF5 virtual dataset 的 API。 |
| 首个 trajectory frame 前的 observable shard | `VdsTrajectoryH5Writer` 要求 observable frame 依附已有 shard；runtime 当前跳过首个 trajectory frame 前的 VDS observable。 |
| force VDS 输出 | force H5 在 T3 中已暂缓，VDS 同步暂缓。 |
| runtime 目录创建/fsync 策略验证 | HighFiveBackend 创建父目录，但未运行实际文件验证。 |
| 21 帧生成两个 shard 验证 | 未运行编译/测试/验证。 |

## 审查点

- `output_h5_trajectory_vds=true` 不新增任何路径 key；仍只使用 `output_h5_trajectory_path` 和内部派生 shard root。
- `output_h5_trajectory_chunk_size` 的语义是 trajectory frame 数，不是 MD step 数。
- wrapper 当前不是完整 VDS materialized file，而是 manifest-bearing wrapper；真实 VDS dataset API 需要 backend 后续扩展。
- wrapper finalize 前会拒绝 incomplete 或不连续 manifest；这不等价于真实 VDS dataset validation。
- VDS observable 与 trajectory shard 绑定，当前不会为首个 trajectory frame 前的 mdout frame 单独创建 observable-only shard。

## 本轮增量：particle VDS 物化

本轮将 T6 从纯 manifest/wrapper 控制流推进到核心 particle 数据集的真实 HDF5 VDS 物化：

- `WriterBackend`/`H5MDWriter` 增加 `Create_Virtual_Dataset` 抽象接口。
- `HighFiveBackend` 通过 HDF5 C API (`H5Pset_virtual`) 创建 virtual dataset。
- `VdsTrajectoryH5Writer::Finalize()` 在 manifest 完整性检查通过后，先物化 wrapper 内的 particle VDS，再写 manifest 并 finalize。
- 当前已覆盖 `/particles/all/step`、`/particles/all/time`、`/particles/all/position/value`、`/particles/all/box/edges/value`，以及显式启用时的 `velocity/value` 和 `force/value`。
- `position/step`、`position/time`、`box/edges/step`、`box/edges/time`、`velocity/*`、`force/*` 的时间轴采用 wrapper 内 hard link 复用主 `/particles/all/{step,time}`。

剩余边界：observable/module stream 的跨分片 VDS 尚未物化；当前 wrapper 仍通过 manifest 记录分片范围，后续应按各 observable 的 shape registry 逐项补齐 VDS。

## 本轮增量：ordinary observable VDS 物化

在 particle VDS 的基础上，本轮继续补齐普通 `mdout` observable stream 的跨分片 VDS：

- `VdsShardManifestEntry` 增加 `observable_frame_count`，用于记录每个 shard 内实际写入的普通 observable 帧数。
- `Append_Observable_Frame()` 会同步更新当前 shard 的 observable frame count，因此不再假设 observable 帧数必然等于 trajectory 帧数。
- `Finalize()` 在 particle VDS 之后继续物化 `/observables/all/step`、`/observables/all/time` 和 `/observables/all/<name>/value`。
- `/observables/all/<name>/step` 与 `/observables/all/<name>/time` 通过 hard link 指向 shared observable time axis。
- wrapper 会重新写入 `/parameters/sponge/mdout/columns/{original_name,hdf5_name}`，使 VDS wrapper 本身具备完整的普通 observable column metadata。

剩余边界：module-specific observable streams 尚未跨分片 VDS 物化，包括 NHC、SITS、metadynamics scalar、QC、ReaxFF terms 等。它们当前仍由 shard 文件承载，wrapper 主要依赖 manifest 和已实现的普通 observable/particle VDS。

## 本轮增量：module-specific stream VDS 物化

本轮继续补齐 VDS wrapper 中 module-specific observable streams 的跨分片物化：

- manifest entry 增加 NHC、SITS nk、metadynamics scalar、QC、ReaxFF 的 per-shard frame count。
- NHC：物化 `/observables/all/nose_hoover_chain/{step,time}`、`coordinate/value`、`velocity/value`，并为 coordinate/velocity 建立 step/time hard link。
- SITS：物化 `SITS/<module>/nk/{step,time,value}`，value shape 为 `float32[frame,k_count]`。
- Metadynamics scalar：物化 metadynamics 独立 step/time 轴，以及 `meta/rbias/rct` 三个 scalar value stream。
- QC：物化 QC 独立 step/time 轴、`energy/value`，并在启用时物化 `spin_square/value`。
- ReaxFF：物化 ReaxFF 独立 step/time 轴，以及每个 energy term 的 scalar value stream。

当前 T6 的 VDS wrapper 已覆盖 particle、ordinary observable 和已接入的 module-specific stream。剩余不再是路径覆盖，而是工程可靠性事项：编译验证、真实 HDF5 VDS 文件检查、失败恢复、truncate/repair/resume 策略。

## 本轮增量：VDS source path 可搬迁语义

VDS wrapper 现在区分两类 shard path：

- manifest 中的 shard path 保持原始 `entry.path`，用于 provenance、调试和分片审计。
- HDF5 virtual dataset mapping 使用相对 wrapper 文件所在目录的 source path；若无法计算相对路径，则回退到原始 shard path。

这样在 `*.spg.h5md` wrapper 与其 shard 目录整体移动时，HDF5 VDS source reference 更可能保持有效，不会默认固化为构建时绝对路径。

## 本轮增量：VDS complete-prefix repair hook

VDS writer 增加显式 `Finalize_With_Repair()`，可在 wrapper finalize 时只保留 complete contiguous shard prefix。该接口不改变默认 strict finalize 行为，也不暴露新的 mdin key。它用于后续 resolver/runtime 在显式 repair policy 下重新生成可用 wrapper。
