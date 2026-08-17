# SPONGE H5 Bundle Track T9 报告

## Track

T9 failure/resume/finalize。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T9 细化章节。契约要求 frame 和 shard 只有在完整后才可见；finalization 必须拒绝不一致 manifest，除非显式 repair。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 首次 append 前写 file-level status | 已完成 | `H5MDWriter::Initialize_Common_Layout` 写入 `open` status，并初始化 `/parameters/sponge/output/{frame_count,last_complete_step,last_complete_time}`。 |
| 仅在完整 frame 后推进 `frame_count` | 已完成 | `TrajectoryH5Writer::Append_Particle_Frame`、`ObservableH5Writer::Append_Observable_Frame` 和 `RestartH5Writer::Write_Structural_State` 均在主数据写完后调用 `Write_Output_Completion`。 |
| 存储 last complete step/time | 已完成 | `H5MDWriter::Write_Output_Completion` 写 `/parameters/sponge/output/last_complete_step` 和 `/parameters/sponge/output/last_complete_time`。 |
| 维护 shard status/frame counts | 部分完成 | VDS wrapper 在每个 complete particle frame 后更新总 `frame_count`/last complete；`VdsShardManifestEntry` 维护 shard frame range/status。真实 VDS dataset materialization 仍未完成。 |
| writer-level failed status | 已完成 | `H5MDWriter::Mark_Failed` 写 `status=failed` 和 `/parameters/sponge/output/error`；trajectory、observable-only、restart 的主 frame 写入失败路径会调用该方法。 |
| 检测 incomplete trailing frame | control 层已完成 | `Has_Incomplete_Frame` 暴露 pending incomplete frame 状态。 |
| truncate 或 ignore incomplete data | 推迟 | 需要具体 backend 支持 dataset truncation 或 manifest repair。 |
| 从 complete manifest entries 重新生成 VDS wrapper | 已准备 | `Validate_Complete_Manifest` 返回 complete shard 和 frame count，供 wrapper regeneration 逻辑使用。 |
| 仅在 complete state 后标记 finalized | 已完成 | `Mark_Finalized` 在存在 incomplete frame 时拒绝 finalize；`VdsTrajectoryH5Writer::Finalize` 写 wrapper manifest 前校验 complete manifest。 |
| manifest 不一致时 hard error，除非显式 repair | 部分完成 | `Validate_Complete_Manifest(manifest, allow_repair=false)` 和 VDS writer 内部 manifest validation 对 incomplete 或不连续 shard 失败；repair mode 尚未接入 runtime。 |

## 架构说明

`SPONGE/utils/h5md/completion_tracker.hpp` 引入：

```text
OutputCompletionTracker
ManifestValidationReport
Validate_Complete_Manifest
```

tracker 强制核心状态机：

```text
open -> closing -> finalized
open -> failed
```

一个 frame 必须遵循：

```text
Begin_Frame(step,time) -> all required writes -> Complete_Frame()
```

只有 `Complete_Frame()` 会递增 `frame_count`。

manifest validator 检查：

```text
all visible shards have status == complete
shard indices are contiguous
frame ranges are contiguous
frame_count > 0
```

VDS wrapper finalize 现在也执行同类检查：

```text
Complete_Current_Shard()
  -> Validate_Complete_Manifest()
  -> Write_Manifest_To_Wrapper()
  -> wrapper_writer_->Finalize()
```

如果 manifest 不完整或不连续，wrapper 会通过 `Mark_Failed` 写入 failed status 和
error reason，并停止 finalize。

当前 runtime writer 还直接写入 completion metadata：

```text
/parameters/sponge/output/status
/parameters/sponge/output/frame_count
/parameters/sponge/output/last_complete_step
/parameters/sponge/output/last_complete_time
```

写入规则：

```text
H5MDWriter::Open
  -> status = open
  -> frame_count = 0
  -> last_complete_step = -1
  -> last_complete_time = 0

TrajectoryH5Writer::Append_Particle_Frame
  -> write position/box/optional velocity/optional force
  -> increment particle frame_count metadata

ObservableH5Writer::Append_Observable_Frame
  -> write all scalar observable values
  -> increment observable frame_count metadata

RestartH5Writer::Write_Structural_State
  -> write coordinate/box/velocity and run metadata
  -> set frame_count = 1

writer-level frame write failure
  -> status = failed
  -> /parameters/sponge/output/error = reason
```

这保证 frame count 不会在主 payload 写入前推进。当前 metadata 采用 append-style
record，后续 reader 应读取最后一条作为最新完成状态。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| resume 时 dataset truncation | 需要具体 HDF5 backend API。 |
| 忽略磁盘上的 trailing incomplete shard | 需要 filesystem/backend 集成。 |
| 重新生成真实 VDS wrapper dataset | 需要具体 HDF5 VDS backend。 |
| shard-local completion metadata 与 wrapper VDS 的严格一致性验证 | 需要真实 VDS dataset materialization 后统一处理。 |
| runtime repair mode | 当前只有 validation，尚未实现显式 repair policy。 |
| runtime kill/resume 验证 | 需要 backend 和 output manager 集成。 |

## 审查点

- `OutputCompletionTracker` 可附带 writer，也可不附带 writer。当前实际 writer 已直接写 completion metadata；tracker 仍可作为后续 resume/repair 控制层使用。
- completion metadata 当前追加 scalar metadata record。具体 backend 后续可按需要物化为 replace/update scalar dataset。
- manifest validation 中 `allow_repair=true` 会在第一个 incomplete shard 停止，并把之前 complete 的 shard 视为可用。

## 本轮增量：VDS finalize gating

VDS wrapper finalize 现在在 manifest 完整性检查之后、写出 wrapper manifest 之前执行 particle virtual dataset 物化。若 VDS 物化失败，wrapper 会进入 `failed` 状态并写入失败原因。这使 T9 的完成性语义从“manifest gate only”扩展为“manifest gate + core particle VDS materialization gate”。

## 本轮增量：observable VDS finalize gate

VDS wrapper finalize 的 gating 范围进一步扩展：manifest 必须完整，particle VDS 必须成功物化，ordinary observable VDS 也必须成功物化。任一阶段失败都会将 wrapper 标记为 `failed` 并记录错误原因。

## 本轮增量：module VDS finalize gate

VDS finalize gate 进一步扩展到 module-specific streams。现在 wrapper finalize 必须完成 manifest 校验、particle VDS、ordinary observable VDS 和 module-specific VDS 物化；任一失败都会进入 `failed` 状态。当前尚未实现的是对已失败 wrapper 的自动 repair/truncate/resume。

## 本轮增量：显式 VDS repair finalize API

本轮补入不改变默认行为的显式 repair finalize 路径：

- `VdsTrajectoryH5Writer::Finalize()` 仍保持 strict mode；遇到 incomplete/non-contiguous manifest 会 hard error 并标记 wrapper failed。
- 新增 `VdsTrajectoryH5Writer::Finalize_With_Repair()`，使用 `allow_complete_prefix` 策略。
- repair mode 会丢弃当前生命周期内未能 finalize 的 trailing shard，或将 manifest 截断到 complete contiguous prefix。
- repair 后会重算 wrapper 的 trajectory/observable frame totals，并追加更新 completion metadata。
- wrapper 写入 `/parameters/sponge/output/repair_policy`、`repair_status` 和 `repaired_shard_count`。

当前 repair 边界：

- 不读取已有 wrapper 或 shard 文件，因此不是跨进程 resume。
- 不对磁盘上的 shard dataset 做 HDF5 truncation。
- 不删除被丢弃的 shard 文件，只保证 wrapper manifest/VDS 不引用它们。
- 默认 runtime 仍调用 strict `Finalize()`；repair 需要上层显式选择接入。

## 本轮增量：runtime repair policy 接入

`output_h5_trajectory_repair_policy` 已接入 runtime finalize：

- 默认 `strict` 调用 `VdsTrajectoryH5Writer::Finalize()`。
- 显式 `complete_prefix` 调用 `VdsTrajectoryH5Writer::Finalize_With_Repair()`。
- `complete_prefix` 仅在 `output_h5_trajectory_path` 存在且 `output_h5_trajectory_vds=true` 时合法。
- 该 key 通过 `[output.h5.trajectory] repair_policy` flatten 到现有 parser-visible command 风格，不引入新的 `[state]` 或内部 HDF5 path 配置。
