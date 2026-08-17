# SPONGE H5 Bundle Track T2 报告

## Track

T2 输出 resolver。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T2 细化章节。本 track 的目标是生成 immutable runtime output plan，集中处理 H5 artifact 启用、VDS 默认值、shard-root 派生、`chunk_size` 校验和 legacy sidecar gating。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 添加 resolver 入口 | 已完成 | `SPONGE/utils/h5md/output_plan.hpp` 中的 `SpongeH5OutputPlan::Resolve_Output_Plan(CONTROLLER*, bool)`。 |
| 读取分组 TOML flatten 后的 H5 键 | 已完成 | resolver 消费 `output_h5_trajectory_path`、`output_h5_trajectory_vds`、`output_h5_trajectory_chunk_size`、`output_h5_restart_path` 和 `output_h5_observable_path`。 |
| 应用默认 `trajectory_vds=false` | 已完成 | `TrajectoryH5OutputPlan::vds` 默认使用 `SpongeH5OutputContract::kDefaultTrajectoryVds`。 |
| 应用默认 `trajectory_chunk_size=20` | 已完成 | resolver 使用 `SpongeH5OutputContract::Trajectory_Chunk_Size`。 |
| 对 `chunk_size <= 0` hard error | 已完成 | resolver 将 plan 标记为 invalid，并在 `throw_on_error=true` 时调用 `Throw_SPONGE_Error`。 |
| 内部派生 shard root | 已完成 | `Derive_Shards_Root` 从 `*.spg.h5md` 派生 `.spg.shards`，否则追加 `.shards`。 |
| 不接受 shard dir 用户键 | 已完成 | resolver 不读取任何 user-facing shard directory key。 |
| 计算 legacy 默认 gating | 已完成 | `LegacyOutputPlan::default_enabled = !Any_H5_Output_Enabled(controller)`。 |
| 只在默认启用或显式路径时激活 legacy | 已完成 | `Resolve_Legacy_Output_Plan` 跟踪 `mdout`、`mdinfo`、`crd`、`box`、`vel`、`frc`、`rst` 和 `qc_scf_output` 的 `enabled` 与 `explicit_path`。 |

## 架构说明

resolver 是 header-only 且依赖较轻。它不构造 H5 writer、不打开文件，也不直接改变现有输出行为。后续 track 可使用 `ResolvedOutputPlan` 驱动 trajectory、observable、restart、VDS 和 legacy writer 构造。

plan 对象结构为：

```text
ResolvedOutputPlan
  trajectory: TrajectoryH5OutputPlan
  restart: RestartH5OutputPlan
  observable: ObservableH5OutputPlan
  legacy: LegacyOutputPlan
```

## 审查点

- 单独设置 `output_h5_trajectory_vds` 但未设置 `output_h5_trajectory_path` 时，不启用 trajectory H5 输出。
- 后缀检查仅记录为 `has_recommended_suffix`，不阻止用户路径。
- 真实 runtime 集成留给后续 track。现有输出 call site 尚未切换到 `ResolvedOutputPlan`。
- 现在可以通过 `Resolve_Output_Plan` 做 parser dry-run 覆盖，但本轮未添加测试框架。

## 本轮增量：repair policy resolver

`ResolvedOutputPlan::trajectory` 增加 `repair_policy` 和 `allow_complete_prefix_repair`。resolver 解析 `output_h5_trajectory_repair_policy`，默认 `strict`；显式 `complete_prefix` 会打开 complete-prefix repair finalize，并且要求 VDS trajectory 输出已启用。
