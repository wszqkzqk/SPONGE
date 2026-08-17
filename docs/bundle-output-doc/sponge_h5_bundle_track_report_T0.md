# SPONGE H5 Bundle Track T0 报告

## Track

T0 契约常量。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T0 细化章节。本轮实现遵循 T0 范围：集中定义 parser 可见的 H5 输出键、默认值、推荐后缀和 legacy sidecar 策略，不引入 HDF5 writer 逻辑。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 为 5 个 H5 输出键添加常量 | 已完成 | `SPONGE/utils/control/h5_output_contract.hpp` 定义 `output_h5_trajectory_path`、`output_h5_trajectory_vds`、`output_h5_trajectory_chunk_size`、`output_h5_restart_path` 和 `output_h5_observable_path`。 |
| 添加推荐文件后缀常量 | 已完成 | 同一头文件定义 `.spg.h5md`、`.spgr.h5` 和 `.obs.spg.h5md`。 |
| 添加默认 trajectory chunk size | 已完成 | 同一头文件定义 `kDefaultTrajectoryChunkSize = 20`。 |
| 添加 H5 输出是否启用的辅助函数 | 已完成 | `Any_H5_Output_Enabled(CONTROLLER*)`。 |
| 添加 legacy sidecar 显式请求检测 | 已完成 | `Legacy_Sidecar_Requested(CONTROLLER*, const char*)`。 |
| 添加 legacy sidecar gating 辅助函数 | 已完成 | `Legacy_Sidecars_Default_Enabled` 和 `Legacy_Sidecar_Enabled`。 |
| 添加 parser 可见分组来源说明 | 已完成 | 头文件注释说明 TOML 分组来源和 flatten 后的 parser 可见键。 |
| 添加 parser dry-run 用例 | 推迟 | 本轮未添加 parser 测试框架；应在 T2 resolver 入口稳定后覆盖。 |

## 架构说明

新增头文件由 `SPONGE/control.h` 在 `CONTROLLER` 定义之后包含。这样 helper 层依赖较轻，同时后续 output resolver 和 writer 可以通过现有 controller include 路径使用这些定义。

该 helper 层不打开文件、不写 HDF5、不直接改变 legacy 输出行为，也不验证运行时输出计划。这些工作属于后续 track。

## 审查点

- `Any_H5_Output_Enabled` 目前只把 H5 路径键视为 canonical H5 输出启用条件。单独设置 `output_h5_trajectory_vds` 但未设置 trajectory path 时，不启用 H5 输出。
- `Trajectory_Chunk_Size` 返回解析到的整数或默认值 `20`；对 `<= 0` 的 hard error 校验留给 T2 resolver。
- 推荐后缀 helper 是非阻塞工具函数。T0 不强制路径后缀。

## 本轮增量：trajectory repair policy key

新增 parser-visible H5 output key：

- `output_h5_trajectory_repair_policy`

该 key 来自 `[output.h5.trajectory] repair_policy`，默认值为 `strict`。当前合法值为 `strict` 和 `complete_prefix`；`complete_prefix` 只允许在 `output_h5_trajectory_vds=true` 且设置 trajectory path 时使用。
