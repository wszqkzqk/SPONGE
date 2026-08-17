# SPONGE H5 Bundle Track T4 报告

## Track

T4 observable-only H5。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T4 细化章节。契约要求 observable-only H5MD artifact 包含 `/h5md`、`/observables` 和 `/parameters`，但不包含 trajectory particle 字段。本轮继续推进 runtime 接入，使 `output_h5_observable_path` 可以独立跟随 mdout cadence 写出。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 添加 observable-only writer construction | 已完成 | `MD_INFORMATION::trajectory_output::Initial_H5_Observable` 在 `output_h5_observable_path` 启用时创建 `HighFiveBackend` 和 `ObservableH5Writer`。 |
| 复用 observable schema 和 mdout metadata | 已完成 | `Initial_H5_Observable` 读取 `controller->outputs_key`，写入原始列名和 HDF5-safe 列名。 |
| 复用 `/parameters/sponge` log/provenance writer | facade 层已具备 | `ObservableH5Writer` 已提供 `Write_Mdinfo_Text` 和 `Write_Provenance_String`；runtime mdinfo/provenance 还未接入。 |
| 防止 `/particles` 创建 | 已完成 | `ObservableH5Writer::Open` 使用 `WriterOptions::observable_only = true`，T1 common layout 在该模式下跳过 `/particles`。 |
| 支持无 trajectory H5 的 observable-only 输出 | 已完成 | `Initial_H5_Observable` 只依赖 `output_h5_observable_path`，可独立于 `output_h5_trajectory_path` 启用。 |
| observable cadence 独立于 trajectory cadence | 已完成 | `Main_Print` 在 `Check_Mdout_Step()` 分支中调用 `Append_H5_Observable_Only_Frame`，不依赖 `Check_Trajectory_Step()`。 |
| 在 mdout reset 前写入 observable values | 已完成 | `Append_H5_Observable_Only_Frame` 位于 `Print_To_Screen_And_Mdout` 之前，读取尚未 reset 的 `outputs_content`。 |
| 独立 finalization | 已完成 | `Main_Clear` 调用 `Finalize_H5_Observable`，独立 finalize/close observable-only writer。 |

## 架构说明

T4 runtime 接入后的数据流为：

```text
Main_Initial
  -> controller.Print_First_Line_To_Mdout()
  -> md_info.output.Initial_H5_Trajectory()
  -> md_info.output.Initial_H5_Observable()

Main_Print mdout cadence
  -> modules Step_Print(...)
  -> md_info.output.Append_H5_Observable_Frame()       # trajectory H5, if enabled
  -> md_info.output.Append_H5_Observable_Only_Frame()  # observable-only H5, if enabled
  -> controller.Print_To_Screen_And_Mdout()

Main_Clear
  -> md_info.output.Finalize_H5_Trajectory()
  -> md_info.output.Finalize_H5_Observable()
```

writer 结构为：

```text
trajectory_output
  -> HighFiveBackend
  -> ObservableH5Writer
  -> H5MDWriter(observable_only=true)
```

observable-only H5 与 trajectory H5 使用相同的 mdout 数据源和列名 sanitization 逻辑，但持有独立 backend 和文件生命周期。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| mdinfo text runtime 写入 | 当前 `controller.printf` 直接写 `mdinfo` FILE；需要单独设计 log capture 或 finalize-time 汇总。 |
| provenance runtime 写入 | 需要统一 provenance owner 和 launch metadata。 |
| HDF5 inspection 验证 `/particles` 不存在 | 未运行编译/测试/验证。 |
| scalar value comparison against `mdout` | 未运行 runtime 验证。 |

## 审查点

- `output_h5_observable_path` 可独立启用，不要求 `output_h5_trajectory_path`。
- 如果同时启用 trajectory H5 和 observable-only H5，同一批 mdout scalar 会分别写入两个独立 H5MD 文件。
- 缺失或无法转换为 double 的 mdout 值会触发 hard error，避免生成不完整 observable frame。
- T4 没有引入新的 mdin key；仍使用契约中的 `output_h5_observable_path`。
