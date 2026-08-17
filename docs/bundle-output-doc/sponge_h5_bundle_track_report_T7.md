# SPONGE H5 Bundle Track T7 报告

## Track

T7 module-specific mappings。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T7 细化章节。目标是在不改变 module 语义的前提下，将 NHC、SITS、metadynamics、QC 和 ReaxFF 输出映射到 canonical H5 path。`eeq_charges.txt` 继续推迟。当前已完成 NHC、SITS `nk`、metad scalar/diagnostic、metad text snapshot restart、QC scalar/SCF log 和 ReaxFF energy terms 的 runtime 接入；metad grid/scatter/hill 结构化 restart schema 继续作为后续子项。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| NHC observable path mapping | 已完成 | `ModuleH5MappingWriter::Ensure_Nose_Hoover_Chain_Observables` 定义 `/observables/all/thermostat/nose_hoover_chain/{step,time,coordinate/value,velocity/value}`。 |
| NHC 使用独立 step/time 轴 | 已完成 | NHC 不再复用 mdout 的 `/observables/all/step,time`，避免 NHC trajectory cadence 与 mdout cadence 混轴。 |
| NHC runtime trajectory H5 接入 | 已完成 | `Main_Print` 在 `nhc.Save_Trajectory_File()` 后调用 `Append_H5_Nose_Hoover_Chain_Frame`，写入 trajectory H5 或 VDS shard。 |
| NHC runtime observable-only H5 接入 | 已完成 | `Append_H5_Nose_Hoover_Chain_Frame` 同时写入 observable-only H5 writer。 |
| NHC restart path | 已完成 | `Export_H5_Restart_File` 接收 NHC coordinate/velocity，并调用 `RestartH5Writer::Write_Nose_Hoover_Chain_State` 写 `/parameters/restart/thermostat/nose_hoover_chain`。 |
| SITS nk observable path mapping | 已完成 | `Ensure_Sits_Nk_Observable` 和 `Append_Sits_Nk_Frame` 定义 `/observables/all/sits/<module>/nk/{step,time,value}`，使用独立 step/time 轴。 |
| SITS restart path | 已完成 | H5 restart export 调用 `RestartH5Writer::Write_Sits_State(module, "nk", ...)` 写 `/parameters/restart/bias/sits/<module>/nk`。 |
| SITS nk runtime trajectory/observable-only 接入 | 已完成 | `CLASSIC_SITS_INFORMATION::SITS_Write_Nk_Norm` 设置 pending 标记，main 在同一步调用 `Append_H5_Sits_Nk_Frame` 写入 trajectory H5/VDS 或 observable-only H5。 |
| Metad scalar path mapping | 已完成 | `Ensure_Metadynamics_Scalars` 和 `Append_Metadynamics_Scalar_Frame` 指向 `/observables/all/metadynamics/{step,time,meta/value,rbias/value,rct/value}`。 |
| Metad scalar 使用独立 step/time 轴 | 已完成 | metad scalar 不再复用 mdout 的 `/observables/all/step,time`，避免 metad canonical scalar 与 generic mdout stream 混轴。 |
| Metad scalar runtime trajectory/observable-only 接入 | 已完成 | `Main_Print` 在 `meta.Step_Print(&controller)` 后调用 `Append_H5_Metadynamics_Scalar_Frame`，写入 trajectory H5/VDS 或 observable-only H5。 |
| Metad potential export runtime 接入 | 已完成 | `Main_Print` 在 `meta.Write_Potential()` 后同步 `write_potential_file_name` 到 `/parameters/sponge/metadynamics/<name>/potential_export`。 |
| Metad hills runtime 接入 | 已完成 | 初始化后 copy-if-present，并在 `potential_update_interval` 命中的 `Do_Metadynamics` 步后同步 `myhill.log` 到 `/parameters/sponge/metadynamics/<name>/hills`。 |
| Metad history/edge runtime 接入 | 已完成 | H5 writer 初始化后 copy-if-present，将 `history.log` 和 `edge_file_name` 同步到 `/parameters/sponge/metadynamics/<name>/{history,edge}`。 |
| Metad direct export | 部分完成 | H5 writer 初始化后 copy-if-present `write_directly_file_name` 到 `/parameters/sponge/metadynamics/<name>/direct_export`；当前源码未发现 `Write_Directly()` 运行时调用点，因此不新增强制写出行为。 |
| Metad restart text snapshots | 部分完成 | H5 restart export 将已存在的 `myhill.log`、`history.log`、`edge_file_name`、`write_potential_file_name` 和 `write_directly_file_name` 写入 `/parameters/restart/bias/meta/<name>/...`。 |
| QC scalar/log mapping | 已完成 | `Ensure_Qc_Observables`、`Append_Qc_Frame` 和 `Write_Qc_Scf_Output` 指向 `/observables/all/qc` 和 `/parameters/sponge/qc/scf_output`。 |
| QC scalar 使用独立 step/time 轴 | 已完成 | QC canonical scalar 使用 `/observables/all/qc/{step,time}`，不复用 generic mdout stream 的 `/observables/all/step,time`。 |
| QC scalar runtime 接入 | 已完成 | `Main_Print` 在 `qc.Step_Print(&controller)` 后调用 `Append_H5_Qc_Frame`，将 `QC` 和可选 `QC_S_sq` mdout 列写入 H5。 |
| QC SCF log runtime 接入 | 已完成 | 显式 `qc_scf_output` 文件在 `Solve_SCF` 后 flush，并同步到 `/parameters/sponge/qc/scf_output`。 |
| ReaxFF energy mapping | 已完成 | `Ensure_Reaxff_Energy_Terms` 和 `Append_Reaxff_Frame` 指向 `/observables/all/reaxff/<term>/value`。 |
| ReaxFF energy runtime 接入 | 已完成 | `Main_Print` 在 `reaxff.Step_Print(&controller, md_info.d_charge)` 后调用 `Append_H5_Reaxff_Frame`，将 `REAXFF` 和 `REAXFF_*` mdout 列写入 H5。 |
| `eeq_charges.txt` deferred | 已完成 | 未为 `eeq_charges.txt` 添加 API 或输入契约。 |

## 架构说明

NHC runtime 接入路径为：

```text
Main_Initial
  -> md_info.output.Initial_H5_Nose_Hoover_Chain(chain_length)
      -> TrajectoryH5Writer / VdsTrajectoryH5Writer / ObservableH5Writer
      -> ModuleH5MappingWriter::Ensure_Nose_Hoover_Chain_Observables

Main_Print trajectory cadence
  -> nhc.Save_Trajectory_File()  # legacy, if explicit/default enabled
  -> md_info.output.Append_H5_Nose_Hoover_Chain_Frame(...)
      -> /observables/all/thermostat/nose_hoover_chain/coordinate/value
      -> /observables/all/thermostat/nose_hoover_chain/velocity/value

Main_Print restart cadence
  -> md_info.output.Export_H5_Restart_File(... nhc state ...)
      -> /parameters/restart/thermostat/nose_hoover_chain
```

NHC 特意使用独立路径：

```text
/observables/all/thermostat/nose_hoover_chain/step
/observables/all/thermostat/nose_hoover_chain/time
/observables/all/thermostat/nose_hoover_chain/coordinate/value
/observables/all/thermostat/nose_hoover_chain/velocity/value
```

原因是 NHC legacy trajectory 跟随 trajectory cadence，而 mdout scalar 跟随 `write_mdout_interval`。如果共用 `/observables/all/step,time`，会破坏不同 observable 的 frame 对齐语义。

SITS runtime 接入路径为：

```text
SITS_Update_Nk
  -> SITS_Write_Nk_Norm()
      -> copy Nk to nk_record_cpu
      -> mark h5_nk_pending

Main_Calculate_Force
  -> if h5_nk_pending:
       md_info.output.Append_H5_Sits_Nk_Frame(...)
       h5_nk_pending = 0

Main_Print restart cadence
  -> Export_H5_Restart_File(... sits nk ...)
      -> /parameters/restart/bias/sits/<module>/nk
```

SITS `nk` 也使用独立路径：

```text
/observables/all/sits/<module>/nk/step
/observables/all/sits/<module>/nk/time
/observables/all/sits/<module>/nk/value
```

Metad scalar runtime 接入路径为：

```text
Main_Print mdout cadence
  -> meta.Step_Print(&controller)
      -> synchronizes potential_local/rbias/rct onto rank 0
      -> writes legacy mdout columns
  -> md_info.output.Append_H5_Metadynamics_Scalar_Frame(...)
      -> /observables/all/metadynamics/step
      -> /observables/all/metadynamics/time
      -> /observables/all/metadynamics/meta/value
      -> /observables/all/metadynamics/rbias/value
      -> /observables/all/metadynamics/rct/value
```

Metad scalar 也使用独立路径：

```text
/observables/all/metadynamics/step
/observables/all/metadynamics/time
/observables/all/metadynamics/meta/value
/observables/all/metadynamics/rbias/value
/observables/all/metadynamics/rct/value
```

原因是 generic mdout H5 stream 仍会按 `controller->outputs_key` 记录所有
mdout 列，而 metad canonical scalar 是 module-specific observable。二者不能
共享同一个 append axis，否则同一步会向 `/observables/all/step,time` 追加两
次不同语义的 frame。

Metad diagnostic runtime 接入路径为：

```text
Main_Initial
  -> md_info.output.Initial_H5_Metadynamics(...)
  -> copy-if-present myhill.log/history.log/edge_file_name/direct_export
       into /parameters/sponge/metadynamics/<name>/

Main_Calculate_Force
  -> meta.Do_Metadynamics(...)
  -> if potential_update_interval is hit:
       copy myhill.log into /parameters/sponge/metadynamics/<name>/hills

Main_Print trajectory cadence
  -> meta.Write_Potential()
  -> copy write_potential_file_name
       into /parameters/sponge/metadynamics/<name>/potential_export
```

`Meta_directly.txt` 当前仅实现 copy-if-present，因为源码中 `META::Write_Directly`
没有运行时调用点。该设计避免为了 H5 bundle 引入新的 legacy sidecar 生成行为。

Metad restart text snapshot 接入路径为：

```text
Main_Print restart cadence
  -> md_info.output.Export_H5_Restart_File(... metad file names ...)
      -> copy-if-present myhill.log
          -> /parameters/restart/bias/meta/<name>/hills
      -> copy-if-present history.log
          -> /parameters/restart/bias/meta/<name>/history
      -> copy-if-present edge_file_name
          -> /parameters/restart/bias/meta/<name>/edge
      -> copy-if-present write_potential_file_name
          -> /parameters/restart/bias/meta/<name>/potential_export
      -> copy-if-present write_directly_file_name
          -> /parameters/restart/bias/meta/<name>/direct_export
```

该接入只处理已存在的 path-dependent text snapshots。它不把 CV 定义、grid 定义、
scatter coordinates 或 protocol-owned 参数复制进 restart，也不在 restart cadence
强制生成新的 legacy text sidecar。结构化二进制 metad restart schema 仍需后续设计。

QC SCF log runtime 接入路径为：

```text
Main_Calculate_Force
  -> qc.Solve_SCF(...)
      -> legacy qc_scf_output file receives SCF iteration diagnostics if explicit
  -> flush qc.scf_output_file
  -> md_info.output.Write_H5_Qc_Scf_Output_File(...)
      -> /parameters/sponge/qc/scf_output
```

该路径只同步显式 `qc_scf_output` 文件。当前 QC 代码在未设置
`qc_scf_output` 时将 SCF iteration diagnostics 写到 stdout；本轮不改变该行为，
也不为 H5 bundle 隐式创建 legacy sidecar。

QC scalar runtime 接入路径为：

```text
Main_Print mdout cadence
  -> qc.Step_Print(&controller)
      -> writes QC
      -> writes QC_S_sq when unrestricted
  -> md_info.output.Append_H5_Qc_Frame(&controller)
      -> parse controller->outputs_content["QC"]
      -> parse controller->outputs_content["QC_S_sq"] when registered
      -> /observables/all/qc/step
      -> /observables/all/qc/time
      -> /observables/all/qc/energy/value
      -> /observables/all/qc/spin_square/value
```

QC scalar 特意使用独立路径：

```text
/observables/all/qc/step
/observables/all/qc/time
/observables/all/qc/energy/value
/observables/all/qc/spin_square/value
```

原因是 generic mdout H5 stream 已经会记录 `QC` 和 `QC_S_sq` 列；module-specific
QC scalar 作为 canonical QC observable 不能再次向 `/observables/all/step,time`
追加 frame，否则会破坏 generic mdout axis。实际是否写入 `spin_square` 由本次
run 是否注册 `QC_S_sq` 决定。

ReaxFF energy runtime 接入路径为：

```text
Main_Print mdout cadence
  -> reaxff.Step_Print(&controller, md_info.d_charge)
      -> ReaxFF 子模块向 mdout 注册并更新 REAXFF/REAXFF_* energy columns
  -> md_info.output.Append_H5_Reaxff_Frame(&controller)
      -> scan controller->outputs_key 中的 REAXFF 和 REAXFF_* key
      -> parse controller->outputs_content
      -> /observables/all/reaxff/step
      -> /observables/all/reaxff/time
      -> /observables/all/reaxff/<term>/value
```

ReaxFF 使用独立路径：

```text
/observables/all/reaxff/step
/observables/all/reaxff/time
/observables/all/reaxff/REAXFF/value
/observables/all/reaxff/REAXFF_BOND/value
/observables/all/reaxff/REAXFF_VDW/value
/observables/all/reaxff/REAXFF_EEQ/value
/observables/all/reaxff/REAXFF_ELP/value
/observables/all/reaxff/REAXFF_OVUN/value
/observables/all/reaxff/REAXFF_ANG/value
/observables/all/reaxff/REAXFF_PEN/value
/observables/all/reaxff/REAXFF_COA/value
/observables/all/reaxff/REAXFF_TOR/value
/observables/all/reaxff/REAXFF_CONJ/value
/observables/all/reaxff/REAXFF_HB/value
```

实际写出的 term 集合由本次 run 中 `controller->outputs_key` 已注册的
ReaxFF mdout 列决定。未初始化的 ReaxFF 子模块不会被强制创建空 dataset。
`eeq_charges.txt` 仍按前述约定不纳入本轮契约。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| Metad grid/scatter/hill 结构化 restart schema | 当前仅写已有 text snapshots；完整二进制 continuation schema 仍需确认最小集合和 reader 行为。 |
| HDF5 inspection 和 module 小 case 验证 | 未运行编译/测试/验证。 |

## 审查点

- NHC 和 SITS `nk` 已形成两个 runtime module mapping 样板：均覆盖 observable trajectory 和 restart state。
- NHC observable 当前会写入启用的 trajectory H5/VDS shard，也会写入 observable-only H5；二者是独立 artifact。
- NHC legacy trajectory/restart 仍由原有文件路径控制；H5 mapping 不改变 NHC 动力学或 legacy 数据源。
- Metad scalar、主要 diagnostic exports 和 restart text snapshots 已 runtime 接入；metad 结构化 restart schema 仍待补。
- `Meta_directly.txt` 当前没有现有运行时调用点，因此只做 copy-if-present。
- QC scalar 和 SCF log 已 runtime 接入；QC scalar 使用独立 `/observables/all/qc` axis。
- ReaxFF energy terms 已 runtime 接入；实际 term 集合跟随 `REAXFF`/`REAXFF_*` mdout 列注册。
