# SPONGE H5 Bundle Tracks 总报告

## 范围

本报告汇总 H5 bundle 实现 track T0 到 T9 的当前状态。

每个 track 在开始实现前均读取了 `docs/sponge_h5_bundle_update_plan.md` 中对应的细化章节，并形成了单项报告：

| Track | 报告 |
|---|---|
| T0 契约常量 | `docs/sponge_h5_bundle_track_report_T0.md` |
| T1 HDF5/H5MD writer 基础层 | `docs/sponge_h5_bundle_track_report_T1.md` |
| T1b Concrete HDF5 backend | `docs/sponge_h5_bundle_track_report_T1b.md` |
| T2 输出 resolver | `docs/sponge_h5_bundle_track_report_T2.md` |
| T3 单文件 trajectory H5 | `docs/sponge_h5_bundle_track_report_T3.md` |
| T4 observable-only H5 | `docs/sponge_h5_bundle_track_report_T4.md` |
| T5 restart H5 | `docs/sponge_h5_bundle_track_report_T5.md` |
| T6 VDS trajectory 输出 | `docs/sponge_h5_bundle_track_report_T6.md` |
| T7 module-specific mappings | `docs/sponge_h5_bundle_track_report_T7.md` |
| T8 legacy compatibility | `docs/sponge_h5_bundle_track_report_T8.md` |
| T9 failure/resume/finalize | `docs/sponge_h5_bundle_track_report_T9.md` |

## 当前实现层级

当前实现是与契约对齐的 C++ facade 加 runtime 初步接入层。HighFive backend、trajectory/observable/restart writer 和主要输出 call site 已接入，但尚未经过编译、运行或 HDF5 inspection 验证，因此不能视为已验证的完整可用 HDF5 输出实现。

已完成的基础能力：

| 区域 | 状态 |
|---|---|
| parser 可见 H5 输出键常量 | 已实现。 |
| 默认 `output_h5_trajectory_chunk_size = 20` | 已实现。 |
| 推荐后缀常量 | 已实现。 |
| H5 output plan resolver | 已实现。 |
| legacy sidecar default-off gating | 已在 `Get_Output_File` 路径上实现。 |
| H5MD writer 抽象 backend API | 已实现。 |
| HighFive concrete HDF5 backend | 已实现，依赖由 pixi 提供，尚未编译验证。 |
| trajectory H5 writer facade | 已实现。 |
| non-VDS trajectory H5 runtime 初步接入 | 已实现 position/box/observable/可选 velocity/可选 force。 |
| observable-only H5 writer facade | 已实现。 |
| observable-only H5 runtime 初步接入 | 已实现，跟随 mdout cadence 写 observables；尚未验证。 |
| restart H5 writer facade | 已实现。 |
| structural restart H5 runtime 初步接入 | 已实现 coordinate/velocity/box 单帧覆盖写出；NHC、SITS `nk` 和 metad text snapshots 已接入。 |
| VDS trajectory writer facade | 已实现。 |
| VDS trajectory runtime 初步接入 | 已实现 wrapper/shard/manifest 控制流，并已接入 particle、ordinary observable 和 module-specific VDS materialization；尚未验证。 |
| module-specific H5 mapping facade | 已实现。 |
| NHC module runtime H5 mapping | 已实现 observable trajectory 和 restart state。 |
| SITS nk runtime H5 mapping | 已实现 nk observable trajectory 和 restart state。 |
| Metad scalar runtime H5 mapping | 已实现 `meta/rbias/rct` module-specific observable，使用独立 step/time 轴。 |
| Metad diagnostic H5 mapping | 已实现 `potential_export`、`hills`、`history`、`edge` 同步；`direct_export` 为 copy-if-present。 |
| Metad restart H5 mapping | 已实现 text snapshot 级别的 `hills/history/edge/potential_export/direct_export` 到 `/parameters/restart/bias/meta`。 |
| QC scalar/log H5 mapping | 已实现 `QC`/`QC_S_sq` 到 `/observables/all/qc`，以及显式 `qc_scf_output` 到 `/parameters/sponge/qc/scf_output` 的同步。 |
| ReaxFF energy H5 mapping | 已实现 `REAXFF`/`REAXFF_*` energy terms 到 `/observables/all/reaxff` 的 runtime 同步。 |
| completion/status/manifest tracker | 已实现。 |
| runtime completion metadata | 已实现 `status/frame_count/last_complete_step/last_complete_time` 的初步写入；主 frame 写入失败会标记 `failed`；resume/repair 未实现。 |

尚未实现的能力：

| 区域 | 原因 |
|---|---|
| 已验证的 runtime 落盘 `*.spg.h5md`、`*.spgr.h5` 或 `*.obs.spg.h5md` | backend 和 call site 已初步接入，但尚未编译、运行或 inspection。 |
| 完整 runtime call-site 集成 | metad grid/scatter/hill 结构化二进制 restart schema 等仍未接入。 |
| 真实 VDS dataset creation | 已在 HighFive backend 中初步实现 HDF5 VDS creation；尚未编译和文件级验证。 |
| 显式 legacy sidecar 的 H5 provenance 记录 | 已实现到 `/parameters/sponge/files/legacy_sidecars/{key,path}`；尚未验证。 |
| 与 legacy 输出的 runtime 验证 | 需要 backend 和 call-site 集成。 |

## Track 状态

| Track | 当前状态 | 说明 |
|---|---|---|
| T0 | contract-helper 层完成 | 定义 key、默认值、后缀和 legacy gating helper。 |
| T1 | 部分完成 | backend-independent writer API 已存在。 |
| T1b | 部分完成 | HighFive backend 已实现并接入构建配置，但尚未编译或 HDF5 inspection 验证。 |
| T2 | resolver 层完成 | 生成 immutable output plan 和确定性 legacy gating。 |
| T3 | 部分完成 | 单文件 trajectory facade 和 non-VDS runtime 初步接入已存在；position/box/observable/可选 velocity/可选 force 已接入，尚未编译或文件验证。 |
| T4 | 部分完成 | observable-only facade 和 runtime 初步接入已存在；尚未编译或文件验证。 |
| T5 | 部分完成 | structural restart runtime 已初步接入并采用单文件覆盖；NHC、SITS `nk` 和 metad text snapshots 已写入，输入 reader 尚未接入。 |
| T6 | 部分完成 | VDS runtime wrapper/shard/manifest、finalize 前 manifest validation、particle/observable/module VDS materialization 和 source path relocation 已接入；尚未编译或文件级验证。 |
| T7 | 部分完成 | NHC、SITS nk、metad scalar、主要 metad diagnostic、metad restart text snapshots、QC scalar/SCF log 和 ReaxFF energy terms runtime mapping 已接入；metad 结构化 restart schema 仍待补。 |
| T8 | 部分完成 | generic `Get_Output_File` legacy gating、restart cadence legacy gating 与显式 legacy sidecar provenance 已实现；尚未 runtime 验证。 |
| T9 | 部分完成 | completion tracker、manifest validation、VDS finalize gate、runtime completion metadata、writer-level failed status 和 explicit complete-prefix repair policy 已接入；跨进程 resume、dataset truncation 与 kill/resume 验证尚未实现。 |

## 主要阻塞依赖

功能性 H5 输出的中心阻塞项已从“缺少 concrete backend”推进为“缺少编译、运行和 HDF5 文件 inspection 验证”。

当前状态：

```text
- 已新增 HighFiveBackend，实现 WriterBackend。
- 已在 pixi/CMake 中加入由 pixi 环境提供的 HighFive/HDF5 依赖与链接。
- 主要 runtime call site 已初步接入，尚未运行 configure/compile，也尚未生成 smoke H5MD 文件。
```

因此 HighFive backend 仍需要编译和文件 inspection 确认；T3 non-VDS trajectory、T4 observable-only、T5 structural restart、T6 VDS trajectory，以及 T7 NHC/SITS nk/metad scalar/metad diagnostic/metad restart text snapshots/QC scalar/QC SCF log/ReaxFF energy module mapping 已有 runtime 初步接入；T7 metad 结构化 restart schema 仍需补齐；T9 已有 explicit repair policy，但跨进程 resume 和 kill/resume 验证仍未完成。

## 推荐下一 track

下一步应从 backend smoke 验证和缺口 runtime output manager 接入开始：

```text
T1c HighFive backend smoke validation
T3r/T4r/T5r runtime call-site integration
```

最小交付项：

| 交付项 | 验收标准 |
|---|---|
| 编译 HighFive backend | configure/compile 能通过 HighFive/HDF5 链接。 |
| 实现最小 smoke writer | 可生成并 inspection 一个小型 `*.spg.h5md` 文件。 |
| 接入 trajectory runtime output | 初步接入已完成；仍需编译、运行和 HDF5 inspection 验证。 |
| 接入 observable-only runtime output | 初步接入已完成；仍需编译、运行和 HDF5 inspection 验证。 |
| 接入 restart runtime output | structural restart 初步接入已完成；module state、reader 和验证仍未完成。 |

## 用户审查清单

建议优先检查以下设计决策：

| 决策 | 当前选择 |
|---|---|
| H5 输出键命名空间 | `[output.h5.*]` flatten 为 `output_h5_*`。 |
| VDS chunk size | 默认 `20` 个 trajectory frame。 |
| shard directory | 由 writer 派生，不允许 mdin 配置。 |
| observable-only H5 | writer option 层保证不创建 `/particles`。 |
| restart H5 | 只保留单个 structural state，不作为 trajectory。 |
| legacy sidecar | 启用任意 H5 输出后，legacy sidecar 仅在显式路径下写出。 |
| `eeq_charges.txt` | 继续推迟。 |

## 验证状态

本轮未运行编译、runtime 或 HDF5 文件验证。

当前实现应视为与设计契约对齐的代码 scaffold 加 HighFive backend 初版。下一步需要编译验证 backend，并接入 runtime 输出路径，之后才能进入真实文件验证阶段。

## 本轮增量：HighFive/HDF5 VDS backend

本轮补入了 backend-level VDS API，并在 HighFive 后端使用 HDF5 virtual dataset 能力实现 `spg.h5md` wrapper 的核心 particle 路径物化。当前 VDS 路径不再只是记录分片 manifest；对坐标、box、可选速度、可选力等 H5MD particle 主数据，wrapper 会在 finalize 时创建虚拟数据集并映射到各 shard。

仍未完成的部分是 observable/module stream 的 VDS 物化，以及基于已失败/半完成 wrapper 的 resume/truncate/repair 策略。

## 本轮增量：ordinary observable VDS

VDS wrapper 现在除核心 particle 数据外，也会物化普通 `mdout` observable stream：`/observables/all/step`、`/observables/all/time` 和每个 `/observables/all/<name>/value`。该实现使用 per-shard observable frame count，不强制假设 observable frame 与 trajectory frame 完全一一对应。

仍未完成：module-specific streams 的 VDS 物化，例如 NHC、SITS、metadynamics scalar、QC 和 ReaxFF terms。

## 本轮增量：module-specific VDS

VDS wrapper 现在覆盖 module-specific streams：NHC、SITS nk、metadynamics scalar、QC 和 ReaxFF energy terms。每类模块使用自己的 per-shard frame count 与独立 step/time 轴，避免把不同输出频率的模块强行绑定到 trajectory 或普通 mdout 轴。

剩余风险集中在未编译验证的 C++/HDF5 API 兼容性，以及失败 wrapper 的修复/恢复策略。

## 本轮增量：VDS source path relocation

VDS source path 已调整为优先相对 wrapper 文件路径记录。manifest 仍保存原始 shard path 供审计，HDF5 VDS mapping 则使用相对路径以支持 wrapper 与 shard 目录整体搬迁。

## 本轮增量：显式 repair finalize

新增 VDS writer 级别的 explicit repair finalize hook。默认输出仍 strict；只有调用 `Finalize_With_Repair()` 时才会丢弃 trailing/inconsistent shard suffix，并使 wrapper manifest/VDS 指向 complete contiguous prefix。该实现尚不是完整跨进程 resume，也不执行 HDF5 dataset truncation。

## 本轮增量：repair policy key

新增并接入 `output_h5_trajectory_repair_policy`。默认 strict；显式 `complete_prefix` 时，VDS trajectory finalize 会调用 complete-prefix repair hook。该设置位于 `[output.h5.trajectory] repair_policy`，与现有 H5 output 分组兼容，不新增 shard path 配置。
