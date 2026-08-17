# SPONGE H5 Bundle 当前轮次审查报告

日期：2026-07-03

## 1. 本轮目标

本轮围绕 H5 bundle 设计与实现 tracklist 逐步推进，要求每个 track 开始前读取细化文档，并在完成当前轮次后形成可审查报告。

本报告覆盖当前已经落入工作树的实现与文档状态。未运行编译、单元测试、runtime smoke test 或 HDF5 文件 inspection，因此本报告不声明实现已验证通过。

## 2. 已形成的 track 报告

| Track | 文件 | 当前状态 |
|---|---|---|
| T0 契约常量 | `docs/sponge_h5_bundle_track_report_T0.md` | 已实现 H5 output key、默认值、后缀、legacy gating helper，并补入 repair policy key。 |
| T1 HDF5/H5MD writer 基础层 | `docs/sponge_h5_bundle_track_report_T1.md` | 已实现 writer facade 和 backend 抽象。 |
| T1b HighFive backend | `docs/sponge_h5_bundle_track_report_T1b.md` | 已接入 HighFive/HDF5 backend；未编译验证。 |
| T2 输出 resolver | `docs/sponge_h5_bundle_track_report_T2.md` | 已实现 H5 output plan resolver、legacy gating 和 repair policy resolver。 |
| T3 单文件 trajectory H5 | `docs/sponge_h5_bundle_track_report_T3.md` | 已初步接入 trajectory H5 runtime 输出；未验证。 |
| T4 observable-only H5 | `docs/sponge_h5_bundle_track_report_T4.md` | 已初步接入 observable-only H5MD；未验证。 |
| T5 restart H5 | `docs/sponge_h5_bundle_track_report_T5.md` | 已初步接入 structural restart、NHC/SITS/metad text state；reader 未接入。 |
| T6 VDS trajectory | `docs/sponge_h5_bundle_track_report_T6.md` | 已实现 wrapper/shard/manifest、particle/observable/module VDS 物化、source path relocation。 |
| T7 module mappings | `docs/sponge_h5_bundle_track_report_T7.md` | 已覆盖 NHC、SITS nk、metad scalar/diagnostics、QC、ReaxFF 等主要输出映射。 |
| T8 legacy compatibility | `docs/sponge_h5_bundle_track_report_T8.md` | 已实现 H5 输出启用时 legacy sidecar 默认关闭、显式 sidecar provenance。 |
| T9 failure/resume/finalize | `docs/sponge_h5_bundle_track_report_T9.md` | 已实现 completion metadata、strict finalize、explicit complete-prefix repair policy；跨进程 resume 未实现。 |

## 3. 关键代码落点

| 区域 | 主要文件 |
|---|---|
| H5 output key contract | `SPONGE/utils/control/h5_output_contract.hpp` |
| H5 output resolver | `SPONGE/utils/h5md/output_plan.hpp` |
| Writer abstraction | `SPONGE/utils/h5md/h5md_writer.hpp` |
| HighFive backend | `SPONGE/utils/h5md/highfive_backend.hpp` |
| 单文件 trajectory writer | `SPONGE/utils/h5md/trajectory_h5_writer.hpp` |
| observable-only writer | `SPONGE/utils/h5md/observable_h5_writer.hpp` |
| restart writer | `SPONGE/utils/h5md/restart_h5_writer.hpp` |
| VDS trajectory writer | `SPONGE/utils/h5md/vds_trajectory_h5_writer.hpp` |
| module-specific mapping | `SPONGE/utils/h5md/module_h5_mappings.hpp` |
| completion tracker | `SPONGE/utils/h5md/completion_tracker.hpp` |
| runtime output 接入 | `SPONGE/MD_core/output.hpp` |
| include 集成 | `SPONGE/MD_core/MD_core.h` |
| build dependency | `pixi.toml`, `cmake/targets/SPONGE.cmake` |

## 4. 当前实现能力

### 4.1 Parser-visible H5 output keys

已接入的 flatten key：

```text
output_h5_trajectory_path
output_h5_trajectory_vds
output_h5_trajectory_chunk_size
output_h5_trajectory_repair_policy
output_h5_restart_path
output_h5_observable_path
```

对应 grouped TOML：

```toml
[output.h5.trajectory]
path = "prod.spg.h5md"
vds = true
chunk_size = 20
repair_policy = "strict"

[output.h5.restart]
path = "prod.spgr.h5"

[output.h5.observable]
path = "prod.obs.spg.h5md"
```

### 4.2 Output artifact suffix

| Artifact | Suffix |
|---|---|
| trajectory H5MD | `*.spg.h5md` |
| restart H5 | `*.spgr.h5` |
| observable-only H5MD | `*.obs.spg.h5md` |

### 4.3 Single-file trajectory

已初步接入：

- `/h5md`
- `/particles/all/{step,time}`
- `/particles/all/position/value`
- `/particles/all/box/edges/value`
- 可选 `/particles/all/velocity/value`
- 可选 `/particles/all/force/value`
- `/observables/all/{step,time}`
- `/observables/all/<name>/value`
- `/parameters/sponge/*`

### 4.4 Observable-only H5MD

已初步接入 observable-only writer，设计上不创建 `/particles`，用于轻量级分析输出。

### 4.5 Restart H5

restart 保持单状态语义，不作为 trajectory history。已覆盖：

- structural coordinates/velocities/box
- step/time metadata
- NHC state
- SITS `nk`
- metadynamics text snapshots
- SPONGE restart extension metadata under parameters

### 4.6 VDS trajectory

VDS 模式已经从单纯 manifest wrapper 推进到真实 HDF5 VDS 物化：

- writer-derived shard root，不允许 mdin 配置 shard directory。
- shard 文件使用 `segment_%06d.spg.h5md`。
- wrapper manifest 记录 shard index/path/frame ranges/status。
- HDF5 VDS source path 优先使用相对 wrapper 路径，支持 wrapper + shard 目录整体搬迁。
- 已物化 particle VDS。
- 已物化 ordinary observable VDS。
- 已物化 module-specific stream VDS。

### 4.7 Failure/finalize/repair

默认 finalize 行为保持 strict：

- incomplete shard hard error。
- non-contiguous manifest hard error。
- VDS materialization 失败会标记 wrapper failed。

显式 repair policy：

```toml
[output.h5.trajectory]
repair_policy = "complete_prefix"
```

该模式只在 VDS trajectory 输出下合法。它会保留 complete contiguous shard prefix，并记录：

```text
/parameters/sponge/output/repair_policy
/parameters/sponge/output/repair_status
/parameters/sponge/output/repaired_shard_count
```

边界：不做跨进程 resume、不读取旧 wrapper、不做 HDF5 dataset truncation、不删除 orphan shard。

## 5. 仍未验证或未完成的边界

| 项目 | 状态 |
|---|---|
| 编译验证 | 未运行。 |
| HighFive/HDF5 API 兼容性 | 未验证。 |
| HDF5 文件 layout inspection | 未运行。 |
| HDF5 VDS 可读性检查 | 未运行。 |
| runtime smoke simulation | 未运行。 |
| restart H5 reader | 未实现。 |
| 跨进程 resume | 未实现。 |
| HDF5 dataset truncation | 未实现。 |
| orphan shard 清理 | 未实现。 |
| metadynamics 结构化二进制 restart schema | 未完全实现。 |
| `eeq_charges.txt` H5 化 | 暂缓。 |

## 6. 当前轮次结论

当前轮次已经完成 H5 bundle 的主要设计到代码 scaffold/runtime 初步接入：

- output key contract 已落地。
- grouped TOML flatten 兼容路径已落地。
- pixi/CMake 已指向 HighFive/HDF5。
- HighFive backend 已实现。
- single-file trajectory、observable-only、restart、VDS trajectory writer 已实现。
- module-specific 输出映射已覆盖主要 runtime 输出。
- VDS wrapper 已具备 particle/observable/module stream 物化能力。
- strict finalize 与 explicit repair policy 已接入 runtime。
- 相关契约文档和 track report 已更新。

本轮不能声明“功能已验证可用”。下一阶段应转入系统化单元测试和 smoke/inspection 测试，先证明新路径的编译、writer 行为、HDF5 layout、VDS mapping 和 failure policy。

## 7. 下一阶段建议

建议将下一个 goal 定为：

```text
对于 h5 bundle 相关的新路径，编写详细且覆盖全面的单元测试。
```

优先测试范围：

| 测试组 | 覆盖目标 |
|---|---|
| contract/resolver tests | key flatten、默认值、非法值、suffix、legacy gating、repair policy。 |
| backend tests | HighFive dataset create/append/string/hardlink/VDS/error path。 |
| writer tests | trajectory、observable-only、restart 的 layout 和 completion metadata。 |
| VDS tests | shard rotation、manifest、relative source path、particle/observable/module VDS。 |
| failure tests | strict failure、complete-prefix repair、failed status/error metadata。 |
| runtime-adjacent tests | output manager 在 mock controller 下选择正确 writer 和 finalize policy。 |
