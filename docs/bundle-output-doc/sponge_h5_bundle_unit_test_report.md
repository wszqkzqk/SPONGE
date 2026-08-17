# SPONGE H5 Bundle 单元测试覆盖报告

日期：2026-07-03

## 1. 本轮目标

为 H5 bundle 相关的新路径编写详细且覆盖全面的单元测试。当前提交的是测试代码和 CMake 入口；本轮未运行 configure、build、ctest 或 HDF5 文件级验证。

## 2. 测试构建入口

新增 opt-in CMake 开关：

```text
SPONGE_BUILD_TESTS=ON
```

新增测试目录：

```text
tests/
tests/h5_bundle/
```

新增 CMake 文件：

```text
tests/CMakeLists.txt
tests/h5_bundle/CMakeLists.txt
```

测试目标默认不参与普通 SPONGE 构建，避免影响现有生产构建路径。

建议后续验证命令：

```bash
pixi run -e dev-cpu cmake -S . -B build-h5-tests -DSPONGE_BUILD_TESTS=ON
pixi run -e dev-cpu cmake --build build-h5-tests --target test_h5_output_plan
pixi run -e dev-cpu ctest --test-dir build-h5-tests --output-on-failure
```

具体环境名应按当前 pixi 环境实际存在情况调整。

## 3. 公共测试基础设施

新增：

```text
tests/h5_bundle/h5_bundle_test_common.hpp
```

提供内容：

| 组件 | 用途 |
|---|---|
| lightweight `CONTROLLER` mock | 允许直接测试 `h5_output_contract.hpp` 和 `output_plan.hpp` 中依赖 controller command 的逻辑。 |
| `MockBackend` | 覆盖 `WriterBackend` API，记录 group/dataset/append/string/hardlink/VDS 操作。 |
| `MockBackendFactory` | 供 `VdsTrajectoryH5Writer` 测试 wrapper/shard backend 创建顺序和失败注入。 |
| `REQUIRE_TRUE`/`REQUIRE_EQ` | 轻量断言，不引入 GTest/Catch/doctest。 |
| `Unique_Temp_Path` | HighFive 文件级测试的临时路径。 |

## 4. 测试目标与覆盖面

### 4.1 `test_h5_output_plan`

文件：

```text
tests/h5_bundle/test_h5_output_plan.cpp
```

覆盖：

| 场景 | 覆盖点 |
|---|---|
| 默认设置 | 无 H5 path 时不启用 H5 output；trajectory VDS 默认 false；chunk size 默认 20；repair policy 默认 strict。 |
| legacy gating | 未启用 H5 output 时 legacy sidecar 默认开启；启用 H5 output 后 legacy 默认关闭但显式 `mdout` 保留。 |
| VDS trajectory plan | `output_h5_trajectory_path`、`output_h5_trajectory_vds`、`output_h5_trajectory_chunk_size`、`output_h5_trajectory_repair_policy` 解析。 |
| shard root derivation | `prod.spg.h5md` 派生为 `prod.spg.shards`。 |
| invalid chunk size | `output_h5_trajectory_chunk_size <= 0` 产生 invalid plan。 |
| invalid repair policy | 非 `strict`/`complete_prefix` 值产生 invalid plan。 |
| repair policy constraint | `complete_prefix` 必须同时启用 trajectory path 和 VDS。 |
| suffix helper | trajectory suffix 推荐值检查。 |

### 4.2 `test_h5md_writers_with_mock_backend`

文件：

```text
tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp
```

覆盖：

| Writer | 覆盖点 |
|---|---|
| `H5MDWriter` | common layout、observable-only 不创建 `/particles`、初始 completion metadata。 |
| `TrajectoryH5Writer` | position/box/velocity/force dataset、ordinary observable stream、hard link、completion metadata、finalize。 |
| `ObservableH5Writer` | observable-only stream、无 `/particles`、frame completion metadata。 |
| `RestartH5Writer` | structural restart 单状态写入、position/velocity/box/run metadata、第二次写 structural state 失败并标记 failed。 |

### 4.3 `test_vds_trajectory_writer_with_mock_backend`

文件：

```text
tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp
```

覆盖：

| 场景 | 覆盖点 |
|---|---|
| shard rotation | `chunk_size=2` 时 3 个 trajectory frame 产生 2 个 manifest entries。 |
| particle VDS | wrapper 生成 position/velocity VDS。 |
| ordinary observable VDS | wrapper 生成 `/observables/all/<name>/value` VDS。 |
| module-specific VDS | NHC、SITS nk、metad scalar、QC、ReaxFF term VDS 均被记录。 |
| relative source path | VDS source path 使用相对 wrapper 目录的 shard path。 |
| strict finalize failure | shard finalize 失败时 strict `Finalize()` 返回 false。 |
| complete-prefix repair | `Finalize_With_Repair()` 丢弃失败 shard，manifest 截断为空，并写 repair metadata。 |

### 4.4 `test_highfive_backend_io`

文件：

```text
tests/h5_bundle/test_highfive_backend_io.cpp
```

覆盖：

| 场景 | 覆盖点 |
|---|---|
| basic HDF5 write/read | `HighFiveBackend` 创建文件、group、dataset、append、hard link、finalize，再用 HighFive 读回 layout 和 shape。 |
| VDS file-level check | source H5 文件写入数据，wrapper 通过 `Create_Virtual_Dataset` 建立 VDS，再用 HighFive 读回 VDS shape 和 values。 |

该测试是最接近真实 HDF5 行为的单元/集成边界测试，后续编译运行时最可能暴露 HighFive/HDF5 API 兼容问题。

## 5. 当前未覆盖或后续建议

| 缺口 | 建议 |
|---|---|
| 完整 MD runtime 输出 | 后续增加小系统 smoke test，而不是放在纯单元测试内。 |
| restart H5 reader | reader 实现后增加 round-trip test。 |
| metad structured binary restart | schema 稳定后增加 shape/value test。 |
| orphan shard 清理 | 若实现文件清理策略，再增加 filesystem side-effect test。 |
| kill/resume 跨进程 | 需要专门 runtime/integration test，不建议混入当前 header-level unit tests。 |
| CUDA/MPI 多 rank 行为 | 需要 runtime test 或 mock MPI rank test；当前单元测试只覆盖 rank-0 writer 组件。 |

## 6. 验证状态

本轮未运行：

```text
cmake configure
cmake build
ctest
HDF5 inspection
SPONGE runtime smoke test
```

因此当前状态是：测试代码已写入，覆盖目标已设计，但尚未证明可编译或通过。

## 7. 本轮增量：失败语义与 VDS 细节断言

本轮继续增强 mock-backend 单元测试覆盖：

### 7.1 Writer failure path

`test_h5md_writers_with_mock_backend.cpp` 新增：

| 测试 | 覆盖点 |
|---|---|
| `Test_Trajectory_Append_Failure_Marks_Failed` | particle frame append 阶段 backend 失败时，trajectory writer 应将 writer status 标记为 `failed`，并写入 `/parameters/sponge/output/error`。 |
| `Test_Observable_Missing_Value_Marks_Failed` | observable-only writer 在定义多个 column 但 append 缺失某个 value 时，应 hard fail，并写入明确错误原因。 |

这两项覆盖了 H5 bundle failure contract 的核心路径：主 payload 未完整写入时不得推进为成功输出，而应进入 failed 状态。

### 7.2 VDS manifest/module counters

`test_vds_trajectory_writer_with_mock_backend.cpp` 增加对 manifest entry 的细粒度断言：

| 字段 | 覆盖点 |
|---|---|
| `observable_frame_count` | ordinary observable frame count 不强行等于全局 trajectory assumption，而是随 shard 记录。 |
| `nhc_frame_count` | NHC module stream per-shard frame count。 |
| `sits_nk_frame_count` | SITS nk module stream per-shard frame count。 |
| `metadynamics_scalar_frame_count` | metad scalar stream per-shard frame count。 |
| `qc_frame_count` | QC stream per-shard frame count。 |
| `reaxff_frame_count` | ReaxFF term stream per-shard frame count。 |

### 7.3 VDS source shape/start

同一测试进一步断言：

- particle VDS source dims 包含 `[frame, atom_count, 3]`。
- ordinary observable VDS source dims/start 按 per-shard observable frame count 拼接。
- NHC vector stream source dims 包含 chain length。
- SITS nk vector stream source dims 包含 k count。
- 第二个 shard 的 virtual start 从第一个 shard 的 frame count 后开始。

这些断言用于保护 VDS wrapper 的核心 invariant：wrapper dataset shape 和每个 shard 的 source hyperslab 必须由 manifest 中的实际 completed frame count 推导，而不是由固定输出频率或 trajectory frame 假设隐式推断。

## 8. 本轮增量：completion tracker 与 module mapping 单元测试

本轮新增两个测试目标，并接入 `tests/h5_bundle/CMakeLists.txt`。

### 8.1 `test_completion_tracker`

文件：

```text
tests/h5_bundle/test_completion_tracker.cpp
```

覆盖：

| 场景 | 覆盖点 |
|---|---|
| tracker 无 writer 状态机 | `open -> begin frame -> complete frame -> closing -> finalized` 状态转移。 |
| incomplete frame hard error | 已有 open frame 时不能再次 begin；有 incomplete frame 时不能 finalize。 |
| tracker 写 metadata | 绑定 `H5MDWriter` 时写 `/parameters/sponge/output/{frame_count,last_complete_step,last_complete_time}`。 |
| failed status | `Mark_Failed` 写 `failed` status 与 error reason。 |
| manifest strict validation | complete manifest 通过；incomplete shard、非连续 index、非连续 frame range、非正 frame_count 均失败并返回具体 error。 |
| manifest repair-prefix validation | `allow_repair=true` 时遇到 trailing incomplete shard 停止，并报告之前 complete prefix 的 shard/frame count。 |

### 8.2 `test_module_h5_mappings_with_mock_backend`

文件：

```text
tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp
```

覆盖：

| 模块 | 覆盖点 |
|---|---|
| NHC | `coordinate/value`、`velocity/value`、step/time hard link、append count。 |
| SITS nk | `SITS/<module>/nk/value` dataset 与 append count。 |
| Metadynamics scalar | `meta/rbias/rct` scalar datasets 和 append。 |
| Metadynamics diagnostics | `hills/history/edge/potential_export/direct_export` 写入 `/parameters/sponge/bias/meta/<name>/...`。 |
| QC | `energy`、可选 `spin_square`、`scf_output`。 |
| ReaxFF | 多 energy term dataset、正常 append、缺失 term 时返回 `missing ReaxFF term: <term>`。 |

这两个测试补齐了 H5 bundle 中与 runtime module 输出最相关的纯单元测试边界，使 T7/T9 的核心路径不只依赖 VDS writer 间接覆盖。

## 9. 本轮增量：contract 边界与 HighFive invalid-path 测试

本轮继续扩展两个既有测试目标。

### 9.1 `test_h5_output_plan` 增量

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Output_Selectors_Do_Not_Enable_H5_Without_Path` | 单独设置 `output_h5_trajectory_vds` 或 `output_h5_trajectory_chunk_size` 不应启用 H5 trajectory output；H5 output 仍由 path key 选择。 |
| `Test_Repair_Policy_Alias` | `allow_complete_prefix` 作为 `complete_prefix` 的兼容别名可被 resolver 接受。 |
| `Test_Restart_And_Observable_Paths` | restart/observable path 启用 H5 output，推荐 suffix 识别正确，显式 `qc_scf_output` legacy sidecar 在 H5 output 下仍保留。 |
| `Test_Suffix_And_Shard_Derivation` 扩展 | 覆盖 trajectory/restart/observable 三类推荐 suffix 和 `Recommended_Suffix_For_Key`。 |

### 9.2 `test_highfive_backend_io` 增量

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_HighFive_Backend_Observable_Only_Layout` | `WriterOptions::observable_only=true` 时真实 HDF5 文件包含 `/h5md`、`/observables`、`/parameters`，且不创建 `/particles`。 |
| `Test_HighFive_Backend_Rejects_Invalid_Operations` | 对缺失 dataset append、空 source/string VDS、缺失 hard-link target 等 invalid operation 返回失败。 |

这些测试补强了两个关键不变量：

- H5 bundle 的选择开关不能替代 path binding；只有 canonical path key 才启用对应 H5 artifact。
- HighFive backend 需要在 invalid operation 上显式失败，不能 silent success。

## 10. 本轮增量：legacy provenance、restart module state 与 VDS wrapper diagnostics

本轮继续补充 runtime-adjacent writer API 覆盖，重点是不会自动由主 trajectory frame 测试覆盖的 side-channel 输出。

### 10.1 `test_h5md_writers_with_mock_backend` 增量

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Restart_Module_State_And_Legacy_Provenance` | restart writer 的 NHC state、SITS state、metad text state，以及 legacy sidecar provenance。 |
| `Test_Legacy_Provenance_On_Trajectory_And_Observable` | trajectory writer 与 observable-only writer 的 `/parameters/sponge/files/legacy_sidecars/{key,path}` 写入。 |

这些测试保护 T5 restart state extension 和 T8 explicit legacy sidecar provenance 的交叉路径。

### 10.2 `test_vds_trajectory_writer_with_mock_backend` 增量

在 VDS wrapper 主测试中新增：

| API | 覆盖点 |
|---|---|
| `Write_Metadynamics_Diagnostic` | VDS wrapper 级 `/parameters/sponge/bias/meta/<name>/<component>` 写入。 |
| `Write_Qc_Scf_Output` | VDS wrapper 级 `/parameters/sponge/qc/scf_output` 写入。 |
| `Write_Legacy_Sidecar_Paths` | VDS wrapper 级 legacy sidecar provenance。 |

这些断言确保 VDS mode 不只覆盖 shard payload 和 virtual datasets，也保留 wrapper 上的 diagnostic/log/provenance 信息。

## 11. 本轮增量：测试目录 README 与覆盖矩阵

新增测试目录入口文档：

```text
tests/h5_bundle/README.md
```

该 README 记录：

- 如何启用 `SPONGE_BUILD_TESTS=ON`。
- 六个 H5 bundle 测试目标的职责边界。
- contract area 到具体测试目标的覆盖矩阵。
- 当前单元测试无法证明的边界，例如完整 MD runtime、CUDA/MPI、跨进程 resume、restart reader、orphan shard cleanup。
- 后续维护规则：新增 H5 key、writer path、module mapping、VDS dataset、HighFive primitive 时应补到哪个测试目标。

该文档用于后续编译修正和审查时快速判断测试覆盖是否与 H5 bundle 契约同步。

## 12. 本轮增量：CTest labels 与分层执行入口

本轮为 H5 bundle 测试目标增加 CTest labels：

| Label | 含义 |
|---|---|
| `h5_bundle` | 全部 H5 bundle 测试。 |
| `contract` | parser-visible key、resolver、suffix、legacy gating。 |
| `mock-backend` | 不依赖真实 HDF5 IO 的 writer facade 测试。 |
| `module` | module-specific H5 path mapping。 |
| `vds` | VDS wrapper/shard/virtual dataset mapping。 |
| `backend-io` | 真实 HighFive/HDF5 文件级测试。 |
| `failure` | failed status、repair、invalid operation 等失败语义。 |

`tests/h5_bundle/README.md` 已补充分层执行命令。推荐后续首次验证时按 `contract -> mock-backend -> module -> vds -> backend-io` 顺序运行，以便快速区分纯逻辑错误和 HighFive/HDF5 backend 错误。

## 13. 本轮增量：resolver 异常路径与 VDS precondition errors

本轮继续补充纯逻辑和 writer precondition 测试。

### 13.1 `test_h5_output_plan` 增量

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Bool_Parsing_Text_Variants` | `Parse_Bool` 支持 `1/true/TRUE/yes/on`，拒绝 `0/false/off/no`，并尊重 null default。 |
| `Test_Throw_On_Error_Uses_Controller_Error_Path` | `Resolve_Output_Plan(..., throw_on_error=true)` 在 invalid plan 时调用 controller error path，并携带具体错误信息。 |

### 13.2 `test_vds_trajectory_writer_with_mock_backend` 增量

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Vds_Precondition_Errors` | particle layout 未定义时 append particle frame 失败；没有已打开 shard 时 append observable frame 失败；SITS layout 未定义时 append SITS frame 失败；SITS module/shape 不匹配时失败。 |

这些测试保护 runtime call-site 接入时最容易产生 silent misuse 的错误路径：writer 必须在 layout 定义和 shard 生命周期满足条件后才接受 frame/module append。

## 14. 本轮增量：HighFive backend 上的 facade 文件级测试

本轮继续扩展 `test_highfive_backend_io`，在真实 HighFive/HDF5 backend 上覆盖 writer facade，而不仅是 backend primitive。

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Trajectory_Writer_With_Real_Backend` | 使用 `TrajectoryH5Writer + HighFiveBackend` 生成真实 `*.spg.h5md`，读回 position、velocity、box、ordinary observable dataset 的存在性和 shape。 |
| `Test_Restart_Writer_With_Real_Backend` | 使用 `RestartH5Writer + HighFiveBackend` 生成真实 `*.spgr.h5`，读回 structural state、run metadata、NHC state、SITS state、metad text state 的存在性和 shape。 |

这两个测试把 `backend-io` label 从“底层 HDF5 primitive 测试”扩展为“writer facade 与真实 backend 的组合测试”，可在后续编译/ctest 阶段更早暴露 facade 与 HighFiveBackend 之间的 shape、append、hard-link 或 layout 不一致问题。

## 15. 本轮增量：真实 HighFive VDS writer 文件级测试

本轮继续扩展 `test_highfive_backend_io`，新增 `VdsTrajectoryH5Writer + HighFiveBackendFactory` 的真实文件级测试。

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Vds_Trajectory_Writer_With_Real_Backend` | 使用真实 HighFive backend 写出 VDS wrapper 和两个 one-frame shard，随后用 HighFive 读回 wrapper。 |

该测试检查：

- wrapper 中 particle VDS dataset 存在且 shape 为 `[2,1,3]`。
- particle VDS 可跨 shard 读回坐标值。
- ordinary observable VDS 可跨 shard 读回 temperature 值。
- shard manifest dataset 存在且 frame counts 为 `[1,1]`。

这补齐了此前缺失的真实 backend 组合路径：mock VDS writer 测试证明 writer 调用逻辑，backend primitive VDS 测试证明 HDF5 VDS API，新增测试证明 writer 与 HighFive backend 组合后能生成可读 wrapper/shard 结构。

## 16. 本轮增量：真实 HighFive observable-only facade 文件级测试

本轮继续扩展 `test_highfive_backend_io`，新增 `ObservableH5Writer + HighFiveBackend` 的真实文件级测试。

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Observable_Writer_With_Real_Backend` | 使用真实 HighFive backend 写出 `*.obs.spg.h5md`，读回 observable-only layout、ordinary observable、QC module observable 和 QC SCF log。 |

该测试检查：

- 文件包含 `/h5md`、`/observables`、`/parameters`。
- 文件不包含 `/particles`。
- `/observables/all/temperature/value` 存在且 frame count 为 1。
- `/observables/all/qc/energy/value` 和 `spin_square/value` 存在。
- `/parameters/sponge/qc/scf_output` 存在。

这补齐了真实 backend facade 组合测试中的 observable-only artifact 路径，使 `backend-io` label 现在覆盖 trajectory、restart、observable-only 和 VDS 四类 H5 bundle artifact。

## 17. 本轮增量：真实 HDF5 completion/status metadata 读回

本轮继续扩展 `test_highfive_backend_io`，新增真实 HDF5 文件中的 output metadata 读回断言。

新增覆盖：

| Artifact | Metadata 覆盖 |
|---|---|
| basic H5 writer file | `/parameters/sponge/output/status = finalized`。 |
| observable-only primitive file | `status = finalized`，且无 `/particles`。 |
| trajectory facade file | `frame_count` 最后一项为 1，`last_complete_step` 为 10，`last_complete_time` 为 0.5，`status = finalized`。 |
| restart facade file | `frame_count` 最后一项为 1，`last_complete_step` 为 20，`status = finalized`。 |
| observable-only facade file | `frame_count` 最后一项为 1，`status = finalized`。 |
| VDS wrapper file | `frame_count` 最后一项为 2，`last_complete_step` 为 20，`status = finalized`，`repair_policy = strict`，`repair_status = not_applied`。 |

这些断言把 backend-io 覆盖从 payload dataset shape/value 扩展到 H5 bundle 的完成性契约字段，保护 finalize/status/completion metadata 不被后续重构破坏。

## 18. 本轮增量：H5 bundle 测试聚合构建目标

本轮为 `tests/h5_bundle` 新增聚合 build target：

```text
sponge_h5_bundle_tests
```

该 target 依赖所有 H5 bundle 测试二进制，可用于先验证测试目标能否全部编译/链接，再按 CTest label 分层运行。

建议后续首次验证顺序：

```bash
cmake -S . -B build-h5-tests -DSPONGE_BUILD_TESTS=ON
cmake --build build-h5-tests --target sponge_h5_bundle_tests
ctest --test-dir build-h5-tests -L contract --output-on-failure
ctest --test-dir build-h5-tests -L mock-backend --output-on-failure
ctest --test-dir build-h5-tests -L module --output-on-failure
ctest --test-dir build-h5-tests -L vds --output-on-failure
ctest --test-dir build-h5-tests -L backend-io --output-on-failure
```

该入口不影响默认构建；只有 `SPONGE_BUILD_TESTS=ON` 时才存在。

## 19. 本轮增量：真实 HighFive VDS module stream 读回

本轮继续扩展 `Test_Vds_Trajectory_Writer_With_Real_Backend`，使真实 `VdsTrajectoryH5Writer + HighFiveBackendFactory` 文件级测试覆盖 module-specific VDS readback。

新增真实 HDF5 wrapper 读回对象：

| Module stream | 覆盖点 |
|---|---|
| NHC coordinate | VDS dataset 存在，shape 为 `[2,2]`，跨 shard 读回数值。 |
| SITS `nk` | VDS dataset 存在，shape 为 `[2,3]`，跨 shard 读回数值。 |
| Metadynamics `meta` | scalar VDS dataset 存在，读回两个 frame 的值。 |
| QC energy | scalar VDS dataset 存在，读回两个 frame 的值。 |
| ReaxFF `bond` term | scalar VDS dataset 存在，读回两个 frame 的值。 |

这使 `backend-io` label 对 VDS 的覆盖从 particle/ordinary observable 扩展到 module-specific streams，形成 mock backend 逻辑断言与真实 HDF5 VDS readback 的双层覆盖。

## 20. 本轮增量：单元测试覆盖审计矩阵

新增覆盖审计文件：

```text
docs/sponge_h5_bundle_unit_test_audit_matrix.md
```

该矩阵逐项映射 H5 bundle 新路径到测试目标，按以下类别组织：

- contract and resolver;
- writer facade and mock backend;
- completion and failure;
- module-specific mapping;
- VDS writer;
- HighFive/HDF5 backend file-level behavior;
- documented integration-test gaps。

矩阵当前所有测试项标记为 `not-run`，因为尚未执行 configure/build/ctest。该文件用于后续完成审计：只有当 `sponge_h5_bundle_tests` 构建成功并且 `h5_bundle` CTest targets 通过后，才能把本 goal 视为完成。

## 21. 本轮增量：分层 build-only targets 与测试目标清单

本轮新增分层 build-only targets：

| Target | 用途 |
|---|---|
| `sponge_h5_bundle_contract_tests` | 编译 contract/completion 纯逻辑测试。 |
| `sponge_h5_bundle_mock_tests` | 编译 mock-backend writer/module/VDS 测试。 |
| `sponge_h5_bundle_backend_io_tests` | 编译真实 HighFive/HDF5 文件级测试。 |

同时新增测试清单：

```text
tests/h5_bundle/TEST_TARGETS.md
```

该清单记录 build-only targets、CTest targets、source files、labels 和推荐分阶段验证顺序。后续新增 H5 bundle 测试时，应同步更新该清单和覆盖审计矩阵。

## 22. 本轮增量：分阶段测试执行脚本

新增脚本：

```text
tests/h5_bundle/run_h5_bundle_tests.sh
```

该脚本封装 H5 bundle 测试的推荐验证顺序：

1. configure `SPONGE_BUILD_TESTS=ON`；
2. build/run contract tests；
3. build/run mock-backend tests；
4. run module and VDS labels；
5. build/run backend-io tests。

脚本支持 `SPONGE_H5_BUNDLE_TEST_BUILD_DIR`、`CMAKE`、`CTEST` 环境变量覆盖，便于通过 pixi 环境调用。当前脚本仅写入，尚未执行。

## 23. 本轮增量：HighFive string-array 与 VDS manifest 字符串数组读回

本轮继续扩展 `test_highfive_backend_io`，覆盖真实 HDF5 string-array 读写路径。

新增覆盖：

| 测试位置 | 覆盖点 |
|---|---|
| basic HighFive writer file | 写入并读回 `/parameters/sponge/test/string_array = ["alpha", "beta"]`。 |
| VDS wrapper file | 读回 `/parameters/sponge/output/shard_manifest/path` 和 `status` string arrays。 |

该测试保护两类实际契约路径：

- legacy sidecar provenance 依赖的 string-array 写入能力；
- VDS shard manifest 中 path/status 的 string-array 写入能力。

## 24. 本轮增量：HighFive VDS source 参数错误覆盖

本轮继续扩展 `Test_HighFive_Backend_Rejects_Invalid_Operations`，新增 VDS source 参数错误覆盖：

| 错误类型 | 覆盖点 |
|---|---|
| source rank mismatch | VDS spec rank 与 `source_dims` / `virtual_start` rank 不一致时，`Create_Virtual_Dataset` 必须失败。 |
| virtual start mismatch | `source_dims` 和 `virtual_start` 长度不一致时，`Create_Virtual_Dataset` 必须失败。 |

这些测试保护 HDF5 VDS mapping 的参数校验逻辑，避免错误 source hyperslab 配置被 silent accept。

## 25. 本轮增量：H5 output contract helper 边界测试

本轮继续扩展 `test_h5_output_plan`，补充 contract helper 的直接单元测试。

新增覆盖：

| 测试 | 覆盖点 |
|---|---|
| `Test_Contract_Helper_Functions` | `Any_H5_Output_Enabled`、`Legacy_Sidecars_Default_Enabled`、`Legacy_Sidecar_Requested`、`Legacy_Sidecar_Enabled` 在 null controller、无 H5 path、有 H5 path、显式 legacy sidecar 下的行为。 |
| `Test_Null_Controller_Resolver` | `Resolve_Output_Plan(nullptr, false)` 返回 invalid plan，并设置 `CONTROLLER is null` 错误信息。 |
| `Test_Suffix_And_Shard_Derivation` 扩展 | unknown output key 的 recommended suffix 返回 null。 |

这些测试保护 H5 output contract helper 本身，而不是只通过 resolver 间接覆盖。

## 26. 真实 HighFive backend 失败状态元数据测试补充

新增 `Test_HighFive_Backend_Failed_Metadata`，覆盖 `H5MDWriter::Mark_Failed` 在真实 HighFive backend 上的文件级落盘结果。

覆盖点：

- 使用真实 `HighFiveBackend` 创建 `failed.spg.h5md`。
- 调用 `Mark_Failed("intentional failure")` 后关闭文件。
- 重新以只读方式打开 HDF5 文件。
- 断言 `/parameters/sponge/output/status == "failed"`。
- 断言 `/parameters/sponge/output/error == "intentional failure"`。

该测试补齐此前 mock backend 已覆盖但真实 HDF5 文件未读回确认的失败元数据路径。

状态：已写入测试代码，尚未执行构建或 CTest。

## 27. 真实 HighFive backend 参数与 legacy 路径读回补充

扩展 `tests/h5_bundle/test_highfive_backend_io.cpp` 中真实 backend 测试，补齐 writer facade 已暴露但此前文件级读回不足的路径。

新增覆盖点：

- trajectory writer:
  - `Define_Particle_Datasets(..., include_force=true)` 创建 `/particles/all/force/value`。
  - `Append_Particle_Frame(..., force_xyz)` 写入 force frame。
  - 真实 HDF5 读回 force dataset shape 与 payload。
  - `Write_Mdinfo_Text` 写入并读回 `path::mdinfo_text`。
  - `Write_Legacy_Sidecar_Paths` 写入并读回 `path::legacy_sidecar_keys` 与 `path::legacy_sidecar_paths`。

- restart writer:
  - 写入并读回 legacy sidecar key/path 数组。
  - 覆盖 restart bundle 对 legacy 输出定位信息的 provenance 记录。

- observable writer:
  - 写入并读回 observable-only 文件中的 `mdinfo` 文本。
  - 写入并读回 legacy sidecar key/path 数组。
  - 写入并读回 `/parameters/sponge/provenance/launch_id`。

这些测试用于确认输出 bundle 在不写 legacy 文件或只保留 legacy provenance 时，相关路径仍能在实际 HDF5 文件中被审计。

状态：已写入测试代码，尚未执行构建或 CTest。

## 28. HighFive backend 状态机测试补充

新增 `Test_HighFive_Backend_Status_State`，覆盖底层 backend 状态枚举和状态数据集写入。

覆盖点：

- 未打开文件时调用 `Append_Int64` 必须失败。
- 失败后 `HighFiveBackend::Status()` 必须进入 `FileStatus::failed`。
- 失败后 `Last_Error()` 必须保留错误信息。
- 打开真实 HDF5 文件后状态为 `FileStatus::open`。
- `Set_Status(FileStatus::closing)` 必须同步更新内存状态并写入 `path::output_status`。
- `Close()` 后内存状态变为 `FileStatus::closed`。
- 重新打开文件读回 `/parameters/sponge/output/status == "closing"`。

该测试覆盖 writer facade 之外的 backend 级状态契约，避免状态写入仅通过 `Finalize` / `Mark_Failed` 被间接测试。

状态：已写入测试代码，尚未执行构建或 CTest。

## 29. 真实 VDS wrapper 参数路径测试补充

扩展 `Test_Vds_Trajectory_Writer_With_Real_Backend`，覆盖分片轨迹 wrapper 文件中的参数与诊断路径。

新增覆盖点：

- `Write_Metadynamics_Diagnostic("meta0", "hills", ...)` 写入 wrapper 的 metad 参数子树。
- `Write_Qc_Scf_Output(...)` 写入 wrapper 的 QC 参数子树。
- `Write_Legacy_Sidecar_Paths(...)` 写入 wrapper 的 legacy sidecar provenance。
- 真实 HighFive 文件读回 `/parameters/sponge/output/trajectory_chunk_size`。
- 真实 HighFive 文件读回 `/parameters/sponge/output/vds_status`。
- 真实 HighFive 文件读回 metad hills 文本、QC SCF 文本、legacy sidecar key/path 数组。

该补充确保 VDS 模式下不只验证 virtual dataset 本身，也验证 wrapper 作为 bundle 入口时必须承载的 SPONGE 参数与诊断记录。

状态：已写入测试代码，尚未执行构建或 CTest。

## 30. output plan 解析层边界测试补充

扩展 `tests/h5_bundle/test_h5_output_plan.cpp`，补充 h5 bundle 入口解析规则的边界覆盖。

新增覆盖点：

- `Recommended_Suffix_For_Key(nullptr)` 与未知 key 返回 `nullptr`。
- `Has_Recommended_Suffix` 对空 path/suffix 返回 false。
- `Trajectory_Chunk_Size` 对未设置 controller 返回默认值 20，对显式 key 返回设置值。
- legacy sidecar 全 key 默认行为：未启用任何 H5 bundle 时，`mdout/mdinfo/crd/box/vel/frc/rst/qc_scf_output` 默认 enabled，但不是 explicit。
- 启用任一 H5 bundle 后，legacy sidecar 默认关闭，仅显式设置的 legacy key 保持 enabled/explicit。
- trajectory/restart/observable 非推荐后缀只设置 `has_recommended_suffix=false`，不使 plan invalid。
- `complete_prefix` repair policy 在未提供 trajectory path 或未启用 VDS 时必须 invalid，并给出 requires 类错误信息。

该补充把 grouped TOML flatten 后的 `output_h5_*` 入口契约、默认 legacy gating 和 VDS repair 前置条件集中成可回归测试。

状态：已写入测试代码，尚未执行构建或 CTest。

## 31. HighFive backend append 数据完整性失败路径补充

扩展 `Test_HighFive_Backend_Rejects_Invalid_Operations`，补充真实 backend 对 append 数据记录的基础校验。

新增覆盖点：

- 创建 appendable scalar observable dataset。
- 使用超出单条 record size 的 `count=2` append，必须失败。
- 使用 `nullptr` 数据指针 append，必须失败。

该测试保护 H5MD frame 写入路径的底层不变量：每次 append 只能追加一个完整 record，且数据指针不可为空。它覆盖 trajectory、observable、module scalar 等共享 backend 逻辑。

状态：已写入测试代码，尚未执行构建或 CTest。

## 32. restart bundle 单帧不变量真实 backend 测试补充

新增 `Test_Restart_Writer_Rejects_Second_State_With_Real_Backend`，覆盖 `.spgr.h5` restart bundle 只保留一帧的核心约束。

覆盖点：

- 使用真实 `HighFiveBackend` 打开 restart writer。
- 写入第一帧 structural state 成功。
- 第二次 `Write_Structural_State` 必须失败。
- 关闭后重新打开 `.spgr.h5` 文件。
- 断言 `/parameters/sponge/output/status == "failed"`。
- 断言 `/parameters/sponge/output/error` 包含 `already contains one structural state`。
- 断言 position dataset 第一维仍为 1，确认未扩展为多帧 restart。

该测试把 restart 仅保留一帧的设计约束从 mock 语义推进到真实 HDF5 文件级审计。

状态：已写入测试代码，尚未执行构建或 CTest。

## 33. ModuleH5MappingWriter 真实 backend metad/ReaxFF 路径补充

新增 `Test_Module_Metad_And_Reaxff_With_Real_Backend`，直接面向 `ModuleH5MappingWriter` 验证模块扩展路径在真实 HDF5 文件中的布局。

覆盖点：

- metadynamics diagnostic convenience API：
  - `Write_Metadynamics_Potential_Export`
  - `Write_Metadynamics_Direct_Export`
  - `Write_Metadynamics_Hills`
  - `Write_Metadynamics_History`
  - `Write_Metadynamics_Edge`
- 真实 HighFive 文件读回 `/parameters/sponge/metadynamics/<name>/potential_export`。
- 真实 HighFive 文件读回 `/parameters/sponge/metadynamics/<name>/direct_export`。
- 真实 HighFive 文件读回 `/parameters/sponge/metadynamics/<name>/hills`。
- 真实 HighFive 文件读回 `/parameters/sponge/metadynamics/<name>/history`。
- 真实 HighFive 文件读回 `/parameters/sponge/metadynamics/<name>/edge`。
- ReaxFF 多 term observable：`bond/angle/over`。
- 真实 HighFive 文件读回三个 ReaxFF term 的 value payload。
- 文件最终状态必须为 `finalized`。

该测试补齐此前 VDS wrapper 只覆盖单个 generic metad diagnostic、ReaxFF 只覆盖单 term 的不足。

状态：已写入测试代码，尚未执行构建或 CTest。

## 34. H5MDWriter common layout 与 detached backend 测试补充

扩展 `tests/h5_bundle/test_highfive_backend_io.cpp`，补充 `H5MDWriter` 自身的基础契约。

新增覆盖点：

- `H5MDWriter(nullptr)` 必须报告未 attached。
- detached writer 的 `Open/Flush/Close/Finalize/Mark_Failed` 均必须安全失败。
- detached writer 的 `Status()` 必须为 `FileStatus::closed`。
- detached writer 的 `Last_Error()` 必须返回明确错误信息。
- 真实 HighFive common layout 初始化后必须写入 schema name/version。
- 真实 HighFive common layout 初始化后必须写入 output completion 初值：
  - `frame_count = 0`
  - `last_complete_step = -1`
  - `last_complete_time = 0.0`

该测试确保 writer facade 的基础失败语义和所有 bundle 文件共享的 common metadata 不只依赖上层 writer 间接覆盖。

状态：已写入测试代码，尚未执行构建或 CTest。

## 35. writer facade open precondition 测试补充

扩展 `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp`，新增 `Test_Writer_Open_Preconditions`。

覆盖点：

- `TrajectoryH5Writer::Open_Single_File` 在 trajectory plan 未 enabled 时必须失败。
- `TrajectoryH5Writer::Open_Single_File` 在 trajectory plan 为 VDS 模式时必须失败。
- `RestartH5Writer::Open` 在 restart plan 未 enabled 时必须失败。
- `ObservableH5Writer::Open` 在 observable plan 未 enabled 时必须失败。
- 每个失败路径均检查 `Last_Error()` 包含对应 precondition 说明。

该测试防止 MD_core 集成层传入错误 plan 时 writer facade 静默打开错误类型的 bundle。

状态：已写入测试代码，尚未执行构建或 CTest。

## 36. ObservableH5Writer 模块输出真实 backend 覆盖补充

扩展 `Test_Observable_Writer_With_Real_Backend`，把 observable-only bundle 入口的模块输出路径纳入真实 HighFive 文件级读回。

新增覆盖点：

- `ObservableH5Writer::Ensure_Nose_Hoover_Chain_Observables` 与 `Append_Nose_Hoover_Chain_Frame`。
- 真实 HighFive 文件读回 NHC coordinate/velocity payload。
- `ObservableH5Writer::Ensure_Sits_Nk_Observable` 与 `Append_Sits_Nk_Frame`。
- 真实 HighFive 文件读回 SITS `nk` payload。
- `ObservableH5Writer::Ensure_Metadynamics_Scalars` 与 `Append_Metadynamics_Scalar_Frame`。
- 真实 HighFive 文件读回 metad `meta/rbias/rct` payload。
- `ObservableH5Writer::Ensure_Reaxff_Energy_Terms` 与 `Append_Reaxff_Frame`。
- 真实 HighFive 文件读回 ReaxFF `bond/angle` payload。

该测试确认 observable-only `.obs.spg.h5md` 不仅能承载普通 mdout/qc 数据，也能承载当前 bundle 设计中归入 observable H5MD 的模块化轻量分析数据。

状态：已写入测试代码，尚未执行构建或 CTest。

## 37. TrajectoryH5Writer 模块输出真实 backend 覆盖补充

扩展 `Test_Trajectory_Writer_With_Real_Backend`，把 single-file trajectory bundle 入口的模块输出路径纳入真实 HighFive 文件级读回。

新增覆盖点：

- `TrajectoryH5Writer::Ensure_Nose_Hoover_Chain_Observables` 与 `Append_Nose_Hoover_Chain_Frame`。
- 真实 HighFive 文件读回 NHC coordinate/velocity payload。
- `TrajectoryH5Writer::Ensure_Sits_Nk_Observable` 与 `Append_Sits_Nk_Frame`。
- 真实 HighFive 文件读回 SITS `nk` payload。
- `TrajectoryH5Writer::Ensure_Metadynamics_Scalars` 与 `Append_Metadynamics_Scalar_Frame`。
- 真实 HighFive 文件读回 metad `meta` payload。
- `TrajectoryH5Writer::Ensure_Qc_Observables` 与 `Append_Qc_Frame`。
- 真实 HighFive 文件读回 QC energy/spin_square payload。
- `TrajectoryH5Writer::Ensure_Reaxff_Energy_Terms` 与 `Append_Reaxff_Frame`。
- 真实 HighFive 文件读回 ReaxFF `bond/angle` payload。

该测试补齐 single-file `.spg.h5md` 在非 VDS 模式下的模块输出覆盖，与 VDS wrapper 和 observable-only bundle 测试形成三路一致性检查。

状态：已写入测试代码，尚未执行构建或 CTest。

## 38. VDS complete-prefix repair 保留完整前缀测试补充

扩展 `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`，新增 `Test_Complete_Prefix_Repair_Retains_Valid_Prefix`。

新增覆盖点：

- 设置 `trajectory.chunk_size = 1`，强制每帧一个 shard。
- 第一个 shard finalize 成功并进入 manifest。
- 第二个 shard finalize 失败。
- `Finalize_With_Repair()` 必须成功，并调用 complete-prefix repair 语义。
- repair 后 manifest 必须只保留第一个完整 shard。
- `Total_Trajectory_Frame_Count()` 必须重算为保留前缀的 frame 数。
- wrapper 中 position virtual dataset 只能引用保留的 `segment_000000.spg.h5md`。
- wrapper 写出：
  - `/parameters/sponge/output/repair_policy = allow_complete_prefix`
  - `/parameters/sponge/output/repair_status = applied`
  - `/parameters/sponge/output/repaired_shard_count`
- shard manifest status 只保留一个 `complete` 项。

该测试补齐此前只覆盖“当前 shard 失败导致 manifest 为空”的 repair 场景，明确验证已有完整前缀不会被丢弃。

状态：已写入测试代码，尚未执行构建或 CTest。

## 39. OutputCompletionTracker 真实 backend 文件级测试补充

扩展 `tests/h5_bundle/test_highfive_backend_io.cpp`，新增 `Test_Output_Completion_Tracker_With_Real_Backend`。

新增覆盖点：

- 使用真实 `HighFiveBackend` 和 `H5MDWriter` 打开 `.spg.h5md` 文件。
- `OutputCompletionTracker::Mark_Open()` 写入 open 状态和初始计数。
- `Begin_Frame/Complete_Frame` 完成一帧。
- `Mark_Closing()` 写入 closing 状态。
- `Mark_Finalized()` 写入 finalized 状态和最终计数。
- 真实 HighFive 文件读回：
  - `/parameters/sponge/output/status = finalized`
  - `frame_count` 最后一项为 1
  - `last_complete_step` 最后一项为 100
  - `last_complete_time` 最后一项为 1.5

该测试把 completion tracker 从纯逻辑状态机推进到真实 HDF5 文件级可审计路径。

状态：已写入测试代码，尚未执行构建或 CTest。

## 40. H5MD step/time hard-link 契约测试补充

扩展 `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp`，补充 H5MD frame axis hard-link 契约。

新增覆盖点：

- `TrajectoryH5Writer`:
  - position `step/time` hard links。
  - velocity `step/time` hard links。
  - force `step/time` hard links。
  - box edges `step/time` hard links。
  - ordinary observable `temperature` 的 `step/time` hard links。
- `ObservableH5Writer`:
  - observable-only `energy` 的 `step/time` hard links。
- `RestartH5Writer`:
  - position `step/time` hard links。
  - velocity `step/time` hard links。
  - box edges `step/time` hard links。

该测试确保 H5MD 中各 value dataset 共享正确 frame axis，而不是各自生成不一致的 step/time 数据集。

状态：已写入测试代码，尚未执行构建或 CTest。

## 41. ModuleH5MappingWriter step/time hard-link 契约测试补充

扩展 `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp`，补充模块 observable 的 frame axis hard-link 契约。

新增覆盖点：

- NHC coordinate/velocity 的 `step/time` hard links。
- metadynamics `meta/rbias/rct` 的 `step/time` hard links。
- QC `energy/spin_square` 的 `step/time` hard links。
- ReaxFF `bond/angle` 的 `step/time` hard links。

注意：当前 SITS `nk` 实现创建独立 step/time 数据集并写入 value，但没有给 `nk/value` 建 `step/time` hard link。本轮测试按当前实现覆盖，不引入会失败的假设。

状态：已写入测试代码，尚未执行构建或 CTest。

## 42. VDS module append precondition 测试补充

扩展 `Test_Vds_Precondition_Errors`，补充 VDS 模块输出路径的运行时 precondition 覆盖。

新增覆盖点：

- NHC layout 未定义时 append 必须失败。
- NHC chain length 与已定义 layout 不一致时必须失败。
- metadynamics scalar layout 未定义时 append 必须失败。
- QC layout 未定义时 append 必须失败。
- ReaxFF layout 未定义时 append 必须失败。

这些测试覆盖 VDS trajectory 中模块 frame 必须在 particle shard 打开后、且对应 layout 已定义后才能写入的约束。

状态：已写入测试代码，尚未执行构建或 CTest。

## 43. restart base layout 与 run metadata 真实 backend 覆盖补充

扩展 `Test_Restart_Writer_With_Real_Backend`，补充 `.spgr.h5` restart bundle 基础分组和 run metadata 读回。

新增覆盖点：

- `/parameters/sponge/restart` root 存在。
- restart thermostat group 存在。
- restart barostat group 存在。
- restart bias group 存在。
- `/run/current_time` dataset 存在并读回最后值。
- `/run/state_type` 写入并读回 `restart`。

该测试补齐 restart base layout 中此前仅被创建但没有真实文件级断言的路径。

状态：已写入测试代码，尚未执行构建或 CTest。

## 44. restart optional velocity 与 output plan helper 边界补充

新增两个小边界覆盖：

- `Test_Restart_Writer_Rejects_Second_State_With_Real_Backend` 中补充 `include_velocity=false` 场景，确认真实 `.spgr.h5` 不创建 `path::velocity_value`。
- `Test_Contract_Helper_Edge_Cases` 中补充 `SpongeH5OutputPlan::Command_String`：
  - null controller 返回空字符串。
  - missing key 返回空字符串。
  - existing key 返回 controller 中的 value。

该补充覆盖 restart optional velocity 语义和 output plan helper 的入口边界。

状态：已写入测试代码，尚未执行构建或 CTest。

## 45. VDS source path 相对化测试补充

扩展 `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`，新增 `Test_Vds_Source_Path_Relativization`。

新增覆盖点：

- wrapper 路径位于 `/tmp/sponge_h5_vds_case/wrappers/prod.spg.h5md`。
- shard root 位于 sibling 目录 `/tmp/sponge_h5_vds_case/shards/prod.spg.shards`。
- `Finalize()` 后 wrapper virtual dataset source path 必须写成相对路径：
  - `../shards/prod.spg.shards/segment_000000.spg.h5md`

该测试覆盖 `Vds_Source_Path` 的相对路径分支，确保 VDS wrapper 和 shard bundle 可整体移动，避免在云端分段同步场景中固化绝对路径。

状态：已写入测试代码，尚未执行构建或 CTest。

## 46. VdsTrajectoryH5Writer open precondition 测试补充

扩展 `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`，新增 `Test_Vds_Open_Precondition_Errors`。

新增覆盖点：

- backend factory 为 null 时 `Open` 必须失败。
- trajectory plan 未 enabled 时 `Open` 必须失败。
- trajectory plan 非 VDS 模式时 `Open` 必须失败。
- 每个失败路径均断言 `Last_Error()` 的明确错误信息。

该测试补齐 VDS writer 在进入 shard/VDS materialization 前的入口保护。

状态：已写入测试代码，尚未执行构建或 CTest。

## 47. VDS 空输出 finalize 边界测试补充

扩展 `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`，新增 `Test_Vds_Finalize_Without_Frames`。

新增覆盖点：

- VDS writer 只 `Open`，不定义 layout、不写 frame，直接 `Finalize()`。
- manifest 必须为空。
- total trajectory/observable frame count 必须为 0。
- wrapper 必须 finalized。
- wrapper 写出：
  - `repair_policy = strict`
  - `repair_status = not_applied`
  - `repaired_shard_count`
  - `vds_status`
- wrapper 不应创建任何 virtual dataset。

该测试覆盖空输出/无帧退出场景，确保 VDS wrapper 的 metadata 仍完整，且不会生成无意义的 VDS datasets。

状态：已写入测试代码，尚未执行构建或 CTest。

## 48. output plan 解析鲁棒性补充

扩展 `tests/h5_bundle/test_h5_output_plan.cpp`，补充 mdin/grouped TOML flatten 后的解析鲁棒性测试。

新增覆盖点：

- `output_h5_trajectory_repair_policy = COMPLETE_PREFIX` 必须被 lower-case 后识别为 `complete_prefix`。
- `output_h5_trajectory_vds = TRUE` 必须被识别为 true。
- `output_h5_trajectory_chunk_size = -3` 必须 invalid。
- `output_h5_trajectory_chunk_size = not_an_integer` 经 `atoi` 解析为 0 后必须 invalid。
- H5 bundle 启用后，explicit legacy sidecar `vel = legacy.vel` 必须保留 path 值。
- legacy sidecar 列表大小必须保持当前契约的 8 项。

该补充覆盖用户可编辑 mdin 输入到 `ResolvedOutputPlan` 的常见错误和大小写边界。

状态：已写入测试代码，尚未执行构建或 CTest。

## 49. VDS source path 无 wrapper parent 边界补充

扩展 `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`，新增 `Test_Vds_Source_Path_Without_Wrapper_Parent`。

新增覆盖点：

- wrapper path 为 `prod.spg.h5md`，没有 parent directory。
- shard root 为 `prod.spg.shards`。
- `Finalize()` 后 wrapper virtual dataset source path 必须保持：
  - `prod.spg.shards/segment_000000.spg.h5md`

该测试覆盖 `Vds_Source_Path` 中 `wrapper_parent.empty()` 的分支，与 sibling directory 相对路径测试互补。

状态：已写入测试代码，尚未执行构建或 CTest。

## 50. H5MDWriter repeated output completion 真实 backend 测试补充

扩展 `tests/h5_bundle/test_highfive_backend_io.cpp`，新增 `Test_H5MD_Writer_Repeated_Output_Completion`。

新增覆盖点：

- `H5MDWriter::Open` 初始化 output completion：`frame_count=0, last_step=-1, last_time=0.0`。
- 连续调用 `Write_Output_Completion(1, 10, 0.5)`。
- 连续调用 `Write_Output_Completion(2, 20, 1.0)`。
- 真实 HighFive 文件读回三条 completion 记录：
  - `frame_count = [0, 1, 2]`
  - `last_complete_step = [-1, 10, 20]`
  - `last_complete_time = [0.0, 0.5, 1.0]`

该测试覆盖 repeated `Create_Dataset` 的 idempotence 和 append history 保留，是 trajectory/restart/observable/VDS 共享 completion 路径的基础回归测试。

状态：已写入测试代码，尚未执行构建或 CTest。

## 51. ObservableH5Writer missing value 真实失败路径补充

扩展 `tests/h5_bundle/test_highfive_backend_io.cpp`，新增 `Test_Observable_Writer_Missing_Value_With_Real_Backend`。

新增覆盖点：

- observable-only `.obs.spg.h5md` 定义 `energy/temperature` 两个 observable。
- append frame 时只提供 `energy`，缺失 `temperature`，必须失败。
- 真实 HighFive 文件读回：
  - `/parameters/sponge/output/status = failed`
  - `/parameters/sponge/output/error = observable value is missing: temperature`
- 失败 frame 不得推进 completion：
  - `frame_count.back() == 0`
  - `last_complete_step.back() == -1`

该测试把 observable missing value 的失败语义从 mock backend 推进到真实文件级契约，防止半写 frame 被误标记为 completed。

状态：已写入测试代码，尚未执行构建或 CTest。

## 52. TrajectoryH5Writer observable missing value 行为补充

扩展 `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp`，新增 `Test_Trajectory_Observable_Missing_Value_Does_Not_Advance`。

新增覆盖点：

- trajectory writer 定义 `energy/temperature` 两个 observable。
- append observable frame 时只提供 `energy`，缺失 `temperature`，必须返回 false。
- `Observable_Frame_Count()` 必须保持 0。
- `Last_Error()` 必须记录 `observable value is missing: temperature`。
- 当前实现下 trajectory writer 不调用 `Mark_Failed`，mock backend status 保持 open，且不写 output error。

该测试明确固定 trajectory writer 与 observable-only writer 在 missing observable value 失败处理上的差异，避免后续误改或误解。

状态：已写入测试代码，尚未执行构建或 CTest。

## 53. HighFive string/string-array overwrite 真实 backend 测试补充

扩展 `tests/h5_bundle/test_highfive_backend_io.cpp`，新增 `Test_HighFive_Backend_String_Overwrite`。

新增覆盖点：

- 同一路径 `path::mdinfo_text` 连续 `Write_String` 两次。
- 真实 HighFive 文件读回时必须只保留新值 `new mdinfo`。
- 同一路径 `path::legacy_sidecar_keys` 连续 `Write_String_Array` 两次。
- 真实 HighFive 文件读回时必须只保留新数组 `new_key/second_key`。

该测试覆盖 HighFive backend 中 `Delete_If_Exists` 后重建 string/string-array dataset 的行为，是 status、mdinfo、legacy provenance、repair metadata 等参数路径重复写入的基础保障。

状态：已写入测试代码，尚未执行构建或 CTest。

## 54. RestartH5Writer custom run metadata 真实 backend 测试补充

扩展 `tests/h5_bundle/test_highfive_backend_io.cpp`，新增 `Test_Restart_Writer_Custom_Run_Metadata_With_Real_Backend`。

新增覆盖点：

- 直接调用 `RestartH5Writer::Write_Run_Metadata(42, 2.5, "checkpoint")`。
- 真实 `.spgr.h5` 文件读回：
  - `/run/current_step = 42`
  - `/run/current_time = 2.5`
  - `/run/state_type = checkpoint`
  - `/parameters/sponge/output/status = finalized`

该测试覆盖 run metadata 的公开写入路径，而不是只通过 structural restart state 间接覆盖默认 `state_type=restart`。

状态：已写入测试代码，尚未执行构建或 CTest。

## 55. VDS shard manifest 持久化字段测试补充

扩展 `Test_Vds_Wrapper_And_Module_Virtual_Datasets`，补充 wrapper shard manifest 数值字段的写出断言。

新增覆盖点：

- `path::shard_manifest_index` dataset 创建并写入。
- `path::shard_manifest_frame_start` dataset 创建并写入。
- `path::shard_manifest_frame_count` dataset 创建并写入。
- `path::shard_manifest_step_start` dataset 创建并写入。
- `path::shard_manifest_step_end` dataset 创建并写入。
- `path::shard_manifest_time_start` dataset 创建并写入。
- `path::shard_manifest_time_end` dataset 创建并写入。
- 对 2 个 shard 的测试场景，以上每个字段 append count 均为 2。

注意：当前 `Write_Manifest_To_Wrapper` 只持久化 path/status/index/frame/step/time 字段；module frame counts 仍只存在内存 manifest，不写入 wrapper 文件。本测试按当前持久化契约覆盖。

状态：已写入测试代码，尚未执行构建或 CTest。

## 56. HighFiveBackend invalid dataset definition 测试补充

扩展 `Test_HighFive_Backend_Rejects_Invalid_Operations`，补充底层 dataset 定义的非法输入覆盖。

新增覆盖点：

- `Create_Dataset` 使用空 path 必须失败。
- `Create_Dataset` 使用空 dims 必须失败。

该测试覆盖所有 H5 bundle writer 共享的 dataset definition 边界，防止空路径或空维度被静默写入 HDF5 backend。

状态：已写入测试代码，尚未执行构建或 CTest。

## 57. HighFiveBackend close 后误写保护测试补充

扩展 `Test_HighFive_Backend_Status_State`，补充文件关闭后的写入保护。

新增覆盖点：

- `HighFiveBackend::Close()` 后内存状态为 `FileStatus::closed`。
- close 后再次 `Write_String` 必须失败。
- 失败后 backend 状态进入 `FileStatus::failed`。
- `Last_Error()` 必须包含文件未打开相关信息。

该测试覆盖 bundle finalization/close 之后的误写保护，避免关闭后的输出路径被静默接受。

状态：已写入测试代码，尚未执行构建或 CTest。

## 58. VDS wrapper parameter write precondition 测试补充

扩展 `Test_Vds_Open_Precondition_Errors`，补充 VDS wrapper 参数写入在 wrapper 未 open 时的失败路径。

新增覆盖点：

- `Write_Metadynamics_Diagnostic` 在 wrapper 未 open 时必须失败。
- `Write_Qc_Scf_Output` 在 wrapper 未 open 时必须失败。
- `Write_Legacy_Sidecar_Paths` 在 wrapper 未 open 时必须失败。
- 每个失败路径均要求 `Last_Error() == "VDS wrapper is not open"`。

该测试覆盖 VDS wrapper 中 `/parameters/sponge/...` 参数和 provenance 写入入口保护，避免在 wrapper 未初始化时静默接受参数输出。

状态：已写入测试代码，尚未执行构建或 CTest。

## 59. output H5 suffix 精确匹配边界测试补充

扩展 `Test_Suffix_And_Shard_Derivation`，补充推荐后缀匹配的边界情况。

新增覆盖点：

- 输入字符串短于 suffix 时 `Ends_With` 必须返回 false。
- 相似但缺少最后字符的 suffix 必须返回 false。
- 推荐 suffix 后还有额外扩展名时必须返回 false。

该测试确保 `*.spg.h5md`、`*.spgr.h5`、`*.obs.spg.h5md` 的推荐后缀识别不会误判临时文件或不完整文件名。

状态：已写入测试代码，尚未执行构建或 CTest。

## 60. OutputCompletionTracker 状态机错误路径测试补充

扩展 `test_completion_tracker.cpp`，新增 `Test_State_Machine_Error_Paths_Without_Writer`。

新增覆盖点：

- 未打开 frame 时调用 `Complete_Frame` 必须失败。
- 未打开 frame 的 completion 失败不得推进 `frame_count` 或 `last_complete_step`。
- frame 未完成时调用 `Mark_Closing` 必须失败。
- close 失败后状态仍保持 `open`，且 incomplete frame 标志仍为 true。
- 显式 `Mark_Failed` 必须设置 `FileStatus::failed` 并记录 failure reason。

该测试补齐 completion tracker 的无 writer 纯状态机错误分支，避免只覆盖成功流和 metadata 写入流。

状态：已写入测试代码，尚未执行构建或 CTest。

## 61. grouped output_h5 parser-visible key 测试补充

扩展 `test_h5_output_plan.cpp`，新增 `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating`。

新增覆盖点：

- 固定 `[output.h5.trajectory]` flatten 后的 key 为 `output_h5_trajectory_path`。
- 固定 VDS 开关 key 为 `output_h5_trajectory_vds`。
- 固定 VDS chunk size key 为 `output_h5_trajectory_chunk_size`。
- 固定 VDS repair policy key 为 `output_h5_trajectory_repair_policy`。
- 固定 restart bundle key 为 `output_h5_restart_path`。
- 固定 observable-only bundle key 为 `output_h5_observable_path`。
- 旧命名 `output_trajectory_h5_path`、`output_traj_h5`、`output_restart_h5`、`output_observable_h5` 不得误触发新 H5 bundle 输出。
- trajectory、restart、observable 任一 H5 path key 显式存在时，legacy sidecar 默认关闭。
- H5 bundle 开启后，只有显式 legacy path 仍允许写 legacy sidecar。

该测试保护本轮约定的 grouped TOML flatten 命名，避免后续实现把早期讨论中的旧 key 重新引入为有效入口。

状态：已写入测试代码，尚未执行构建或 CTest。

## 62. RestartH5Writer base layout 路径测试补充

扩展 `test_h5md_writers_with_mock_backend.cpp`，新增 `Test_Restart_Writer_Base_Layout_Paths`。

新增覆盖点：

- restart bundle 打开时必须创建 `/run`。
- restart bundle 打开时必须创建 H5MD particle 根路径 `/particles/all`。
- restart bundle 打开时必须创建 `/particles/all/position`。
- restart bundle 打开时必须创建 `/particles/all/velocity`。
- restart bundle 打开时必须创建 `/particles/all/box`。
- restart bundle 打开时必须创建 `/particles/all/box/edges`。
- restart bundle 打开时必须创建 `/parameters/sponge/restart`。
- restart bundle 打开时必须创建 thermostat restart 子树。
- restart bundle 打开时必须创建 barostat restart 子树。
- restart bundle 打开时必须创建 bias restart 子树。
- restart bundle 打开时必须创建 SITS restart 子树。
- restart bundle 打开时必须创建 metadynamics restart 子树。
- restart bundle schema name 必须写为 `sponge.restart.h5`。
- restart bundle schema version 必须来自 writer open 参数。

该测试把 restart H5 bundle 的基础路径结构固定为显式单元测试契约，避免后续只在具体 state 写入路径中间接覆盖。

状态：已写入测试代码，尚未执行构建或 CTest。

## 63. Trajectory/Observable H5 writer base layout 路径测试补充

扩展 `test_h5md_writers_with_mock_backend.cpp`，新增 `Test_Trajectory_And_Observable_Base_Layout_Paths`。

新增覆盖点：

- trajectory writer 打开单文件输出时必须创建 `/particles/all`。
- trajectory writer 打开时必须创建 position、velocity、force、box、box/edges 子树。
- trajectory writer 打开时必须创建 `/observables/all`。
- trajectory writer 打开时必须创建 `/parameters/sponge/mdout`。
- trajectory writer 打开时必须创建 mdout columns 子树。
- trajectory writer 打开时必须创建 `/parameters/sponge/log`。
- trajectory writer schema name 必须为 `sponge.output.h5md`。
- observable-only writer 打开时不得创建 `/particles`。
- observable-only writer 打开时必须创建 `/observables/all`。
- observable-only writer 打开时必须创建 mdout/log 参数子树。
- observable-only writer schema name 必须为 `sponge.output.h5md`。
- 两类 writer 的 schema version 必须来自 open 参数。

该测试把 trajectory 与 observable-only bundle 的 open 阶段 layout 契约固定下来，避免仅依赖后续 dataset 定义间接覆盖。

状态：已写入测试代码，尚未执行构建或 CTest。

## 64. Metadynamics diagnostic parameter 路径测试补充

扩展 `test_module_h5_mappings_with_mock_backend.cpp` 的 `Test_Metadynamics_And_Diagnostics`。

新增覆盖点：

- metadynamics diagnostic 必须创建 `module_path::metad_parameter_root`。
- metadynamics diagnostic 必须创建具体 bias 名称子 group。
- `Write_Metadynamics_History` 必须写入 `<metad_parameter_root>/<name>/history`。
- `Write_Metadynamics_Edge` 必须写入 `<metad_parameter_root>/<name>/edge`。
- `Write_Metadynamics_Potential_Export` 必须写入 `<metad_parameter_root>/<name>/potential_export`。
- 既有 hills/direct_export 断言保留，形成 hills/history/edge/potential_export/direct_export 全覆盖。

该测试补齐增强采样输出中 metad diagnostic 参数区的路径级覆盖，避免只检查部分 convenience writer。

状态：已写入测试代码，尚未执行构建或 CTest。

## 65. QC optional spin_square observable 路径测试补充

扩展 `test_module_h5_mappings_with_mock_backend.cpp`，新增 `Test_Qc_Optional_Spin_Square_Path`。

新增覆盖点：

- `Ensure_Qc_Observables(false)` 必须创建 QC observable 根 group。
- `Ensure_Qc_Observables(false)` 必须创建 QC energy value 路径。
- `Ensure_Qc_Observables(false)` 不得创建 QC spin_square value 路径。

该测试固定 QC 输出中可选 `spin_square` observable 的路径创建语义，避免不启用 spin-square 时向 H5MD observable bundle 写入额外字段。

状态：已写入测试代码，尚未执行构建或 CTest。

## 66. Restart dynamic extension group 与 legacy provenance 路径测试补充

扩展 `test_h5md_writers_with_mock_backend.cpp` 的 restart/module 与 legacy provenance 测试。

新增覆盖点：

- restart SITS state 写入时必须创建动态 module group：`<restart_sits>/sits_a`。
- restart metad state 写入时必须创建动态 bias group：`<restart_meta>/meta0`。
- restart legacy sidecar provenance 写入时必须创建 `/parameters/sponge/files`。
- restart legacy sidecar provenance 写入时必须创建 `path::legacy_sidecars`。
- trajectory legacy sidecar provenance 写入时必须创建 `/parameters/sponge/files`。
- trajectory legacy sidecar provenance 写入时必须创建 `path::legacy_sidecars`。
- observable-only legacy sidecar provenance 写入时必须创建 `/parameters/sponge/files`。
- observable-only legacy sidecar provenance 写入时必须创建 `path::legacy_sidecars`。

该测试补齐动态 restart extension group 与三类 writer 共享 legacy sidecar provenance group 的路径级覆盖，防止只验证 string array 数据而漏掉 provenance 容器路径。

状态：已写入测试代码，尚未执行构建或 CTest。

## 67. VDS trajectory chunk size 参数路径测试补充

扩展 `test_vds_trajectory_writer_with_mock_backend.cpp` 的 `Test_Vds_Wrapper_And_Module_Virtual_Datasets`。

新增覆盖点：

- VDS wrapper 打开时必须写入 `/parameters/sponge/output/trajectory_chunk_size`。
- `Make_Vds_Plan()` 中 `trajectory.chunk_size=2` 必须以字符串 `"2"` 写入 wrapper 参数区。

该测试固定 grouped `output_h5_trajectory_chunk_size` 解析结果在 H5MD wrapper 内的持久化位置，避免只通过 shard 数量间接观察 chunk 行为。

状态：已写入测试代码，尚未执行构建或 CTest。

## 68. VDS wrapper shard manifest group 与 schema metadata 路径测试补充

扩展 `test_vds_trajectory_writer_with_mock_backend.cpp` 的 `Test_Vds_Wrapper_And_Module_Virtual_Datasets`。

新增覆盖点：

- VDS wrapper 打开时必须创建 `path::shard_manifest` group。
- VDS wrapper schema name 必须写为 `sponge.output.h5md`。
- VDS wrapper schema version 必须来自 `Open(plan, schema_version)` 的参数。

该测试固定 VDS wrapper 自身的结构元数据，避免只通过最终 manifest datasets 和 virtual datasets 间接覆盖 wrapper 初始化行为。

状态：已写入测试代码，尚未执行构建或 CTest。

## 69. VDS shard file path 与 schema metadata 测试补充

扩展 `test_vds_trajectory_writer_with_mock_backend.cpp` 的 `Test_Vds_Wrapper_And_Module_Virtual_Datasets`。

新增覆盖点：

- 第一个 VDS shard 必须打开为 `<derived_shard_root>/segment_000000.spg.h5md`。
- 第二个 VDS shard 必须打开为 `<derived_shard_root>/segment_000001.spg.h5md`。
- 每个 VDS shard schema name 必须写为 `sponge.output.h5md`。
- 每个 VDS shard schema version 必须来自 wrapper `Open(plan, schema_version)` 参数。

该测试直接约束 VDS shard 文件命名和 shard 自身 H5MD schema metadata，避免只通过 wrapper VDS source path 间接覆盖 shard 创建逻辑。

状态：已写入测试代码，尚未执行构建或 CTest。

## 70. Real HighFive VDS shard 文件读回测试补充

扩展 `test_highfive_backend_io.cpp` 的 `Test_Vds_Trajectory_Writer_With_Real_Backend`。

新增覆盖点：

- finalize 后真实文件系统中必须存在 `segment_000000.spg.h5md`。
- finalize 后真实文件系统中必须存在 `segment_000001.spg.h5md`。
- 第一个真实 shard 文件必须写入 `/parameters/sponge/schema/name = sponge.output.h5md`。
- 第二个真实 shard 文件必须写入 `/parameters/sponge/schema/name = sponge.output.h5md`。
- 第一个真实 shard 文件必须写入 `/parameters/sponge/schema/version = test`。
- 第二个真实 shard 文件必须写入 `/parameters/sponge/schema/version = test`。
- 两个真实 shard 文件都必须包含标准 H5MD position value dataset。

该测试把 mock 层的 shard path/schema 约束提升到真实 HighFive backend 文件读回，覆盖 VDS 分片在文件系统和 HDF5 schema 层面的实际落盘结果。

状态：已写入测试代码，尚未执行构建或 CTest。

## 71. Observable-only provenance 与 mdout columns 路径测试补充

扩展 observable-only writer 的 mock 和 real HighFive 测试。

新增覆盖点：

- mock observable-only writer 调用 `Write_Provenance_String` 时必须创建 `/parameters/sponge/provenance`。
- mock observable-only writer 必须写入 `/parameters/sponge/provenance/launch_id`。
- real HighFive observable-only writer 必须写入 `/parameters/sponge/mdout/columns/original_name`。
- real HighFive observable-only writer 必须写入 `/parameters/sponge/mdout/columns/hdf5_name`。
- real HighFive 读回的 original column name 必须保留 SPONGE 原始列名 `TEMP`。
- real HighFive 读回的 HDF5 column name 必须为规范化名称 `temperature`。

该测试补齐 observable-only bundle 的参数区路径覆盖，尤其是轻量分析输出中的 provenance 与 mdout column mapping。

状态：已写入测试代码，尚未执行构建或 CTest。

## 72. Trajectory single-file schema 与 mdout columns 真实读回测试补充

扩展 `test_highfive_backend_io.cpp` 的 `Test_Trajectory_Writer_With_Real_Backend`。

新增覆盖点：

- trajectory 单文件 H5MD 必须写入 `/parameters/sponge/schema/name = sponge.output.h5md`。
- trajectory 单文件 H5MD 必须写入 `/parameters/sponge/schema/version = test`。
- trajectory 单文件 H5MD 必须创建 `/parameters/sponge/log`。
- trajectory 单文件 H5MD 必须写入 `/parameters/sponge/mdout/columns/original_name`。
- trajectory 单文件 H5MD 必须写入 `/parameters/sponge/mdout/columns/hdf5_name`。
- real HighFive 读回的 original column name 必须保留 `TEMP`。
- real HighFive 读回的 HDF5 column name 必须为 `temperature`。

该测试补齐单文件 trajectory output bundle 的参数区真实读回覆盖，和 observable-only 的 mdout column mapping 测试保持一致。

状态：已写入测试代码，尚未执行构建或 CTest。

## 73. Restart real-backend schema 与 extension group 真实读回测试补充

扩展 `test_highfive_backend_io.cpp` 的 `Test_Restart_Writer_With_Real_Backend`。

新增覆盖点：

- restart `.spgr.h5` 必须写入 `/parameters/sponge/schema/name = sponge.restart.h5`。
- restart `.spgr.h5` 必须写入 `/parameters/sponge/schema/version = test`。
- restart `.spgr.h5` 必须创建 `path::restart_sits` group。
- restart `.spgr.h5` 必须创建 `path::restart_meta` group。
- restart `.spgr.h5` 必须创建 `path::legacy_sidecars` group。

该测试把 restart bundle 的真实 HighFive 文件读回覆盖扩展到 schema metadata、动态 extension 根路径和 legacy sidecar provenance 根路径，补齐 mock 层路径断言对应的真实落盘验证。

状态：已写入测试代码，尚未执行构建或 CTest。

## 74. Real HighFive VDS manifest numeric fields 读回测试补充

扩展 `test_highfive_backend_io.cpp` 的 `Test_Vds_Trajectory_Writer_With_Real_Backend`。

新增覆盖点：

- 真实 VDS wrapper 必须写出 `path::shard_manifest_index`。
- 真实 VDS wrapper 必须写出 `path::shard_manifest_frame_start`。
- 真实 VDS wrapper 必须写出 `path::shard_manifest_step_start`。
- 真实 VDS wrapper 必须写出 `path::shard_manifest_step_end`。
- 真实 VDS wrapper 必须写出 `path::shard_manifest_time_start`。
- 真实 VDS wrapper 必须写出 `path::shard_manifest_time_end`。
- 两个 shard 的 index 必须分别为 0 和 1。
- 两个 shard 的 frame_start 必须分别为 0 和 1。
- 两个 shard 的 step range 必须分别为 10 和 20。
- 两个 shard 的 time range 必须分别为 0.1 和 0.2。

该测试把 mock 层已覆盖的 manifest 数值字段提升到真实 HighFive wrapper 文件读回，确保 VDS 分片索引、frame/step/time 范围完整落盘。

状态：已写入测试代码，尚未执行构建或 CTest。

## 75. Restart real-backend shared step/time link 路径读回测试补充

扩展 `test_highfive_backend_io.cpp` 的 `Test_Restart_Writer_With_Real_Backend`。

新增覆盖点：

- restart `.spgr.h5` 必须暴露 `path::position_step`。
- restart `.spgr.h5` 必须暴露 `path::position_time`。
- restart `.spgr.h5` 必须暴露 `path::velocity_step`。
- restart `.spgr.h5` 必须暴露 `path::velocity_time`。
- restart `.spgr.h5` 必须暴露 `path::box_edges_step`。
- restart `.spgr.h5` 必须暴露 `path::box_edges_time`。
- position/velocity/box 的 step 路径读回值必须为 restart step `20`。
- position/velocity/box 的 time 路径读回值必须为 restart time `1.0`。

该测试补齐 restart bundle 对 H5MD step/time 共享路径的真实读回覆盖，验证标准 particle state 通过各 element 的 step/time path 可访问，而不是只验证底层 `/particles/all/step` 和 `/particles/all/time`。

状态：已写入测试代码，尚未执行构建或 CTest。

## 76. Trajectory real-backend shared step/time link 路径读回测试补充

扩展 `test_highfive_backend_io.cpp` 的 `Test_Trajectory_Writer_With_Real_Backend`。

新增覆盖点：

- trajectory `.spg.h5md` 必须暴露 `path::position_step` 与 `path::position_time`。
- trajectory `.spg.h5md` 必须暴露 `path::velocity_step` 与 `path::velocity_time`。
- trajectory `.spg.h5md` 必须暴露 `path::force_step` 与 `path::force_time`。
- trajectory `.spg.h5md` 必须暴露 `path::box_edges_step` 与 `path::box_edges_time`。
- trajectory `.spg.h5md` 必须暴露 `/observables/all/temperature/step` 与 `/observables/all/temperature/time`。
- position/velocity/force/box/temperature 的 step 读回值必须为 `10`。
- position/velocity/force/box/temperature 的 time 读回值必须为 `0.5`。

该测试补齐 trajectory output bundle 对 H5MD element-level step/time link 的真实读回覆盖，和 restart bundle 的 step/time link 读回测试保持一致。

状态：已写入测试代码，尚未执行构建或 CTest。

## 77. Observable-only real-backend shared step/time link 路径读回测试补充

本轮补充 `tests/h5_bundle/test_highfive_backend_io.cpp` 中的
`Test_Observable_Writer_With_Real_Backend`，目标是把 observable-only
H5MD 的真实 HighFive 文件读回覆盖从 value dataset 扩展到对应的
`step/time` 数据集或 link 路径。

新增覆盖点：

- `/observables/all/temperature/step` 与 `/observables/all/temperature/time`
- `/observables/all/qc/energy/step` 与 `/observables/all/qc/energy/time`
- `/observables/all/qc/spin_square/step` 与 `/observables/all/qc/spin_square/time`
- Nose-Hoover chain coordinate/velocity 的 `step/time`
- SITS `nk` observable 的 `step/time`
- metadynamics `meta`、`rbias`、`rct` scalar 的 `step/time`
- ReaxFF `bond`、`angle` energy observable 的 `step/time`

这些断言用于确认轻量级 `*.obs.spg.h5md` 不含 `/particles` 时，仍然保留
H5MD observable 语义所需的帧索引和时间轴信息；同时覆盖 SPONGE module
mapping 写出的增强采样、QC、thermostat 与 ReaxFF 可观测量。

状态：已写入测试代码，尚未执行构建或 CTest。

## 78. Public H5MD/module path constants 契约锁定测试补充

本轮补充两个不依赖真实 HDF5 后端的路径常量契约测试：

- `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` 新增
  `Test_Public_H5MD_Path_Constants`。
- `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` 新增
  `Test_Module_Path_Constants`。

新增覆盖点包括：

- H5MD common layout path：`/h5md`、`/particles`、`/observables`、
  `/parameters`、`/parameters/sponge/*`。
- output completion/status path：`frame_count`、`last_complete_step`、
  `last_complete_time`、`status`。
- VDS shard 与 shard manifest path：`shard/status`、`frame_start`、
  `frame_count`、`last_complete_*`、`output/shard_manifest/*`。
- legacy sidecar provenance path：`/parameters/sponge/files/legacy_sidecars/*`。
- particle element path：position、velocity、force、box edges 的
  `value/step/time`。
- observable axis path：`/observables/all/step` 与
  `/observables/all/time`。
- restart state path：`/run/*`、`/parameters/restart/*`、thermostat、barostat、
  bias、SITS、metadynamics restart extension root。
- module extension path：NHC、SITS、metadynamics、QC、ReaxFF 的 canonical
  observable/parameter paths。
- SITS 动态路径拼接函数：`Sits_Nk_Root`、`Sits_Nk_Value_Path`、
  `Sits_Nk_Step_Path`、`Sits_Nk_Time_Path`。

该补充用于防止后续重构时路径常量漂移。它不替代真实 HighFive 后端测试，
而是把本轮 h5 bundle 规范中新增的公开路径字符串作为单元测试契约固定下来。

状态：已写入测试代码，尚未执行构建或 CTest。

## 79. Legacy sidecar resolution matrix 单元测试补充

本轮补充 `tests/h5_bundle/test_h5_output_plan.cpp` 中的
`Test_Resolve_Legacy_Output_Plan_Matrix`，用于锁定 H5 bundle 新输出路径启用后
legacy sidecar 的默认关闭和显式恢复规则。

新增覆盖场景：

- `Resolve_Legacy_Output_Plan(nullptr)` 必须返回 8 个 legacy sidecar，默认启用，
  且均非显式请求。
- 空 `CONTROLLER` 必须返回 8 个 legacy sidecar，默认启用，且均非显式请求。
- 设置 `output_h5_trajectory_path` 后，8 个 legacy sidecar 必须全部默认关闭。
- 设置 `output_h5_observable_path` 后，如果 8 个 legacy sidecar key 均显式设置路径，
  它们必须全部恢复启用，并保留对应显式路径。
- 8 个 legacy key 的顺序和集合固定为：`mdout`、`mdinfo`、`crd`、`box`、`vel`、
  `frc`、`rst`、`qc_scf_output`。

该测试补齐从 grouped H5 output key 到 legacy sidecar suppression/provenance
策略的完整矩阵覆盖，防止新 bundle 路径启用后意外继续写出 legacy 文件，或显式
legacy 路径被错误丢弃。

状态：已写入测试代码，尚未执行构建或 CTest。

## 80. VDS repair metadata 与 manifest string path 测试补充

本轮补充 VDS wrapper 相关的新路径测试，覆盖 `Write_Repair_Metadata` 和
`Write_Manifest_To_Wrapper` 已实现的落盘契约。

新增覆盖点：

- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 检查
    `/parameters/sponge/output/repaired_shard_count` dataset 被创建。
  - 同一测试检查 `path::shard_manifest_path` 与
    `path::shard_manifest_status` string array 的长度和 `complete` 状态。
  - `Test_Vds_Finalize_Without_Frames`、`Test_Complete_Prefix_Repair_Finalize`、
    `Test_Complete_Prefix_Repair_Retains_Valid_Prefix` 均补充
    `repaired_shard_count` dataset 存在断言。
- `tests/h5_bundle/test_highfive_backend_io.cpp`
  - `Test_Vds_Trajectory_Writer_With_Real_Backend` 检查真实 HDF5 wrapper 中
    `/parameters/sponge/output/repaired_shard_count` 存在。
  - 同一测试读回 `repaired_shard_count[0] == 0`，与 strict finalize 且未应用
    repair 的路径一致。

该补充避免只覆盖 `repair_policy` 和 `repair_status` 字符串，而遗漏实际用于记录
修复 shard 数量的 numeric metadata path。

状态：已写入测试代码，尚未执行构建或 CTest。

## 81. H5 bundle test manifest 与 README 覆盖矩阵同步

本轮同步 `tests/h5_bundle/README.md` 与 `tests/h5_bundle/TEST_TARGETS.md`，使测试入口文档反映当前已经写入的 H5 bundle 新路径单元测试覆盖。

更新内容包括：

- `test_h5_output_plan` 的职责扩展为 parser-visible key、suffix、VDS policy、完整 legacy sidecar resolution matrix。
- `test_h5md_writers_with_mock_backend` 明确覆盖 public H5MD/restart/VDS path constants。
- `test_module_h5_mappings_with_mock_backend` 明确覆盖 module path constants 与 SITS 动态路径生成。
- `test_vds_trajectory_writer_with_mock_backend` 明确覆盖 VDS repair metadata、manifest path/status arrays、wrapper diagnostics/provenance。
- `test_highfive_backend_io` 明确覆盖真实 HighFive element-level step/time、observable-only、restart、VDS wrapper/shard 和 repaired shard count 读回。
- README coverage matrix 补充 path constants、legacy 8-key matrix、observable-only step/time、trajectory/restart step/time、VDS repair metadata 等条目。
- TEST_TARGETS 增加 expanded path-contract coverage 段落，说明哪些 test target 锁定哪些公开路径契约。

该更新不新增测试断言，但使测试套件的入口文档与当前测试代码一致，避免后续维护者误判新路径覆盖范围。

状态：已写入文档，尚未执行构建或 CTest。

## 82. Output plan bundle path 组合语义与 parser helper 边界测试补充

本轮补充 `tests/h5_bundle/test_h5_output_plan.cpp` 的纯 resolver 单元测试，目标是锁定 grouped `output_h5_*` path key 与 selector key 的边界。

新增覆盖点：

- `Test_Output_Selectors_Do_Not_Enable_H5_Without_Path` 增加
  `output_h5_trajectory_repair_policy = strict`，确认 selector/config key 本身不启用 H5 output。
- 新增 `Test_H5_Output_Path_Keys_Enable_Only_Their_Bundle`：
  - `output_h5_trajectory_path` 只启用 trajectory bundle。
  - `output_h5_restart_path` 只启用 restart bundle。
  - `output_h5_observable_path` 只启用 observable bundle。
  - 任一 bundle path 启用后 legacy sidecar 默认关闭。
- 新增 `Test_All_H5_Output_Bundles_Can_Be_Enabled_Together`：
  - trajectory、restart、observable 三类 H5 bundle path 可同时启用。
  - `output_h5_trajectory_vds = off` 必须解析为 non-VDS。
  - `output_h5_trajectory_chunk_size = 20` 保持默认 chunk size 语义。
  - 三个推荐后缀 flag 均为 true。
  - 未显式设置 legacy path 时 legacy sidecar 保持关闭。
- 扩展 `Test_Bool_Parsing_Text_Variants`：
  - 未知 bool 文本如 `maybe` 和空字符串必须解析为 false。
  - `Lowercase` 必须把 repair policy 文本规范化为小写。

该补充覆盖的是 mdin/TOML flatten 后进入 resolver 的组合语义，避免 selector key 被误当作 output path，也避免未来把三种 bundle 错误处理为互斥输出。

状态：已写入测试代码，尚未执行构建或 CTest。

## 83. VDS observable axis 与 mdout column mapping 测试补充

本轮补充 VDS wrapper 中 observable 侧路径的单元测试覆盖，避免只验证 observable `value` 和 shard manifest，而遗漏 H5MD observable axis 与 SPONGE mdout column metadata。

新增覆盖点：

- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 检查 wrapper 创建 `path::observables_all_step` 与 `path::observables_all_time` virtual datasets。
  - 检查 `/observables/all/temperature/step` hard-link 到 `path::observables_all_step`。
  - 检查 `/observables/all/temperature/time` hard-link 到 `path::observables_all_time`。
  - 检查 `/parameters/sponge/mdout/columns/original_name[0] == TEMP`。
  - 检查 `/parameters/sponge/mdout/columns/hdf5_name[0] == temperature`。
- `tests/h5_bundle/test_highfive_backend_io.cpp`
  - `Test_Vds_Trajectory_Writer_With_Real_Backend` 检查真实 HDF5 wrapper 中 `/observables/all/temperature/step` 与 `/time` 存在。
  - 读回 VDS observable step 为 `10, 20`，time 为 `0.1, 0.2`。
  - 读回 VDS wrapper 的 mdout original/hdf5 column mapping。

该补充使 VDS trajectory wrapper 与 single-file trajectory、observable-only writer 在 observable `step/time` 与 mdout column metadata 上的测试覆盖保持一致。

状态：已写入测试代码，尚未执行构建或 CTest。

## 84. VDS particle axis 与 element step/time 测试补充

本轮补充 VDS wrapper 中 particle 侧 H5MD axis 与 element `step/time` path 的测试覆盖。

新增覆盖点：

- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 检查 wrapper 创建 `path::particles_all_step` 与 `path::particles_all_time` virtual datasets。
  - 检查 wrapper 创建 `path::box_edges_value` virtual dataset。
  - 检查 position、velocity、box edges 的 `step/time` hard links 指向 shared particle axis。
- `tests/h5_bundle/test_highfive_backend_io.cpp`
  - `Test_Vds_Trajectory_Writer_With_Real_Backend` 检查真实 HDF5 wrapper 中 `path::position_step`、`path::position_time`、`path::box_edges_step`、`path::box_edges_time` 存在。
  - 读回 position step 为 `10, 20`，time 为 `0.1, 0.2`。
  - 读回 box edges step 为 `10, 20`，time 为 `0.1, 0.2`。

该补充使 VDS wrapper 的 particle side coverage 与 single-file trajectory/restart 的 element-level `step/time` 覆盖对齐，避免只检查 VDS `value` dataset 而遗漏 H5MD axis 语义。

状态：已写入测试代码，尚未执行构建或 CTest。

## 85. VDS module observable axis 与 element step/time 测试补充

本轮补充 VDS wrapper 中 module observable 的 `step/time` 路径覆盖，和此前 particle/ordinary observable 的 VDS axis 测试保持一致。

新增覆盖点：

- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 检查 NHC shared `step/time` VDS。
  - 检查 NHC coordinate/velocity 的 `step/time` hard links。
  - 检查 SITS `nk` 的 `step/time` VDS。
  - 检查 metadynamics shared `step/time` VDS，及 `meta`、`rbias`、`rct` 的 hard links。
  - 检查 QC shared `step/time` VDS，及 `energy` hard links。
  - 检查 ReaxFF shared `step/time` VDS，及 `bond` hard links。
- `tests/h5_bundle/test_highfive_backend_io.cpp`
  - `Test_Vds_Trajectory_Writer_With_Real_Backend` 检查真实 HDF5 wrapper 中 NHC coordinate/velocity、SITS `nk`、metad `meta`、QC `energy`、ReaxFF `bond` 的 `step/time` 路径存在。
  - 读回上述 module 路径的 step 为 `10, 20`，time 为 `0.1, 0.2`。

该补充覆盖 VDS wrapper 的 module axis 语义，避免 VDS module 测试只读 `value` dataset 而遗漏 H5MD element axis。

状态：已写入测试代码，尚未执行构建或 CTest。

## 86. VDS conditional module branches 测试补充

本轮补充 VDS module materialization 的条件分支覆盖，重点是此前没有进入 VDS 测试的 QC spin-square、ReaxFF 多 term，以及真实 wrapper 中 metadynamics `rbias/rct` 的读回。

新增覆盖点：

- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 改为 `Ensure_Qc_Observables(true)`，覆盖 `qc_spin_square_enabled_` 分支。
  - ReaxFF terms 从 `{bond}` 扩展为 `{bond, angle}`，覆盖多 term materialization。
  - 检查 QC `spin_square/value` VDS 和 `spin_square/step,time` hard links。
  - 检查 ReaxFF `angle/value` VDS 和 `angle/step,time` hard links。
- `tests/h5_bundle/test_highfive_backend_io.cpp`
  - `Test_Vds_Trajectory_Writer_With_Real_Backend` 写入并读回 QC `spin_square` 值 `0.11, 0.22`。
  - 写入并读回 ReaxFF `angle` 值 `3.5, 4.5`。
  - 读回 metadynamics `rbias` 值 `2.0, 5.0` 和 `rct` 值 `3.0, 6.0`。
  - 检查并读回上述扩展 observable 的 `step/time` 路径。

该补充使 VDS module 测试覆盖从基础单变量路径扩展到可选/多项 module 分支，避免 future regression 只影响 VDS wrapper 的 optional module dataset 时被遗漏。

状态：已写入测试代码，尚未执行构建或 CTest。
## 87. Optional trajectory/restart particle path negative tests 补充

本轮补充了 `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` 中对可选 particle 字段的负向路径覆盖：

- `Test_Trajectory_Optional_Velocity_And_Force_Paths`：在 `Define_Particle_Datasets(atom_count, false, false)` 下只启用 position 和 box，确认不会创建 `/particles/all/velocity/*` 与 `/particles/all/force/*` 数据集。
- 同一测试确认 trajectory writer 仍会为 position 和 box 建立标准 H5MD step/time hard-link，但不会为未启用的 velocity/force 建立 hard-link。
- 同一测试确认未启用的 velocity/force 不会产生 append 写入计数，避免 optional branch 意外落盘。
- `Test_Restart_Optional_Velocity_Path`：在 `Define_Structural_State(atom_count, false)` 下只写入 coordinates 和 box，确认 restart bundle 不创建 velocity 数据集。
- 同一测试确认 restart 的 position 和 box 仍保持标准 H5MD step/time 链接，且 `/parameters/run/state_type` 仍记录为 `restart`。

该补充覆盖的是 schema 分支的“禁用路径”契约，防止未来实现把可选 velocity/force 当作总是存在的数据集写出。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 88. Observable-only bundle mock path coverage 补充

本轮继续补充 `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` 中 observable-only bundle 的 mock backend 路径契约：

- 扩展 `Test_Observable_Only_Writer`，在不创建 `/particles` group 的前提下覆盖两个 mdout observable column，确认 `/observables/all/<name>/value`、step hard-link 与 time hard-link 都写入正确位置。
- 同一测试补充 `/parameters/sponge/mdout/columns/original_name` 与 `/parameters/sponge/mdout/columns/hdf5_name` 的数组写入断言。
- 同一测试补充 `/parameters/sponge/log/mdinfo_text`、`/parameters/sponge/files/legacy_sidecars/{key,path}`、`/parameters/sponge/provenance/launch_id` 与 finalized status 的路径断言。
- 新增 `Test_Observable_Only_Module_Proxy_Paths`，通过 `ObservableH5Writer` 的 public module proxy API 覆盖 NHC、SITS nk、metadynamics scalars、QC energy、ReaxFF energy term、metad diagnostics 与 `qc_scf_output` 在 observable-only 文件中的路径落点。
- 对 SITS nk 保持当前实现契约：`/observables/all/sits/<module>/nk/{step,time,value}` 是独立 dataset，而不是 value 子组 hard-link。
- 对 QC optional spin-square 保持负向契约：`Ensure_Qc_Observables(false)` 时只写 `/observables/all/qc/energy/value`，不创建 `/observables/all/qc/spin_square/value`。

该补充提升了 `.obs.spg.h5md` 轻量分析输出的 schema 覆盖，尤其是确认 observable-only 文件可以承载 module observable 与 parameters/sponge 下的日志、legacy 和 provenance 元数据。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 89. HighFive optional particle-field disabled branches 补充

本轮补充 `tests/h5_bundle/test_highfive_backend_io.cpp` 中真实 HighFive/HDF5 文件层的 optional particle-field 负向路径覆盖：

- 新增 `Test_Trajectory_Optional_Particle_Fields_With_Real_Backend`，使用 `TrajectoryH5Writer + HighFiveBackend` 创建单文件 `.spg.h5md`，在 `Define_Particle_Datasets(atom_count, false, false)` 下写入 position 和 box。
- 该测试读回真实 HDF5 文件，确认 position 与 box 的 `value/step/time` 路径存在，且 `/particles/all/velocity/{value,step,time}` 与 `/particles/all/force/{value,step,time}` 均不存在。
- 同一测试读回 position dataset shape、position value、position/box step-time、completion frame count、last complete step 和 finalized status。
- 新增 `Test_Restart_Optional_Velocity_With_Real_Backend`，使用 `RestartH5Writer + HighFiveBackend` 创建 `.spgr.h5`，在 `Define_Structural_State(atom_count, false)` 下写入单帧 structural state。
- 该测试读回真实 HDF5 restart 文件，确认 position 与 box 的 `value/step/time` 路径存在，velocity 的 `value/step/time` 均不存在。
- 同一测试读回 `/run/current_step`、`/run/current_time`、`/run/state_type`、completion metadata、position/box step-time 和 finalized status。

该补充把此前 mock backend 层的 optional disabled 契约推进到真实 HighFive 后端层，重点防止实际 HDF5 文件中误创建禁用的 velocity/force dataset 或 hard-link。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 90. VDS optional particle-field disabled branches 补充

本轮补充 VDS trajectory wrapper/shard 层的 optional particle-field 负向路径覆盖：

- 扩展 `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` 的 `Test_Vds_Wrapper_And_Module_Virtual_Datasets`，在 `include_velocity=true, include_force=false` 下确认 wrapper 不 materialize `/particles/all/force/value`，也不建立 force step/time hard-link。
- 同一测试确认每个 mock shard 不创建 force value dataset，避免 shard 与 wrapper schema 不一致。
- 新增 `Test_Vds_Optional_Particle_Fields_Disabled`，在 `Define_Particle_Datasets(atom_count, false, false)` 下写入两个 shard，并确认 wrapper 只 materialize particle step/time、position value 与 box edges value。
- 新测试确认 wrapper 不 materialize velocity/force value virtual dataset，也不建立 velocity/force step/time hard-link。
- 新测试确认每个 shard 只包含 position/box 数据集，不包含 velocity/force 数据集，并确认 manifest status 仍为 complete。
- 扩展 `tests/h5_bundle/test_highfive_backend_io.cpp` 的 `Test_Vds_Trajectory_Writer_With_Real_Backend`，在真实 HighFive VDS wrapper 中确认 velocity/force value/step/time 路径不存在。
- 同一真实后端测试继续打开两个 shard 文件，确认 shard 侧 position/box 存在而 velocity/force value/step/time 不存在。

该补充覆盖了 VDS writer 的 `include_velocity_` 与 `include_force_` finalize materialization 分支，防止分片输出场景中禁用字段被错误写入 wrapper virtual dataset 或 shard dataset。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 91. H5 output route helper contract tests 补充

本轮补充 `MD_core/output.hpp` 集成入口中可单元测试的 H5 输出路由 helper 覆盖：

- 新增 `SPONGE/utils/h5md/output_route_helpers.hpp`，把原先位于 `output.hpp` 匿名 namespace 中的纯 helper 逻辑抽出为可测试 inline 函数。
- `output.hpp` 保留原 helper 名称 wrapper，内部委托到 `SpongeH5OutputRoute`，以降低对现有调用点的扰动。
- 新增 `tests/h5_bundle/test_output_route_helpers.cpp`，覆盖 SPONGE mdout key 到 H5MD-safe observable name 的转换规则。
- 同一测试覆盖 sanitized name 发生 collision 后的 `_1`、`_2`、`_3` suffix 规则。
- 同一测试覆盖 mdout scalar string 到 double 的严格解析规则，包括前后空白、科学计数法、空字符串、非数字、尾随非空白字符和 null output pointer。
- 同一测试覆盖 ReaxFF 输出 key 识别规则：`REAXFF` 和 `REAXFF_*` 进入 H5 ReaxFF module routing，大小写不做兼容放宽。
- 同一测试覆盖 output key lookup 语义与 optional text sidecar read 语义，包括空路径、缺失文件、null text pointer 和 binary/text 内容保留。
- `tests/h5_bundle/CMakeLists.txt` 新增 `test_output_route_helpers`，并纳入 `sponge_h5_bundle_contract_tests`。

该补充锁住了数据进入 writer 前的路径命名和路由语义，尤其是 observable HDF5 dataset name 的稳定性。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 92. HighFive backend factory and repeated dataset contract 补充

本轮补充 `tests/h5_bundle/test_highfive_backend_io.cpp` 中真实 HighFive backend 底层行为覆盖：

- 新增 `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics`，通过 `HighFiveBackendFactory::Create_Backend()` 创建 backend，并用 `H5MDWriter` 打开嵌套路径 `nested/output/factory.spg.h5md`。
- 该测试确认 factory 创建的 backend 可以完整执行 common H5MD layout 初始化、schema metadata 写入、append、finalize 和 close。
- 该测试覆盖 appendable dataset 的维度增长：连续 append 三条 width=2 的记录后，真实 HDF5 dataset shape 应为 `[3,2]`。
- 该测试覆盖重复 `Create_Dataset` 调用语义：在同一路径已有数据后再次定义同一 dataset，不应删除或截断已 append 的数据，后续 append 应继续增长。
- 该测试读回真实 HDF5 文件中的 dataset 值 `[1,2,3,4,5,6]`，确认记录顺序和内容保留。
- 该测试确认嵌套输出路径实际创建，且 `/parameters/sponge/output/status` 最终为 `finalized`。

该补充锁住 writer facade 之下的 backend 基础契约，避免后续 HighFive backend 重构破坏 VDS shard 创建、append shape 或重复 dataset definition 的幂等行为。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 93. HighFive string-array metadata overwrite coverage 补充

本轮继续扩展 `tests/h5_bundle/test_highfive_backend_io.cpp` 中真实 HighFive backend 的 string-array metadata 覆盖：

- 扩展 `Test_HighFive_Backend_String_Overwrite`，新增嵌套 metadata path `/parameters/sponge/test/nested/string_array`。
- 该测试先写入包含空字符串元素、带空格字符串和相对路径字符串的 string array，再在同一路径写入 replacement array。
- 读回真实 HDF5 文件时确认 replacement array 生效，旧数组内容被删除，说明 `Write_String_Array` 的 `Delete_If_Exists` 路径没有留下旧 dataset。
- 该测试确认 replacement array 可以保留空字符串元素，覆盖 legacy sidecar、manifest path/status、mdout columns 等 metadata 未来可能遇到的空值边界。
- 该测试同时覆盖 nested parent group 自动创建路径，避免 metadata 子树新增时需要预先手动创建 group。

该补充锁住 SPONGE H5 metadata string-array 的覆盖写入语义，避免 HighFive backend 后续重构导致旧数组残留、嵌套路径失败或空字符串元素丢失。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 94. HighFive VDS zero-frame finalize coverage 补充

本轮补充 `tests/h5_bundle/test_highfive_backend_io.cpp` 中真实 HighFive backend 的 VDS zero-frame finalize 覆盖：

- 新增 `Test_Vds_Finalize_Without_Frames_With_Real_Backend`，使用 `VdsTrajectoryH5Writer + HighFiveBackendFactory` 打开 VDS wrapper 后不写入任何 frame，直接 finalize。
- 测试确认 writer manifest size、total trajectory frame count、total observable frame count 均为 0。
- 读回真实 wrapper 文件，确认 common H5MD roots `/h5md`、`/particles`、`/observables` 存在。
- 测试确认 zero-frame wrapper 不创建 particle step/time、position、box、observable step/time，也不创建 shard manifest index/path dataset。
- 测试确认未创建 shard root 目录，避免空输出段产生空 shard 传输负担。
- 测试确认 wrapper 仍写入 trajectory chunk size、repair policy/status、repaired shard count、VDS status、completion metadata 和 finalized output status。

该补充锁住分段运行中“本段无 trajectory frame”的合法输出行为，防止 zero-frame finalize 误创建空 shard、空 VDS dataset 或遗漏必要 metadata。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 95. Restart metad text component coverage 补充

本轮补充 restart bundle 中 metadynamics text state component 的路径覆盖：

- 扩展 `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` 的 `Test_Restart_Module_State_And_Legacy_Provenance`。
- mock backend 层现在覆盖 `/parameters/restart/bias/meta/meta0/{hills,history,edge,potential_export,direct_export}` 五类 text component。
- 扩展 `tests/h5_bundle/test_highfive_backend_io.cpp` 的 `Test_Restart_Writer_With_Real_Backend`。
- 真实 HighFive restart 文件现在检查上述五类 metad component path 均存在，并读回对应文本内容。
- 该覆盖对应 `MD_core/output.hpp` 中 `Write_H5_Restart_Text_File_If_Present` 会写出的 restart metad state component，而不仅是此前覆盖的 `hills`。

该补充锁住 `.spgr.h5` 中 metad restart text state 的完整 component path，避免后续只保留 myhill/hills 而遗漏 history、edge、potential/direct export 等重启相关状态。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 96. Restart SITS dynamic state component coverage 补充

本轮补充 restart bundle 中 SITS dynamic state component 的路径覆盖：

- 扩展 `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` 的 `Test_Restart_Module_State_And_Legacy_Provenance`。
- mock backend 层现在覆盖同一 module 下的 `/parameters/restart/bias/sits/sits_a/nk` 与 `/parameters/restart/bias/sits/sits_a/weight` 两个 state component。
- 测试确认两个 component dataset 都被创建，且 append count 分别匹配输入长度。
- 扩展 `tests/h5_bundle/test_highfive_backend_io.cpp` 的 `Test_Restart_Writer_With_Real_Backend`。
- 真实 `.spgr.h5` 文件现在检查 `sits_a/nk` 与 `sits_a/weight` 都存在，并读回 `weight` payload 的 shape 与数值。

该补充锁住 `RestartH5Writer::Write_Sits_State(module_name, state_name, ...)` 的动态 state_name 路径语义，避免 restart bundle 只能覆盖固定 `nk` 数据集而遗漏未来 SITS/adaptive state component。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 97. HighFive repeated hard-link coverage 补充

本轮补充真实 HighFive backend 的 hard-link 幂等语义覆盖：

- 扩展 `tests/h5_bundle/test_highfive_backend_io.cpp` 的 `Test_HighFive_Backend_Basic_File_Layout`。
- 同一路径 `/observables/all/temperature_alias` 连续调用两次 `Create_Hard_Link`，确认第二次调用成功。
- 读回真实 HDF5 文件时检查 alias dataset 存在、shape 与原 value dataset 一致，并可读回原始 temperature 值。

该补充锁住 H5MD step/time link 机制依赖的 backend 幂等行为，避免重复建链时误报失败或破坏已有 link。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 98. HighFive VDS complete-prefix repair coverage 补充

本轮补充真实 HighFive backend 的 VDS complete-prefix repair 覆盖：

- 新增 `MaybeFailFinalizeHighFiveBackend` 与 `SelectiveFailHighFiveBackendFactory`，用于在真实 HighFive backend 上按 backend 创建序号注入 finalize 失败。
- 新增 `Test_Vds_Complete_Prefix_Repair_With_Real_Backend`，构造 wrapper、完整 shard0 和 finalize 失败的 shard1。
- 测试调用 `VdsTrajectoryH5Writer::Finalize_With_Repair()`，确认 writer manifest 被修剪为 1 个完整 shard，trajectory frame count 回退为 1。
- 读回真实 wrapper `.spg.h5md`，确认 virtual position/box dataset 只包含 shard0 的一帧，step/time 回退到 step 10/time 0.1。
- 读回 shard manifest，确认只保留 `segment_000000.spg.h5md`，status 为 `complete`，不会把失败 shard 写入 wrapper manifest。
- 读回 output completion metadata，确认最终 `frame_count/last_complete_step/last_complete_time` 与完整前缀一致。
- 读回 repair metadata，确认 `repair_policy = complete_prefix`、`repair_status = applied`、`repaired_shard_count = 1`。
- 读回 shard0 文件，确认完整 shard 仍是可读的真实 HDF5 shard，并保持 velocity/force disabled 分支。

该补充把 VDS repair 从 mock backend 行为推进到真实 HDF5 wrapper/shard 读回层，锁住云端分段输出场景中“失败尾段不进入最终 VDS”的关键契约。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 99. H5MD common layout roots/schema/output metadata coverage 补充

本轮补充 mock backend 层 `H5MDWriter` common layout 覆盖：

- 将 `Test_Common_Layout_Observable_Only` 扩展并重命名为 `Test_Common_Layout_Roots_And_Output_Metadata`。
- 普通 trajectory-style H5MD 打开路径现在检查 `/h5md`、`/particles`、`/observables`、`/parameters`、`/parameters/sponge`、`/parameters/sponge/schema`、`/parameters/sponge/output` 都由 common layout 创建。
- observable-only 打开路径现在检查同一 common metadata root，同时确认不会创建 `/particles`。
- 两种模式都检查 schema name/version 写入，以及 output status 初始为 `open`。
- 两种模式都检查 `frame_count`、`last_complete_step`、`last_complete_time` 三个 output completion dataset 的初始化 append。
- 测试表达尽量使用 `path::*` public constants，减少 literal path 漂移风险。

该补充锁住 writer facade 的最底层 common layout 契约，避免后续只在具体 trajectory/restart/observable writer 中间接覆盖 root/schema/output metadata。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 100. H5 output plan helper null/unknown key coverage 补充

本轮补充 parser-visible H5 output plan helper 的边界覆盖：

- 扩展 `Test_Contract_Helper_Edge_Cases`，直接覆盖 `Command_String(&controller, nullptr)` 返回空字符串。
- 扩展 `Test_Resolve_Legacy_Output_Plan_Matrix`，在 legacy default enabled、空 controller、有 H5 bundle 禁用 legacy default、全部 explicit legacy sidecar 四种场景中都检查 `Enabled(nullptr)` 与 `Explicitly_Requested(nullptr)`。
- 同一矩阵同时检查 unknown legacy key 对 `Enabled` 与 `Explicitly_Requested` 均返回 false。
- 该覆盖对应 `LegacyOutputPlan` 的公共 helper 分支，避免后续为未知 legacy sidecar 或 null key 引入错误的 permissive 行为。

该补充锁住 legacy sidecar gating 的边界语义，尤其是启用 H5 bundle 后只允许显式 legacy sidecar 回写的规则不会被 unknown/null key 绕过。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 101. H5 output route name global uniqueness coverage 补充

本轮补充并修正 H5 output routing helper 的 observable name uniqueness 契约：

- 修正 `Make_Unique_Output_Names`，从按 sanitized base 计数改为维护全局 `used_names`，避免已存在 suffix 名称被后续 collision 重复生成。
- 扩展 `Test_Unique_Output_Name_Collision_Resolution`，覆盖 `A`、`A_1`、`A`、`A-1`、`A` 这类用户输出名已经包含 suffix 的情况。
- 新测试期望输出为 `A`、`A_1`、`A_2`、`A_1_1`、`A_3`，锁住所有 HDF5 observable dataset name 全局唯一。
- 扩展 scalar parsing 测试，确认纯空白字符串不被解析为数值，且不会改写旧 value。
- 扩展 ReaxFF key 测试，确认 `REAXFFBOND` 不会被误识别为 ReaxFF module term，只有 `REAXFF` 或 `REAXFF_` 前缀进入该路由。
- 扩展 output key lookup 与 optional text sidecar read 测试，覆盖 empty key vector 与 null filename 边界。

该补充修复了一个会导致 mdout column 映射到重复 HDF5 observable path 的实际风险，并用单元测试锁住路径命名全局唯一性。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 102. Explicit legacy sidecar provenance collection coverage 补充

本轮补充 H5 bundle legacy sidecar provenance routing 覆盖：

- 新增 `SpongeH5OutputPlan::Collect_Explicit_Legacy_Sidecars`，将四个 `Write_H5_Legacy_Sidecar_Provenance` overload 中重复的 explicit sidecar 筛选逻辑抽成公共 helper。
- 四个 runtime writer 路径现在共用同一 helper，避免 trajectory、VDS、observable-only、restart 之间的 legacy provenance 筛选规则漂移。
- 新增 `Test_Explicit_Legacy_Sidecar_Collection`，确认 legacy default enabled 时不会把默认 legacy sidecar 写入 H5 provenance。
- 同一测试确认启用 H5 bundle 后，只有显式路径的 legacy sidecar 会被收集。
- 同一测试确认输出顺序保持 canonical legacy sidecar 顺序，而不是用户设置顺序：`mdout`、`rst`、`qc_scf_output`。
- 同一测试覆盖 null output vector pointer 时 helper 不改写已有 paths。

该补充锁住“bundle 后默认不写 legacy，除非显式指定；显式指定则写入 provenance”的核心契约，并让四类 H5 writer 使用同一套可测试逻辑。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 103. Writer dataset type/shape/chunk contract coverage 补充

本轮补充 mock backend 层 DatasetSpec 级别的 type/shape/chunk 覆盖：

- 新增 `Require_Dataset_Spec` helper，用于断言 dataset path、`DataType`、初始 shape、max shape、chunk shape 和 appendable flag。
- 扩展 single-file trajectory writer 测试，覆盖 particle step/time、position、velocity、force、box edges 以及 observable scalar dataset 的 type/shape/chunk。
- 扩展 observable-only writer 测试，覆盖 observable step/time 与 energy/temperature scalar dataset 的 type/shape/chunk。
- 扩展 restart writer 测试，覆盖 single-state structural step/time、position、velocity、box edges、`/run/current_step`、`/run/current_time` 的 type/shape/chunk。
- 扩展 restart module state 测试，覆盖 NHC state 与 SITS dynamic `nk/weight` state 的 type/shape/chunk。

该补充把测试从“路径存在”推进到“路径上的 dataset schema 正确”，能更早发现 H5 bundle schema 维度、chunk 或类型漂移。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 104. VDS wrapper virtual dataset schema coverage 补充

本轮补充 mock backend 层 VDS wrapper virtual dataset 的 DatasetSpec 覆盖：

- 扩展 `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`，为 VDS 专用测试加入 `Require_Dataset_Spec` helper。
- `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 现在覆盖 wrapper 中 particle step/time、position、velocity、box edges 的 virtual dataset type/shape/chunk/appendable 契约。
- 同一测试覆盖 observable stream 在 VDS wrapper 中的 step/time 与 scalar value virtual dataset 契约。
- 同一测试覆盖 NHC、SITS `nk`、metadynamics scalar、QC scalar、ReaxFF term scalar 的 module VDS wrapper dataset 契约。
- `Test_Vds_Optional_Particle_Fields_Disabled` 现在覆盖 velocity/force disabled 分支下 wrapper 仍正确创建 position/box virtual dataset specs，且不创建被禁用字段。
- 所有 VDS wrapper dataset 都固定为 total-frame shape、max shape、chunk shape 相同，并设置 `appendable=false`，与 shard 内部 appendable stream 明确区分。

该补充锁住分段 H5MD wrapper 的 schema 语义，避免后续 VDS materialization 只保留路径存在但悄悄改变 type、shape、chunk 或 appendable flag。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 105. Module mapping dataset schema coverage 补充

本轮补充 module mapping mock backend 层的 DatasetSpec 覆盖：

- 扩展 `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp`，为 module mapping 测试加入 `Require_Dataset_Spec` helper。
- `Test_Nhc_And_Sits_Mappings` 现在覆盖 NHC step/time、coordinate、velocity，以及 SITS `nk` step/time/value 的 type/shape/chunk/appendable 契约。
- `Test_Metadynamics_And_Diagnostics` 现在覆盖 metadynamics step/time 与 `meta/rbias/rct` scalar value dataset 契约。
- `Test_Qc_And_Reaxff_Mappings` 现在覆盖 QC step/time、energy、spin_square，以及 ReaxFF step/time、bond、angle scalar value dataset 契约。
- `Test_Qc_Optional_Spin_Square_Path` 现在覆盖 spin_square disabled 分支下 QC energy schema 正确且不创建 spin_square value dataset。

该补充把 module extension 测试从“路径和 hard link 存在”推进到“路径上的 dataset schema 正确”，与 trajectory、restart、observable-only 和 VDS wrapper 的 schema 级测试保持一致。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 106. Output error metadata path constant coverage 补充

本轮补充 H5 output failure metadata path 的 public constant 覆盖：

- 在 `SPONGE/utils/h5md/h5md_writer.hpp` 的 `SpongeH5MD::path` namespace 中新增 `output_error = "/parameters/sponge/output/error"`。
- `OutputCompletionTracker::Mark_Failed` 现在通过 `path::output_error` 写入失败原因，避免 tracker 内部 literal path 漂移。
- `Test_Public_H5MD_Path_Constants` 现在显式锁住 `path::output_error` 的路径值。
- `test_completion_tracker`、`test_h5md_writers_with_mock_backend` 和 `test_highfive_backend_io` 中的失败原因断言改为使用 `path::output_error`。

该补充把失败状态元数据路径纳入与 `output_status/frame_count/last_complete_*` 相同的 public path constant 契约。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 107. HighFive backend max-dims/chunk layout coverage 补充

本轮补充真实 HighFive backend 的 dataset layout 读回覆盖：

- 扩展 `tests/h5_bundle/test_highfive_backend_io.cpp`，新增 `Dataset_Max_Dims` 和 `Dataset_Chunk_Dims` helper，通过 HDF5 C API 读回真实 dataset 的 max dims 与 chunk dims。
- 扩展 `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics`。
- appendable dataset `{{0,2},{0,0},{0,2}}` 写入三帧后，测试确认真实 HDF5 当前 shape 为 `[3,2]`。
- 同一 appendable dataset 现在确认 frame 维 max dim 为 `H5S_UNLIMITED`，第二维 max dim 为 `2`。
- 同一 appendable dataset 现在确认 zero chunk 经过 backend normalization 后落盘 chunk 为 `[1,2]`。
- 新增 fixed dataset `{{2,2},{2,2},{2,2}}`，确认真实 HDF5 当前 shape、max dims、chunk dims 都为 `[2,2]`。

该补充把 mock 层 DatasetSpec 断言推进到真实 HDF5 layout 层，锁住 `HighFiveBackend::Normalize_Spec`、appendable unlimited frame 维和 fixed dataset layout 的关键语义。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 108. Empty H5 output path contract coverage 补充

本轮补充 H5 output path 空值语义覆盖：

- 新增 `SpongeH5OutputContract::Command_Has_Non_Empty_Value`，将 bundle 启用条件从“path key 存在”收紧为“path key 存在且 value 非空”。
- `Any_H5_Output_Enabled` 现在使用该 helper，避免空 `output_h5_*_path` 错误关闭默认 legacy sidecar 输出。
- `Resolve_Output_Plan` 的 trajectory、restart、observable 三个 bundle enabled 判断同步改为非空 path 判断。
- 新增 `Test_Empty_H5_Output_Paths_Do_Not_Enable_Bundles`，覆盖三个 `output_h5_*_path = ""` 同时存在时不启用任何 H5 bundle。
- 同一测试确认 selector 类参数仍可被解析：`output_h5_trajectory_vds=true` 与 `output_h5_trajectory_chunk_size=20` 保留在 plan 中，但不单独启用 trajectory H5。
- 同一测试确认空 H5 path 不会禁用默认 legacy sidecars，保持“显式有效 bundle path 才切换 canonical H5 输出”的契约。

该补充消除了 mdin/TOML flattening 后空字符串路径造成的歧义：空 path 不是有效 bundle 选择，也不能作为关闭 legacy 输出的信号。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 109. SPONGE parameter path constant coverage 补充

本轮补充 H5 bundle 参数区路径的 public constant 与测试覆盖：

- 在 `SpongeH5MD::path` 中新增 `sponge_mdout`、`sponge_log`、`sponge_files`、`sponge_provenance` 四个参数区 root constant。
- 在 `SpongeH5MD::path` 中新增 `mdout_columns_original_name` 与 `mdout_columns_hdf5_name` 两个 mdout column leaf constant。
- `TrajectoryH5Writer` 和 `ObservableH5Writer` 的 base layout 与 mdout column metadata 写入改为使用这些 constants。
- `RestartH5Writer` 和 `VdsTrajectoryH5Writer` 的 legacy sidecar/provenance 相关路径同步改为 public constants。
- `Test_Public_H5MD_Path_Constants` 现在显式锁住新增 root/leaf path literal。
- mock writer、VDS wrapper 与 HighFive real-backend 测试中的 mdout/log/files/provenance 断言改为使用 public constants，减少路径 literal 漂移。

该补充把 `/parameters/sponge/mdout`、`/parameters/sponge/log`、`/parameters/sponge/files`、`/parameters/sponge/provenance` 以及 mdout column leaf paths 纳入与其他 H5 bundle 新路径相同的 contract 测试层。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 110. VDS output metadata path constant coverage 补充

本轮补充 VDS output metadata 路径的 public constant 与测试覆盖：

- 在 `SpongeH5MD::path` 中新增 `output_trajectory_chunk_size`、`output_vds_status`、`output_repair_policy`、`output_repair_status`、`output_repaired_shard_count`。
- `VdsTrajectoryH5Writer` 写入 trajectory chunk size、VDS status、repair policy/status 和 repaired shard count 时改为使用这些 public constants。
- `Test_Public_H5MD_Path_Constants` 显式锁住上述五个路径 literal。
- VDS mock backend 测试中的 wrapper metadata 断言改为通过 public constants 访问。
- HighFive real-backend VDS 测试中的 metadata existence/readback 断言同步改为通过 public constants 访问。

该补充将 VDS wrapper 的控制/恢复元数据路径纳入与 completion metadata、manifest metadata 相同的 path contract 层，避免后续 VDS repair 或 chunking 相关路径漂移。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 111. Schema leaf and particle group path coverage 补充

本轮补充 schema leaf 与 particle component group 路径覆盖：

- 在 `SpongeH5MD::path` 中新增 `sponge_schema_name` 与 `sponge_schema_version`。
- 在 `SpongeH5MD::path` 中新增 `particles_all_position`、`particles_all_velocity`、`particles_all_force`、`particles_all_box`、`particles_all_box_edges`。
- `TrajectoryH5Writer` 与 `RestartH5Writer` 的 base layout group 创建改为使用 particle component group constants。
- `Test_Public_H5MD_Path_Constants` 显式锁住新增 schema leaf 与 particle group path literal。
- mock writer 与 HighFive backend 测试中的 schema/path group 断言同步使用 public constants。

该补充覆盖了 previously 只通过 value/step/time dataset constants 间接约束的 group-level layout，避免 base layout group 名称和 dataset leaf 名称在后续重构中漂移。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 112. VDS particle virtual dataset path constant coverage 补充

本轮补充 VDS particle virtual dataset 物化路径的 constant 使用覆盖：

- `VdsTrajectoryH5Writer::Materialize_Particle_Virtual_Datasets` 中 particle step/time、position、velocity、force、box edges 的 virtual dataset path 改为使用 `SpongeH5MD::path` public constants。
- 同一函数中的 step/time hard link target 也改为使用 `path::position_step/time`、`path::velocity_step/time`、`path::force_step/time`、`path::box_edges_step/time`。
- 现有 `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 已覆盖这些 wrapper virtual dataset paths、hard links、type/shape/chunk 和 source dims。
- 现有 `Test_Vds_Optional_Particle_Fields_Disabled` 已覆盖 velocity/force disabled 分支下 public particle path constants 不被创建。
- 窄范围文本检查确认 `vds_trajectory_h5_writer.hpp` 中不再残留 `/particles/all/*` particle virtual dataset literal。

该补充把 VDS wrapper particle 路径与 single-file trajectory/restart writer 的 public path constants 对齐，避免同一 H5MD particle layout 在不同输出模式中漂移。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 113. Writer open precondition coverage 补充

本轮补充 H5 writer 入口绑定前置条件覆盖：

- 新增 `Test_Writer_Open_Preconditions_Reject_Unbound_Bundles`。
- `TrajectoryH5Writer::Open_Single_File` 现在有负向测试覆盖：未启用 trajectory plan 时必须拒绝打开。
- 同一测试覆盖：trajectory plan 处于 VDS mode 时，single-file trajectory writer 必须拒绝打开。
- `ObservableH5Writer::Open` 现在有负向测试覆盖：即使 trajectory path 已启用，只要 observable bundle 未启用就必须拒绝打开。
- `RestartH5Writer::Open` 现在有负向测试覆盖：即使 observable path 已启用，只要 restart bundle 未启用就必须拒绝打开。
- 每个分支同时检查 `Last_Error()` 的具体错误文本，确保错误路径清晰且不会 silent fallback 到错误 writer。

该补充锁住 bundle path plan 到 writer facade 的入口边界，避免后续新增路径时 trajectory、observable、restart writer 误接收其他 bundle 的 path plan。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 114. Bare H5MD writer no-backend coverage 补充

本轮补充裸 `H5MDWriter` 未绑定 backend 时的安全行为覆盖：

- 新增 `Test_H5MD_Writer_No_Backend_Is_Safe`。
- 测试确认 `Is_Attached()` 为 false，`Status()` 为 `FileStatus::closed`。
- 测试覆盖 `Open`、`Flush`、`Close`、`Finalize` 在未绑定 backend 时均返回 false。
- 测试覆盖 group、dataset、virtual dataset、hard link、int64 append、float64 append、string、string array、status、failure metadata、output completion 等写入入口均返回 false。
- 测试确认 `Last_Error()` 返回 `H5MD writer backend is not attached`。

该补充锁住 H5 bundle writer facade 的最底层防御边界，避免后续新增 H5 path 写入入口在 backend 缺失时 silent succeed 或崩溃。本文档仅记录测试代码更新；本轮未执行构建、单元测试或 CTest。

## 115. Module dynamic leaf path builder coverage 补充

本轮补充 module extension 动态 leaf path builder 覆盖：

- 在 `module_h5_mappings.hpp` 中新增 metadynamics scalar path builders：`Metadynamics_Scalar_Root`、`Metadynamics_Scalar_Value_Path`、`Metadynamics_Scalar_Step_Path`、`Metadynamics_Scalar_Time_Path`。
- 新增 metadynamics diagnostic path builders：`Metadynamics_Diagnostic_Root`、`Metadynamics_Diagnostic_Path`。
- 新增 QC observable path builders：`Qc_Observable_Root`、`Qc_Observable_Value_Path`、`Qc_Observable_Step_Path`、`Qc_Observable_Time_Path`、`Qc_Scf_Output_Path`。
- 新增 ReaxFF term path builders：`Reaxff_Term_Root`、`Reaxff_Term_Value_Path`、`Reaxff_Term_Step_Path`、`Reaxff_Term_Time_Path`。
- module writer 内部 metadynamics、QC、ReaxFF 的 value/step/time/scf/diagnostic 路径改为使用这些 helpers。
- `Test_Module_Path_Constants` 现在直接锁住这些动态 helper 的 literal 输出。
- module mapping 行为测试中的 metadynamics/QC/ReaxFF dataset、hard-link、string metadata 断言改为通过 helper 访问。

该补充把原先分散在生产代码和测试代码中的字符串拼接集中为可测试 path builder，降低 module extension 新路径漂移风险。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 116. VDS module dynamic leaf helper alignment 补充

本轮补充 VDS wrapper 对 module dynamic leaf path helper 的复用覆盖：

- `VdsTrajectoryH5Writer::Materialize_Metadynamics_Virtual_Datasets` 改为使用 `Metadynamics_Scalar_Root`。
- `VdsTrajectoryH5Writer::Materialize_Qc_Virtual_Datasets` 改为使用 `Qc_Observable_Root`。
- `VdsTrajectoryH5Writer::Materialize_Reaxff_Virtual_Datasets` 改为使用 `Reaxff_Term_Root`。
- `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 中 metadynamics/QC/ReaxFF wrapper virtual dataset existence 和 DatasetSpec 断言改为通过 helper 访问。
- 这样 single-file writer、observable-only proxy、module mock test 和 VDS wrapper 对同一组 module leaf path 使用同一个 path-builder 契约。

该补充避免 VDS wrapper 在 module extension 路径上继续保留独立字符串拼接，降低分段输出与单文件输出之间的 schema drift 风险。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 117. HighFive module dynamic helper readback coverage 补充

本轮补充 HighFive real-backend 测试对 module dynamic leaf path helper 的 readback 覆盖：

- `test_highfive_backend_io` 中 metadynamics scalar value/step/time readback 改为使用 `Metadynamics_Scalar_*_Path` helpers。
- metadynamics diagnostic group/string readback 改为使用 `Metadynamics_Diagnostic_Root` 与 `Metadynamics_Diagnostic_Path`。
- QC energy、spin_square value/step/time readback 改为使用 `Qc_Observable_*_Path` helpers。
- QC SCF text readback 改为使用 `Qc_Scf_Output_Path`。
- ReaxFF bond、angle、over 等 term readback 改为使用 `Reaxff_Term_*_Path` helpers。
- 窄范围文本检查确认 `test_highfive_backend_io.cpp` 中不再残留 metadynamics/QC/ReaxFF dynamic root 的旧式字符串拼接。

该补充使 dynamic helper 契约同时被 mock mapping、VDS wrapper 和真实 HDF5 文件 readback 覆盖，避免生产代码和测试代码分别维护独立路径拼接。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 118. Restart dynamic state path builder coverage 补充

本轮补充 restart bundle 中动态 state 路径的 helper 与测试覆盖：

- 在 `h5md_writer.hpp` 中新增 `Restart_Sits_State_Root` 与 `Restart_Sits_State_Path`，用于 `/parameters/restart/bias/sits/<module>/<state>`。
- 在 `h5md_writer.hpp` 中新增 `Restart_Metad_State_Root` 与 `Restart_Metad_State_Path`，用于 `/parameters/restart/bias/meta/<name>/<component>`。
- `RestartH5Writer::Write_Sits_State` 改为通过 SITS restart helper 创建 module group 和 state dataset。
- `RestartH5Writer::Write_Metad_State_Text` 改为通过 metad restart helper 创建 group 和 text component。
- `Test_Public_H5MD_Path_Constants` 新增 restart dynamic helper literal 断言。
- `Test_Restart_Module_State_And_Legacy_Provenance` 中 SITS/metad restart extension 的 dataset/group/string/append-count 断言改为通过 helper 访问。
- `test_highfive_backend_io` 中真实 `.spgr.h5` readback 的 SITS state dataset 与 metad text component 路径改为通过 helper 访问。

该补充把 restart bundle 的动态扩展路径与 metadynamics/QC/ReaxFF observable dynamic helper 放在同一层 path contract 下，避免 writer 和测试继续分散维护字符串拼接。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 119. Ordinary observable dynamic path builder coverage 补充

本轮补充普通 observable 动态路径的 helper 与测试覆盖：

- 在 `h5md_writer.hpp` 中新增 `Observable_Root`、`Observable_Value_Path`、`Observable_Step_Path`、`Observable_Time_Path`，用于 `/observables/all/<name>/value|step|time`。
- `TrajectoryH5Writer::Define_Observable_Stream` 与 `Append_Observable_Frame` 改为通过普通 observable helper 创建 group、value dataset、step/time hard link 和 append path。
- `ObservableH5Writer::Define_Observable_Stream` 与 `Append_Observable_Frame` 改为通过同一组 helper 写 observable-only 文件。
- `VdsTrajectoryH5Writer` wrapper 普通 observable virtual dataset 和 step/time hard link materialization 改为通过同一组 helper。
- `Test_Public_H5MD_Path_Constants` 新增 ordinary observable helper literal 断言。
- mock writer、VDS mock writer 和 HighFive real-backend 测试中的 `temperature`、`energy`、`pressure` 普通 observable value/step/time 断言改为通过 helper 访问。
- 窄范围文本检查确认生产 writer 中普通 observable schema 路径不再残留独立 `path::observables_all + name` 拼接，除 helper 自身外。

该补充把普通 mdout/observable scalar 的动态 H5MD 路径纳入与 module dynamic path 相同的 path-builder 契约，降低 single-file、observable-only 与 VDS wrapper 之间的 schema drift 风险。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 120. SPONGE provenance dynamic path builder coverage 补充

本轮补充 SPONGE provenance 动态路径的 helper 与测试覆盖：

- 在 `h5md_writer.hpp` 中新增 `Sponge_Provenance_Path`，用于 `/parameters/sponge/provenance/<name>`。
- `ObservableH5Writer::Write_Provenance_String` 改为通过 `Sponge_Provenance_Path` 写入 provenance string。
- `Test_Public_H5MD_Path_Constants` 新增 `Sponge_Provenance_Path("launch_id")` literal 断言。
- `Test_Observable_Only_Writer` 中 launch provenance string 断言改为通过 helper 访问。
- `test_highfive_backend_io` 中真实 observable-only `.obs.spg.h5md` 文件的 launch provenance existence/readback 改为通过 helper 访问。
- 同步清理 mock writer 与 VDS mock writer 中残留的 metadynamics/QC/ReaxFF module dynamic leaf 旧式字符串拼接，使这些断言统一使用已有 module helper。

该补充把 `/parameters/sponge/provenance/<name>` 纳入 path-builder contract，避免后续 provenance key 扩展继续手写 HDF5 path。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 121. NHC/SITS dynamic root path builder coverage 补充

本轮补充 module extension 中 NHC 与 SITS dynamic root 路径的 helper 与测试覆盖：

- 在 `module_h5_mappings.hpp` 中新增 `Nose_Hoover_Chain_Coordinate_Root` 与 `Nose_Hoover_Chain_Velocity_Root`，用于 `/observables/all/thermostat/nose_hoover_chain/coordinate|velocity`。
- 在 `module_h5_mappings.hpp` 中新增 `Sits_Module_Root`，用于 `/observables/all/sits/<module>`。
- `Sits_Nk_Root` 改为基于 `Sits_Module_Root` 组合，避免重复拼接 SITS module root。
- `ModuleH5MappingWriter::Ensure_Nose_Hoover_Chain_Observables` 改为通过 NHC root helpers 创建 coordinate/velocity groups。
- `ModuleH5MappingWriter::Ensure_Sits_Nk_Observable` 改为通过 `Sits_Module_Root` 创建 module group。
- `VdsTrajectoryH5Writer` 的 NHC/SITS wrapper group materialization 改为通过同一组 helpers。
- `Test_Module_Path_Constants` 新增 NHC coordinate/velocity root 与 SITS module root literal 断言。
- module mock 与 VDS mock 测试新增 helper 对应 group existence 断言。

该补充覆盖此前只通过 value/step/time dataset 间接约束的 module dynamic root group，降低 NHC/SITS wrapper 与 single-file module mapping 之间的路径漂移风险。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 122. Generic module scalar observable leaf path builder coverage 补充

本轮补充 module scalar observable 通用 leaf 路径 helper 与测试覆盖：

- 在 `module_h5_mappings.hpp` 中新增 `Scalar_Observable_Value_Path`、`Scalar_Observable_Step_Path`、`Scalar_Observable_Time_Path`，用于 `<root>/value|step|time`。
- `Sits_Nk_*_Path`、`Metadynamics_Scalar_*_Path`、`Qc_Observable_*_Path`、`Reaxff_Term_*_Path` 改为复用通用 scalar observable leaf helpers。
- `ModuleH5MappingWriter::Create_Scalar_Observable_With_Axis` 改为通过通用 leaf helpers 创建 value dataset 和 step/time hard links。
- `VdsTrajectoryH5Writer::Create_Module_Scalar_With_Axis` 改为通过同一组 helpers 创建 wrapper virtual dataset 和 hard links。
- `Test_Module_Path_Constants` 新增通用 scalar observable leaf helper literal 断言。
- 窄范围文本检查确认生产代码中 module scalar `group/root + "/value|step|time"` 旧式拼接只剩 helper 实现本身。

该补充把 metadynamics、QC、ReaxFF、SITS `nk` 共享的 scalar observable leaf schema 收敛到同一套 path-builder 契约，避免 VDS wrapper 和 single-file module writer 分别维护 leaf path 组合逻辑。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 123. VDS source dataset path mapping coverage 补充

本轮补充 VDS virtual dataset source mapping 的路径覆盖：

- `Test_Vds_Wrapper_And_Module_Virtual_Datasets` 现在显式检查 particle VDS source 的 `dataset_path` 等于 `path::position_value`。
- 同一测试检查 ordinary observable VDS source 的 `dataset_path` 等于 `Observable_Value_Path("temperature")`。
- 同一测试检查 NHC、SITS、metadynamics、QC、ReaxFF 代表性 module VDS source 的 `dataset_path` 分别等于对应 helper path。
- 该覆盖与已有的 source file path 相对化、source dims、virtual start 断言互补，锁住 wrapper VDS 指向 shard 内 dataset 的 HDF5 path 契约。
- 同步清理 HighFive module diagnostic readback 中残留的局部 `metad_root + "/component"` 拼接，改为 `Metadynamics_Diagnostic_Path`。

该补充确保 VDS wrapper 不仅创建正确的 virtual dataset path，也把 source mapping 指向 shard 文件内的同名 canonical dataset path。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 124. HighFive nested group idempotence coverage 补充

本轮补充真实 HighFive backend 对嵌套 group 创建的覆盖：

- 新增 `Test_HighFive_Backend_Nested_Group_Idempotence`。
- 测试通过 `H5MDWriter::Ensure_Group` 两次创建 `/parameters/sponge/test/deep/group`，锁住重复 group 创建必须幂等。
- 测试随后写入 `/parameters/sponge/test/deep/group/value` 并在真实 HDF5 文件中读回，确认父级 group 链被递归创建。
- 测试已加入 `test_highfive_backend_io.cpp` 的 `Run_Test` 列表。

该补充覆盖 HighFive backend 的 path materialization 基础能力，避免新 H5 bundle schema 增加深层 `/parameters/sponge/...` 路径时只在 mock backend 中通过。本文档仅记录代码与测试更新；本轮未执行构建、单元测试或 CTest。

## 125. H5 bundle test registration coverage audit 补充

本轮补充 h5 bundle 测试注册链审计：

- 顶层 `CMakeLists.txt` 在 `SPONGE_BUILD_TESTS=ON` 时调用 `enable_testing()` 并 `add_subdirectory(tests)`。
- `tests/CMakeLists.txt` 调用 `add_subdirectory(h5_bundle)`。
- `tests/h5_bundle/CMakeLists.txt` 中当前 7 个 `test_*.cpp` 均通过 `add_sponge_h5_bundle_test` 注册为 executable 和 CTest test。
- `tests/h5_bundle/TEST_TARGETS.md` 中的 test executable/file/label 表与当前 `test_*.cpp` 文件集合一致。
- aggregate targets 覆盖 contract、mock、backend-io 三类测试：`sponge_h5_bundle_contract_tests`、`sponge_h5_bundle_mock_tests`、`sponge_h5_bundle_backend_io_tests`，以及总目标 `sponge_h5_bundle_tests`。
- HighFive/HDF5 dependency 在 `tests/h5_bundle/CMakeLists.txt` 中通过 `find_package(HighFive CONFIG REQUIRED)` 和 `find_package(HDF5 REQUIRED COMPONENTS C)` 接入，并链接到公共测试 interface target。

该补充确认新增 h5 bundle 路径测试不是孤立源码，而是处于 `SPONGE_BUILD_TESTS=ON` 的 CMake/CTest 注册链中。本文档仅记录源码级审计；本轮未执行 CMake configure、构建、单元测试或 CTest。

## 126. Test schema literal cleanup coverage 补充

本轮对 `tests/h5_bundle` 中残留的公开 H5 schema 字面量做了窄范围审计，目标是避免测试在普通断言中绕过 `SpongeH5MD::path::*` 或 path helper，从而在 schema 路径调整时形成隐性漂移。

处理原则如下：

- 保留专门用于锁定 public path constant/helper literal 的测试，这类测试的职责就是确认 helper 对应的规范路径没有意外变化。
- 保留 HighFive backend 低层能力测试中的任意 HDF5 路径，这些路径用于验证 backend 的通用 group/dataset/link 行为，不代表 SPONGE H5 bundle schema owner。
- 替换普通 writer 行为断言中仍硬编码的公开 schema 路径。

本轮代码更新：

- `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` 中 observable-only writer 的 `frame_count` 追加计数断言改为使用 `path::output_frame_count`。

该补充覆盖的风险点是：observable-only 输出路径已经由 public constants 管理，测试不应在 writer 行为断言处重新复制 `/parameters/sponge/output/frame_count` 字面量。

未运行 CMake、构建或 CTest；本节只记录契约和测试源码层面的收尾。

## 127. Restart writer failure contract coverage 补充

本轮继续补齐 wrapper 级失败路径测试，重点是 `restart.spgr.h5` 的单帧 restart 契约。

新增覆盖：

- `Test_Restart_Writer_Is_Single_State` 现在不仅断言二次 `Write_Structural_State` 会使 writer 进入 `FileStatus::failed`，还断言：
  - `/parameters/sponge/output/status` 写入 `failed`；
  - `/parameters/sponge/output/error` 写入 `restart H5 already contains one structural state`。

该测试与 trajectory/observable writer 的 failure-path 断言保持一致，确保 restart bundle 的失败原因不会只停留在内存 `last_error_`，而是落入 H5 bundle 的标准 output 状态路径。

未运行 CMake、构建或 CTest；本节只记录单元测试源码补充。

## 128. VDS wrapper shard-finalize failure coverage 补充

本轮审计 `VdsTrajectoryH5Writer::Finalize_Internal` 时发现：当当前 shard 在 `Complete_Current_Shard` 阶段 finalize 失败时，函数会直接返回 `false`，但此前没有把失败写回 wrapper H5MD 的标准 output 状态字段。这会使 VDS wrapper 在 shard finalize failure 场景下缺少 `/parameters/sponge/output/status` 和 `/parameters/sponge/output/error` 记录。

本轮实现和测试更新：

- `SPONGE/utils/h5md/vds_trajectory_h5_writer.hpp`
  - 在 `Complete_Current_Shard` 失败分支调用 `wrapper_writer_->Mark_Failed(last_error_)`。
- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - `Test_Strict_Finalize_Fails_On_Shard_Finalize_Error` 现在断言 wrapper backend：
    - `FileStatus::failed`；
    - `/parameters/sponge/output/status == failed`；
    - `/parameters/sponge/output/error == mock finalize failure`。

该补充保证 VDS/chunked H5MD 输出的 wrapper 文件在 shard finalize 失败时仍携带标准 bundle 失败元数据，便于后续 resolver、修复工具或运行系统判断输出是否完整。

未运行 CMake、构建或 CTest；本节记录源码和单元测试契约补充。

## 129. VDS materialization failure coverage 补充

本轮继续补齐 VDS/chunked H5MD wrapper 的失败路径覆盖。此前已有 shard finalize failure 覆盖，但 materialization 阶段的 wrapper virtual dataset 创建失败没有专门测试。

本轮实现和测试更新：

- `tests/h5_bundle/h5_bundle_test_common.hpp`
  - `BackendLog` 增加 `fail_next_virtual_dataset` 测试注入开关。
  - `MockBackend::Create_Virtual_Dataset` 在该开关置位时返回失败，并设置 `last_error = "mock virtual dataset failure: <path>"`。
- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - 新增 `Test_Vds_Materialize_Failure_Marks_Wrapper_Failed`。
  - 测试使用最小 VDS trajectory，触发 wrapper 第一个 virtual dataset `/particles/all/step` 创建失败。
  - 断言 `writer.Last_Error()`、wrapper `FileStatus::failed`、`/parameters/sponge/output/status` 和 `/parameters/sponge/output/error` 一致。

该补充覆盖的契约是：VDS wrapper 在 materialization 阶段失败时，不能只返回 `false`，必须把失败原因写入标准 H5 bundle output 状态路径。

未运行 CMake、构建或 CTest；本节记录单元测试源码补充。

## 130. VDS finalize tail failure metadata coverage 补充

本轮继续审计 VDS/chunked H5MD `Finalize_Internal` 的后半段失败路径。此前 shard finalize、manifest validation、VDS materialization 等分支已经会写 wrapper failed metadata，但 `Write_Manifest_To_Wrapper`、`Write_Repair_Metadata`、`output_vds_status` 写入和最终 wrapper `Finalize` 失败分支仍可能只返回 `false`。

本轮实现更新：

- `SPONGE/utils/h5md/vds_trajectory_h5_writer.hpp`
  - `Write_Manifest_To_Wrapper` 失败时现在会同步 `last_error_ = wrapper_writer_->Last_Error()`。
  - `Finalize_Internal` 在 manifest 写入、repair metadata 写入、`output_vds_status` 写入、最终 wrapper finalize 失败时，都会调用 `wrapper_writer_->Mark_Failed(last_error_)`。

本轮测试更新：

- `tests/h5_bundle/h5_bundle_test_common.hpp`
  - `BackendLog` 增加 `fail_next_string_array` 测试注入开关。
  - `MockBackend::Write_String_Array` 可模拟 manifest string-array 写入失败。
- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - 新增 `Test_Vds_Manifest_Write_Failure_Marks_Wrapper_Failed`。
  - 测试在 VDS materialization 成功后，让 manifest path array 写入失败，并断言：
    - `writer.Last_Error()` 精确记录 `mock string-array failure: /parameters/sponge/shards/path`；
    - wrapper `FileStatus::failed`；
    - `/parameters/sponge/output/status == failed`；
    - `/parameters/sponge/output/error` 记录相同失败原因。

该补充覆盖的契约是：VDS wrapper 在 finalize 后半段失败时，仍必须写出标准 H5 bundle failure metadata，而不是只依赖调用栈返回值。

未运行 CMake、构建或 CTest；本节记录源码和单元测试契约补充。

## 131. VDS string metadata failure coverage 补充

本轮继续补齐 VDS/chunked H5MD finalize tail 中的 string metadata 失败路径。此前 manifest string-array 写入失败已有覆盖，但 repair metadata 和 final VDS status 这两个 string path 仍缺少精确失败注入测试。

本轮测试基础设施更新：

- `tests/h5_bundle/h5_bundle_test_common.hpp`
  - `BackendLog` 增加 `fail_string_path`。
  - `MockBackend::Write_String` 在写入指定 path 时返回失败，并设置 `last_error = "mock string failure: <path>"`。

本轮新增测试：

- `Test_Vds_Repair_Metadata_Failure_Marks_Wrapper_Failed`
  - 触发 `path::output_repair_policy` 写入失败。
  - 断言 `writer.Last_Error()`、wrapper `FileStatus::failed`、`path::output_status` 和 `path::output_error` 一致。
- `Test_Vds_Status_Write_Failure_Marks_Wrapper_Failed`
  - 触发 `path::output_vds_status` 写入失败。
  - 断言同样的 wrapper failed metadata 契约。

该补充覆盖的契约是：VDS wrapper 的 repair metadata 和 VDS materialized status 属于 H5 bundle 新路径，写入失败时必须通过标准 output failed/error path 对外可见。

未运行 CMake、构建或 CTest；本节记录单元测试源码补充。

## 132. VDS wrapper final finalize failure coverage 补充

本轮补齐 VDS/chunked H5MD finalize tail 的最后一个独立失败场景：shard finalize 已成功、manifest 和 metadata 已写完，但 wrapper 自身最终 `Finalize()` 失败。

新增测试：

- `Test_Vds_Wrapper_Finalize_Failure_Marks_Wrapper_Failed`
  - 使用 `factory.fail_finalize = {true, false}`，其中第一个 backend 是 wrapper，第二个 backend 是 shard。
  - shard finalize 成功后，最终 wrapper finalize 返回失败。
  - 断言 `writer.Last_Error() == "mock finalize failure"`。
  - 断言 wrapper `FileStatus::failed`、`path::output_status == failed`、`path::output_error == mock finalize failure`。

该测试与 `Test_Strict_Finalize_Fails_On_Shard_Finalize_Error` 区分开：后者覆盖 shard finalize failure，本轮新增测试覆盖 wrapper final finalize failure。

未运行 CMake、构建或 CTest；本节记录单元测试源码补充。

## 133. Strict trajectory chunk-size parsing coverage 补充

本轮审计 `output_h5_trajectory_chunk_size` 时发现旧实现通过 `atoi` 解析，因此 `12frames` 这类带尾随文本的值会被部分解析为 `12`，存在绕过整数类型约束的风险。

本轮实现更新：

- `SPONGE/utils/control/h5_output_contract.hpp`
  - `Trajectory_Chunk_Size` 改为基于 `strtol` 的完整整数解析。
  - 拒绝空值、无数字、溢出、以及尾随非空白字符。

本轮测试更新：

- `tests/h5_bundle/test_h5_output_plan.cpp`
  - `Test_Invalid_Values` 增加 `output_h5_trajectory_chunk_size = "12frames"`。
  - 断言 resolver 返回 invalid，并包含 `chunk_size` 错误信息。

该补充覆盖的契约是：`output_h5_trajectory_chunk_size` 是整数型 run policy，不允许被 `atoi` 式部分解析接受。

未运行 CMake、构建或 CTest；本节记录源码和单元测试契约补充。

## 134. Strict trajectory VDS bool parsing coverage 补充

本轮审计 `output_h5_trajectory_vds` 时发现：底层 `Parse_Bool` 会把未识别文本返回为 `false`，因此 resolver 层若不额外校验，`output_h5_trajectory_vds = "maybe"` 会被静默当作非 VDS 输出。

本轮实现更新：

- `SPONGE/utils/h5md/output_plan.hpp`
  - 新增 `Is_Bool_Text`，用于识别 parser-visible bool 文本。
  - `Resolve_Output_Plan` 在读取 `output_h5_trajectory_vds` 时先校验文本是否合法。
  - 非法值返回 invalid plan，错误信息为 `output_h5_trajectory_vds must be boolean`。
  - `Parse_Bool` 保持原有转换 helper 语义。

本轮测试更新：

- `tests/h5_bundle/test_h5_output_plan.cpp`
  - `Test_Bool_Parsing_Text_Variants` 增加 `Is_Bool_Text` 合法/非法值覆盖。
  - `Test_Invalid_Values` 增加 `output_h5_trajectory_vds = "maybe"`，断言 plan invalid，并确认错误信息包含 `output_h5_trajectory_vds`。

该补充覆盖的契约是：`output_h5_trajectory_vds` 是 bool 型 run policy，不允许未识别文本静默退化为 `false`。

未运行 CMake、构建或 CTest；本节记录源码和单元测试契约补充。

## 135. Trajectory repair-policy error contract coverage 补充

本轮继续收紧 `output_h5_trajectory_repair_policy` 的 plan 层错误契约。原测试已经覆盖非法 repair policy 和 `complete_prefix` 需要 VDS trajectory 的 invalid 状态，但没有钉住错误信息指向具体 parser-visible key/条件。

本轮测试更新：

- `tests/h5_bundle/test_h5_output_plan.cpp`
  - `repair_policy = "truncate"` 时，除 invalid 外，现在断言错误信息包含 `output_h5_trajectory_repair_policy`。
  - `repair_policy = "complete_prefix"` 但未设置 `output_h5_trajectory_vds=true` 时，除 invalid 外，现在断言错误信息包含 `output_h5_trajectory_vds=true`。

该补充覆盖的契约是：repair policy 相关配置错误必须能被用户定位到对应的 H5 bundle parser-visible key 或缺失条件。

未运行 CMake、构建或 CTest；本节记录单元测试源码补充。

## 136. Directory-aware VDS shard-root derivation coverage 补充

本轮审计 trajectory H5MD suffix 与 VDS shard root 派生时，已有 `x.spg.h5md -> x.spg.shards` 和非推荐 suffix fallback 的基础覆盖，但缺少带目录路径的断言。

本轮测试更新：

- `tests/h5_bundle/test_h5_output_plan.cpp`
  - `Test_Suffix_And_Shard_Derivation` 增加：
    - `Derive_Shards_Root("runs/prod.spg.h5md") == "runs/prod.spg.shards"`。

该补充覆盖的契约是：VDS/chunked H5MD 的派生 shard root 必须保留 trajectory bundle 所在目录，不能在路径派生时丢失 parent directory。

未运行 CMake、构建或 CTest；本节记录单元测试源码补充。

## 137. VDS shard filename sequence coverage 补充

本轮审计 VDS shard path 生成时，已有 `segment_000000.spg.h5md` 与 `segment_000001.spg.h5md` 的 opened path/source path 覆盖，但缺少对连续编号和六位 zero-padding 规则的直接测试。

本轮测试更新：

- `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp`
  - 新增 `Test_Vds_Shard_Filename_Sequence_Uses_Six_Digit_Padding`。
  - 使用 `chunk_size = 1` 写入三帧，强制生成三个 shard。
  - 断言第三个 shard manifest path 和实际 opened path 均为：
    - `/tmp/sponge_h5_vds_case/prod.spg.shards/segment_000002.spg.h5md`。

该补充覆盖的契约是：VDS shard 文件命名使用 `segment_%06d.spg.h5md`，并且 shard 序列递增时路径保持稳定。

未运行 CMake、构建或 CTest；本节记录单元测试源码补充。
