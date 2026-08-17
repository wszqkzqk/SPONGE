# SPONGE H5 Bundle Unit Test Audit Matrix

日期：2026-07-03

## 1. Scope

This matrix maps each new H5 bundle path or behavior to the unit test target
that is intended to cover it.

Status meanings:

| Status | Meaning |
|---|---|
| `unit-covered` | Covered by a unit or lightweight file-level test target in `tests/h5_bundle`. |
| `backend-covered` | Covered by a real HighFive/HDF5 file-level test in `test_highfive_backend_io`. |
| `documented-gap` | Intentionally not covered by unit tests; requires runtime/integration validation. |
| `not-run` | Test code exists but has not been built or executed in this goal yet. |

All test targets are currently `not-run` until configure/build/ctest is executed.

## 2. Contract and resolver

| Item | Test target | Status |
|---|---|---|
| `output_h5_trajectory_path` enables trajectory H5 output | `test_h5_output_plan` | unit-covered, not-run |
| `output_h5_trajectory_vds` does not enable output without trajectory path | `test_h5_output_plan` | unit-covered, not-run |
| `output_h5_trajectory_chunk_size` default is 20 | `test_h5_output_plan` | unit-covered, not-run |
| `output_h5_trajectory_chunk_size <= 0` is invalid | `test_h5_output_plan` | unit-covered, not-run |
| `output_h5_trajectory_repair_policy = strict` default | `test_h5_output_plan` | unit-covered, not-run |
| `complete_prefix` repair policy requires VDS trajectory path | `test_h5_output_plan` | unit-covered, not-run |
| `allow_complete_prefix` compatibility alias | `test_h5_output_plan` | unit-covered, not-run |
| bool parser accepts `1/true/TRUE/yes/on` | `test_h5_output_plan` | unit-covered, not-run |
| bool parser rejects `0/false/off/no` | `test_h5_output_plan` | unit-covered, not-run |
| `throw_on_error=true` uses controller error path | `test_h5_output_plan` | unit-covered, not-run |
| trajectory suffix `*.spg.h5md` | `test_h5_output_plan` | unit-covered, not-run |
| restart suffix `*.spgr.h5` | `test_h5_output_plan` | unit-covered, not-run |
| observable suffix `*.obs.spg.h5md` | `test_h5_output_plan` | unit-covered, not-run |
| H5 output disables implicit legacy sidecars | `test_h5_output_plan` | unit-covered, not-run |
| explicit legacy sidecars remain enabled | `test_h5_output_plan` | unit-covered, not-run |

## 3. Writer facade and mock backend

| Item | Test target | Status |
|---|---|---|
| H5MD common layout creates `/h5md`, `/observables`, `/parameters` | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| observable-only common layout omits `/particles` | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| trajectory position/box/velocity/force datasets | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| ordinary observable dataset and hard links | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| trajectory append failure marks failed | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| observable missing value marks failed | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| restart structural state is single-state | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| restart NHC state | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| restart SITS state | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| restart metad text state | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |
| legacy sidecar provenance in trajectory/observable/restart | `test_h5md_writers_with_mock_backend` | unit-covered, not-run |

## 4. Completion and failure

| Item | Test target | Status |
|---|---|---|
| completion state machine open/begin/complete/closing/finalized | `test_completion_tracker` | unit-covered, not-run |
| cannot begin a second frame while one is open | `test_completion_tracker` | unit-covered, not-run |
| cannot finalize incomplete frame | `test_completion_tracker` | unit-covered, not-run |
| `Mark_Failed` writes failed status and reason | `test_completion_tracker` | unit-covered, not-run |
| manifest strict validation accepts contiguous complete entries | `test_completion_tracker` | unit-covered, not-run |
| manifest rejects incomplete shard in strict mode | `test_completion_tracker` | unit-covered, not-run |
| manifest rejects non-contiguous shard indices | `test_completion_tracker` | unit-covered, not-run |
| manifest rejects non-contiguous frame ranges | `test_completion_tracker` | unit-covered, not-run |
| manifest rejects non-positive frame count | `test_completion_tracker` | unit-covered, not-run |
| repair-prefix manifest validation stops at trailing incomplete shard | `test_completion_tracker` | unit-covered, not-run |

## 5. Module-specific mapping

| Item | Test target | Status |
|---|---|---|
| NHC coordinate/velocity observables | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |
| SITS `nk` observable | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |
| metadynamics `meta/rbias/rct` scalar observables | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |
| metadynamics diagnostics under parameters | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |
| QC energy and optional spin square | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |
| QC SCF log under parameters | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |
| ReaxFF energy terms | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |
| ReaxFF missing term error | `test_module_h5_mappings_with_mock_backend` | unit-covered, not-run |

## 6. VDS writer

| Item | Test target | Status |
|---|---|---|
| shard rotation by trajectory frame count | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| manifest frame ranges and status | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| ordinary observable per-shard count | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| NHC/SITS/metad/QC/ReaxFF per-shard counts | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| particle VDS source dims and virtual starts | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| ordinary observable VDS source dims and virtual starts | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| module VDS source dims | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| relative VDS source path | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| strict finalize fails on shard finalize error | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| complete-prefix repair finalize | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| VDS wrapper metad diagnostics | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| VDS wrapper QC SCF log | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |
| VDS writer precondition errors | `test_vds_trajectory_writer_with_mock_backend` | unit-covered, not-run |

## 7. HighFive/HDF5 backend file-level behavior

| Item | Test target | Status |
|---|---|---|
| real HDF5 group/dataset/string/hardlink creation | `test_highfive_backend_io` | backend-covered, not-run |
| observable-only file omits `/particles` | `test_highfive_backend_io` | backend-covered, not-run |
| invalid backend operations fail | `test_highfive_backend_io` | backend-covered, not-run |
| raw VDS readback | `test_highfive_backend_io` | backend-covered, not-run |
| trajectory facade + real backend | `test_highfive_backend_io` | backend-covered, not-run |
| restart facade + real backend | `test_highfive_backend_io` | backend-covered, not-run |
| observable-only facade + real backend | `test_highfive_backend_io` | backend-covered, not-run |
| VDS trajectory writer + real backend | `test_highfive_backend_io` | backend-covered, not-run |
| particle VDS readback through wrapper | `test_highfive_backend_io` | backend-covered, not-run |
| ordinary observable VDS readback through wrapper | `test_highfive_backend_io` | backend-covered, not-run |
| module-specific VDS readback through wrapper | `test_highfive_backend_io` | backend-covered, not-run |
| real output status/completion metadata readback | `test_highfive_backend_io` | backend-covered, not-run |
| VDS repair metadata readback | `test_highfive_backend_io` | backend-covered, not-run |

## 8. Documented gaps requiring integration tests

| Item | Status | Reason |
|---|---|---|
| Full SPONGE MD run with H5 outputs enabled | documented-gap | Requires runtime simulation input and executable validation. |
| CUDA/MPI rank behavior | documented-gap | Requires runtime or rank-aware integration harness. |
| Restart H5 reader round-trip | documented-gap | Reader is not implemented in this unit-test scope. |
| Cross-process kill/resume | documented-gap | Requires process-level integration test. |
| HDF5 dataset truncation | documented-gap | Not implemented as a feature. |
| Orphan shard cleanup | documented-gap | Not implemented as a feature. |
| Metadynamics structured binary restart | documented-gap | Schema and runtime data path remain future work. |

## 9. Completion condition for this goal

The unit test writing portion is complete when:

1. all test files listed in this matrix exist;
2. the tests are registered in CMake under `SPONGE_BUILD_TESTS=ON`;
3. `sponge_h5_bundle_tests` can build;
4. all `h5_bundle` CTest targets pass, or failures are documented and fixed;
5. this matrix is updated from `not-run` to validated evidence.

Current status: test code and audit matrix are written, but build/CTest evidence
is still missing.

## 10. Additional backend string-array coverage

| Item | Test target | Status |
|---|---|---|
| real HDF5 string-array write/read | `test_highfive_backend_io` | backend-covered, not-run |
| VDS shard manifest `path/status` string arrays | `test_highfive_backend_io` | backend-covered, not-run |

## 11. Additional VDS negative backend coverage

| Item | Test target | Status |
|---|---|---|
| VDS source rank mismatch is rejected | `test_highfive_backend_io` | backend-covered, not-run |
| VDS source virtual-start rank mismatch is rejected | `test_highfive_backend_io` | backend-covered, not-run |

## 12. Additional contract helper coverage

| Item | Test target | Status |
|---|---|---|
| `Any_H5_Output_Enabled(nullptr)` returns false | `test_h5_output_plan` | unit-covered, not-run |
| legacy sidecar default/explicit helper behavior | `test_h5_output_plan` | unit-covered, not-run |
| null controller resolver error | `test_h5_output_plan` | unit-covered, not-run |
| unknown output key has no recommended suffix | `test_h5_output_plan` | unit-covered, not-run |

## 13. 真实 backend 失败元数据审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| 失败输出必须写出 `failed` 状态 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Failed_Metadata` | 真实 HighFive 文件读回 `/parameters/sponge/output/status` | 已覆盖，未运行 |
| 失败原因必须写入 sponge output parameters | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Failed_Metadata` | 真实 HighFive 文件读回 `/parameters/sponge/output/error` | 已覆盖，未运行 |

该补充用于确保失败路径不仅在 mock backend 中被记录，也在实际 `.spg.h5md` 文件中具备可审计状态。

## 14. 真实 backend 参数与 legacy provenance 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory bundle 可选 force field 必须落在 H5MD particle force 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `path::force_value` shape 与 payload | 已覆盖，未运行 |
| `mdinfo` 文本应进入 `/parameters/sponge` 子树 | `tests/h5_bundle/test_highfive_backend_io.cpp` / trajectory 与 observable real backend tests | 真实 HighFive 文件读回 `path::mdinfo_text` | 已覆盖，未运行 |
| legacy sidecar 路径只作为 bundle provenance 记录 | `tests/h5_bundle/test_highfive_backend_io.cpp` / trajectory、restart、observable real backend tests | 真实 HighFive 文件读回 `path::legacy_sidecar_keys` 与 `path::legacy_sidecar_paths` | 已覆盖，未运行 |
| observable-only bundle 可记录 launch provenance | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `/parameters/sponge/provenance/launch_id` | 已覆盖，未运行 |

## 15. HighFive backend 状态机审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| 未打开 HDF5 文件时 backend 写操作必须失败 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Status_State` | 直接调用 `Append_Int64` 并检查 `FileStatus::failed` 与 `Last_Error()` | 已覆盖，未运行 |
| backend 显式状态更新必须写入 output status dataset | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Status_State` | `Set_Status(FileStatus::closing)` 后真实 HighFive 文件读回 `path::output_status` | 已覆盖，未运行 |
| backend close 后内存状态必须为 closed | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Status_State` | 调用 `Close()` 后检查 `HighFiveBackend::Status()` | 已覆盖，未运行 |

## 16. VDS wrapper 参数路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 必须记录 chunk size 语义 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `/parameters/sponge/output/trajectory_chunk_size` | 已覆盖，未运行 |
| VDS wrapper 必须记录 materialization 状态 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `/parameters/sponge/output/vds_status` | 已覆盖，未运行 |
| VDS wrapper 可承载 metad diagnostic 文本 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 metad parameter `hills` 文本 | 已覆盖，未运行 |
| VDS wrapper 可承载 QC SCF legacy 文本 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 QC parameter `scf_output` 文本 | 已覆盖，未运行 |
| VDS wrapper 可承载 legacy sidecar provenance | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 legacy sidecar key/path 数组 | 已覆盖，未运行 |

## 17. output plan 解析边界审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| grouped TOML flatten 后的 trajectory chunk size 默认值为 20 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Contract_Helper_Edge_Cases` | 直接调用 `Trajectory_Chunk_Size` 验证默认与显式值 | 已覆盖，未运行 |
| 未启用 H5 bundle 时 legacy sidecar 默认全部开启 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Legacy_Output_Plan_All_Keys` | 遍历 `mdout/mdinfo/crd/box/vel/frc/rst/qc_scf_output` | 已覆盖，未运行 |
| 启用 H5 bundle 后 legacy sidecar 默认关闭，仅显式 key 开启 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Legacy_Output_Plan_All_Keys` | observable bundle + explicit `vel` 场景 | 已覆盖，未运行 |
| 非推荐文件后缀只产生标记，不阻止 plan 有效 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Output_Path_Suffix_Flags_Are_Non_Fatal` | trajectory/restart/observable 错后缀仍 `valid=true` | 已覆盖，未运行 |
| `complete_prefix` repair 必须绑定 VDS trajectory | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Complete_Prefix_Repair_Requires_Vds_Trajectory` | 无 path 与 path 非 VDS 两个 invalid 场景 | 已覆盖，未运行 |

## 18. HighFive backend append 数据完整性审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| append count 必须匹配一个完整 record | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Rejects_Invalid_Operations` | 对 scalar observable dataset 使用 `count=2` append 并要求失败 | 已覆盖，未运行 |
| append 数据指针不可为空 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Rejects_Invalid_Operations` | 对真实 HighFive dataset 使用 `nullptr` append 并要求失败 | 已覆盖，未运行 |

## 19. restart 单帧不变量审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart bundle 只允许一个 structural state | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_Rejects_Second_State_With_Real_Backend` | 第二次 `Write_Structural_State` 必须失败 | 已覆盖，未运行 |
| restart 单帧违规必须写出 failed 状态和错误原因 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_Rejects_Second_State_With_Real_Backend` | 真实 `.spgr.h5` 读回 output status/error | 已覆盖，未运行 |
| restart 单帧违规不能把 position dataset 扩展为多帧 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_Rejects_Second_State_With_Real_Backend` | 真实 `.spgr.h5` 读回 position 第一维仍为 1 | 已覆盖，未运行 |

## 20. ModuleH5MappingWriter metad/ReaxFF 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| metad potential export 应写入 `/parameters/sponge/metadynamics` | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Module_Metad_And_Reaxff_With_Real_Backend` | 真实 HighFive 文件读回 `potential_export` 文本 | 已覆盖，未运行 |
| metad direct export 应写入 `/parameters/sponge/metadynamics` | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Module_Metad_And_Reaxff_With_Real_Backend` | 真实 HighFive 文件读回 `direct_export` 文本 | 已覆盖，未运行 |
| metad hills/history/edge 应作为参数文本保留 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Module_Metad_And_Reaxff_With_Real_Backend` | 真实 HighFive 文件读回 `hills/history/edge` 文本 | 已覆盖，未运行 |
| ReaxFF 多 energy term 应各自形成 observable value | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Module_Metad_And_Reaxff_With_Real_Backend` | 真实 HighFive 文件读回 `bond/angle/over/value` payload | 已覆盖，未运行 |

## 21. H5MDWriter common layout 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| detached H5MD writer 必须安全失败 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_H5MD_Writer_Detached_Backend_Semantics` | 直接调用 `Open/Flush/Close/Finalize/Mark_Failed` 并要求 false | 已覆盖，未运行 |
| detached H5MD writer 必须报告明确错误 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_H5MD_Writer_Detached_Backend_Semantics` | 检查 `Status()==closed` 与 `Last_Error()` | 已覆盖，未运行 |
| common layout 必须写入 schema name/version | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Basic_File_Layout` | 真实 HighFive 文件读回 `/parameters/sponge/schema/name` 与 `/version` | 已覆盖，未运行 |
| common layout 必须初始化 output completion | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Basic_File_Layout` | 真实 HighFive 文件读回 frame count/last step/last time 初值 | 已覆盖，未运行 |

## 22. writer facade open precondition 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| single-file trajectory writer 不得接受未 enabled trajectory plan | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Writer_Open_Preconditions` | mock backend 上调用 `Open_Single_File` 并检查失败与错误消息 | 已覆盖，未运行 |
| single-file trajectory writer 不得接受 VDS trajectory plan | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Writer_Open_Preconditions` | 设置 `plan.trajectory.vds=true` 后要求失败 | 已覆盖，未运行 |
| restart writer 不得接受未 enabled restart plan | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Writer_Open_Preconditions` | mock backend 上调用 `RestartH5Writer::Open` 并检查错误消息 | 已覆盖，未运行 |
| observable writer 不得接受未 enabled observable plan | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Writer_Open_Preconditions` | mock backend 上调用 `ObservableH5Writer::Open` 并检查错误消息 | 已覆盖，未运行 |

## 23. ObservableH5Writer 模块输出审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| observable-only bundle 可承载 NHC observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_With_Real_Backend` | 真实 HighFive 文件读回 NHC coordinate/velocity payload | 已覆盖，未运行 |
| observable-only bundle 可承载 SITS nk observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `Sits_Nk_Value_Path("obs_sits")` payload | 已覆盖，未运行 |
| observable-only bundle 可承载 metad scalar observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `meta/rbias/rct` payload | 已覆盖，未运行 |
| observable-only bundle 可承载 ReaxFF 多 term observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `bond/angle/value` payload | 已覆盖，未运行 |

## 24. TrajectoryH5Writer 模块输出审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| single-file trajectory bundle 可承载 NHC observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 NHC coordinate/velocity payload | 已覆盖，未运行 |
| single-file trajectory bundle 可承载 SITS nk observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `Sits_Nk_Value_Path("traj_sits")` payload | 已覆盖，未运行 |
| single-file trajectory bundle 可承载 metad scalar observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `meta/value` payload | 已覆盖，未运行 |
| single-file trajectory bundle 可承载 QC observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `energy/spin_square` payload | 已覆盖，未运行 |
| single-file trajectory bundle 可承载 ReaxFF 多 term observables | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HighFive 文件读回 `bond/angle/value` payload | 已覆盖，未运行 |

## 25. VDS complete-prefix repair 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| complete-prefix repair 必须保留已完成 shard 前缀 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Complete_Prefix_Repair_Retains_Valid_Prefix` | 两个 shard 中第二个 finalize 失败后，manifest 只保留第一个 complete shard | 已覆盖，未运行 |
| repair 后 total trajectory frame count 必须按保留前缀重算 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Complete_Prefix_Repair_Retains_Valid_Prefix` | 检查 `Total_Trajectory_Frame_Count()==1` | 已覆盖，未运行 |
| repair 后 VDS source 只能引用保留 shard | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Complete_Prefix_Repair_Retains_Valid_Prefix` | mock wrapper virtual dataset source 只包含 `segment_000000.spg.h5md` | 已覆盖，未运行 |
| repair metadata 必须记录 applied 和 repaired shard count | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Complete_Prefix_Repair_Retains_Valid_Prefix` | 检查 repair policy/status 和 `repaired_shard_count` append | 已覆盖，未运行 |

## 26. OutputCompletionTracker 真实 backend 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| completion tracker 状态必须落盘到 output status | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Output_Completion_Tracker_With_Real_Backend` | 真实 HighFive 文件读回 finalized status | 已覆盖，未运行 |
| completion tracker frame count 必须落盘 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Output_Completion_Tracker_With_Real_Backend` | 真实 HighFive 文件读回 `frame_count.back()==1` | 已覆盖，未运行 |
| completion tracker last step/time 必须落盘 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Output_Completion_Tracker_With_Real_Backend` | 真实 HighFive 文件读回 last complete step/time | 已覆盖，未运行 |

## 27. H5MD step/time hard-link 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory particle datasets 必须共享 `/particles/all` step/time axis | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_Writer_Paths_And_Completion` | mock backend 检查 position/velocity/force/box hard-link 创建 | 已覆盖，未运行 |
| trajectory observables 必须共享 `/observables/all` step/time axis | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_Writer_Paths_And_Completion` | mock backend 检查 temperature step/time hard-link 创建 | 已覆盖，未运行 |
| observable-only observables 必须共享 `/observables/all` step/time axis | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Observable_Only_Writer` | mock backend 检查 energy step/time hard-link 创建 | 已覆盖，未运行 |
| restart structural datasets 必须共享 `/particles/all` step/time axis | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Writer_Is_Single_State` | mock backend 检查 position/velocity/box hard-link 创建 | 已覆盖，未运行 |

## 28. ModuleH5MappingWriter step/time hard-link 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| NHC coordinate/velocity 必须链接到 NHC step/time axis | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Nhc_And_Sits_Mappings` | mock backend 检查 NHC coordinate/velocity hard-link 创建 | 已覆盖，未运行 |
| metad scalar observables 必须链接到 metad step/time axis | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | mock backend 检查 `meta/rbias/rct` hard-link 创建 | 已覆盖，未运行 |
| QC observables 必须链接到 QC step/time axis | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Qc_And_Reaxff_Mappings` | mock backend 检查 `energy/spin_square` hard-link 创建 | 已覆盖，未运行 |
| ReaxFF term observables 必须链接到 ReaxFF step/time axis | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Qc_And_Reaxff_Mappings` | mock backend 检查 `bond/angle` hard-link 创建 | 已覆盖，未运行 |
| SITS nk 当前未建立 value step/time hard-link | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Nhc_And_Sits_Mappings` | 当前仅覆盖 dataset 与 append 行为，保留实现现状 | 记录，未运行 |

## 29. VDS module append precondition 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS NHC frame 必须先定义 NHC layout | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Precondition_Errors` | 未定义 layout 时 append NHC 要求失败 | 已覆盖，未运行 |
| VDS NHC chain length 不得改变 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Precondition_Errors` | 定义 chain length 后用不同长度 append 要求失败 | 已覆盖，未运行 |
| VDS metad scalar frame 必须先定义 metad layout | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Precondition_Errors` | 未定义 layout 时 append metad scalar 要求失败 | 已覆盖，未运行 |
| VDS QC frame 必须先定义 QC layout | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Precondition_Errors` | 未定义 layout 时 append QC 要求失败 | 已覆盖，未运行 |
| VDS ReaxFF frame 必须先定义 ReaxFF layout | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Precondition_Errors` | 未定义 layout 时 append ReaxFF 要求失败 | 已覆盖，未运行 |

## 30. restart base layout 与 run metadata 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart bundle 必须包含 restart parameter root | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 真实 `.spgr.h5` 检查 `path::parameters_restart` | 已覆盖，未运行 |
| restart bundle 必须预留 thermostat/barostat/bias state 分组 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 真实 `.spgr.h5` 检查 `restart_thermostat/restart_barostat/restart_bias` | 已覆盖，未运行 |
| restart bundle 必须记录 current time 和 state type | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 真实 `.spgr.h5` 读回 `run_current_time` 与 `run_state_type` | 已覆盖，未运行 |

## 31. optional velocity 与 output helper 边界审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart structural state 不请求 velocity 时不得创建 velocity value dataset | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_Rejects_Second_State_With_Real_Backend` | 真实 `.spgr.h5` 检查 `!file.exist(path::velocity_value)` | 已覆盖，未运行 |
| output plan `Command_String` 对 null/missing key 必须安全返回空字符串 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Contract_Helper_Edge_Cases` | 直接调用 helper 验证 null 和 missing 场景 | 已覆盖，未运行 |
| output plan `Command_String` 对 existing key 必须返回 value | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Contract_Helper_Edge_Cases` | 直接调用 helper 验证 `present=value` | 已覆盖，未运行 |

## 32. VDS source path 相对化审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 应使用相对 source path 引用 sibling shard root | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Source_Path_Relativization` | mock wrapper virtual dataset source 读回 `../shards/prod.spg.shards/segment_000000.spg.h5md` | 已覆盖，未运行 |
| VDS shard path 格式必须使用 `segment_%06d.spg.h5md` | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Source_Path_Relativization` | source path 后缀断言 `segment_000000.spg.h5md` | 已覆盖，未运行 |

## 33. VDS writer open precondition 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS writer 必须要求 backend factory | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Open_Precondition_Errors` | null factory 调用 `Open` 并检查错误消息 | 已覆盖，未运行 |
| VDS writer 必须要求 trajectory plan enabled | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Open_Precondition_Errors` | `trajectory.enabled=false` 调用 `Open` 并检查错误消息 | 已覆盖，未运行 |
| VDS writer 必须要求 trajectory plan 为 VDS 模式 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Open_Precondition_Errors` | `trajectory.vds=false` 调用 `Open` 并检查错误消息 | 已覆盖，未运行 |

## 34. VDS 空输出 finalize 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 无帧 finalize 时 manifest 必须为空 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Finalize_Without_Frames` | 只 open/finalize 后检查 manifest size 为 0 | 已覆盖，未运行 |
| VDS wrapper 无帧 finalize 时总 frame count 必须为 0 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Finalize_Without_Frames` | 检查 trajectory/observable frame count 均为 0 | 已覆盖，未运行 |
| VDS wrapper 无帧 finalize 仍需写出 repair/vds metadata | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Finalize_Without_Frames` | mock wrapper 检查 repair policy/status/count 和 vds_status | 已覆盖，未运行 |
| VDS wrapper 无帧 finalize 不应创建 virtual datasets | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Finalize_Without_Frames` | mock wrapper 检查 `virtual_datasets.empty()` | 已覆盖，未运行 |

## 35. output plan 解析鲁棒性审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| repair policy 解析应大小写不敏感 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Repair_Policy_Is_Case_Insensitive` | `COMPLETE_PREFIX` 解析为 `complete_prefix` | 已覆盖，未运行 |
| trajectory VDS bool 解析应接受大写 TRUE | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Repair_Policy_Is_Case_Insensitive` | `output_h5_trajectory_vds=TRUE` 后 plan valid 且 repair enabled | 已覆盖，未运行 |
| trajectory chunk size 负数必须 invalid | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Invalid_Values` | `chunk_size=-3` 要求 invalid 且错误包含 chunk_size | 已覆盖，未运行 |
| trajectory chunk size 非数字必须 invalid | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Invalid_Values` | `chunk_size=not_an_integer` 要求 invalid | 已覆盖，未运行 |
| explicit legacy sidecar path 必须保留 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Legacy_Output_Plan_All_Keys` | H5 enabled + `vel=legacy.vel` 后检查 sidecar path | 已覆盖，未运行 |

## 36. VDS source path 无 wrapper parent 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| wrapper 无 parent directory 时 VDS source path 应保持 shard root 相对形式 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Source_Path_Without_Wrapper_Parent` | mock wrapper virtual dataset source 读回 `prod.spg.shards/segment_000000.spg.h5md` | 已覆盖，未运行 |
| wrapper 无 parent directory 时 shard 文件名仍使用固定 segment 格式 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Source_Path_Without_Wrapper_Parent` | source path 后缀断言 `segment_000000.spg.h5md` | 已覆盖，未运行 |

## 37. H5MDWriter repeated output completion 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| output completion 初始化记录必须保留 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_H5MD_Writer_Repeated_Output_Completion` | 真实 HighFive 文件读回 `frame_count[0]=0` 和 `last_step[0]=-1` | 已覆盖，未运行 |
| repeated output completion 不得覆盖历史记录 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_H5MD_Writer_Repeated_Output_Completion` | 真实 HighFive 文件读回三条 completion 记录 | 已覆盖，未运行 |
| repeated `Create_Dataset` 后 append 必须保持可用 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_H5MD_Writer_Repeated_Output_Completion` | 连续两次 `Write_Output_Completion` 后读回 `[0,1,2]` | 已覆盖，未运行 |

## 38. observable-only missing value 失败路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| observable-only 缺失 observable value 必须写出 failed 状态 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_Missing_Value_With_Real_Backend` | 真实 HighFive 文件读回 output status | 已覆盖，未运行 |
| observable-only 缺失 observable value 必须写出明确错误原因 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_Missing_Value_With_Real_Backend` | 真实 HighFive 文件读回 output error | 已覆盖，未运行 |
| observable-only 失败 frame 不得推进 completion | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_Missing_Value_With_Real_Backend` | 真实 HighFive 文件读回 `frame_count=0` 与 `last_step=-1` | 已覆盖，未运行 |

## 39. trajectory observable missing value 行为审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory observable 缺失 value 不得推进 observable frame count | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_Observable_Missing_Value_Does_Not_Advance` | mock backend 下检查 `Observable_Frame_Count()==0` | 已覆盖，未运行 |
| trajectory observable 缺失 value 必须记录 Last_Error | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_Observable_Missing_Value_Does_Not_Advance` | 检查 `Last_Error()` 文本 | 已覆盖，未运行 |
| 当前 trajectory writer 缺失 observable value 不标记 bundle failed | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_Observable_Missing_Value_Does_Not_Advance` | mock backend status 保持 open 且无 output error | 已覆盖，未运行 |

## 40. HighFive string overwrite 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| string dataset 重复写入必须覆盖旧值 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_String_Overwrite` | 真实 HighFive 文件读回 `path::mdinfo_text == new mdinfo` | 已覆盖，未运行 |
| string-array dataset 重复写入必须覆盖旧数组 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_String_Overwrite` | 真实 HighFive 文件读回 legacy key 数组为两项新值 | 已覆盖，未运行 |

## 41. RestartH5Writer custom run metadata 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart bundle run metadata 应支持自定义 state type | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_Custom_Run_Metadata_With_Real_Backend` | 真实 `.spgr.h5` 读回 `/run/state_type=checkpoint` | 已覆盖，未运行 |
| restart bundle run metadata 应独立记录 step/time | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_Custom_Run_Metadata_With_Real_Backend` | 真实 `.spgr.h5` 读回 `/run/current_step` 与 `/run/current_time` | 已覆盖，未运行 |

## 42. VDS shard manifest 字段审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper manifest 必须写出 shard index | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `path::shard_manifest_index` dataset 与 append count | 已覆盖，未运行 |
| VDS wrapper manifest 必须写出 frame range | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 frame start/count datasets 与 append count | 已覆盖，未运行 |
| VDS wrapper manifest 必须写出 step range | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 step start/end datasets 与 append count | 已覆盖，未运行 |
| VDS wrapper manifest 必须写出 time range | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 time start/end datasets 与 append count | 已覆盖，未运行 |
| module frame counts 当前不属于 wrapper manifest 持久化字段 | `SPONGE/utils/h5md/vds_trajectory_h5_writer.hpp` / `Write_Manifest_To_Wrapper` | 记录当前实现只写 path/status/index/frame/step/time | 记录 |

## 43. HighFive invalid dataset definition 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| HDF5 dataset path 不能为空 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Rejects_Invalid_Operations` | 真实 HighFive backend 调用空 path `Create_Dataset` 并要求失败 | 已覆盖，未运行 |
| HDF5 dataset dims 不能为空 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Rejects_Invalid_Operations` | 真实 HighFive backend 调用空 dims `Create_Dataset` 并要求失败 | 已覆盖，未运行 |

## 44. HighFive close 后误写保护审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| HDF5 backend close 后不得接受写入 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Status_State` | close 后调用 `Write_String` 并要求失败 | 已覆盖，未运行 |
| close 后误写必须设置 failed 状态和错误信息 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_HighFive_Backend_Status_State` | 检查 `Status()==failed` 与 `Last_Error()` | 已覆盖，未运行 |

## 45. VDS wrapper parameter write precondition 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 未 open 时不得写 metad diagnostic | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Open_Precondition_Errors` | 调用 `Write_Metadynamics_Diagnostic` 并检查错误消息 | 已覆盖，未运行 |
| VDS wrapper 未 open 时不得写 QC SCF text | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Open_Precondition_Errors` | 调用 `Write_Qc_Scf_Output` 并检查错误消息 | 已覆盖，未运行 |
| VDS wrapper 未 open 时不得写 legacy sidecar provenance | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Open_Precondition_Errors` | 调用 `Write_Legacy_Sidecar_Paths` 并检查错误消息 | 已覆盖，未运行 |

## 46. output H5 suffix 精确匹配审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| 推荐 suffix 匹配不得接受短字符串 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Suffix_And_Shard_Derivation` | `Ends_With("h5md", ".spg.h5md") == false` | 已覆盖，未运行 |
| 推荐 suffix 匹配不得接受不完整 suffix | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Suffix_And_Shard_Derivation` | `Ends_With("x.spg.h5m", ".spg.h5md") == false` | 已覆盖，未运行 |
| 推荐 suffix 匹配不得接受 suffix 后追加扩展名 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Suffix_And_Shard_Derivation` | `Ends_With("x.spg.h5md.tmp", ".spg.h5md") == false` | 已覆盖，未运行 |

## 47. OutputCompletionTracker 状态机错误路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| 未打开 frame 时不得完成 frame | `tests/h5_bundle/test_completion_tracker.cpp` / `Test_State_Machine_Error_Paths_Without_Writer` | 调用 `Complete_Frame` 并检查错误消息 | 已覆盖，未运行 |
| 未打开 frame 的完成失败不得推进 completion 计数 | `tests/h5_bundle/test_completion_tracker.cpp` / `Test_State_Machine_Error_Paths_Without_Writer` | 检查 `frame_count=0` 与 `last_complete_step=-1` | 已覆盖，未运行 |
| frame 未完成时不得关闭输出 | `tests/h5_bundle/test_completion_tracker.cpp` / `Test_State_Machine_Error_Paths_Without_Writer` | 调用 `Mark_Closing` 并检查错误消息 | 已覆盖，未运行 |
| close 失败后状态机应保持 open 且 frame incomplete | `tests/h5_bundle/test_completion_tracker.cpp` / `Test_State_Machine_Error_Paths_Without_Writer` | 检查 `status=open` 与 `Has_Incomplete_Frame()` | 已覆盖，未运行 |
| 显式 failed 状态必须记录 failure reason | `tests/h5_bundle/test_completion_tracker.cpp` / `Test_State_Machine_Error_Paths_Without_Writer` | 调用 `Mark_Failed` 并检查 `status=failed` 与 `Last_Error()` | 已覆盖，未运行 |

## 48. grouped output_h5 parser-visible key 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory H5 path key 必须为 grouped flatten 名 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating` | 检查 `kTrajectoryPathKey == output_h5_trajectory_path` | 已覆盖，未运行 |
| trajectory VDS 相关 key 必须为 grouped flatten 名 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating` | 检查 vds/chunk_size/repair_policy 三个常量 | 已覆盖，未运行 |
| restart H5 path key 必须为 grouped flatten 名 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating` | 检查 `kRestartPathKey == output_h5_restart_path` | 已覆盖，未运行 |
| observable H5 path key 必须为 grouped flatten 名 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating` | 检查 `kObservablePathKey == output_h5_observable_path` | 已覆盖，未运行 |
| 旧 output_trajectory_h5/output_traj_h5/output_restart_h5/output_observable_h5 命名不得启用新 bundle | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating` | 用旧 key 设置 controller 后检查 `any_h5_output_enabled=false` | 已覆盖，未运行 |
| 任一新 H5 path key 启用后 legacy sidecar 默认关闭 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating` | 分别设置 trajectory/restart/observable path 并检查 legacy gating | 已覆盖，未运行 |
| H5 bundle 启用后显式 legacy path 仍保留 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Grouped_Output_H5_Key_Names_And_Legacy_Gating` | 设置 `mdout` 后检查 requested/enabled 均为 true | 已覆盖，未运行 |

## 49. RestartH5Writer base layout 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart bundle 必须创建 `/run` | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Writer_Base_Layout_Paths` | mock backend 检查 group 创建记录 | 已覆盖，未运行 |
| restart bundle 必须创建 H5MD particle 根路径 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Writer_Base_Layout_Paths` | 检查 `/particles/all` group | 已覆盖，未运行 |
| restart bundle 必须创建 position/velocity/box/edges 子树 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Writer_Base_Layout_Paths` | 检查四个 particle 子 group | 已覆盖，未运行 |
| restart bundle 必须创建 restart parameter 根路径 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Writer_Base_Layout_Paths` | 检查 `path::parameters_restart` group | 已覆盖，未运行 |
| restart bundle 必须创建 thermostat/barostat/bias/SITS/metad 子树 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Writer_Base_Layout_Paths` | 检查对应 `path::restart_*` groups | 已覆盖，未运行 |
| restart bundle schema metadata 必须标识 restart container | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Writer_Base_Layout_Paths` | 检查 schema name/version string | 已覆盖，未运行 |

## 50. Trajectory/Observable H5 writer base layout 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory bundle 打开时必须创建 H5MD particle 根路径 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_And_Observable_Base_Layout_Paths` | mock backend 检查 `/particles/all` group | 已覆盖，未运行 |
| trajectory bundle 打开时必须创建 position/velocity/force/box/edges 子树 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_And_Observable_Base_Layout_Paths` | 检查五个 particle 子 group | 已覆盖，未运行 |
| trajectory bundle 打开时必须创建 observable/mdout/log 子树 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_And_Observable_Base_Layout_Paths` | 检查 `/observables/all`、mdout columns 和 log group | 已覆盖，未运行 |
| trajectory bundle schema metadata 必须标识 output H5MD | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_And_Observable_Base_Layout_Paths` | 检查 schema name/version string | 已覆盖，未运行 |
| observable-only bundle 不得创建 `/particles` | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_And_Observable_Base_Layout_Paths` | mock backend 检查 `/particles` group 不存在 | 已覆盖，未运行 |
| observable-only bundle 打开时必须创建 observable/mdout/log 子树 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_And_Observable_Base_Layout_Paths` | 检查 `/observables/all`、mdout columns 和 log group | 已覆盖，未运行 |
| observable-only bundle schema metadata 必须标识 output H5MD | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Trajectory_And_Observable_Base_Layout_Paths` | 检查 schema name/version string | 已覆盖，未运行 |

## 51. Metadynamics diagnostic parameter 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| metad diagnostic 必须创建参数根 group | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | mock backend 检查 `module_path::metad_parameter_root` group | 已覆盖，未运行 |
| metad diagnostic 必须创建 bias 名称子 group | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | mock backend 检查 `<metad_parameter_root>/meta0` group | 已覆盖，未运行 |
| metad hills 输出必须进入参数区 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | 检查 `<metad_parameter_root>/meta0/hills` string | 已覆盖，未运行 |
| metad history 输出必须进入参数区 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | 检查 `<metad_parameter_root>/meta0/history` string | 已覆盖，未运行 |
| metad edge 输出必须进入参数区 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | 检查 `<metad_parameter_root>/meta0/edge` string | 已覆盖，未运行 |
| metad potential export 输出必须进入参数区 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | 检查 `<metad_parameter_root>/meta0/potential_export` string | 已覆盖，未运行 |
| metad direct export 输出必须进入参数区 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Metadynamics_And_Diagnostics` | 检查 `<metad_parameter_root>/meta0/direct_export` string | 已覆盖，未运行 |

## 52. QC optional spin_square observable 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| QC observables 必须创建 QC 根 group | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Qc_Optional_Spin_Square_Path` | mock backend 检查 `module_path::qc_root` group | 已覆盖，未运行 |
| QC observables 默认必须创建 energy value 路径 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Qc_Optional_Spin_Square_Path` | 检查 `<qc_root>/energy/value` dataset | 已覆盖，未运行 |
| 未请求 spin_square 时不得创建 spin_square value 路径 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Qc_Optional_Spin_Square_Path` | 检查 `<qc_root>/spin_square/value` dataset 不存在 | 已覆盖，未运行 |

## 53. Restart dynamic extension group 与 legacy provenance 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart SITS state 必须创建动态 module group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Module_State_And_Legacy_Provenance` | mock backend 检查 `<restart_sits>/sits_a` group | 已覆盖，未运行 |
| restart metad state 必须创建动态 bias group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Module_State_And_Legacy_Provenance` | mock backend 检查 `<restart_meta>/meta0` group | 已覆盖，未运行 |
| restart legacy sidecar provenance 必须创建 files group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Module_State_And_Legacy_Provenance` | 检查 `/parameters/sponge/files` group | 已覆盖，未运行 |
| restart legacy sidecar provenance 必须创建 legacy_sidecars group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Restart_Module_State_And_Legacy_Provenance` | 检查 `path::legacy_sidecars` group | 已覆盖，未运行 |
| trajectory legacy sidecar provenance 必须创建 files group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Legacy_Provenance_On_Trajectory_And_Observable` | 检查 `/parameters/sponge/files` group | 已覆盖，未运行 |
| trajectory legacy sidecar provenance 必须创建 legacy_sidecars group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Legacy_Provenance_On_Trajectory_And_Observable` | 检查 `path::legacy_sidecars` group | 已覆盖，未运行 |
| observable-only legacy sidecar provenance 必须创建 files group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Legacy_Provenance_On_Trajectory_And_Observable` | 检查 `/parameters/sponge/files` group | 已覆盖，未运行 |
| observable-only legacy sidecar provenance 必须创建 legacy_sidecars group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Legacy_Provenance_On_Trajectory_And_Observable` | 检查 `path::legacy_sidecars` group | 已覆盖，未运行 |

## 54. VDS trajectory chunk size 参数路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 必须持久化 trajectory chunk size | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `/parameters/sponge/output/trajectory_chunk_size` string | 已覆盖，未运行 |
| chunk size 参数值必须来自 resolved output plan | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | `trajectory.chunk_size=2` 时检查写入字符串 `"2"` | 已覆盖，未运行 |

## 55. VDS wrapper shard manifest group 与 schema metadata 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 必须创建 shard manifest group | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `path::shard_manifest` group | 已覆盖，未运行 |
| VDS wrapper schema name 必须标识 output H5MD | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 检查 `/parameters/sponge/schema/name` string | 已覆盖，未运行 |
| VDS wrapper schema version 必须来自 open 参数 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 检查 `/parameters/sponge/schema/version == test` | 已覆盖，未运行 |

## 56. VDS shard file path 与 schema metadata 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| 第一个 VDS shard 必须使用固定 segment 文件名 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock shard0 检查 `opened_path` 为 `segment_000000.spg.h5md` | 已覆盖，未运行 |
| 第二个 VDS shard 必须使用固定 segment 文件名 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock shard1 检查 `opened_path` 为 `segment_000001.spg.h5md` | 已覆盖，未运行 |
| VDS shard schema name 必须标识 output H5MD | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 检查两个 shard 的 `/parameters/sponge/schema/name` | 已覆盖，未运行 |
| VDS shard schema version 必须来自 open 参数 | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 检查两个 shard 的 `/parameters/sponge/schema/version == test` | 已覆盖，未运行 |

## 57. Real HighFive VDS shard 文件读回审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| 真实 VDS 第一个 shard 文件必须存在 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | `std::filesystem::exists(segment_000000.spg.h5md)` | 已覆盖，未运行 |
| 真实 VDS 第二个 shard 文件必须存在 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | `std::filesystem::exists(segment_000001.spg.h5md)` | 已覆盖，未运行 |
| 真实 VDS shard schema name 必须标识 output H5MD | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 直接打开两个 shard HDF5 文件读回 schema name | 已覆盖，未运行 |
| 真实 VDS shard schema version 必须来自 open 参数 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 直接打开两个 shard HDF5 文件读回 schema version | 已覆盖，未运行 |
| 真实 VDS shard 必须包含标准 position value dataset | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 两个 shard 均检查 `path::position_value` 存在 | 已覆盖，未运行 |

## 58. Observable-only provenance 与 mdout columns 路径审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| observable-only provenance 写入必须创建 provenance group | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Observable_Only_Writer` | mock backend 检查 `/parameters/sponge/provenance` group | 已覆盖，未运行 |
| observable-only provenance 写入必须落到 launch_id path | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Observable_Only_Writer` | mock backend 检查 `/parameters/sponge/provenance/launch_id` string | 已覆盖，未运行 |
| real observable-only 必须写出 mdout original_name column mapping | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/mdout/columns/original_name` | 已覆盖，未运行 |
| real observable-only 必须写出 mdout hdf5_name column mapping | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Observable_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/mdout/columns/hdf5_name` | 已覆盖，未运行 |

## 59. Trajectory single-file schema 与 mdout columns 真实读回审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory 单文件 schema name 必须标识 output H5MD | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/schema/name` | 已覆盖，未运行 |
| trajectory 单文件 schema version 必须来自 open 参数 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/schema/version == test` | 已覆盖，未运行 |
| trajectory 单文件必须创建 log 参数 group | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 检查 `/parameters/sponge/log` 存在 | 已覆盖，未运行 |
| trajectory 单文件必须写出 mdout original_name column mapping | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/mdout/columns/original_name` | 已覆盖，未运行 |
| trajectory 单文件必须写出 mdout hdf5_name column mapping | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/mdout/columns/hdf5_name` | 已覆盖，未运行 |

## 60. Restart real-backend schema 与 extension group 真实读回审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart bundle schema name 必须标识 restart H5 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/schema/name` | 已覆盖，未运行 |
| restart bundle schema version 必须来自 open 参数 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 真实 HDF5 读回 `/parameters/sponge/schema/version == test` | 已覆盖，未运行 |
| restart bundle 必须创建 SITS extension root | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 检查 `path::restart_sits` group 存在 | 已覆盖，未运行 |
| restart bundle 必须创建 metad extension root | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 检查 `path::restart_meta` group 存在 | 已覆盖，未运行 |
| restart bundle 必须创建 legacy sidecar provenance root | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 检查 `path::legacy_sidecars` group 存在 | 已覆盖，未运行 |

## 61. Real HighFive VDS manifest numeric fields 读回审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| 真实 VDS wrapper 必须写出 shard index manifest | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `path::shard_manifest_index` | 已覆盖，未运行 |
| 真实 VDS wrapper 必须写出 frame start manifest | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `path::shard_manifest_frame_start` | 已覆盖，未运行 |
| 真实 VDS wrapper 必须写出 step range manifest | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 step start/end datasets | 已覆盖，未运行 |
| 真实 VDS wrapper 必须写出 time range manifest | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 time start/end datasets | 已覆盖，未运行 |
| manifest numeric values 必须对应两个单帧 shard | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 检查 index/frame_start/step/time 数值 | 已覆盖，未运行 |

## 62. Restart real-backend shared step/time link 路径读回审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| restart position element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 检查 `path::position_step` 和 `path::position_time` 存在并读回 | 已覆盖，未运行 |
| restart velocity element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 检查 `path::velocity_step` 和 `path::velocity_time` 存在并读回 | 已覆盖，未运行 |
| restart box edges element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | 检查 `path::box_edges_step` 和 `path::box_edges_time` 存在并读回 | 已覆盖，未运行 |
| restart element step/time link 读回必须匹配 structural state | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Restart_Writer_With_Real_Backend` | position/velocity/box step 为 20，time 为 1.0 | 已覆盖，未运行 |

## 63. Trajectory real-backend shared step/time link 路径读回审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory position element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 检查 `path::position_step` 和 `path::position_time` 存在并读回 | 已覆盖，未运行 |
| trajectory velocity element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 检查 `path::velocity_step` 和 `path::velocity_time` 存在并读回 | 已覆盖，未运行 |
| trajectory force element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 检查 `path::force_step` 和 `path::force_time` 存在并读回 | 已覆盖，未运行 |
| trajectory box edges element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 检查 `path::box_edges_step` 和 `path::box_edges_time` 存在并读回 | 已覆盖，未运行 |
| trajectory observable element 必须暴露 step/time 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | 检查 `/observables/all/temperature/step` 和 `/observables/all/temperature/time` 存在并读回 | 已覆盖，未运行 |
| trajectory element step/time link 读回必须匹配 frame state | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Trajectory_Writer_With_Real_Backend` | position/velocity/force/box/temperature step 为 10，time 为 0.5 | 已覆盖，未运行 |

## 64. Observable-only real-backend shared step/time link 路径读回审计补充

| 审计项 | 覆盖文件 | 覆盖方式 | 状态 |
|---|---|---|---|
| 普通 observable 的元素级 `step/time` 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` | 真实 HighFive 读回 `/observables/all/temperature/step` 和 `/time` | 已写入，未执行 |
| QC energy 与 spin-square 的 `step/time` 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` | 真实 HighFive 读回 `/observables/all/qc/{energy,spin_square}/step,time` | 已写入，未执行 |
| Nose-Hoover chain observable 的 `step/time` 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` | 真实 HighFive 读回 coordinate/velocity 的 module mapping 路径 | 已写入，未执行 |
| SITS `nk` observable 的 `step/time` 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` | 真实 HighFive 读回 `Sits_Nk_Step_Path` 与 `Sits_Nk_Time_Path` | 已写入，未执行 |
| metadynamics scalar 的 `step/time` 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` | 真实 HighFive 读回 `meta`、`rbias`、`rct` 的 `step/time` | 已写入，未执行 |
| ReaxFF observable 的 `step/time` 路径 | `tests/h5_bundle/test_highfive_backend_io.cpp` | 真实 HighFive 读回 `bond`、`angle` 的 `step/time` | 已写入，未执行 |

说明：该审计补充面向 `*.obs.spg.h5md`，用于确认 observable-only bundle
虽然不写 `/particles` 轨迹字段，但仍提供可独立分析的 H5MD observable
时间轴契约。

## 65. Public H5MD/module path constants 契约锁定审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| H5MD common layout path 常量必须固定 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Public_H5MD_Path_Constants` | 直接断言 `/h5md`、`/particles`、`/observables`、`/parameters`、`/parameters/sponge/*` 字符串 | 已写入，未执行 |
| output completion/status path 常量必须固定 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Public_H5MD_Path_Constants` | 直接断言 `frame_count`、`last_complete_step`、`last_complete_time`、`status` 路径 | 已写入，未执行 |
| VDS shard/shard manifest path 常量必须固定 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Public_H5MD_Path_Constants` | 直接断言 shard metadata 与 `output/shard_manifest/*` 路径 | 已写入，未执行 |
| particle element `value/step/time` path 常量必须固定 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Public_H5MD_Path_Constants` | 直接断言 position、velocity、force、box edges 路径 | 已写入，未执行 |
| restart state/extension path 常量必须固定 | `tests/h5_bundle/test_h5md_writers_with_mock_backend.cpp` / `Test_Public_H5MD_Path_Constants` | 直接断言 `/run/*`、`/parameters/restart/*`、thermostat、barostat、bias、SITS、metadynamics 路径 | 已写入，未执行 |
| module extension path 常量必须固定 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Module_Path_Constants` | 直接断言 NHC、SITS、metadynamics、QC、ReaxFF canonical paths | 已写入，未执行 |
| SITS 动态 nk path 生成规则必须固定 | `tests/h5_bundle/test_module_h5_mappings_with_mock_backend.cpp` / `Test_Module_Path_Constants` | 直接断言 `Sits_Nk_Root/Value/Step/Time_Path("sits_a")` 结果 | 已写入，未执行 |

说明：该审计补充锁定的是公开路径字符串契约，用于覆盖“新路径是否稳定”这一类
风险；真实 HDF5 写入、hard link、VDS manifest 和 observable-only 读回仍由对应
HighFive 后端测试覆盖。

## 66. Legacy sidecar resolution matrix 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| null controller 下 legacy sidecar 默认启用 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Resolve_Legacy_Output_Plan_Matrix` | 直接调用 `Resolve_Legacy_Output_Plan(nullptr)` 并检查 8 个 key | 已写入，未执行 |
| 空 controller 下 legacy sidecar 默认启用 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Resolve_Legacy_Output_Plan_Matrix` | 检查 8 个 key 均 enabled 且非 explicit | 已写入，未执行 |
| H5 bundle trajectory path 启用后 legacy sidecar 默认关闭 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Resolve_Legacy_Output_Plan_Matrix` | 设置 `output_h5_trajectory_path` 后检查 8 个 key 均 disabled | 已写入，未执行 |
| H5 bundle observable path 启用后显式 legacy sidecar 可恢复 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Resolve_Legacy_Output_Plan_Matrix` | 设置 `output_h5_observable_path` 与全部 8 个 legacy path，检查 enabled/explicit/path | 已写入，未执行 |
| legacy sidecar key 集合和顺序必须固定 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Resolve_Legacy_Output_Plan_Matrix` | 检查 sidecars size 为 8，并逐项比较 key/path | 已写入，未执行 |

说明：该审计补充覆盖的是 H5 bundle 输出路径启用后的 legacy 输出抑制策略，和
writer 层的 legacy sidecar provenance 写入测试互补。

## 67. VDS repair metadata 与 manifest string path 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 必须创建 `repaired_shard_count` dataset | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock backend 检查 dataset path 存在并 append 一次 | 已写入，未执行 |
| VDS 无帧 finalize 仍必须写 repair metadata | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Finalize_Without_Frames` | 检查 `repair_policy`、`repair_status` 和 `repaired_shard_count` path | 已写入，未执行 |
| VDS complete-prefix repair 必须写 repair metadata | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Complete_Prefix_Repair_Finalize` | 检查 repair applied 时 `repaired_shard_count` path 存在 | 已写入，未执行 |
| VDS repaired prefix manifest 必须保留 status string array | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Complete_Prefix_Repair_Retains_Valid_Prefix` | 检查 `path::shard_manifest_status` 长度和 `complete` 值 | 已写入，未执行 |
| VDS normal manifest 必须写 path/status string arrays | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 检查 manifest path/status array 长度和 complete 状态 | 已写入，未执行 |
| 真实 HighFive VDS wrapper 必须读回 repaired shard count | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 检查 path 存在并读回 `0` | 已写入，未执行 |

说明：该审计补充覆盖 VDS wrapper 的 repair metadata numeric path，与此前已覆盖的
`repair_policy`、`repair_status`、`vds_status` 和 manifest numeric fields 互补。

## 68. H5 bundle test manifest 与 README 覆盖矩阵同步审计

| 契约项 | 文档位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| README test target scope 必须反映新增路径测试职责 | `tests/h5_bundle/README.md` | 更新 target scope 表，加入 path constants、legacy matrix、VDS repair metadata 和 real step/time 覆盖 | 已写入，未执行 |
| README coverage matrix 必须列出新增路径契约 | `tests/h5_bundle/README.md` | 增加 public path constants、legacy 8-key matrix、observable-only/trajectory/restart step-time、VDS repair metadata 条目 | 已写入，未执行 |
| VDS backend-IO coverage 描述必须包含 repair metadata | `tests/h5_bundle/README.md` | 在 Backend-IO VDS writer coverage 中加入 chunk size、vds status、repair policy/status、repaired shard count | 已写入，未执行 |
| Observable-only backend-IO coverage 描述必须包含 module step/time | `tests/h5_bundle/README.md` | 在 observable-only facade coverage 中加入 NHC/SITS/metad/QC/ReaxFF `step/time` 读回 | 已写入，未执行 |
| TEST_TARGETS build target scope 必须反映新增路径测试职责 | `tests/h5_bundle/TEST_TARGETS.md` | 更新 aggregate target purpose 描述 | 已写入，未执行 |
| TEST_TARGETS 必须说明 expanded path-contract coverage | `tests/h5_bundle/TEST_TARGETS.md` | 追加 `Expanded path-contract coverage` 段落 | 已写入，未执行 |

说明：该审计补充面向测试套件可维护性，不作为构建或运行通过的证据；实际通过状态仍需后续配置、编译和 CTest 验证。

## 69. Output plan bundle path 组合语义与 parser helper 边界审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| trajectory selector/config key 不应单独启用 H5 output | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Output_Selectors_Do_Not_Enable_H5_Without_Path` | 设置 vds、chunk_size、repair_policy 但不设置 path，检查 `any_h5_output_enabled == false` | 已写入，未执行 |
| trajectory path key 只启用 trajectory bundle | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_H5_Output_Path_Keys_Enable_Only_Their_Bundle` | 设置 `output_h5_trajectory_path` 并检查 restart/observable disabled | 已写入，未执行 |
| restart path key 只启用 restart bundle | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_H5_Output_Path_Keys_Enable_Only_Their_Bundle` | 设置 `output_h5_restart_path` 并检查 trajectory/observable disabled | 已写入，未执行 |
| observable path key 只启用 observable bundle | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_H5_Output_Path_Keys_Enable_Only_Their_Bundle` | 设置 `output_h5_observable_path` 并检查 trajectory/restart disabled | 已写入，未执行 |
| 三类 H5 bundle path 可同时启用 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_All_H5_Output_Bundles_Can_Be_Enabled_Together` | 同时设置 trajectory/restart/observable path 并检查三者 enabled | 已写入，未执行 |
| 三类 H5 bundle 推荐后缀 flag 必须同时计算 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_All_H5_Output_Bundles_Can_Be_Enabled_Together` | 检查 trajectory/restart/observable `has_recommended_suffix` | 已写入，未执行 |
| 未知 bool 文本必须解析为 false | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Bool_Parsing_Text_Variants` | 检查 `maybe` 和空字符串 | 已写入，未执行 |
| repair policy 大小写规范化 helper 必须稳定 | `tests/h5_bundle/test_h5_output_plan.cpp` / `Test_Bool_Parsing_Text_Variants` | 检查 `Lowercase("CoMpLeTe_PrEfIx") == complete_prefix` | 已写入，未执行 |

说明：该审计补充面向 parser-visible grouped `output_h5_*` key 的组合边界，不涉及真实 HDF5 IO。

## 70. VDS observable axis 与 mdout column mapping 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 必须创建 observable shared step/time VDS | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `path::observables_all_step` 与 `path::observables_all_time` virtual datasets | 已写入，未执行 |
| VDS observable element 必须暴露 step/time hard links | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 检查 temperature `step/time` hard link 指向 shared observable axis | 已写入，未执行 |
| VDS wrapper 必须写 mdout original/hdf5 column mapping | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `original_name == TEMP`、`hdf5_name == temperature` | 已写入，未执行 |
| 真实 HighFive VDS observable step/time 必须可读 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 temperature step `10,20` 和 time `0.1,0.2` | 已写入，未执行 |
| 真实 HighFive VDS wrapper 必须读回 mdout column mapping | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 original/hdf5 column arrays | 已写入，未执行 |

说明：该审计补充覆盖 VDS wrapper 的 observable metadata 路径，与此前 single-file trajectory 和 observable-only 的 step/time、mdout metadata 覆盖对齐。

## 71. VDS particle axis 与 element step/time 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS wrapper 必须创建 particle shared step/time VDS | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `path::particles_all_step` 与 `path::particles_all_time` virtual datasets | 已写入，未执行 |
| VDS wrapper 必须创建 box edges VDS | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `path::box_edges_value` virtual dataset | 已写入，未执行 |
| VDS particle elements 必须暴露 step/time hard links | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 检查 position、velocity、box edges 的 hard links 指向 shared particle axis | 已写入，未执行 |
| 真实 HighFive VDS position step/time 必须可读 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 position step `10,20` 和 time `0.1,0.2` | 已写入，未执行 |
| 真实 HighFive VDS box edges step/time 必须可读 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 box step `10,20` 和 time `0.1,0.2` | 已写入，未执行 |

说明：该审计补充覆盖 VDS wrapper 的 particle axis 语义，与 VDS observable axis、single-file trajectory 和 restart bundle 的 element-level `step/time` 覆盖对齐。

## 72. VDS module observable axis 与 element step/time 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS NHC module 必须创建 shared step/time VDS | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `module_path::nhc_step/time` virtual datasets | 已写入，未执行 |
| VDS NHC coordinate/velocity 必须暴露 step/time hard links | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 coordinate/velocity hard links | 已写入，未执行 |
| VDS SITS nk 必须创建 step/time VDS | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `Sits_Nk_Step_Path` 与 `Sits_Nk_Time_Path` | 已写入，未执行 |
| VDS metadynamics scalar 必须暴露 step/time hard links | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `meta`、`rbias`、`rct` hard links | 已写入，未执行 |
| VDS QC energy 必须暴露 step/time hard links | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `energy` hard links | 已写入，未执行 |
| VDS ReaxFF term 必须暴露 step/time hard links | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | mock wrapper 检查 `bond` hard links | 已写入，未执行 |
| 真实 HighFive VDS module step/time 必须可读 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 NHC、SITS、metad、QC、ReaxFF step `10,20` 和 time `0.1,0.2` | 已写入，未执行 |

说明：该审计补充覆盖 VDS wrapper 的 module observable axis，与 single-file trajectory、observable-only 和 VDS ordinary observable 的 `step/time` 覆盖对齐。

## 73. VDS conditional module branches 审计补充

| 契约项 | 测试位置 | 覆盖方式 | 状态 |
|---|---|---|---|
| VDS QC spin-square 分支必须 materialize | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | `Ensure_Qc_Observables(true)` 后检查 `spin_square/value` VDS 与 `step/time` hard links | 已写入，未执行 |
| 真实 VDS QC spin-square 必须可读 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `spin_square` 值 `0.11,0.22` 和 step/time | 已写入，未执行 |
| VDS ReaxFF 多 term 必须 materialize | `tests/h5_bundle/test_vds_trajectory_writer_with_mock_backend.cpp` / `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | terms `{bond, angle}` 后检查 `angle/value` VDS 与 `step/time` hard links | 已写入，未执行 |
| 真实 VDS ReaxFF angle term 必须可读 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `angle` 值 `3.5,4.5` 和 step/time | 已写入，未执行 |
| 真实 VDS metadynamics rbias/rct 必须可读 | `tests/h5_bundle/test_highfive_backend_io.cpp` / `Test_Vds_Trajectory_Writer_With_Real_Backend` | 真实 HDF5 读回 `rbias`、`rct` 值和 step/time | 已写入，未执行 |

说明：该审计补充覆盖 VDS module materialization 的 optional/multi-term 分支，与此前基础 module path 和 axis 覆盖互补。
## 74. Optional trajectory/restart particle path negative tests 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| trajectory 关闭 velocity/force 时不创建对应 H5MD 数据集 | `Test_Trajectory_Optional_Velocity_And_Force_Paths` | 已覆盖 |
| trajectory 关闭 velocity/force 时不创建对应 step/time hard-link | `Test_Trajectory_Optional_Velocity_And_Force_Paths` | 已覆盖 |
| trajectory 关闭 velocity/force 时不发生 append 写入 | `Test_Trajectory_Optional_Velocity_And_Force_Paths` | 已覆盖 |
| trajectory 关闭 optional 字段后 position/box 仍保持标准 H5MD path 与 hard-link | `Test_Trajectory_Optional_Velocity_And_Force_Paths` | 已覆盖 |
| restart 关闭 velocity 时不创建 velocity 数据集 | `Test_Restart_Optional_Velocity_Path` | 已覆盖 |
| restart 关闭 velocity 后 position/box 仍保持标准 H5MD path 与 hard-link | `Test_Restart_Optional_Velocity_Path` | 已覆盖 |
| restart 关闭 velocity 后仍记录 restart state type | `Test_Restart_Optional_Velocity_Path` | 已覆盖 |

剩余风险：本轮仅补充 mock backend 层的 schema/路径契约测试，尚未执行编译或运行验证；真实 HighFive 后端的 optional disabled 分支是否完全通过，需要后续显式授权后再运行对应测试目标确认。

## 75. Observable-only bundle mock path coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| observable-only 文件不创建 `/particles` | `Test_Observable_Only_Writer`, `Test_Observable_Only_Module_Proxy_Paths` | 已覆盖 |
| 多个 mdout observable 的 value/step/time 路径 | `Test_Observable_Only_Writer` | 已覆盖 |
| mdout original/hdf5 column name 数组 | `Test_Observable_Only_Writer` | 已覆盖 |
| observable-only completion frame_count 更新 | `Test_Observable_Only_Writer` | 已覆盖 |
| observable-only mdinfo text 路径 | `Test_Observable_Only_Writer` | 已覆盖 |
| observable-only legacy sidecar key/path 路径 | `Test_Observable_Only_Writer` | 已覆盖 |
| observable-only provenance launch_id 路径 | `Test_Observable_Only_Writer` | 已覆盖 |
| observable-only finalized status | `Test_Observable_Only_Writer` | 已覆盖 |
| observable-only NHC coordinate/velocity module observable 路径 | `Test_Observable_Only_Module_Proxy_Paths` | 已覆盖 |
| observable-only SITS nk step/time/value 独立 dataset 路径 | `Test_Observable_Only_Module_Proxy_Paths` | 已覆盖 |
| observable-only metad meta/rbias/rct 路径与 diagnostic parameters 路径 | `Test_Observable_Only_Module_Proxy_Paths` | 已覆盖 |
| observable-only QC energy 路径与 `qc_scf_output` legacy parameters 路径 | `Test_Observable_Only_Module_Proxy_Paths` | 已覆盖 |
| observable-only QC spin_square disabled 负向路径 | `Test_Observable_Only_Module_Proxy_Paths` | 已覆盖 |
| observable-only ReaxFF energy term 路径 | `Test_Observable_Only_Module_Proxy_Paths` | 已覆盖 |

剩余风险：本轮仅补充 mock backend 层测试代码，尚未执行编译或运行验证；真实 HighFive backend 的 observable-only module proxy 组合路径仍需后续授权后用对应测试目标确认。

## 76. HighFive optional particle-field disabled branches 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| 真实 `.spg.h5md` trajectory 关闭 velocity/force 时不创建 velocity value/step/time | `Test_Trajectory_Optional_Particle_Fields_With_Real_Backend` | 已覆盖 |
| 真实 `.spg.h5md` trajectory 关闭 velocity/force 时不创建 force value/step/time | `Test_Trajectory_Optional_Particle_Fields_With_Real_Backend` | 已覆盖 |
| 真实 `.spg.h5md` trajectory 关闭 optional 字段后 position/box value/step/time 仍存在 | `Test_Trajectory_Optional_Particle_Fields_With_Real_Backend` | 已覆盖 |
| 真实 `.spg.h5md` trajectory position dataset shape 与 value 读回 | `Test_Trajectory_Optional_Particle_Fields_With_Real_Backend` | 已覆盖 |
| 真实 `.spg.h5md` trajectory completion metadata 与 finalized status | `Test_Trajectory_Optional_Particle_Fields_With_Real_Backend` | 已覆盖 |
| 真实 `.spgr.h5` restart 关闭 velocity 时不创建 velocity value/step/time | `Test_Restart_Optional_Velocity_With_Real_Backend` | 已覆盖 |
| 真实 `.spgr.h5` restart 关闭 velocity 后 position/box value/step/time 仍存在 | `Test_Restart_Optional_Velocity_With_Real_Backend` | 已覆盖 |
| 真实 `.spgr.h5` restart position dataset shape 与 value 读回 | `Test_Restart_Optional_Velocity_With_Real_Backend` | 已覆盖 |
| 真实 `.spgr.h5` restart `/run/current_step`、`/run/current_time`、`/run/state_type` | `Test_Restart_Optional_Velocity_With_Real_Backend` | 已覆盖 |
| 真实 `.spgr.h5` restart completion metadata 与 finalized status | `Test_Restart_Optional_Velocity_With_Real_Backend` | 已覆盖 |

剩余风险：本轮新增真实 HighFive 后端测试代码，但尚未执行编译或 CTest；是否存在 HighFive API、link existence 或类型读回层面的实际失败，需要后续显式授权后验证。

## 77. VDS optional particle-field disabled branches 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| VDS wrapper 在 `include_force=false` 时不 materialize force value virtual dataset | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper 在 `include_force=false` 时不建立 force step/time hard-link | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS shard 在 `include_force=false` 时不创建 force value dataset | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper 在 `include_velocity=false, include_force=false` 时仍 materialize particle step/time、position、box | `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| VDS wrapper 在 `include_velocity=false, include_force=false` 时不 materialize velocity/force value virtual dataset | `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| VDS wrapper 在 `include_velocity=false, include_force=false` 时不建立 velocity/force step/time hard-link | `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| VDS shard 在 `include_velocity=false, include_force=false` 时只创建 position/box，不创建 velocity/force | `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| VDS manifest status 在 optional disabled 分支中仍为 complete | `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| 真实 HighFive VDS wrapper 在 `include_velocity=false, include_force=false` 时不包含 velocity/force value/step/time | `Test_Vds_Trajectory_Writer_With_Real_Backend` | 已覆盖 |
| 真实 HighFive VDS shard 在 `include_velocity=false, include_force=false` 时不包含 velocity/force value/step/time | `Test_Vds_Trajectory_Writer_With_Real_Backend` | 已覆盖 |

剩余风险：本轮新增 mock 与真实 HighFive 后端测试代码，但尚未执行编译或 CTest；VDS virtual dataset absence 与 shard absence 的实际行为仍需后续显式授权后验证。

## 78. H5 output route helper contract tests 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| SPONGE mdout key 到 H5MD-safe observable name 的字符替换规则 | `Test_Output_Name_Sanitization` | 已覆盖 |
| 空 key 与数字开头 key 的 `_` 前缀规则 | `Test_Output_Name_Sanitization` | 已覆盖 |
| sanitized observable name collision suffix 规则 | `Test_Unique_Output_Name_Collision_Resolution` | 已覆盖 |
| mdout scalar string 到 double 的合法解析路径 | `Test_Output_Double_Parsing` | 已覆盖 |
| mdout scalar string 非法输入不覆盖旧值 | `Test_Output_Double_Parsing` | 已覆盖 |
| null output pointer 解析失败语义 | `Test_Output_Double_Parsing` | 已覆盖 |
| ReaxFF H5 routing key 识别 | `Test_Reaxff_Output_Key_Recognition` | 已覆盖 |
| output key lookup 与 null key 语义 | `Test_Output_Key_Exists` | 已覆盖 |
| optional text sidecar 存在时完整读入 | `Test_Text_File_Read_If_Present` | 已覆盖 |
| optional text sidecar 空路径、缺失文件、null text pointer 失败语义 | `Test_Text_File_Read_If_Present` | 已覆盖 |
| `output.hpp` helper wrapper 继续保留原调用名 | `SPONGE/MD_core/output.hpp` wrapper 委托 | 已覆盖源码结构 |
| contract aggregate target 纳入 output route helper 测试 | `tests/h5_bundle/CMakeLists.txt` | 已覆盖 |

剩余风险：本轮包含一个小型 production helper 抽取与新测试目标，但尚未执行编译或 CTest；`output.hpp` wrapper 与新 header 的实际编译状态仍需后续显式授权后验证。

## 79. HighFive backend factory and repeated dataset contract 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `HighFiveBackendFactory::Create_Backend()` 返回可用 backend | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| factory-created backend 可完成 H5MD common layout 初始化 | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| nested output file path 可被创建 | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| appendable dataset append 后第一维增长 | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| 已存在 dataset 上重复 `Create_Dataset` 不截断旧数据 | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| 重复定义后继续 append 保持记录顺序和值 | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| factory-created backend finalize 写入 finalized status | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |

剩余风险：本轮新增真实 HighFive backend 测试代码，但尚未执行编译或 CTest；nested path 创建与 repeated dataset definition 的实际行为仍需后续显式授权后验证。

## 80. HighFive string-array metadata overwrite coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| string-array metadata 支持 nested parent group 自动创建 | `Test_HighFive_Backend_String_Overwrite` | 已覆盖 |
| string-array metadata 可写入空字符串元素 | `Test_HighFive_Backend_String_Overwrite` | 已覆盖 |
| string-array metadata 可写入带空格字符串 | `Test_HighFive_Backend_String_Overwrite` | 已覆盖 |
| string-array metadata 可写入相对路径字符串 | `Test_HighFive_Backend_String_Overwrite` | 已覆盖 |
| 同一路径重复 `Write_String_Array` 会覆盖旧 dataset | `Test_HighFive_Backend_String_Overwrite` | 已覆盖 |
| 覆盖后旧 string-array 内容不会残留 | `Test_HighFive_Backend_String_Overwrite` | 已覆盖 |
| 覆盖后的空字符串元素可被真实 HDF5 readback 保留 | `Test_HighFive_Backend_String_Overwrite` | 已覆盖 |

剩余风险：本轮新增真实 HighFive backend string-array 测试代码，但尚未执行编译或 CTest；空字符串元素在当前 HighFive/HDF5 版本中的实际读写行为仍需后续显式授权后验证。

## 81. HighFive VDS zero-frame finalize coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| VDS writer 无 frame finalize 时 manifest size 为 0 | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| VDS writer 无 frame finalize 时 trajectory/observable frame count 为 0 | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS wrapper 创建 common H5MD roots | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS wrapper 不创建 particle step/time 与 position/box datasets | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS wrapper 不创建 observable step/time datasets | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS wrapper 不创建 shard manifest index/path datasets | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS finalize 不创建 shard root 目录 | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS wrapper 写入 chunk size metadata | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS wrapper 写入 repair policy/status/repaired count | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |
| zero-frame VDS wrapper 写入 vds_status、completion metadata、finalized status | `Test_Vds_Finalize_Without_Frames_With_Real_Backend` | 已覆盖 |

剩余风险：本轮新增真实 HighFive VDS zero-frame 测试代码，但尚未执行编译或 CTest；zero-frame wrapper metadata 与 shard root absence 的实际行为仍需后续显式授权后验证。

## 82. Restart metad text component coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| restart metad `hills` text component 路径 | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| restart metad `history` text component 路径 | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| restart metad `edge` text component 路径 | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| restart metad `potential_export` text component 路径 | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| restart metad `direct_export` text component 路径 | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| 真实 `.spgr.h5` restart 文件读回 metad component 文本内容 | `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |

剩余风险：本轮新增 mock 与真实 HighFive restart metad component 测试代码，但尚未执行编译或 CTest；实际 `.spgr.h5` 读写行为仍需后续显式授权后验证。

## 83. Restart SITS dynamic state component coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| restart SITS dynamic module group `/parameters/restart/bias/sits/sits_a` | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| restart SITS `nk` state component 路径 | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| restart SITS `weight` state component 路径 | `Test_Restart_Module_State_And_Legacy_Provenance`, `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |
| restart SITS dynamic state append count 与输入长度匹配 | `Test_Restart_Module_State_And_Legacy_Provenance` | 已覆盖 |
| 真实 `.spgr.h5` restart 文件读回 SITS `weight` shape 与 payload | `Test_Restart_Writer_With_Real_Backend` | 已覆盖 |

剩余风险：本轮新增 mock 与真实 HighFive restart SITS dynamic state 测试代码，但尚未执行编译或 CTest；实际 `.spgr.h5` 读写行为仍需后续显式授权后验证。

## 84. HighFive repeated hard-link coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| 同一 target/link path 重复 `Create_Hard_Link` 返回成功 | `Test_HighFive_Backend_Basic_File_Layout` | 已覆盖 |
| 重复 hard-link 后 alias dataset 仍存在 | `Test_HighFive_Backend_Basic_File_Layout` | 已覆盖 |
| 重复 hard-link 后 alias dataset shape 与原 dataset 一致 | `Test_HighFive_Backend_Basic_File_Layout` | 已覆盖 |
| 重复 hard-link 后 alias dataset 可读回原始 payload | `Test_HighFive_Backend_Basic_File_Layout` | 已覆盖 |

剩余风险：本轮新增真实 HighFive hard-link 测试代码，但尚未执行编译或 CTest；实际 HDF5 hard-link 幂等行为仍需后续显式授权后验证。

## 85. HighFive VDS complete-prefix repair coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| 真实 HighFive backend 可按创建序号注入 shard finalize 失败 | `MaybeFailFinalizeHighFiveBackend`, `SelectiveFailHighFiveBackendFactory` | 已覆盖 |
| `Finalize_With_Repair()` 在最后一个 shard 失败时返回成功 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| writer manifest 修剪为完整前缀 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| writer trajectory frame count 回退到完整前缀帧数 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| 真实 wrapper VDS position dataset 只引用完整 shard0 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| 真实 wrapper position step/time 回退到 shard0 的最后完成帧 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| wrapper shard manifest 不包含失败 shard1 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| wrapper output completion metadata 回退到完整前缀 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| repair metadata 写入 `complete_prefix/applied/1` | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| 完整 shard0 仍为可读真实 HDF5 shard | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |
| repair 场景保持 velocity/force disabled 分支 | `Test_Vds_Complete_Prefix_Repair_With_Real_Backend` | 已覆盖 |

剩余风险：本轮新增真实 HighFive VDS repair 测试代码，但尚未执行编译或 CTest；`Finalize_With_Repair()` 在当前 HighFive/HDF5 版本下的实际文件关闭和 VDS materialization 行为仍需后续显式授权后验证。

## 86. H5MD common layout roots/schema/output metadata coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| 普通 H5MD common layout 创建 `/h5md` root | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| 普通 H5MD common layout 创建 `/particles` root | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| 普通 H5MD common layout 创建 `/observables` root | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| 普通 H5MD common layout 创建 `/parameters` root | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| common layout 创建 `/parameters/sponge` root | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| common layout 创建 schema/output metadata roots | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| common layout 写入 schema name/version | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| common layout 初始化 output status 为 `open` | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| common layout 初始化 output frame_count/last_step/last_time datasets | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |
| observable-only common layout 不创建 `/particles` root | `Test_Common_Layout_Roots_And_Output_Metadata` | 已覆盖 |

剩余风险：本轮新增 mock backend common layout 测试代码，但尚未执行编译或 CTest；真实 HighFive common layout 已由 `test_highfive_backend_io` 间接读回覆盖，仍需后续显式授权后统一验证。

## 87. H5 output plan helper null/unknown key coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `Command_String` 接收 null key 时返回空字符串 | `Test_Contract_Helper_Edge_Cases` | 已覆盖 |
| default legacy enabled 场景下 `Enabled(nullptr)` 为 false | `Test_Resolve_Legacy_Output_Plan_Matrix` | 已覆盖 |
| default legacy enabled 场景下 `Explicitly_Requested(nullptr)` 为 false | `Test_Resolve_Legacy_Output_Plan_Matrix` | 已覆盖 |
| default legacy enabled 场景下 unknown key 不被启用 | `Test_Resolve_Legacy_Output_Plan_Matrix` | 已覆盖 |
| H5 bundle 禁用 legacy default 场景下 null key 不被启用 | `Test_Resolve_Legacy_Output_Plan_Matrix` | 已覆盖 |
| H5 bundle 禁用 legacy default 场景下 unknown key 不被启用 | `Test_Resolve_Legacy_Output_Plan_Matrix` | 已覆盖 |
| 全部 explicit legacy sidecar 场景下 null key 不被误判 explicit | `Test_Resolve_Legacy_Output_Plan_Matrix` | 已覆盖 |
| 全部 explicit legacy sidecar 场景下 unknown key 不被误判 explicit | `Test_Resolve_Legacy_Output_Plan_Matrix` | 已覆盖 |

剩余风险：本轮新增 parser helper 边界测试代码，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 88. H5 output route name global uniqueness coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| sanitized output name collision 使用全局唯一结果 | `Test_Unique_Output_Name_Collision_Resolution` | 已覆盖 |
| 用户已有 suffix 名称 `A_1` 时后续 `A` 不重复生成 `A_1` | `Test_Unique_Output_Name_Collision_Resolution` | 已覆盖 |
| 同时存在 sanitized collision `A-1 -> A_1` 时生成 `A_1_1` | `Test_Unique_Output_Name_Collision_Resolution` | 已覆盖 |
| 多次重复 base 名称继续生成递增 suffix | `Test_Unique_Output_Name_Collision_Resolution` | 已覆盖 |
| 纯空白 mdout scalar 解析失败且不改写旧值 | `Test_Output_Double_Parsing` | 已覆盖 |
| `REAXFFBOND` 不进入 ReaxFF module routing | `Test_Reaxff_Output_Key_Recognition` | 已覆盖 |
| empty output key vector lookup 返回 false | `Test_Output_Key_Exists` | 已覆盖 |
| null sidecar filename 不读取且不改写旧文本 | `Test_Text_File_Read_If_Present` | 已覆盖 |

剩余风险：本轮修正了 routing helper 并新增单元测试代码，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 89. Explicit legacy sidecar provenance collection coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| legacy default enabled 时默认 sidecar 不进入 H5 provenance collection | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| 启用 H5 bundle 后只收集 explicit legacy sidecar | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| explicit sidecar collection 会清除 stale output vectors | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| explicit sidecar collection 保持 canonical legacy key 顺序 | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| `mdout` explicit path 被收集为 provenance key/path | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| `rst` explicit path 被收集为 provenance key/path | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| `qc_scf_output` explicit path 被收集为 provenance key/path | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| null output vector pointer 不改写已有 collection | `Test_Explicit_Legacy_Sidecar_Collection` | 已覆盖 |
| trajectory/VDS/observable/restart provenance overload 共用同一 collection helper | `MD_core/output.hpp` | 已覆盖 |

剩余风险：本轮抽取了 runtime provenance helper 并新增单元测试代码，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 90. Writer dataset type/shape/chunk contract coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| trajectory particle step dataset 为 int64 appendable scalar stream | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| trajectory particle time dataset 为 float64 appendable scalar stream | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| trajectory position dataset 为 float32 `[frame,N,3]` stream | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| trajectory velocity dataset 为 float32 `[frame,N,3]` stream | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| trajectory force dataset 为 float32 `[frame,N,3]` stream | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| trajectory box edges dataset 为 float32 `[frame,3,3]` stream | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| trajectory observable scalar dataset 为 float64 stream | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| observable-only step/time dataset type/shape | `Test_Observable_Only_Writer` | 已覆盖 |
| observable-only scalar value dataset type/shape | `Test_Observable_Only_Writer` | 已覆盖 |
| restart structural step/time dataset 为 single-state stream | `Test_Restart_Writer_Is_Single_State` | 已覆盖 |
| restart position/velocity dataset 为 float32 single-state `[1,N,3]` contract | `Test_Restart_Writer_Is_Single_State` | 已覆盖 |
| restart box edges dataset 为 float32 single-state `[1,3,3]` contract | `Test_Restart_Writer_Is_Single_State` | 已覆盖 |
| restart run current step/time dataset type/shape | `Test_Restart_Writer_Is_Single_State` | 已覆盖 |
| restart NHC state dataset type/shape | `Test_Restart_Module_State_And_Legacy_Provenance` | 已覆盖 |
| restart SITS dynamic state dataset type/shape | `Test_Restart_Module_State_And_Legacy_Provenance` | 已覆盖 |

剩余风险：本轮新增 mock backend DatasetSpec 测试代码，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 91. VDS wrapper virtual dataset schema coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| VDS wrapper particle step dataset 为 int64 fixed-frame virtual dataset | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper particle time dataset 为 float64 fixed-frame virtual dataset | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper position dataset 为 float32 `[frame,N,3]` fixed-frame virtual dataset | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper velocity dataset 为 float32 `[frame,N,3]` fixed-frame virtual dataset | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper box edges dataset 为 float32 `[frame,3,3]` fixed-frame virtual dataset | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper observable step/time dataset type/shape/chunk | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper observable scalar value dataset type/shape/chunk | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper NHC coordinate/velocity virtual dataset type/shape/chunk | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper SITS `nk` virtual dataset type/shape/chunk | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper metadynamics scalar virtual dataset type/shape/chunk | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper QC scalar virtual dataset type/shape/chunk | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper ReaxFF term scalar virtual dataset type/shape/chunk | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper disabled velocity/force 分支不创建对应 dataset spec | `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| VDS wrapper virtual dataset 均为 fixed shape 且 `appendable=false` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets`, `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |

剩余风险：本轮新增 VDS wrapper DatasetSpec 测试代码，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 92. Module mapping dataset schema coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| NHC step/time dataset type/shape/chunk | `Test_Nhc_And_Sits_Mappings` | 已覆盖 |
| NHC coordinate/velocity vector dataset type/shape/chunk | `Test_Nhc_And_Sits_Mappings` | 已覆盖 |
| SITS `nk` step/time/value dataset type/shape/chunk | `Test_Nhc_And_Sits_Mappings` | 已覆盖 |
| metadynamics step/time dataset type/shape/chunk | `Test_Metadynamics_And_Diagnostics` | 已覆盖 |
| metadynamics `meta/rbias/rct` scalar value dataset type/shape/chunk | `Test_Metadynamics_And_Diagnostics` | 已覆盖 |
| QC step/time dataset type/shape/chunk | `Test_Qc_And_Reaxff_Mappings`, `Test_Qc_Optional_Spin_Square_Path` | 已覆盖 |
| QC energy/spin_square scalar value dataset type/shape/chunk | `Test_Qc_And_Reaxff_Mappings` | 已覆盖 |
| QC spin_square disabled 分支不创建 spin_square value dataset | `Test_Qc_Optional_Spin_Square_Path` | 已覆盖 |
| ReaxFF step/time dataset type/shape/chunk | `Test_Qc_And_Reaxff_Mappings` | 已覆盖 |
| ReaxFF term scalar value dataset type/shape/chunk | `Test_Qc_And_Reaxff_Mappings` | 已覆盖 |
| module mapping datasets 均保持 appendable stream 语义 | `Test_Nhc_And_Sits_Mappings`, `Test_Metadynamics_And_Diagnostics`, `Test_Qc_And_Reaxff_Mappings` | 已覆盖 |

剩余风险：本轮新增 module mapping DatasetSpec 测试代码，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 93. Output error metadata path constant coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `path::output_error` public constant 映射到 `/parameters/sponge/output/error` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| completion tracker failure metadata 使用 `path::output_error` | `Test_State_Machine_Writes_Metadata` | 已覆盖 |
| trajectory append failure 写入 `path::output_error` | `Test_Trajectory_Append_Failure_Marks_Failed` | 已覆盖 |
| observable-only missing value failure 写入 `path::output_error` | `Test_Observable_Missing_Value_Marks_Failed`, `test_highfive_backend_io` observable failure case | 已覆盖 |
| restart duplicate structural state failure 写入 `path::output_error` | `test_highfive_backend_io` restart single-state failure case | 已覆盖 |
| real HighFive failed metadata readback 使用 `path::output_error` | `Test_HighFive_Backend_Failed_Metadata` | 已覆盖 |

剩余风险：本轮新增 output error path constant 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 94. HighFive backend max-dims/chunk layout coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| appendable dataset 追加后真实 shape 增长为 `[frame,width]` | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| appendable dataset frame 维 max dim 为 `H5S_UNLIMITED` | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| appendable dataset 非 frame 维 max dim 保持固定宽度 | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| appendable dataset zero chunk 归一化为合法 chunk `[1,width]` | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| fixed dataset 当前 shape 保持 fixed spec | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| fixed dataset max dims 保持 fixed spec | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |
| fixed dataset chunk dims 保持 fixed spec | `Test_HighFive_Backend_Factory_And_Dataset_Reopen_Semantics` | 已覆盖 |

剩余风险：本轮新增真实 HighFive layout 读回测试代码，但尚未执行编译或 CTest；HDF5 C API helper 与当前 HighFive/HDF5 版本的实际通过情况仍需后续显式授权后验证。

## 95. Empty H5 output path contract 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `Command_Has_Non_Empty_Value` 对缺失 H5 path 返回 false | `Test_Contract_Helper_Functions` | 已覆盖 |
| `Command_Has_Non_Empty_Value` 对非空 trajectory path 返回 true | `Test_Contract_Helper_Functions` | 已覆盖 |
| `output_h5_trajectory_path = ""` 不启用 trajectory bundle | `Test_Empty_H5_Output_Paths_Do_Not_Enable_Bundles` | 已覆盖 |
| `output_h5_restart_path = ""` 不启用 restart bundle | `Test_Empty_H5_Output_Paths_Do_Not_Enable_Bundles` | 已覆盖 |
| `output_h5_observable_path = ""` 不启用 observable bundle | `Test_Empty_H5_Output_Paths_Do_Not_Enable_Bundles` | 已覆盖 |
| 三个空 H5 path 不设置 `any_h5_output_enabled` | `Test_Empty_H5_Output_Paths_Do_Not_Enable_Bundles` | 已覆盖 |
| 空 H5 path 不关闭默认 legacy sidecar 输出 | `Test_Empty_H5_Output_Paths_Do_Not_Enable_Bundles` | 已覆盖 |
| selector 参数可解析但不单独启用 H5 bundle | `Test_Empty_H5_Output_Paths_Do_Not_Enable_Bundles` | 已覆盖 |

剩余风险：本轮新增 contract helper 与 output plan 空 path 测试，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 96. SPONGE parameter path constant coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `path::sponge_mdout` 映射 `/parameters/sponge/mdout` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::sponge_log` 映射 `/parameters/sponge/log` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::sponge_files` 映射 `/parameters/sponge/files` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::sponge_provenance` 映射 `/parameters/sponge/provenance` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::mdout_columns_original_name` 映射 mdout original column leaf | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::mdout_columns_hdf5_name` 映射 mdout HDF5 column leaf | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| trajectory/observable base layout 使用参数区 root constants | `Test_H5_Output_Writer_Base_Layout_Paths` | 已覆盖 |
| observable-only mdout column metadata 使用 leaf constants | `Test_Observable_Only_Writer` | 已覆盖 |
| VDS wrapper mdout column metadata 使用 leaf constants | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| real HighFive trajectory/observable mdout column metadata 使用 leaf constants 读回 | `test_highfive_backend_io` trajectory/observable cases | 已覆盖 |
| provenance launch id 路径以 `path::sponge_provenance` 为 root | `Test_Observable_Only_Writer`, `test_highfive_backend_io` observable case | 已覆盖 |
| legacy sidecar root 以 `path::sponge_files` 为 root | `Test_Restart_Module_State_And_Legacy_Provenance`, writer legacy provenance cases | 已覆盖 |

剩余风险：本轮新增参数区 path constants 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 97. VDS output metadata path constant coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `path::output_trajectory_chunk_size` 映射 `/parameters/sponge/output/trajectory_chunk_size` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::output_vds_status` 映射 `/parameters/sponge/output/vds_status` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::output_repair_policy` 映射 `/parameters/sponge/output/repair_policy` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::output_repair_status` 映射 `/parameters/sponge/output/repair_status` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::output_repaired_shard_count` 映射 `/parameters/sponge/output/repaired_shard_count` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| VDS wrapper 写入 chunk size 使用 `path::output_trajectory_chunk_size` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets`, `test_highfive_backend_io` VDS cases | 已覆盖 |
| VDS wrapper 写入 VDS status 使用 `path::output_vds_status` | `Test_Vds_Finalize_Without_Frames`, `test_highfive_backend_io` VDS cases | 已覆盖 |
| VDS wrapper 写入 repair policy/status 使用 public constants | `Test_Vds_Wrapper_And_Module_Virtual_Datasets`, repair tests, `test_highfive_backend_io` VDS cases | 已覆盖 |
| VDS wrapper 写入 repaired shard count 使用 public constant | `Test_Vds_Wrapper_And_Module_Virtual_Datasets`, repair tests, `test_highfive_backend_io` VDS cases | 已覆盖 |

剩余风险：本轮新增 VDS output metadata path constants 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 98. Schema leaf and particle group path coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `path::sponge_schema_name` 映射 `/parameters/sponge/schema/name` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::sponge_schema_version` 映射 `/parameters/sponge/schema/version` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::particles_all_position` 映射 `/particles/all/position` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::particles_all_velocity` 映射 `/particles/all/velocity` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::particles_all_force` 映射 `/particles/all/force` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::particles_all_box` 映射 `/particles/all/box` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `path::particles_all_box_edges` 映射 `/particles/all/box/edges` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| trajectory base layout 使用 particle component group constants | `Test_H5_Output_Writer_Base_Layout_Paths` | 已覆盖 |
| restart base layout 使用 particle component group constants | `Test_Restart_Writer_Base_Layout_Paths` | 已覆盖 |
| real HighFive tests 通过 schema leaf constants 读回 schema metadata | `test_highfive_backend_io` | 已覆盖 |

剩余风险：本轮新增 schema/group path constants 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 99. VDS particle virtual dataset path constant coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| VDS wrapper particle step/time virtual datasets 使用 `path::particles_all_step/time` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper position virtual dataset 使用 `path::position_value` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper velocity virtual dataset 使用 `path::velocity_value` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper force disabled 分支不创建 `path::force_value` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets`, `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| VDS wrapper box edges virtual dataset 使用 `path::box_edges_value` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper particle hard links 使用 public step/time target constants | `Test_Vds_Wrapper_And_Module_Virtual_Datasets`, `Test_Vds_Optional_Particle_Fields_Disabled` | 已覆盖 |
| VDS writer 生产代码中不再残留 `/particles/all/*` particle virtual dataset literal | 窄范围文本检查 | 已覆盖 |

剩余风险：本轮新增 VDS particle virtual dataset path constant 使用更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 100. Writer open precondition coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| single-file trajectory writer 拒绝未启用 trajectory plan | `Test_Writer_Open_Preconditions_Reject_Unbound_Bundles` | 已覆盖 |
| single-file trajectory writer 拒绝 VDS trajectory plan | `Test_Writer_Open_Preconditions_Reject_Unbound_Bundles` | 已覆盖 |
| observable-only writer 拒绝未启用 observable plan | `Test_Writer_Open_Preconditions_Reject_Unbound_Bundles` | 已覆盖 |
| restart writer 拒绝未启用 restart plan | `Test_Writer_Open_Preconditions_Reject_Unbound_Bundles` | 已覆盖 |
| writer open failure 暴露具体 `Last_Error()` 文本 | `Test_Writer_Open_Preconditions_Reject_Unbound_Bundles` | 已覆盖 |

剩余风险：本轮新增 writer open precondition 测试，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 101. Bare H5MD writer no-backend coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| 未绑定 backend 时 `Is_Attached()` 为 false | `Test_H5MD_Writer_No_Backend_Is_Safe` | 已覆盖 |
| 未绑定 backend 时 open/flush/close/finalize 均失败 | `Test_H5MD_Writer_No_Backend_Is_Safe` | 已覆盖 |
| 未绑定 backend 时 group/dataset/VDS/hard link 写入均失败 | `Test_H5MD_Writer_No_Backend_Is_Safe` | 已覆盖 |
| 未绑定 backend 时 scalar/string/string-array/status/failure/completion 写入均失败 | `Test_H5MD_Writer_No_Backend_Is_Safe` | 已覆盖 |
| 未绑定 backend 时状态保持 `FileStatus::closed` | `Test_H5MD_Writer_No_Backend_Is_Safe` | 已覆盖 |
| 未绑定 backend 时 `Last_Error()` 返回明确错误文本 | `Test_H5MD_Writer_No_Backend_Is_Safe` | 已覆盖 |

剩余风险：本轮新增 bare writer no-backend 测试，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 102. Module dynamic leaf path builder coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| metadynamics scalar root/value/step/time helper literal | `Test_Module_Path_Constants` | 已覆盖 |
| metadynamics diagnostic root/component helper literal | `Test_Module_Path_Constants` | 已覆盖 |
| QC observable root/value/step/time helper literal | `Test_Module_Path_Constants` | 已覆盖 |
| QC SCF output helper literal | `Test_Module_Path_Constants` | 已覆盖 |
| ReaxFF term root/value/step/time helper literal | `Test_Module_Path_Constants` | 已覆盖 |
| metadynamics mapping behavior使用 dynamic helper 断言 dataset/hard-link/string path | `Test_Metadynamics_And_Diagnostics` | 已覆盖 |
| QC mapping behavior使用 dynamic helper 断言 dataset/hard-link/scf path | `Test_Qc_And_Reaxff_Mappings`, `Test_Qc_Optional_Spin_Square_Path` | 已覆盖 |
| ReaxFF mapping behavior使用 dynamic helper 断言 term dataset/hard-link path | `Test_Qc_And_Reaxff_Mappings` | 已覆盖 |

剩余风险：本轮新增 module dynamic leaf path builders 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 103. VDS module dynamic leaf helper alignment 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| VDS metadynamics wrapper scalar roots 使用 `Metadynamics_Scalar_Root` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS metadynamics wrapper value datasets 使用 `Metadynamics_Scalar_Value_Path` 断言 | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS QC wrapper scalar roots 使用 `Qc_Observable_Root` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS QC wrapper value datasets 使用 `Qc_Observable_Value_Path` 断言 | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS ReaxFF wrapper scalar roots 使用 `Reaxff_Term_Root` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS ReaxFF wrapper value datasets 使用 `Reaxff_Term_Value_Path` 断言 | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |

剩余风险：本轮新增 VDS module dynamic helper 对齐更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 104. HighFive module dynamic helper readback coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| metadynamics scalar value/step/time 真实 HDF5 readback 使用 dynamic helpers | `test_highfive_backend_io` | 已覆盖 |
| metadynamics diagnostic root/component 真实 HDF5 readback 使用 dynamic helpers | `test_highfive_backend_io` | 已覆盖 |
| QC energy 与 spin_square value/step/time 真实 HDF5 readback 使用 dynamic helpers | `test_highfive_backend_io` | 已覆盖 |
| QC SCF text 真实 HDF5 readback 使用 `Qc_Scf_Output_Path` | `test_highfive_backend_io` | 已覆盖 |
| ReaxFF bond/angle/over term 真实 HDF5 readback 使用 dynamic helpers | `test_highfive_backend_io` | 已覆盖 |
| HighFive real-backend 测试不再残留 metadynamics/QC/ReaxFF dynamic root 旧式拼接 | 窄范围文本检查 | 已覆盖 |

剩余风险：本轮新增 HighFive real-backend helper readback 对齐，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 105. Restart dynamic state path builder coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `Restart_Sits_State_Root` literal 映射 `/parameters/restart/bias/sits/<module>` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `Restart_Sits_State_Path` literal 映射 `/parameters/restart/bias/sits/<module>/<state>` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `Restart_Metad_State_Root` literal 映射 `/parameters/restart/bias/meta/<name>` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `Restart_Metad_State_Path` literal 映射 `/parameters/restart/bias/meta/<name>/<component>` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| restart writer 创建 SITS dynamic state group/dataset 时使用 helper 契约 | `Test_Restart_Module_State_And_Legacy_Provenance` | 已覆盖 |
| restart writer 创建 metad dynamic text component 时使用 helper 契约 | `Test_Restart_Module_State_And_Legacy_Provenance` | 已覆盖 |
| 真实 `.spgr.h5` readback 使用 restart SITS/metad dynamic helper 访问 dataset/string | `test_highfive_backend_io` | 已覆盖 |

剩余风险：本轮新增 restart dynamic path helper 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 106. Ordinary observable dynamic path builder coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `Observable_Root` literal 映射 `/observables/all/<name>` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `Observable_Value_Path` literal 映射 `/observables/all/<name>/value` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `Observable_Step_Path` literal 映射 `/observables/all/<name>/step` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| `Observable_Time_Path` literal 映射 `/observables/all/<name>/time` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| single-file trajectory ordinary observable dataset/hard-link/append path 使用 helper 契约 | `Test_Trajectory_Writer_Paths_And_Completion` | 已覆盖 |
| observable-only ordinary observable dataset/hard-link/append path 使用 helper 契约 | `Test_Observable_Only_Writer` | 已覆盖 |
| VDS wrapper ordinary observable virtual dataset 与 step/time link 使用 helper 契约 | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| 真实 HDF5 readback 使用 ordinary observable helper 访问 value/step/time | `test_highfive_backend_io` | 已覆盖 |
| 生产 writer 中普通 observable schema 路径不再残留独立 `path::observables_all + name` 拼接 | 窄范围文本检查 | 已覆盖 |

剩余风险：本轮新增 ordinary observable dynamic path helper 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 107. SPONGE provenance dynamic path builder coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `Sponge_Provenance_Path` literal 映射 `/parameters/sponge/provenance/<name>` | `Test_Public_H5MD_Path_Constants` | 已覆盖 |
| observable-only writer 写 launch provenance 时使用 `Sponge_Provenance_Path` | `Test_Observable_Only_Writer` | 已覆盖 |
| 真实 observable-only HDF5 readback 使用 `Sponge_Provenance_Path("launch_id")` | `test_highfive_backend_io` | 已覆盖 |
| mock writer 与 VDS mock writer 不再残留 metadynamics/QC/ReaxFF module dynamic root 旧式拼接 | 窄范围文本检查 | 已覆盖 |
| SPONGE provenance 旧式 `path::sponge_provenance + "/launch_id"` 拼接不再残留于测试断言 | 窄范围文本检查 | 已覆盖 |

剩余风险：本轮新增 SPONGE provenance dynamic path helper 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 108. NHC/SITS dynamic root path builder coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `Nose_Hoover_Chain_Coordinate_Root` literal 映射 NHC coordinate group | `Test_Module_Path_Constants` | 已覆盖 |
| `Nose_Hoover_Chain_Velocity_Root` literal 映射 NHC velocity group | `Test_Module_Path_Constants` | 已覆盖 |
| `Sits_Module_Root` literal 映射 `/observables/all/sits/<module>` | `Test_Module_Path_Constants` | 已覆盖 |
| module writer 创建 NHC coordinate/velocity groups 时使用 helper 契约 | `Test_Nhc_And_Sits_Mappings` | 已覆盖 |
| module writer 创建 SITS module group 时使用 helper 契约 | `Test_Nhc_And_Sits_Mappings` | 已覆盖 |
| VDS wrapper 创建 NHC coordinate/velocity groups 时使用 helper 契约 | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| VDS wrapper 创建 SITS module group 时使用 helper 契约 | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| 生产代码中 NHC coordinate/velocity 与 SITS module root 旧式拼接只剩 helper 实现本身 | 窄范围文本检查 | 已覆盖 |

剩余风险：本轮新增 NHC/SITS dynamic root path helper 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 109. Generic module scalar observable leaf path builder coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| `Scalar_Observable_Value_Path` literal 映射 `<root>/value` | `Test_Module_Path_Constants` | 已覆盖 |
| `Scalar_Observable_Step_Path` literal 映射 `<root>/step` | `Test_Module_Path_Constants` | 已覆盖 |
| `Scalar_Observable_Time_Path` literal 映射 `<root>/time` | `Test_Module_Path_Constants` | 已覆盖 |
| SITS `nk` path helpers 复用 generic scalar observable leaf helpers | `Test_Module_Path_Constants`, `Test_Nhc_And_Sits_Mappings` | 已覆盖 |
| metadynamics/QC/ReaxFF path helpers 复用 generic scalar observable leaf helpers | `Test_Module_Path_Constants`, module mapping tests | 已覆盖 |
| module writer 创建 scalar observable value dataset 和 step/time hard links 时使用 generic leaf helpers | `Test_Metadynamics_And_Diagnostics`, `Test_Qc_And_Reaxff_Mappings` | 已覆盖 |
| VDS wrapper 创建 module scalar virtual dataset 和 step/time hard links 时使用 generic leaf helpers | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| 生产代码中 module scalar `group/root + "/value|step|time"` 旧式拼接只剩 helper 实现本身 | 窄范围文本检查 | 已覆盖 |

剩余风险：本轮新增 generic module scalar observable leaf path helpers 与测试断言更新，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 110. VDS source dataset path mapping coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| particle VDS source `dataset_path` 等于 `path::position_value` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| ordinary observable VDS source `dataset_path` 等于 `Observable_Value_Path("temperature")` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| NHC VDS source `dataset_path` 等于 `module_path::nhc_coordinate_value` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| SITS VDS source `dataset_path` 等于 `Sits_Nk_Value_Path("sits_a")` | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| metadynamics/QC/ReaxFF VDS source `dataset_path` 等于对应 dynamic helper path | `Test_Vds_Wrapper_And_Module_Virtual_Datasets` | 已覆盖 |
| HighFive module diagnostic readback 不再残留局部 `metad_root + "/component"` 拼接 | 窄范围文本检查 | 已覆盖 |

剩余风险：本轮新增 VDS source dataset path mapping 断言，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 111. HighFive nested group idempotence coverage 审计补充

| 覆盖项 | 测试位置 | 审计结论 |
|---|---|---|
| HighFive backend 可递归创建多层 group | `Test_HighFive_Backend_Nested_Group_Idempotence` | 已覆盖 |
| HighFive backend 重复创建同一 group 必须幂等 | `Test_HighFive_Backend_Nested_Group_Idempotence` | 已覆盖 |
| 深层 group 下的 string dataset 可真实写入并读回 | `Test_HighFive_Backend_Nested_Group_Idempotence` | 已覆盖 |
| 新增 HighFive nested group 测试已加入 `Run_Test` 执行列表 | `test_highfive_backend_io.cpp` | 已覆盖 |

剩余风险：本轮新增真实 HighFive nested group 测试，但尚未执行编译或 CTest；实际通过情况仍需后续显式授权后验证。

## 112. H5 bundle test registration coverage audit 补充

| 覆盖项 | 证据位置 | 审计结论 |
|---|---|---|
| 顶层 CMake 可通过 `SPONGE_BUILD_TESTS=ON` 进入测试目录 | `CMakeLists.txt` | 已覆盖 |
| `tests/` 目录会进入 `tests/h5_bundle` | `tests/CMakeLists.txt` | 已覆盖 |
| 当前 7 个 h5 bundle `test_*.cpp` 均有 CMake executable 和 CTest 注册 | `tests/h5_bundle/CMakeLists.txt` | 已覆盖 |
| TEST_TARGETS executable/file/label 表与当前测试文件集合一致 | `tests/h5_bundle/TEST_TARGETS.md` | 已覆盖 |
| aggregate build targets 覆盖 contract、mock、backend-io 与总测试目标 | `tests/h5_bundle/CMakeLists.txt` | 已覆盖 |
| HighFive/HDF5 backend-io 测试依赖通过 CMake interface target 接入 | `tests/h5_bundle/CMakeLists.txt` | 已覆盖 |

剩余风险：本轮只做 CMake/manifest 源码审计，尚未执行 CMake configure、构建、单元测试或 CTest；实际通过情况仍需后续显式授权后验证。

## 113. Test schema literal cleanup coverage 审计补充

| 范围 | 结论 | 处理 |
|---|---|---|
| `test_h5md_writers_with_mock_backend.cpp` 普通 writer 行为断言 | 发现 observable-only writer 中 `frame_count` 仍硬编码公开 schema path | 已改为 `path::output_frame_count` |
| `test_h5md_writers_with_mock_backend.cpp` public path constant/helper literal 测试 | 字面量是测试对象本身，用于锁定规范路径 | 保留 |
| `test_highfive_backend_io.cpp` backend primitive/group/link 测试 | 部分 `/observables`、`/particles` 字面量是任意 HDF5 backend 能力测试路径，不作为 bundle schema owner | 保留 |
| H5 bundle writer/proxy/VDS/restart/module helper 测试 | 公开 schema 路径已经主要通过 `path::*`、`Observable_*_Path`、`Restart_*_State_Path`、module mapping helper 访问 | 维持 |

审计结论：本轮没有发现还需要新增测试目标的公开路径缺口；剩余字面量属于路径常量锁定测试或低层 backend 通用行为测试。后续如果新增 H5 bundle schema path，应同时新增 public helper/constant，并在普通 writer/VDS/restart 测试中优先使用 helper，而不是复制字符串。

未运行 CMake、构建或 CTest。

## 114. Restart writer failure contract coverage 审计补充

| 契约点 | 当前测试证据 | 本轮处理 |
|---|---|---|
| `restart.spgr.h5` 只允许一个 structural state | `Test_Restart_Writer_Is_Single_State` 已二次调用 `Write_Structural_State` 并期望失败 | 保持 |
| 二次写入失败必须进入 H5 bundle failed 状态 | 同一测试已断言 `FileStatus::failed` | 保持 |
| 二次写入失败必须写入标准 output status path | 之前未显式断言 `/parameters/sponge/output/status` | 新增断言 `failed` |
| 二次写入失败必须写入标准 output error path | 之前未显式断言 `/parameters/sponge/output/error` | 新增精确错误文本断言 |

审计结论：restart writer 的 wrapper-level failure path 现在与 trajectory/observable writer 一样，覆盖了状态枚举、标准状态字符串和标准错误路径。

未运行 CMake、构建或 CTest。

## 115. VDS wrapper shard-finalize failure coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| shard finalize 失败时 VDS strict finalize 返回失败 | `Test_Strict_Finalize_Fails_On_Shard_Finalize_Error` 已断言 `Finalize()` 返回 `false` | 保持 |
| shard finalize 失败时 wrapper 进入 failed 状态 | 旧实现未在 `Complete_Current_Shard` 失败分支写 wrapper failed 状态 | 在 `Finalize_Internal` 中调用 `wrapper_writer_->Mark_Failed(last_error_)` |
| wrapper 标准 output status path 记录失败 | 旧测试未覆盖 | 新增 `path::output_status == failed` 断言 |
| wrapper 标准 output error path 记录失败原因 | 旧测试未覆盖 | 新增 `path::output_error == mock finalize failure` 断言 |

审计结论：VDS/chunked H5MD 输出在 shard finalize failure 下的 wrapper-level failure metadata 现在有实现和 mock 单元测试共同覆盖。

未运行 CMake、构建或 CTest。

## 116. VDS materialization failure coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| wrapper virtual dataset 创建失败可被 mock 注入 | 旧 `MockBackend` 没有 VDS 创建失败开关 | 新增 `fail_next_virtual_dataset` |
| VDS materialization 失败会返回 `false` | 由生产路径逻辑间接支持，但缺少专门测试 | 新增 `Test_Vds_Materialize_Failure_Marks_Wrapper_Failed` |
| 失败原因进入 `writer.Last_Error()` | 旧测试未覆盖 | 新增精确错误文本断言 |
| wrapper 标准 output status path 记录失败 | 旧测试未覆盖 materialization failure 分支 | 新增 `path::output_status == failed` 断言 |
| wrapper 标准 output error path 记录失败原因 | 旧测试未覆盖 materialization failure 分支 | 新增 `path::output_error == mock virtual dataset failure: /particles/all/step` 断言 |

审计结论：VDS/chunked H5MD 输出在 shard finalize failure 与 wrapper VDS materialization failure 两类 finalize 阶段失败下，均已有 wrapper-level status/error metadata 覆盖。

未运行 CMake、构建或 CTest。

## 117. VDS finalize tail failure metadata coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| manifest 写入失败会传播 `Last_Error()` | `Write_Manifest_To_Wrapper` 旧分支未统一设置 `last_error_` | 为 group、string-array、dataset create、append 失败补充 `last_error_` 传播 |
| manifest 写入失败会写 wrapper failed metadata | 旧 `Finalize_Internal` 对 `Write_Manifest_To_Wrapper` 失败直接返回 | 新增 `wrapper_writer_->Mark_Failed(last_error_)` |
| repair metadata 写入失败会写 wrapper failed metadata | 旧 `Finalize_Internal` 对 `Write_Repair_Metadata` 失败直接返回 | 新增 `wrapper_writer_->Mark_Failed(last_error_)` |
| `output_vds_status` 写入失败会写 wrapper failed metadata | 旧分支直接返回 | 新增 `last_error_` 传播和 `Mark_Failed` |
| wrapper 最终 finalize 失败会写 wrapper failed metadata | 旧实现直接返回 `wrapper_writer_->Finalize()` | 改为显式检查失败并 `Mark_Failed` |
| manifest string-array 写入失败测试注入 | 旧 mock backend 无注入点 | 新增 `fail_next_string_array` |
| manifest path array 写入失败单元测试 | 旧测试未覆盖 | 新增 `Test_Vds_Manifest_Write_Failure_Marks_Wrapper_Failed` |

审计结论：VDS/chunked H5MD wrapper 的 finalize tail 失败路径现在具备实现层 failed metadata 落盘和 mock 单元测试覆盖。

未运行 CMake、构建或 CTest。

## 118. VDS string metadata failure coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| mock backend 可按指定 string path 注入失败 | 旧 mock 只能失败 append、virtual dataset、string-array | 新增 `fail_string_path` |
| `output_repair_policy` 写入失败可见 | 旧测试未覆盖 | 新增 `Test_Vds_Repair_Metadata_Failure_Marks_Wrapper_Failed` |
| `output_vds_status` 写入失败可见 | 旧测试未覆盖 | 新增 `Test_Vds_Status_Write_Failure_Marks_Wrapper_Failed` |
| VDS string metadata 写入失败进入 wrapper failed 状态 | 旧覆盖仅到 manifest string-array | 新增 `FileStatus::failed` 与 `path::output_status == failed` 断言 |
| VDS string metadata 写入失败原因进入标准 error path | 旧覆盖不足 | 新增 `path::output_error` 精确错误文本断言 |

审计结论：VDS/chunked H5MD wrapper finalize tail 的 string metadata 失败路径现在覆盖了 repair metadata 和 final VDS status 两类新路径。

未运行 CMake、构建或 CTest。

## 119. VDS wrapper final finalize failure coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| shard finalize failure 写 wrapper failed metadata | 已由 `Test_Strict_Finalize_Fails_On_Shard_Finalize_Error` 覆盖 | 保持 |
| wrapper final finalize failure 独立于 shard failure | 旧测试未覆盖 `{wrapper fails, shard succeeds}` 顺序 | 新增 `Test_Vds_Wrapper_Finalize_Failure_Marks_Wrapper_Failed` |
| wrapper final finalize failure 传播到 `Last_Error()` | 旧测试未覆盖 | 新增精确错误文本断言 |
| wrapper final finalize failure 写标准 output status/error | 旧测试未覆盖 | 新增 `path::output_status` 与 `path::output_error` 断言 |

审计结论：VDS/chunked H5MD wrapper 在最后 `Finalize()` 阶段失败时，现在有独立 mock 单元测试覆盖其标准 failed metadata 落盘行为。

未运行 CMake、构建或 CTest。

## 120. Strict trajectory chunk-size parsing coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| `output_h5_trajectory_chunk_size` 必须为完整整数 | 旧实现使用 `atoi`，可能接受 `12frames` | 改为 `strtol` 完整解析 |
| 非数字 chunk size 无效 | 已有 `not_an_integer` 测试 | 保持 |
| 部分数字加尾随文本无效 | 旧测试未覆盖 | 新增 `12frames` invalid 测试 |
| 错误信息指向 chunk size | 已有错误信息约束 | 新增测试复用 `chunk_size` substring 断言 |

审计结论：trajectory VDS chunk-size 新输入路径现在覆盖了 0、负数、非数字和部分数字文本四类无效输入。

未运行 CMake、构建或 CTest。

## 121. Strict trajectory VDS bool parsing coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| `Parse_Bool` 转换 true/false 文本 | 已有 `Test_Bool_Parsing_Text_Variants` | 保持 |
| resolver 拒绝非法 bool 文本 | 旧实现会把 `maybe` 静默解析为 false | 新增 `Is_Bool_Text` 和 resolver invalid plan |
| `output_h5_trajectory_vds=maybe` 无效 | 旧测试未覆盖 | 新增 invalid plan 测试 |
| 错误信息指向 VDS key | 旧测试未覆盖 | 新增 `output_h5_trajectory_vds` substring 断言 |

审计结论：trajectory VDS bool 新输入路径现在覆盖了合法 true/false 文本、非法文本识别和 resolver 层错误传播。

未运行 CMake、构建或 CTest。

## 122. Trajectory repair-policy error contract coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| 非法 repair policy 返回 invalid | 已有 `truncate` invalid 测试 | 保持 |
| 非法 repair policy 错误信息可定位 key | 旧测试未覆盖 | 新增 `output_h5_trajectory_repair_policy` substring 断言 |
| `complete_prefix` 需要 VDS trajectory | 已有 invalid 测试 | 保持 |
| `complete_prefix` 缺少 VDS true 时错误信息可定位条件 | 旧测试未覆盖 | 新增 `output_h5_trajectory_vds=true` substring 断言 |

审计结论：trajectory repair-policy 新输入路径现在覆盖合法值、alias、大小写、非法值、缺失 VDS 条件，以及相应错误信息契约。

未运行 CMake、构建或 CTest。

## 123. Directory-aware VDS shard-root derivation coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| 推荐 suffix trajectory path 派生 `.spg.shards` | 已有 `x.spg.h5md -> x.spg.shards` 测试 | 保持 |
| 非推荐 suffix fallback 派生 `.shards` | 已有 `x.h5 -> x.h5.shards` 测试 | 保持 |
| 带目录 trajectory path 派生时保留目录 | 旧测试未覆盖 | 新增 `runs/prod.spg.h5md -> runs/prod.spg.shards` |

审计结论：VDS shard root 派生现在覆盖了 basename 与 directory-preserving 两类路径形态。

未运行 CMake、构建或 CTest。

## 124. VDS shard filename sequence coverage 审计补充

| 契约点 | 旧证据 | 本轮处理 |
|---|---|---|
| 第一个 shard 命名为 `segment_000000.spg.h5md` | 已有 opened/source path 断言 | 保持 |
| 第二个 shard 命名为 `segment_000001.spg.h5md` | 已有 opened path 断言 | 保持 |
| 连续编号和六位 zero-padding 覆盖到第三个 shard | 旧测试未覆盖 | 新增 `Test_Vds_Shard_Filename_Sequence_Uses_Six_Digit_Padding` |
| manifest path 与实际 opened shard path 一致 | 旧测试局部覆盖前两个 shard | 新增第三个 shard manifest/opened path 双断言 |

审计结论：VDS/chunked H5MD shard filename 规则现在覆盖连续编号、六位 zero-padding 和 manifest/opened path 一致性。

未运行 CMake、构建或 CTest。
