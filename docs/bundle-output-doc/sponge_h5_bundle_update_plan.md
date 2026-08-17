# SPONGE H5 Bundle Update Development Plan

## 1. Objective

This plan turns the H5 input/output contracts into an implementation tracklist
for SPONGE.

The target canonical H5 outputs are:

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

The parser-visible flattened keys are:

```text
output_h5_trajectory_path
output_h5_trajectory_vds
output_h5_trajectory_chunk_size
output_h5_trajectory_repair_policy
output_h5_restart_path
output_h5_observable_path
```

The target artifacts are:

| Artifact | Role |
|---|---|
| `*.spg.h5md` | Trajectory H5MD output containing `/particles`, `/observables`, and `/parameters`. |
| `*.spgr.h5` | Single-state restart H5 output. |
| `*.obs.spg.h5md` | Observable-only H5MD output containing `/observables` and `/parameters`, but no trajectory particle fields. |

The migration goal is to make these H5 outputs the canonical owners while
keeping legacy sidecars available only as explicit compatibility outputs.

## 2. Contract Inputs and Outputs

### Input keys

| Key | Type | Default | Meaning |
|---|---|---|---|
| `output_h5_trajectory_path` | path | unset | Canonical trajectory H5MD output path. |
| `output_h5_trajectory_vds` | bool | `false` | Use chunked H5MD shards plus HDF5 VDS wrapper. |
| `output_h5_trajectory_chunk_size` | int | `20` | VDS file-level shard size in trajectory frames. |
| `output_h5_trajectory_repair_policy` | string | `strict` | VDS finalize policy: `strict` or explicit `complete_prefix`. |
| `output_h5_restart_path` | path | unset | Canonical restart H5 output path. |
| `output_h5_observable_path` | path | unset | Optional observable-only H5MD output path. |

`output_h5_trajectory_chunk_size` is not an HDF5 dataset chunk shape. It does
not change striping, compression, `buffer_frame`, or
`write_trajectory_interval`.

### Legacy output behavior

Legacy output paths remain compatibility sidecars:

```text
mdout
mdinfo
crd
box
vel
frc
rst
qc_scf_output
```

If any canonical H5 output is enabled, legacy sidecars are disabled by default.
Legacy sidecars are written only when their legacy path keys are explicitly set.

`default_out_file_prefix` remains a legacy prefix mechanism and is not the
recommended interface for H5 bundle naming.

### Output ownership

| Current output | Canonical H5 owner |
|---|---|
| `crd` | `*.spg.h5md:/particles/all/position/value` |
| `box` | `*.spg.h5md:/particles/all/box/edges/value` |
| `vel` | `*.spg.h5md:/particles/all/velocity/value` |
| `frc` | `*.spg.h5md:/particles/all/force/value` |
| `mdout` columns | `*.spg.h5md:/observables/all/<column>/value` |
| `mdinfo` text | `*.spg.h5md:/parameters/sponge/log/mdinfo_text` |
| restart coordinate/velocity/box | `*.spgr.h5:/particles/all/...` |
| restart internal state | `*.spgr.h5:/parameters/restart/...` |
| observable-only scalars | `*.obs.spg.h5md:/observables/...` |

## 3. Tracklist

| Track | Goal | Deliverable | Verification |
|---|---|---|---|
| T0 Contract constants | 固化 key、默认值、文件后缀、legacy gating 规则 | 输出 key 常量、默认 `chunk_size=20`、后缀检查辅助 | parser dry-run 能识别五个 `output_h5_*` key |
| T1 HDF5/H5MD writer foundation | 提供统一 H5 写入基础设施 | H5 file lifecycle、group/dataset append、step/time hard-link helper | writer API 能覆盖后续 H5 输出路径 |
| T1b Concrete HDF5 backend | 使用 HighFive 实现真实 HDF5 backend | HighFiveBackend、HighFiveBackendFactory、HighFive 构建依赖 | 最小 `*.spg.h5md` 可真实落盘并可用 HDF5 工具读取 |
| T2 Output resolver | 将 mdin 输出配置解析为 runtime output plan | `TrajectoryH5OutputPlan`, `RestartH5OutputPlan`, `ObservableH5OutputPlan`, legacy sidecar plan | bundle 启用时 legacy 默认关闭；显式 legacy path 时 sidecar 开启 |
| T3 Trajectory H5 single-file | 实现 `output_h5_trajectory_path` 且 `vds=false` | `/h5md`, `/particles`, `/observables`, `/parameters/sponge` 单文件输出 | 与 legacy `crd/box/vel/frc/mdout` 对照帧数和关键数值 |
| T4 Observable-only H5 | 实现 `output_h5_observable_path` | `*.obs.spg.h5md`，只含 `/h5md`, `/observables`, `/parameters` | 文件中不存在 `/particles`，`mdout` 标量完整 |
| T5 Restart H5 | 实现 `output_h5_restart_path` | 单帧 `*.spgr.h5`，H5MD-compatible `/particles` + `/parameters/restart` | 用生成的 restart 启动下一段 structural continuation |
| T6 VDS trajectory output | 实现 `output_h5_trajectory_vds=true` | VDS wrapper + writer-derived shard directory，`chunk_size=20` frames | 21 个 trajectory frame 产生 2 个 shard，wrapper VDS 可读 |
| T7 Module-specific mappings | 接入非核心输出 | NHC、SITS、metad、QC、ReaxFF scalar/log 映射 | 各模块小 case 能在 H5 中找到契约路径 |
| T8 Legacy compatibility | 保持旧输出可显式开启 | `mdout/mdinfo/crd/box/vel/frc/rst/qc_scf_output` sidecar | bundle 开启且未显式 legacy path 时不生成 sidecar |
| T9 Failure/resume/finalize | 处理 incomplete frame/shard | status、frame_count、manifest、finalized marker | 中断后重启忽略或截断 incomplete trailing frame |

### T0 Contract constants

Architecture:

T0 introduces the stable contract surface used by every later track. It should
be a small, dependency-light layer that centralizes parser-visible key names,
default values, recommended suffixes, and compatibility policy. No HDF5 writer
logic should live here.

Implementation scheme:

| Component | Responsibility |
|---|---|
| output key constants | Define `output_h5_trajectory_path`, `output_h5_trajectory_vds`, `output_h5_trajectory_chunk_size`, `output_h5_restart_path`, `output_h5_observable_path`. |
| output defaults | Define `trajectory_vds=false` and `trajectory_chunk_size=20`. |
| suffix helpers | Validate or warn on non-recommended suffixes without blocking explicitly requested paths unless strict mode is enabled. |
| legacy policy helper | Determine whether legacy sidecars are default-enabled or default-disabled. |

TODO list:

| Status | Task |
|---|---|
| TODO | Add constants for all five H5 output keys. |
| TODO | Add constants for suffixes: `.spg.h5md`, `.spgr.h5`, `.obs.spg.h5md`. |
| TODO | Add default `output_h5_trajectory_chunk_size = 20`. |
| TODO | Add helper for `any_h5_output_enabled`. |
| TODO | Add helper for `legacy_sidecar_requested(key)`. |
| TODO | Add parser-visible documentation comments or equivalent command metadata. |
| TODO | Add dry-run parser case for grouped TOML flattening. |

Completion criteria:

The resolver can query one central contract layer and does not duplicate key
strings or defaults in writer code.

### T1 HDF5/H5MD writer foundation

Architecture:

T1 provides the reusable HDF5/H5MD writer substrate. It should be independent of
MD physics modules and accept already-resolved data buffers, steps, times,
labels, and metadata. Higher tracks should not call raw HDF5 APIs directly
except through this layer.

Implementation scheme:

| Component | Responsibility |
|---|---|
| H5 file handle | Open, close, flush, and finalize one HDF5 file. |
| H5MD layout initializer | Create `/h5md`, `/particles`, `/observables`, `/parameters`, and `/parameters/sponge`. |
| appendable dataset writer | Append scalar, vector, matrix, and particle-frame values. |
| step/time helper | Create shared step/time datasets and hard links for element-local paths. |
| string/table writer | Store mdout columns, provenance, legacy filename records, and logs. |
| status writer | Maintain `open`, `closing`, `finalized`, and `failed` markers. |

TODO list:

| Status | Task |
|---|---|
| TODO | Select the internal HDF5 wrapper style and isolate it behind one writer API. |
| TODO | Implement group creation with idempotent `ensure_group` behavior. |
| TODO | Implement extendable dataset creation for frame-major datasets. |
| TODO | Implement append for `float32[N,3]`, `float32[3,3]`, `float64`, `int64`, and string data. |
| TODO | Implement shared step/time dataset creation and hard-link helper. |
| TODO | Implement `/h5md` metadata writer with creator and version fields. |
| TODO | Implement `/parameters/sponge/schema` and `/parameters/sponge/output` metadata writer. |
| TODO | Add a tiny standalone writer smoke path that writes a minimal readable file. |

Completion criteria:

A minimal `*.spg.h5md` can be written without involving MD output call sites,
and HDF5 inspection tools can read its top-level groups and datasets.

### T1b Concrete HDF5 backend

Architecture:

T1b turns the T1 writer abstraction into a real HDF5-capable backend using
HighFive. Higher-level trajectory, observable, restart, and VDS writers must
still depend only on `WriterBackend`; HighFive types remain isolated to the
backend adapter.

Implementation scheme:

| Component | Responsibility |
|---|---|
| `HighFiveBackend` | Implement `WriterBackend` using HighFive and HDF5. |
| file lifecycle | Create parent directories, open output file, flush, finalize, and close. |
| group creation | Provide idempotent recursive `Ensure_Group`. |
| extendable datasets | Map `DatasetSpec` to HighFive dataspace, max dims, and chunking. |
| append | Extend the first dimension by one record and write frame-major payload. |
| hard links | Use HDF5 hard links for shared step/time paths. |
| strings | Store scalar string and string-array datasets under parameter paths. |
| factory | Provide `HighFiveBackendFactory` for VDS wrapper/shard construction. |

Dependency policy:

HDF5 and HighFive are supplied by pixi environments, not by vendored source or
ad-hoc system installations. SPONGE configure tasks pass the pixi
`$CONDA_PREFIX` as `CMAKE_PREFIX_PATH`, so `find_package(HighFive CONFIG
REQUIRED)` and `find_package(HDF5 REQUIRED COMPONENTS C)` resolve against the
active pixi environment.

TODO list:

| Status | Task |
|---|---|
| DONE | Add pixi-provided `highfive`/`hdf5` dependencies to environments that build SPONGE. |
| DONE | Link the SPONGE target against pixi-provided HighFive and HDF5 CMake targets. |
| DONE | Implement `HighFiveBackend` behind the existing `WriterBackend` API. |
| DONE | Implement frame-major append semantics for numeric datasets. |
| DONE | Implement HDF5 hard-link creation for shared step/time paths. |
| DONE | Implement string and string-array dataset writing. |
| TODO | Add a minimal smoke writer executable or test harness. |
| TODO | Run configure/compile and inspect a generated minimal H5MD file. |

Completion criteria:

The backend is available to higher-level H5 writers through `WriterBackend`, and
a future smoke path can generate a real HDF5 file without involving MD runtime
call sites. Full completion still requires compile and HDF5 inspection
validation.

### T2 Output resolver

Architecture:

T2 converts raw controller commands into an immutable runtime output plan. It is
the only layer that decides which canonical H5 artifacts and which legacy
sidecars are active. Writers consume the plan and do not inspect mdin commands
directly.

Implementation scheme:

| Plan object | Fields |
|---|---|
| `TrajectoryH5OutputPlan` | `enabled`, `path`, `vds`, `chunk_size`, `derived_shard_root`. |
| `RestartH5OutputPlan` | `enabled`, `path`, `write_interval`, `retention_count`. |
| `ObservableH5OutputPlan` | `enabled`, `path`. |
| `LegacyOutputPlan` | explicit path flags for `mdout`, `mdinfo`, `crd`, `box`, `vel`, `frc`, `rst`, `qc_scf_output`. |
| `ResolvedOutputPlan` | all plans plus common provenance and cadence values. |

TODO list:

| Status | Task |
|---|---|
| TODO | Add resolver entrypoint after mdin parsing and before output file creation. |
| TODO | Read grouped TOML flattened H5 keys. |
| TODO | Apply default `trajectory_vds=false`. |
| TODO | Apply default `trajectory_chunk_size=20`. |
| TODO | Hard error on `trajectory_chunk_size <= 0`. |
| TODO | Derive shard root internally from `output_h5_trajectory_path`. |
| TODO | Do not accept user-facing shard directory keys. |
| TODO | Compute `legacy_default_enabled = !any_h5_output_enabled`. |
| TODO | Mark each legacy output active only if legacy default is enabled or the path key is explicit. |

Completion criteria:

All output call sites can be driven by `ResolvedOutputPlan`, and legacy output
gating is deterministic and testable.

### T3 Trajectory H5 single-file

Architecture:

T3 connects the existing core trajectory and scalar output flow to the H5 writer
for the non-VDS case. It should preserve existing cadence semantics:
`write_trajectory_interval` controls particle trajectory frames, while
`write_mdout_interval` controls scalar observables.

Implementation scheme:

| Current data source | H5 target |
|---|---|
| coordinate trajectory | `/particles/all/position/value` |
| box trajectory | `/particles/all/box/edges/value` |
| velocity trajectory | `/particles/all/velocity/value` |
| force trajectory | `/particles/all/force/value` |
| `Step_Print` columns | `/observables/all/<column>/value` |
| mdout column labels | `/parameters/sponge/mdout/columns` |
| mdinfo content | `/parameters/sponge/log/mdinfo_text` |

TODO list:

| Status | Task |
|---|---|
| TODO | Add trajectory H5 writer construction when `output_h5_trajectory_path` is set and `vds=false`. |
| TODO | Write initial H5MD metadata and SPONGE schema metadata. |
| TODO | Append position and box frames at trajectory cadence. |
| TODO | Append velocity frames when velocity output is active or required by output policy. |
| TODO | Append force frames when force output is active. |
| TODO | Capture `Step_Print` registration into H5 observable metadata. |
| TODO | Append scalar observables at mdout cadence. |
| TODO | Store original mdout labels and sanitized HDF5 dataset names. |
| TODO | Store provenance hashes when available. |
| TODO | Finalize status at normal termination. |

Completion criteria:

For a short run, H5 frame counts and selected numeric values match legacy
`crd`, `box`, `vel`, `frc`, and `mdout` outputs when those sidecars are
explicitly enabled for comparison.

### T4 Observable-only H5

Architecture:

T4 reuses the observable and parameter writer paths from T3 but disables all
particle groups. It is intended for lightweight analysis and cloud sync, so it
must not accidentally create trajectory datasets.

Implementation scheme:

| Included | Excluded |
|---|---|
| `/h5md` | `/particles/all/position` |
| `/observables` | `/particles/all/velocity` |
| `/parameters` | `/particles/all/force` |
| `/parameters/sponge` | `/particles/all/box` |

TODO list:

| Status | Task |
|---|---|
| TODO | Add observable-only writer construction when `output_h5_observable_path` is set. |
| TODO | Reuse observable schema and mdout column metadata from trajectory H5. |
| TODO | Reuse `/parameters/sponge` log/provenance writers where applicable. |
| TODO | Prevent creation of `/particles` in observable-only files. |
| TODO | Support observable-only output without `output_h5_trajectory_path`. |
| TODO | Ensure observable-only cadence follows observable cadence, not trajectory cadence. |
| TODO | Finalize observable-only status independently from trajectory H5 status. |

Completion criteria:

A run with only `output_h5_observable_path` creates a valid lightweight H5MD
file with scalar observables and no `/particles` group.

### T5 Restart H5

Architecture:

T5 implements canonical restart output as a single launchable H5 state, not a
trajectory. It writes standard coordinate, velocity, and box state through an
H5MD-compatible one-frame particles layout, while non-standard continuation
state lives under `/parameters/restart`.

Implementation scheme:

| Restart domain | H5 target |
|---|---|
| step/time | `/particles/all/step`, `/particles/all/time` |
| coordinates | `/particles/all/position/value` |
| velocities | `/particles/all/velocity/value` |
| box | `/particles/all/box/edges/value` |
| NHC state | `/parameters/restart/thermostat/nose_hoover_chain` |
| SITS state | `/parameters/restart/bias/sits/<module>/...` |
| metad state | `/parameters/restart/bias/meta/<name>/...` |

TODO list:

| Status | Task |
|---|---|
| TODO | Add restart H5 writer construction when `output_h5_restart_path` is set. |
| TODO | Write exactly one launchable state per restart file. |
| TODO | Write coordinate and box state unconditionally when available. |
| TODO | Write velocity state when present or required by continuation mode. |
| TODO | Write restart metadata under `/run` and `/parameters/restart`. |
| TODO | Route NHC restart state into restart H5. |
| TODO | Route SITS `nk` restart state into restart H5. |
| TODO | Route metad path-dependent accumulated state into restart H5. |
| TODO | Do not copy topology-owned or protocol-owned definitions into restart as owners. |
| TODO | Implement replacement/retention behavior based on `max_restart_export_count`. |

Completion criteria:

A generated `*.spgr.h5` can be bound as restart input for the next run and
provides the expected structural continuation state.

### T6 VDS trajectory output

Architecture:

T6 implements segmented trajectory output. The user-facing file is
`output_h5_trajectory_path`; shard files are writer-derived internals. VDS
chunking uses trajectory frame count, not MD step count and not HDF5 internal
chunk shape.

Implementation scheme:

| Concept | Rule |
|---|---|
| shard size | `output_h5_trajectory_chunk_size` trajectory frames. |
| default shard size | `20` trajectory frames. |
| shard root | Derived internally from `output_h5_trajectory_path`. |
| wrapper | User-facing `*.spg.h5md` file. |
| manifest | Stored under `/parameters/sponge/output/shard_manifest`. |

TODO list:

| Status | Task |
|---|---|
| TODO | Derive internal shard directory from trajectory H5 path. |
| TODO | Open shard 0 on first trajectory frame. |
| TODO | Rotate shard when `trajectory_frame_index % chunk_size == 0` after the first shard. |
| TODO | Write shard-local H5MD-like datasets. |
| TODO | Mark shard `complete` only after all required frame data is flushed. |
| TODO | Maintain manifest with shard index, path, frame range, step/time range, and status. |
| TODO | Generate or update VDS wrapper datasets for particle and observable streams. |
| TODO | Support observable frame counts that differ from trajectory frame counts. |
| TODO | Finalize wrapper without copying all shard payload into one file. |

Completion criteria:

With `chunk_size=20`, a run with 21 trajectory frames produces two shards and a
readable VDS wrapper spanning both shards.

### T7 Module-specific mappings

Architecture:

T7 connects non-core modules to the canonical H5 paths without changing their
scientific semantics. Each module should keep producing the same runtime values,
but canonical storage moves to `/observables`, `/parameters/sponge`, or
`/parameters/restart` according to ownership.

Implementation scheme:

| Module | Canonical output |
|---|---|
| NHC | Observable chain coordinate/velocity trajectories plus restart chain state. |
| SITS | Observable `nk` trajectory plus restart `nk` state. |
| Metad | Observable `meta/rbias/rct`, diagnostic exports, hill/history/edge logs, restart accumulated state. |
| QC | Observable QC energy/spin plus SCF diagnostic log. |
| ReaxFF | Observable energy terms. |

TODO list:

| Status | Task |
|---|---|
| DONE | Route NHC coordinate/velocity trajectories to `/observables/all/thermostat/nose_hoover_chain`. |
| DONE | Route NHC restart output to restart H5. |
| DONE | Route SITS `nk_traj_file` data to `/observables/all/sits/<module>/nk/value`. |
| DONE | Route SITS `nk_rest_file` state to restart H5. |
| DONE | Route metad scalar columns to `/observables/all/metadynamics`. |
| DONE | Route metad `Meta_Potential.txt` to `/parameters/sponge/metadynamics/<name>/potential_export`. |
| PARTIAL | Route metad `Meta_directly.txt` to `/parameters/sponge/metadynamics/<name>/direct_export`. |
| DONE | Route metad `myhill.log` to `/parameters/sponge/metadynamics/<name>/hills`. |
| DONE | Route metad `history.log` to `/parameters/sponge/metadynamics/<name>/history`. |
| DONE | Route metad edge/sumhill data to `/parameters/sponge/metadynamics/<name>/edge`. |
| PARTIAL | Route metad path-dependent text snapshots to restart H5. |
| DONE | Route QC energy/spin scalar columns to `/observables/all/qc`. |
| DONE | Route QC SCF log to `/parameters/sponge/qc/scf_output`. |
| DONE | Route ReaxFF energy terms to `/observables/all/reaxff`. |
| DONE | Leave `eeq_charges.txt` handling deferred. |

Completion criteria:

Each enabled module can be run in a small case and its documented H5 path is
present with expected shape and frame count.

### T8 Legacy compatibility

Architecture:

T8 preserves legacy output behavior for existing workflows while allowing H5
output to become canonical. The rule is explicit and centralized: when H5 is
enabled, legacy defaults turn off; explicit legacy paths still write sidecars.

Implementation scheme:

| Case | Behavior |
|---|---|
| no H5 output keys | Current legacy behavior remains. |
| H5 output key set, no legacy path | Do not write that legacy sidecar. |
| H5 output key set, explicit legacy path | Write H5 and requested legacy sidecar. |
| `default_out_file_prefix` only | Applies to legacy outputs, not H5 bundle names. |

TODO list:

| Status | Task |
|---|---|
| DONE | Centralize explicit legacy path detection. |
| DONE | Prevent implicit `mdout` and `mdinfo` sidecars when H5 output is enabled unless explicitly requested. |
| DONE | Prevent implicit `crd`, `box`, `vel`, and `frc` sidecars when H5 trajectory output is enabled unless explicitly requested. |
| DONE | Prevent implicit `rst` sidecars when H5 restart output is enabled unless explicitly requested. |
| DONE | Keep `qc_scf_output` as an explicit legacy diagnostic path. |
| DONE | Record explicit legacy sidecar paths under `/parameters/sponge/files` or `/parameters/sponge/restart_exports` where relevant. |

Completion criteria:

Existing legacy-only workflows behave as before, and H5-enabled workflows do not
emit unexpected sidecars.

### T9 Failure/resume/finalize

Architecture:

T9 makes H5 output robust to interrupted runs. A frame or shard must be either
complete or ignored. Finalization should expose only complete data to downstream
readers.

Implementation scheme:

| State | Meaning |
|---|---|
| `open` | File or shard is being written. |
| `closing` | Writer is finalizing metadata and links. |
| `finalized` | File or wrapper is complete and reader-facing. |
| `failed` | Writer detected an unrecoverable output error. |

TODO list:

| Status | Task |
|---|---|
| DONE | Write file-level status before first frame append. |
| DONE | Advance `frame_count` only after all required datasets for that frame are written. |
| DONE | Store `last_complete_step` and `last_complete_time`. |
| PARTIAL | In VDS mode, maintain shard-level status and frame counts. |
| DONE | Mark failed status on writer-level frame write failure. |
| TODO | On resume, detect incomplete trailing frame or shard. |
| TODO | Truncate or ignore incomplete trailing data before appending new data. |
| TODO | Regenerate VDS wrapper from complete manifest entries. |
| DONE | Mark wrapper `finalized` only after complete manifest validation. |
| PARTIAL | Provide hard error for inconsistent shard manifest unless repair mode is explicitly requested. |

Completion criteria:

Killing a run during output does not expose partial frames as complete data, and
resume/finalization produces a consistent readable H5 artifact.

## 4. Step-by-step Implementation Plan

### Step 1: Add output contract constants and resolver

Implement a small output configuration layer that consumes flattened TOML keys:

```text
output_h5_trajectory_path
output_h5_trajectory_vds
output_h5_trajectory_chunk_size
output_h5_restart_path
output_h5_observable_path
```

The resolver should produce immutable runtime plans:

```text
TrajectoryH5OutputPlan
RestartH5OutputPlan
ObservableH5OutputPlan
LegacyOutputPlan
```

Required behavior:

| Condition | Behavior |
|---|---|
| `output_h5_trajectory_chunk_size` unset | Use `20`. |
| `output_h5_trajectory_chunk_size <= 0` | Hard error. |
| H5 output enabled and legacy path unset | Do not write that legacy sidecar. |
| H5 output enabled and legacy path explicitly set | Write both canonical H5 and legacy sidecar. |
| no H5 output enabled | Preserve current legacy behavior. |

### Step 2: Build HDF5/H5MD writer primitives

Add writer primitives before changing MD output call sites.

Required primitives:

| Primitive | Purpose |
|---|---|
| file open/close/finalize | Manage HDF5 file lifecycle and status. |
| ensure group | Create `/h5md`, `/particles`, `/observables`, `/parameters`, `/parameters/sponge`. |
| append dataset frame | Append scalar, vector, and particle-frame datasets. |
| shared step/time | Write one stream step/time dataset and link element-local step/time paths. |
| string/table writer | Store labels, logs, provenance, and legacy filename records. |
| status marker | Mark `open`, `closing`, `finalized`, or `failed`. |

### Step 3: Implement single-file trajectory H5

Route existing trajectory-producing data into `output_h5_trajectory_path` when
`output_h5_trajectory_vds = false`.

Minimum required datasets:

```text
/h5md
/particles/all/position/value
/particles/all/position/step
/particles/all/position/time
/particles/all/box/edges/value
/observables/all/<column>/value
/parameters/sponge/schema/name
/parameters/sponge/output/status
/parameters/sponge/mdout/columns
```

Optional datasets are written only when their legacy equivalents are active or
the runtime produces the data:

```text
/particles/all/velocity/value
/particles/all/force/value
```

### Step 4: Implement observable-only H5

Implement `output_h5_observable_path` by reusing observable and parameter
writers while disabling all particle output.

Required top-level groups:

```text
/h5md
/observables
/parameters
/parameters/sponge
```

Forbidden groups:

```text
/particles/all/position
/particles/all/velocity
/particles/all/force
/particles/all/box
```

Observable-only output may be produced with or without trajectory H5 output.

### Step 5: Implement restart H5

Implement `output_h5_restart_path` as a single launchable restart state.

Minimum required datasets:

```text
/h5md
/particles/all/step
/particles/all/time
/particles/all/position/value
/particles/all/box/edges/value
/parameters/restart
/run/current_step
/run/current_time
```

Conditional datasets:

| Dataset | Condition |
|---|---|
| `/particles/all/velocity/value` | Velocity exists or launch policy requires continuation velocities. |
| `/parameters/restart/thermostat/nose_hoover_chain` | Nose-Hoover chain state exists. |
| `/parameters/restart/bias/sits/<module>/nk` | SITS adaptive state exists. |
| `/parameters/restart/bias/meta/<name>/...` | Metad path-dependent restart state exists. |

The restart H5 must not copy topology owner data or protocol owner data.

### Step 6: Implement VDS trajectory output

When `output_h5_trajectory_vds = true`, write completed frame segments to
writer-derived shards and expose one user-facing VDS wrapper at
`output_h5_trajectory_path`.

Required behavior:

| Rule | Meaning |
|---|---|
| shard size | At most `output_h5_trajectory_chunk_size` trajectory frames per shard. |
| default shard size | `20` trajectory frames. |
| shard directory | Derived from trajectory H5 path; not user configurable. |
| wrapper | User-facing `*.spg.h5md` file with VDS datasets. |
| manifest | Store shard index, path, frame range, step/time range, status. |

Observable streams may have a different frame count from trajectory streams in a
shard because `write_mdout_interval` and `write_trajectory_interval` are
separate.

### Step 7: Add module-specific mappings

Route module outputs according to the contract.

| Module | H5 mapping |
|---|---|
| NHC | Trajectory diagnostics under `/observables/all/thermostat/nose_hoover_chain`; restart state under `/parameters/restart/thermostat/nose_hoover_chain`. |
| SITS | `nk_traj_file` under `/observables/all/sits/<module>/nk/value`; restart `nk` under `/parameters/restart/bias/sits/<module>/nk`. |
| Metad | `meta/rbias/rct` under `/observables/all/metadynamics`; `Meta_Potential.txt`, `Meta_directly.txt`, `myhill.log`, `history.log`, and edge/sumhill diagnostics under `/parameters/sponge/metadynamics/<name>`. |
| QC | `QC` and `QC_S_sq` under `/observables/all/qc`; `qc_scf_output` under `/parameters/sponge/qc/scf_output`. |
| ReaxFF | Energy terms under `/observables/all/reaxff/<term>/value`. |

`eeq_charges.txt` remains deferred and should not be added to the input contract
in this implementation round.

### Step 8: Add failure, resume, and finalize handling

Implement markers that prevent incomplete frames or shards from becoming visible
as complete output.

Required markers:

```text
/parameters/sponge/output/status
/parameters/sponge/output/frame_count
/parameters/sponge/output/last_complete_step
/parameters/sponge/output/last_complete_time
```

VDS shards additionally require:

```text
/parameters/sponge/shard/status
/parameters/sponge/shard/frame_start
/parameters/sponge/shard/frame_count
/parameters/sponge/shard/last_complete_step
/parameters/sponge/shard/last_complete_time
```

## 5. Verification Matrix

| Gate | Scenario | Acceptance |
|---|---|---|
| G1 parser/config | Grouped TOML and flattened keys | Five `output_h5_*` keys are consumed; invalid `chunk_size <= 0` errors. |
| G2 single-file trajectory | Short NVT/NVE run with `output_h5_trajectory_vds=false` | `*.spg.h5md` contains particles, observables, and parameters. |
| G3 observable-only | Run with only `output_h5_observable_path` | `*.obs.spg.h5md` has no `/particles` group and contains `mdout` observables. |
| G4 restart continuation | Run A writes `*.spgr.h5`; run B binds it | Run B starts from the expected coordinate/box state. |
| G5 VDS chunking | `chunk_size=20`, at least 21 trajectory frames | Frame 0-19 in shard 0; frame 20 starts shard 1; wrapper VDS reads both. |
| G6 legacy gating | H5 output enabled, no legacy paths | No `mdout`, `mdinfo`, `crd`, `box`, `vel`, `frc`, or `rst` sidecars are created by default. |
| G7 explicit legacy | H5 output enabled plus explicit `mdout` or `crd` | H5 output and requested sidecar are both written. |
| G8 metad | Metadynamics run | `meta/rbias/rct`, potential export, direct export, hills, history, and edge diagnostics land in documented H5 paths. |
| G9 SITS/NHC/QC/ReaxFF | Module-specific small cases | Module outputs land in documented H5 observable, parameter, or restart paths. |
| G10 failure/resume | Kill during frame/shard write | Incomplete trailing frame/shard is ignored or repaired before finalization. |

## 6. Migration and Compatibility Rules

The first implementation should preserve current behavior when no canonical H5
output key is set.

When any canonical H5 output is set:

```text
legacy sidecars default off
explicit legacy path writes sidecar
H5 output remains canonical source of truth
```

Compatibility outputs must never become the canonical owner of data that is
available in H5.

The restart file remains separate from trajectory output. Starting a new run
from a trajectory frame requires an explicit conversion from `*.spg.h5md` to
`*.spgr.h5`.

The VDS shard directory is implementation-internal. Users configure only
`output_h5_trajectory_path`, `output_h5_trajectory_vds`, and
`output_h5_trajectory_chunk_size`.

## 7. Open Non-blocking Follow-ups

- `eeq_charges.txt` H5 handling remains deferred.
- Full bitwise continuation restart is deferred unless integrator and RNG state
  are added.
- Portable VDS materialization into one standalone H5MD archive is an optional
  repack operation.
- Additional compression and HDF5 internal dataset chunk-shape tuning can be
  added after the canonical writer path is stable.

## 实施补充：HighFive/HDF5 VDS backend

已补入 `WriterBackend::Create_Virtual_Dataset` 和 `HighFiveBackend` 的 HDF5 VDS 实现，并将 VDS trajectory wrapper 的 finalize 流程推进为：完成当前 shard、校验完整 manifest、物化核心 particle virtual datasets、写入 manifest、finalize wrapper。

当前完成范围：

- `/particles/all/step`
- `/particles/all/time`
- `/particles/all/position/value`
- `/particles/all/box/edges/value`
- 显式启用时的 `/particles/all/velocity/value`
- 显式启用时的 `/particles/all/force/value`

后续仍需补齐：

- observable stream 的跨分片 VDS 物化。
- module stream 的跨分片 VDS 物化，例如 NHC、SITS、metadynamics scalar、QC、ReaxFF terms。
- VDS wrapper 的失败恢复、截断和 repair 策略。

## 实施补充：ordinary observable VDS

已在 VDS trajectory wrapper 中补入普通 `mdout` observable stream 的 VDS 物化：

- per-shard `observable_frame_count` 记录。
- `/observables/all/step` 和 `/observables/all/time` VDS。
- `/observables/all/<name>/value` VDS。
- `/observables/all/<name>/{step,time}` hard link。
- wrapper 级别的 mdout column metadata。

后续未完成项收敛为 module-specific stream VDS 与失败恢复策略。

## 实施补充：module-specific VDS

已在 VDS trajectory wrapper 中补入 module-specific stream 的 VDS 物化：

- NHC chain coordinate/velocity streams。
- SITS `nk` vector stream。
- Metadynamics `meta/rbias/rct` scalar streams。
- QC `energy` 与可选 `spin_square` scalar streams。
- ReaxFF energy term scalar streams。

这些 stream 使用各自 per-shard frame count，不与 trajectory frame count 或普通 mdout observable frame count 混用。

后续重点转为：编译/API 验证、HDF5 VDS 文件级检查、failed wrapper repair/truncate/resume 策略。

## 实施补充：VDS source path relocation

VDS wrapper 的 source path 语义已细化：manifest 记录原始 shard path，HDF5 VDS mapping 优先使用相对 wrapper 所在目录的 source path。该策略避免把 VDS wrapper 固化到单一机器路径，适合云端分片同步后整体下载/搬迁的场景。

## 实施补充：T9 explicit repair finalize

T9 现在具备 writer-local 的显式 repair finalize 能力：`Finalize()` 为 strict mode，`Finalize_With_Repair()` 为 complete-prefix repair mode。repair mode 会重算 completion metadata，记录 repair policy/status/count，并避免 wrapper 引用不完整 shard。

仍需后续接入：

- mdin/parser 层面的显式 repair policy key。
- 跨进程读取已有 shard manifest 后重建 wrapper。
- HDF5 dataset truncation 或 orphan shard 清理策略。
