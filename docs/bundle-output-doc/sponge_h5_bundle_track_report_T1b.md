# SPONGE H5 Bundle Track T1b 报告

## Track

T1b Concrete HDF5 backend，使用 HighFive 提供 HDF5 backend。

## 已读取的细化文档

本 track 继承 `docs/sponge_h5_bundle_update_plan.md` 中 T1 的 writer foundation 细化章节，并落实总报告中提出的 `T1b Concrete HDF5 backend` 后续项。实现前已读取当前 `H5MDWriter`/`WriterBackend` 抽象、trajectory/observable/restart writer 对 `DatasetSpec` 的实际用法，以及当前 CMake/pixi 构建入口。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 使用 HighFive 实现 concrete backend | 已完成 | `SPONGE/utils/h5md/highfive_backend.hpp` 新增 `SpongeH5MD::HighFiveBackend`。 |
| 保持 backend 隔离在 writer API 后面 | 已完成 | `HighFiveBackend` 实现 `WriterBackend`，高层 writer 不暴露 HighFive 类型。 |
| 实现文件生命周期 | 已完成 | `Open`、`Flush`、`Close`、`Finalize` 均由 `HighFiveBackend` 实现。 |
| 实现 idempotent group creation | 已完成 | `Ensure_Group` 逐级创建 group，已存在时返回成功。 |
| 实现 extendable dataset creation | 已完成 | `Create_Dataset` 使用 HighFive `DataSpace`、max dims 和 chunking 创建可扩展 dataset。 |
| 实现 frame-major append | 已完成 | `Append_Int64`、`Append_Float32`、`Append_Float64` 扩展第一维并写入一条 record。 |
| 实现 hard link | 已完成 | `Create_Hard_Link` 使用 HDF5 C API `H5Lcreate_hard`，在 HighFive file handle 上创建 link。 |
| 实现 string 和 string-array 写入 | 已完成 | `Write_String` 和 `Write_String_Array` 写入 HDF5 string dataset，并支持覆盖已有路径。 |
| 实现 status marker | 已完成 | `Set_Status` 写 `/parameters/sponge/output/status`。 |
| 添加 backend factory | 已完成 | `SpongeH5MD::HighFiveBackendFactory` 实现 `WriterBackendFactory`。 |
| 添加构建依赖 | 已完成 | `pixi.toml` 为 CPU/CUDA/HIP 环境加入 `highfive` 和 `hdf5`；`cmake/targets/SPONGE.cmake` 增加 `find_package(HighFive CONFIG REQUIRED)` 与 `find_package(HDF5 REQUIRED COMPONENTS C)`，并链接 pixi 环境中导出的 HighFive/HDF5 target。 |
| 固定依赖解析来源 | 已完成 | pixi configure task 显式传入 `-DCMAKE_PREFIX_PATH="$CONDA_PREFIX"`，避免优先解析到非 pixi 的系统 HighFive/HDF5。 |

## 依赖来源策略

HDF5 和 HighFive 由 pixi 环境提供。工程不 vendoring HighFive，也不要求用户在系统路径安装 HDF5。CMake 层只声明：

```cmake
find_package(HighFive CONFIG REQUIRED)
find_package(HDF5 REQUIRED COMPONENTS C)
```

pixi task 负责把当前环境的 `$CONDA_PREFIX` 传给 `CMAKE_PREFIX_PATH`，使 CMake package resolution 指向 pixi 环境中的 HighFive/HDF5。这样 CPU、CUDA、HIP 的构建入口共享同一套依赖策略。

## append 语义

当前 facade 中的 `DatasetSpec` 统一按以下语义解释：

```text
第一维 = frame/record 维度
一次 Append_* = 扩展第一维一条 record
count = 单条 record 的元素数 = dims[1] * dims[2] * ...
rank-1 dataset 的单条 record 大小为 1
```

对 `max_dims` 和 `chunk_dims` 的处理规则：

```text
max_dims[0] == 0 且 appendable=true -> 第一维 unlimited
非第一维 max_dims[i] == 0 -> 使用初始 dims[i] 作为固定维度
chunk_dims[i] == 0 -> 使用 1，避免非法 chunk size
```

该规则与现有 trajectory、observable、restart 和 module mapping facade 的 dataset shape 用法对齐。

## 架构说明

T1b 后的 writer 结构为：

```text
TrajectoryH5Writer / ObservableH5Writer / RestartH5Writer / VdsTrajectoryH5Writer
  -> H5MDWriter
    -> WriterBackend
      -> HighFiveBackend
        -> HighFive / HDF5
```

HighFive backend 负责真实 HDF5 文件、group、dataset、link 和 string 写入；高层 writer 仍只依赖 `WriterBackend` 接口。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| 剩余 runtime output manager 缺口 | trajectory、observable-only、restart 和部分 module-specific call site 已初步接入；force H5、QC canonical scalar、metad accumulated restart 等仍待补。 |
| VDS virtual dataset 真实 materialization | 当前 `WriterBackend` 还没有 VDS-specific API；T6 runtime/backend 子阶段需要扩展。 |
| smoke writer 和 HDF5 inspection | 本轮未运行编译或文件验证，遵循当前不验证的工作约束。 |
| backend ABI/API 编译确认 | 未运行 configure/compile。 |

## 审查点

- `HighFiveBackend` 当前使用 overwrite 打开输出文件；append/resume 语义应在 failure/resume track 接入时进一步细化。
- string 写入采用删除后重建 dataset 的方式支持覆盖 status/provenance path。
- hard link 通过 HDF5 C API 创建，因为 HighFive 的高层接口不直接暴露该契约所需的 hard-link helper。
- 该 backend 使 T3/T4/T5 的真实落盘成为可能；当前 runtime call site 已初步接入，但尚未编译、运行或 HDF5 inspection，因此当前仍不是已验证的完整 bundle 输出实现。
