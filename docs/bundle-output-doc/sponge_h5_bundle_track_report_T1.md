# SPONGE H5 Bundle Track T1 报告

## Track

T1 HDF5/H5MD writer 基础层。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T1 细化章节。本 track 的目标是提供可复用的 H5MD writer substrate，包括文件生命周期、layout 初始化、可追加 dataset API、step/time link 接口、字符串/表写入接口和状态标记。

## 当前状态发现

| 发现 | 证据 |
|---|---|
| SPONGE 当前没有已有 HDF5 writer 层。 | 源码扫描未发现已有 `H5F`、`H5D`、`H5S` 或 H5MD writer 实现。 |
| 当前 build/dependency 文件未声明 HDF5。 | 已检查的 `CMakeLists.txt` 和 `pixi.toml` 依赖面中未发现 HDF5 包。 |

由于当前尚未声明 HDF5 依赖，本轮只添加 backend-independent writer foundation，避免直接引入当前构建无法链接的 HDF5 调用。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 定义可复用 writer API | 已完成 | `SPONGE/utils/h5md/h5md_writer.hpp` 定义 `WriterBackend` 和 `H5MDWriter`。 |
| 文件生命周期接口 | 已完成 | backend API 包含 `Open`、`Flush`、`Close`、`Finalize`、`Set_Status` 和 `Status`。 |
| H5MD layout 初始化器 | 已完成 | `H5MDWriter::Initialize_Common_Layout` 创建 `/h5md`、必要时创建 `/particles`、以及 `/observables`、`/parameters`、`/parameters/sponge`。 |
| 可追加 dataset API | 已完成 | backend API 暴露 `Create_Dataset`、`Append_Int64`、`Append_Float32` 和 `Append_Float64`。 |
| step/time hard-link 接口 | 已完成 | backend API 暴露 `Create_Hard_Link`。 |
| 字符串/表写入接口 | 已完成 | backend API 暴露 `Write_String` 和 `Write_String_Array`。 |
| 状态标记接口 | 已完成 | backend API 暴露 `Set_Status(FileStatus)`。 |
| observable-only layout 行为 | 已完成 | `WriterOptions::observable_only` 会阻止创建 `/particles`。 |
| 公共路径常量 | 已完成 | `SpongeH5MD::path` 集中定义常见 H5MD 和 SPONGE 扩展路径。 |
| include 路径暴露 | 已完成 | `SPONGE/control.h` 包含 `utils/h5md/h5md_writer.hpp`。 |

## 推迟的 T1 项

| 项目 | 原因 |
|---|---|
| 真实 HDF5 backend | 需要先在 CMake/pixi 中添加并链接 HDF5 依赖。 |
| 最小可读的落盘 `*.spg.h5md` smoke 文件 | 依赖真实 HDF5 backend。 |
| HDF5 inspection 验证 | 本轮未运行；当前还没有 backend。 |

## 架构说明

基础层采用依赖反转：

```text
MD/output code -> H5MDWriter -> WriterBackend -> concrete HDF5 backend
```

这样未来的输出 call site 不需要依赖原始 HDF5 API。引入 HDF5 后，应在不修改高层 trajectory、observable、restart 或 VDS 逻辑的前提下实现具体 `WriterBackend`。

## 审查点

- 这是 foundation step，不是可工作的 HDF5 backend。
- API 有意保持窄接口，并以 frame append 为中心。
- `observable_only=true` 在 writer layout 层表达，符合输出契约。
- T3/T4/T5 要写出真实 H5 文件前，需要先补上第一个 backend 实现。
