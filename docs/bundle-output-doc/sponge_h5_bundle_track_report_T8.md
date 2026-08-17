# SPONGE H5 Bundle Track T8 报告

## Track

T8 legacy compatibility。

## 已读取的细化文档

实现前已读取 `docs/sponge_h5_bundle_update_plan.md` 中的 T8 细化章节。契约要求：未启用 H5 输出时保持当前 legacy 行为；启用 canonical H5 输出时禁用 implicit legacy sidecar；显式 legacy path key 仍必须写 sidecar。

## 已完成变更

| 项目 | 状态 | 证据 |
|---|---|---|
| 集中显式 legacy path 检测 | 已完成 | T2 `LegacyOutputPlan` 跟踪显式 legacy path；T8 使用 T0 `Legacy_Sidecars_Default_Enabled`。 |
| 启用 H5 输出时阻止隐式 `mdout` 和 `mdinfo` sidecar | 已完成 | `CONTROLLER::Get_Output_File` 在任意 H5 output path 设置时不再 fallback 到默认 suffix/default filename。 |
| 启用 H5 trajectory 输出时阻止隐式 `crd`、`box`、`vel` 和 `frc` sidecar | 通过共享 helper 完成 | 这些输出使用 `Get_Output_File`；显式 path 仍生效，隐式默认值被 gating。 |
| 启用 H5 restart 输出时阻止隐式 `rst` sidecar | 已完成 | `Main_Print` 的 restart cadence 使用 `Should_Write_Legacy_Restart` gating；启用任意 H5 输出后，只有显式 `rst` 才写 legacy restart。 |
| 保持 `qc_scf_output` 显式 | 按策略已完成 | `qc_scf_output` 在 generic helper 中没有隐式默认值；显式 key 行为不变。 |
| 将显式 legacy sidecar path 记录到 H5 metadata | 已完成 | H5 trajectory、VDS wrapper、observable-only 和 restart writer 均写 `/parameters/sponge/files/legacy_sidecars/{key,path}`。 |

## 架构说明

核心 gating 点是 `CONTROLLER::Get_Output_File`：

```text
explicit legacy key -> 写 requested sidecar
no explicit key + no H5 output -> 保持 legacy 默认行为
no explicit key + any H5 output -> 不创建 implicit sidecar
```

由于没有 canonical H5 output path 时 `Legacy_Sidecars_Default_Enabled(controller)` 返回 true，现有 legacy-only workflow 会保持不变。

该规则适用于所有通过 `Get_Output_File` 创建的输出文件，包括 module-specific `prefix_key` output helper 等带前缀形式。

显式 legacy sidecar provenance 写入路径为：

```text
/parameters/sponge/files/legacy_sidecars/key
/parameters/sponge/files/legacy_sidecars/path
```

该记录只包含显式设置的 legacy path key，例如：

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

未显式设置且因 H5 bundle 启用而被关闭的 implicit sidecar 不会被记录为输出。
VDS 模式下 provenance 写入 user-facing wrapper，不重复写入每个 shard。

## 推迟的集成项

| 项目 | 原因 |
|---|---|
| runtime 证明不再产生 sidecar | 需要 backend 集成后运行 H5-enabled case。 |

## 审查点

- 显式 legacy path 仍被保留，因为显式 `Command_Exist` 分支仍在最前面。
- `default_out_file_prefix` 保留为 legacy naming 机制；启用 canonical H5 输出时，它不会用于隐式 sidecar 创建。
- 该变化也可能抑制使用 `Get_Output_File` 的 module sidecar 隐式默认输出。这符合“H5 输出启用后 legacy sidecar 必须显式请求”的全局规则。
- 显式 legacy sidecar provenance 只记录 key/path，不复制 sidecar 内容；sidecar 内容仍由对应 legacy 文件拥有。
