# Bundled I/O A/B 下一轮补强计划（PR 59–74）

## 1. 目标与声明边界

本轮针对 2026-07-14 复核确认的行为证据缺口补强测试与门禁。目标不是用更多
case 数量替代行为证据，而是让每个 supported contract 都能回答以下问题：

1. 哪个 legacy surface 是参考；
2. 哪个 bundled payload 或 writer 被实际消费；
3. 独立 oracle 如何排除默认值、错误字段映射和其他模块贡献；
4. 正常演化、restart continuation、输出组合和执行环境分别覆盖到什么程度；
5. promotion 报告允许声明什么，以及必须排除什么。

本轮**不实现 VDS cross-process reopen-and-append/resume**。契约
`output.vds.cross_process_append_resume` 必须继续保持 `unsupported`，最终声明必须
明确排除该能力。不得使用“全功能等价”“完全替代 legacy I/O”等未限定表述。

Artifact `scope_statement` 固定为：

> Equivalence is limited to explicitly enumerated supported contracts and runtime scenarios and excludes cross-process VDS reopen-and-append/resume.

允许的阶段性结论为：

> 在注册表明确列出的 supported contract、已执行的 CPU/GPU/MPI 场景、给定容差
> 以及不包含跨进程 VDS append/resume 的边界内，legacy 与 bundled I/O 的行为
> 证据满足声明的 E3/E4/F1 门槛。

## 2. 强制 PR/commit 工作方式

本计划继承总计划的规则：每个 PR 恰好一个 commit。

- 从上一 PR 已合并的 HEAD 创建新分支；不得提前混入后续 PR 的文件。
- 完成实现、fixture、测试、注册表、文档和 completion-log 更新后，先运行该 PR
  的全部验收命令，再创建唯一 commit。
- 不 amend 已完成 PR。后续发现的问题必须新开单提交 PR。
- 不提交 `.codegraph/`、`.tmp-*`、build 目录、production 输出或本地 evidence。
- GPU/clean-source evidence 必须基于不可变 commit hash，因此 PR 72/74 的远端证明
  在唯一 commit 创建后运行；若失败，该 PR 不得标记完成，修复另开新 PR。
- 每个 PR 的 completion-log 先写 `This commit`；提交完成后在交付说明中报告 hash。

所有 commit 使用作者：

```text
Xiaoxuan Yu <xiaoxuan_yu@pku.edu.cn>
```

## 3. 共用验证门槛

除各 PR 的专项验证外，每个 PR 至少运行：

```bash
pixi run -e dev-cpu smoke-bundled-io-contract
pixi run -e dev-cpu python -m benchmarks.bundled_io.execution_matrix
pixi run -e dev-cpu ruff check benchmarks/bundled_io
git diff --check
```

涉及 C++ reader/writer/runtime 的 PR 还要运行相关 CTest label，优先使用：

```bash
pixi run -e dev-cpu smoke-bundled-io-staged
pixi run -e dev-cpu smoke-bundled-io-runtime
```

涉及真实 A/B case 时，先运行 medium profile 的目标 case，再运行受影响 case
集合。最终 promotion 前运行完整 production profile。所有命令必须设置独立的
`SPONGE_BUNDLED_IO_AB_RUN_ROOT`，不得复用旧 evidence。

每个 PR 都必须新增至少一个 negative/mutation check，证明新增 assertion 会在目标
缺陷重新出现时失败。

## 4. PR 59：Structural restart 来源闭环

### 目标

修复 `_compare_restart_continuation` 的来源错误，直接证明 bundled-input 生产运行生成
的 H5 restart 能继续出与 legacy text restart 等价的行为。

### 依赖

- 以 PR 58 为基线；无其他新 PR 依赖。

### 主要触点

- `benchmarks/bundled_io/tests/test_bundled_io_ab_production.py`
- `benchmarks/bundled_io/contracts/ab_contracts.json`
- `benchmarks/bundled_io/tests/test_bundled_io_production_gate_manifest.py`
- 本计划和总计划 completion log

### 实施与测试设计

- 删除 H5 continuation 前从 `run.legacy_dir` 覆盖 bundled restart 的复制。
- legacy continuation 只绑定 legacy producer 的 coordinate/velocity text restart。
- bundled continuation 只绑定 bundled producer 自己的 `output/ab.spgr.h5`。
- 在 continuation 启动前记录两个 producer restart 的来源目录、step、time、结构
  payload 摘要和 SHA-256；证据中不得出现 `legacy_restart_copied_to_h5` 路径。
- 比较 continuation 的全部 mdout、position、velocity、force、box 和结束状态。
- mutation test 将 bundled continuation 改指 legacy H5 restart，必须由来源断言拒绝，
  即使数值碰巧相同也不能通过。

### 验收

- `input.restart_load.structural` 的 E4 记录明确标记
  `producer_branch=bundled`。
- 删除 bundled producer restart 或交换来源时，测试在启动 continuation 前失败。
- 真实两步 CPU continuation 通过，且没有放宽既有容差。

### 唯一 commit

```text
test(bundled-io): load bundled structural restart in E4
```

## 5. PR 60：Dynamic restart 契约完整枚举

### 目标

把动态状态从一个泛化 contract 拆成可审计的 module/state contract，避免 NHC 的 E4
证据覆盖 Bussi、pressure barostat 或明确不支持的随机状态。

### 依赖

- PR 59。

### 契约变更

- NHC：`supported/E4`，继续由现有 continuation case 证明。
- Bussi RNG/lambda：暂列 `deferred/E4`，由 PR 61 提升。
- pressure-based barostat RNG/g：暂列 `deferred/E4`，由 PR 62 提升。
- integrator mode/step/time：独立登记并映射到实际 continuation assertion。
- Middle Langevin Philox RNG：`unsupported`，写明无法恢复设备 RNG state。
- Andersen thermostat Philox RNG：`unsupported`，理由同上。
- Monte Carlo barostat C rand state：`unsupported`，写明不可移植恢复限制。

### 主要触点

- `benchmarks/bundled_io/contracts/ab_contracts.json`
- `benchmarks/bundled_io/ab_contracts.py`
- `SPONGE/utils/h5md/input_validation.hpp`
- `tests/h5_bundle/test_h5_input_validation.cpp`
- `tests/control/test_h5_dynamic_state_modules.cpp`

### 验收

- 注册表不能再用一个 NHC case 满足所有 dynamic restart state。
- 删除任一已知 supported/deferred/unsupported module 条目时 inventory gate 失败。
- unsupported payload 的错误类别和稳定诊断 token 有真实 reader/validation 测试。
- coverage 报告分别显示 supported、deferred、unsupported，不把 deferred 计入通过率。

### 唯一 commit

```text
test(bundled-io): enumerate dynamic restart contracts
```

## 6. PR 61：Bussi dynamic continuation E4

### 目标

证明 Bussi RNG engine 和 lambda 的 H5 导出/导入能够复现 continuation，而不只是在
单元测试中序列化 roundtrip。

### 依赖

- PR 60 的独立 Bussi contract。

### 测试设计

- 使用 legacy inputs、固定 seed、非零 `dt` 和非平凡初速度运行一个 Bussi NVT
  reference 到 `checkpoint_step + continuation_steps`，并在 checkpoint/final step 输出
  可审计的 H5 状态作为参考快照。
- 使用同初态的 bundled inputs 另运行 producer 到 checkpoint，输出 structural +
  Bussi H5 state；在 continuation 前先比较 legacy reference 与 bundled producer 的
  prefix、checkpoint 和动态状态，避免 restart 测试掩盖输入分支分歧。
- bundled continuation 从 bundled-input producer 的 H5 checkpoint 恢复并运行相同
  步数。
- 比较 reference suffix 与 restarted suffix 的 step/time、mdout、temperature、lambda、
  position、velocity、force、最终 RNG 文本和 restart state。
- mutation 分别损坏 RNG 文本、删除 lambda、替换 seed；每种 mutation 必须被拒绝或
  产生可检测的 continuation 差异。

### 验收

- Bussi contract 从 `deferred` 提升为 `supported/E4`。
- 至少两步 continuation，`dt > 0`，状态和输出全部有限且非平凡。
- 证据不能由 reader/writer roundtrip 单独满足。

### 唯一 commit

```text
test(bundled-io): close Bussi restart continuation
```

## 7. PR 62：Pressure-based barostat dynamic continuation E4

### 目标

证明 pressure-based barostat 的 RNG 和六分量 `g` 状态能在 NPT continuation 中恢复
box/pressure 行为。

### 依赖

- PR 60；建议在 PR 61 后复用 continuation harness。

### 测试设计

- 使用 legacy inputs、固定 seed、非零 `dt`、非零 barostat 更新次数和可检测的 box
  演化生成连续 reference。
- 使用同初态 bundled inputs 生成 H5 checkpoint；先比较 checkpoint 前缀和动态状态，
  再把 reference 连续后缀与 bundled checkpoint/restart 后缀比较。
- 比较 `g` 六分量、box matrix/volume、pressure/Pxx/Pyy/Pzz、粒子状态、mdout、最终
  RNG 和 restart state。
- mutation 删除 `g`、交换分量、损坏 RNG、让 barostat owner 不匹配，分别要求稳定
  F1 或 E4 差异失败。

### 验收

- pressure-based barostat contract 从 `deferred` 提升为 `supported/E4`。
- box 和 `g` 都必须发生非零变化；静态 `dt=0` 不能满足该 contract。
- 至少覆盖一种实际使用 RNG 的 pressure barostat algorithm。

### 唯一 commit

```text
test(bundled-io): close pressure barostat continuation
```

## 8. PR 63：Bonded/listed 独立 payload-sensitivity

### 目标

为只依赖 full-contract rerun 的 bond、angle、Urey-Bradley、dihedral 和
`custom.listed` 建立独立 payload fingerprint。

### 测试设计

- 每个 contract 使用最小、隔离、确定性 fixture；其他 module contribution 为零或有
  独立解析值。
- correct payload 必须具有解析或预计算的能量与完整 force oracle。
- 每个 case 只修改一个关键字段，例如 `k`、平衡几何、periodicity/phase 或 listed
  参数；要求目标 observable 和 force 按预测改变。
- bundled branch 必须是纯 typed owner，不得保留对应 legacy key/sidecar。
- full-contract rerun 继续作为组合回归，但不再是这些 contract 的唯一 E3 证据。

### 验收

- 五个 contract 各有独立 case、独立 mutation 和 contract-owned oracle details。
- 默认值、零贡献、错误字段映射和其他模块贡献均不能通过。
- 注册表 case_ids 不再只包含 `rerun_full_contract_*`。

### 唯一 commit

```text
test(bundled-io): add bonded payload sensitivity oracles
```

## 9. PR 64：其余 broad-only 与 EDIP payload-sensitivity

### 目标

补齐 sidecar positional restraint、sidecar soft wall、CV、QC energy、ReaxFF 的独立
控制，并给 focused EDIP 增加真正的 payload mutation。

### 测试设计

- sidecar restraint/wall：只修改 sidecar payload，保持 typed/其他协议不变。
- CV：修改 CV 定义或输入几何，验证 CV observable 的解析值和依赖它的消费者。
- QC energy：修改 multiplicity/charge 之外的独立 QC payload，并检查 QC energy 与
  SCF 结果，而不是仅检查列存在。
- ReaxFF：修改一个隔离参数或 atom type，要求目标 ReaxFF 分项与 force 改变。
- EDIP：修改一个 pair/triple parameter，要求 `EDIP` 和完整 force fingerprint 改变。
- `input.full_contract.inventory` 仅保留 E0 结构职责；pure/sidecar aggregate 只汇总
  已有独立 oracle，不得替代模块 contract。

### 验收

- 原 13 个 full-contract-only contract 的集合只剩三个 aggregate/inventory contract。
- 所有功能 contract 均拥有 focused E3 case；EDIP 不再只有“非零且相近”断言。
- mutation tests 对每个新增 oracle 至少拒绝一个错误 payload。

### 唯一 commit

```text
test(bundled-io): add remaining payload sensitivity oracles
```

## 10. PR 65：输入契约非零 dt/多步演化门禁

### 目标

消除 supported 输入 contract 仅由一步或 `dt=0` 证据支撑的情况，覆盖积分、PBC、
邻居表、约束及状态更新。

### 实施设计

- 从注册表自动生成“仅 input_behavior_only / 仅一步 / 仅 dt=0”审计结果，禁止手工
  维护计数。
- 为拓扑/PBC、manybody/custom、protocol/stateful 三个 cohort 增加非零 dt 的短程
  deterministic evolution 或多帧 rerun mutation。
- 每个动态 case 至少三个可比较 frame，并强制一次相关状态变化：邻居表重建、PBC
  穿越、约束投影、bias/state update 或 force-driven coordinate change。
- 若某 contract 物理上只能单点消费，注册表必须写明 `single_point_justification`，并
  由独立 payload mutation 补足；不得静默例外。

### 验收

- 自动审计不再报告无理由的 supported contract 仅一步或仅 `dt=0`。
- 25 个原 input-only contract 均具有动态证据或显式、可审查的单点理由。
- 修改 frame schedule、冻结坐标或跳过状态更新的 mutation 会失败。

### 唯一 commit

```text
test(bundled-io): add nonzero-dt input evolution gates
```

## 11. PR 66：输出族进程级组合矩阵

### 目标

独立证明 trajectory、observable、restart 的全部七种非空组合。

### Case 矩阵

- trajectory-only
- observable-only
- restart-only
- trajectory + observable
- trajectory + restart
- observable + restart
- trajectory + observable + restart

每个组合至少覆盖：默认 legacy suppression、显式 legacy coexistence、正确 provenance、
只生成预期 H5 文件、每个启用 writer 正确 finalize。

### 验收

- 禁用的 writer 不得创建目标文件、临时文件或错误 provenance。
- 单族启用时不依赖其他 writer 初始化。
- 三个单族和三个双族中任一组合缺失都会使 manifest gate 失败。
- restart-only case 仍能产生可读取、可 continuation 的 structural state。

### 唯一 commit

```text
test(bundled-io): cover process output family combinations
```

## 12. PR 67：Writer finalize 故障隔离

### 目标

证明一个 output family 的 append/finalize 失败不会静默破坏其他 output family。

### 实施与测试设计

- 为测试构建提供窄作用域 fault-injection backend/factory；生产默认路径不启用注入。
- 分别注入 trajectory、observable、restart 的 append 和 finalize 失败。
- 定义并测试一致的失败语义：进程非零退出、失败 family 标记 `failed` 和稳定错误；
  已完成的其他 family 保持可打开、completion metadata 一致。
- 若生产 orchestration 当前在首个错误后跳过其他 finalize，本 PR 同时修复该清理顺序，
  但不改变正常路径输出。

### 验收

- 至少六个故障点（3 family × append/finalize）被进程级测试覆盖。
- mutation 若吞掉错误、遗漏其他 writer finalize 或留下 `open` 状态，测试失败。
- 正常七组合矩阵全部回归通过。

### 唯一 commit

```text
test(bundled-io): isolate output writer finalize failures
```

## 13. PR 68：Rerun velocity/force/optional-field 直接闭环

### 目标

把 rerun 的直接语义桥从 position/box 扩展到 velocity、force 和可选字段调度。

### 测试设计

- legacy rerun 使用与输入文件不同的显式输出路径，避免把输入误当输出。
- 将 legacy 输出 sidecar 直接与 bundled H5 trajectory 对应 dataset 比较。
- 覆盖 velocity present/absent、force enabled/disabled、VDS on/off、非默认 particle
  stream、start/strip/frame-limit 边界。
- 对每个字段比较 dataset presence、shape、frame count、step/time schedule 和值。
- mutation 删除一个 frame、错移 step、交换 velocity/force 或错误保留 optional dataset，
  必须失败。

### 验收

- `h5_rerun_semantic_equivalence` details 列出四类粒子 payload 和可选字段决策。
- 不再通过注释中的传递关系替代 velocity/force 的直接比较。
- 现有 position/box 和 mdout 语义保持不变。

### 唯一 commit

```text
test(bundled-io): close rerun optional field semantics
```

## 14. PR 69：正常成功路径全量有限值门禁

### 目标

禁止 legacy 与 bundled 同时产生相同 NaN/Inf 时仍被正常成功 case 判为等价。

### 实施设计

- 在比较器入口统一分类 case：`normal_success`、`failure_semantics`、
  `nonfinite_propagation`。
- `normal_success` 对全部 mdout 数值、trajectory/observable/restart 数值 dataset 和所有
  focused oracle 输入先执行 finite scan，再进行等价比较。
- 非有限值模式比较只允许显式的错误传播 case，且不能产生 supported functionality
  evidence。
- evidence details 记录扫描文件数、dataset 数和值数量。

### 验收

- 对任一正常 case 注入 matching NaN、+Inf 或 -Inf 都必须失败。
- 非有限错误传播测试仍能验证 sign/kind pattern，但 coverage 不计作功能可用。
- 所有正常 production case 的 finite scan 结果可追溯。

### 唯一 commit

```text
test(bundled-io): reject nonfinite successful outputs
```

## 15. PR 70：Oracle-aware evidence 与分层 coverage

### 目标

使 `supported coverage=100%` 不再等价于“只出现过一个 assertion ID”，并明确区分
E0 结构覆盖、E3 行为覆盖、E4 continuation 和 F1 失败语义。

### 实施设计

- coverage 输出至少包含：inventory、conversion、runtime_behavior、continuation、
  failure_semantics 五个分母/分子。
- `input.full_contract.inventory` 的 E0 只进入 inventory coverage。
- 为 assertion details 定义类型化 schema；输入行为 assertion 必须包含
  `oracle_contract_id`、control mutation、预期/实际 delta 和被比较 payload。
- 一个 assertion 可服务多个 contract 时，必须为每个 contract 提供独立 oracle
  details；不能只依赖集合相交。
- promotion evaluator 要求所有 supported runtime contract 达到 E3，restart/dynamic
  达到 E4，invalid contract 达到 F1。

### 验收

- 伪造 assertion ID、空 details、错 contract oracle 或 E0 冒充 E3 都被 mutation test
  拒绝。
- 报告同时展示总体注册表状态和各 evidence class coverage，名称不含歧义。
- 现有真实 evidence 可迁移，不接受手工补写通过记录。

### 唯一 commit

```text
test(bundled-io): make behavior coverage oracle-aware
```

## 16. PR 71：功能 × CPU/MPI 执行矩阵

### 目标

把执行环境轴与高风险功能轴交叉，而不是只运行 TIP3P core。

### 矩阵设计

- 新增 `feature_family` 轴和 `vds` 轴。
- 风险驱动覆盖至少包括：EDIP/ReaxFF、QC、SITS、metadynamics、typed positional/CV
  restraint、soft wall、virtual atom、constraint。
- CPU rank-1 覆盖全部 feature family 的最小 executable case。
- CPU rank-2 覆盖至少 manybody、stateful restart、virtual atom/constraint 三类，并检查
  rank-0 output ownership。
- VDS on/off 分配到不同 feature family，避免矩阵始终固定 VDS-off。
- required combinations 必须声明 environment × feature 组合；只有轴 presence 不算覆盖。

### 验收

- 删除任一 required feature/environment combination 时 matrix validation 失败。
- evidence metadata 同时证明 backend、OMP、MPI、feature family、VDS 和 case IDs。
- CPU rank-1/rank-2 真实 medium matrix 通过，performance metadata 完整有限。

### 唯一 commit

```text
test(bundled-io): cross features with CPU MPI matrix
```

## 17. PR 72：正式 GPU shadow evidence workflow

### 目标

让已枚举 GPU rank-1/rank-2 场景进入可追溯、clean-source、同 run-id 的正式 shadow
workflow，并运行适用的 feature-family 交叉场景。

### 实施设计

- 新增 GPU rank-1 和 GPU rank-2 jobs，使用明确的 self-hosted runner labels。
- CPU、CPU-MPI、GPU、GPU-MPI、contract 和 comparator artifact 共享 GitHub run-id。
- 每个 GPU job 记录 device map、driver/runtime、GPU 型号、MPI rank、OMP 和 rank-0
  ownership。
- 对不适用于 GPU 的 feature/backend 组合在注册表明确 `unsupported` 或
  `not_applicable`；不得用 skip 伪装 evidence。
- artifact merge 只接受相同 source commit/run-id，重复 scenario payload 必须相同。

### 验收

- workflow manifest 测试证明四类 matrix job 都存在并上传 evidence。
- 缺少 GPU artifact、run-id 不一致、环境 metadata 伪造或 scenario skip 均阻止
  promotion derivation。
- 至少完成一次非 promotion 的 clean-source GPU shadow run。

### 唯一 commit

```text
ci(bundled-io): collect GPU matrix evidence
```

### PR 72a：固定 shadow converter 工具链

PR 72 的首次 immutable clean-source GPU run 暴露出 Pixi 中已发布的 Xponge 不含
`legacy-to-bundle`，且本地开发 commit 尚未推送，不能作为 CI pin。单独使用一个修正
PR：checkout 远端可取得的 XPONGE commit，并应用由本仓库审计、SHA-256 固定的最小
runtime patch。所有 shadow jobs 使用同一 converter 根目录；manifest 测试必须锁定
base commit、patch 路径、应用次数和 patch digest。该修正不改变比较容差，也不把
GPU core 统计失败改写成通过。

```text
ci(bundled-io): pin shadow converter toolchain
```

## 18. PR 73：VDS cross-process resume 支持边界门禁

### 目标

本轮不实现 append/resume；本 PR 防止 unsupported 能力被覆盖统计或最终文案误宣称。

### 实施设计

- 保持 `output.vds.cross_process_append_resume` 为 `unsupported/E4` 且无 executable
  case；reason 必须明确 writer 没有 append/resume open mode。
- evidence 和 promotion artifact 新增 `scope_exclusions`，必须包含该 contract。
- coverage 分母分开显示 supported 与 unsupported，unsupported 不能增加或降低
  supported behavior coverage，但必须出现在声明边界。
- 文档、README、promotion summary 和 artifact schema 使用统一限定语。
- mutation tests 删除 exclusion、把 contract 静默改成 supported、或输出“full
  equivalence”时必须失败。

### 验收

- 所有发布候选摘要均包含“不含 cross-process VDS reopen-and-append/resume”。
- 当前 same-process shard、tail repair 和 complete-prefix no-op 证据不能映射到该
  unsupported contract。
- 本 PR 不增加 writer open mode、append API 或跨进程 resume runtime case。

### 唯一 commit

```text
test(bundled-io): enforce VDS resume scope boundary
```

## 19. PR 73b：独立 replica 统计推断与 GPU core 阻断修正

### 触发原因

PR 72 的完整 GPU selector 在 Middle NVT/SETTLE 和 NHC NVT 场景暴露出两个统计
阻断。最初增加单 replica 步数仍不能形成更多独立样本，并且扁平化 H5 分量会把
坐标分量或同一轨迹内 block 误当成独立观测。因此本修正必须先于 PR 74，且不能
amend 已完成的 PR 72/73。

### 实施设计

- production statistical matrix 固定使用至少 48 个独立 seed 的 replica；每个
  replica 保持 10,000 steps、100-step write interval、20-frame burn-in 和 10-frame
  block，共至少 384 个原始 paired blocks。
- 均值等价的 SEM 以每个 replica 的 post-burn block mean 为独立推断单位；block
  数只用于单条轨迹内降噪，不再扩大独立样本数。
- promotion metadata 必须证明 `profile=production`、`fast_mode=false`、至少 48 个
  replica 和 384 个 paired blocks；伪造任一字段都阻止 promotion。
- statistical H5 matrix 只比较有时间演化的 trajectory/observable family；单帧
  restart snapshot 不做伪统计，restart 继续由确定性和 E4 continuation gate 证明。
- position、velocity、force 和 box 不比较扁平原始分量的伪样本，而比较物理时间序列。
  replica 级系综门禁使用全局/分布特征；逐原子 mutation 灵敏度仍由 block-mode
  mutation gate 保留。
- position 保留确定性均匀抽取的 4096 对 PBC 距离分位数；PBC、velocity、force 和
  quantile 计算使用向量化实现，不降低 oracle 分辨率。
- fluctuation 的数值零判定使用机器尺度，不得把均值 absolute margin 当作 std 的
  零阈值；实际波动仍受原 `maximum_std_ratio` 约束。

### 验收

- mutation test 证明大量相关 blocks 不能冒充独立 replicas。
- production/fast sample-plan、restart-family 排除、逐原子 feature 边界和 std 零阈值
  悬崖均有静态或单元回归。
- RTX 4090 上 `middle_nvt_settle_gpu_omp4_rank1` 与
  `nhc_nvt_unconstrained_gpu_omp1_rank1` 均以 48 replicas、原 3σ/容差通过 mdout、
  trajectory 和 observable 门禁。
- 本 PR 只修正统计证据设计；不改变 promotion state，不计入三次 clean production
  evidence。

### 唯一 commit

```text
test(bundled-io): power matrix statistics by replica
```

## 20. PR 74：受限范围 promotion 与三次 clean production evidence

### 目标

在 PR 59–73 全部完成后，用真实 evidence 决定是否从 `shadow` 提升；promotion 只适用
于明确排除 VDS cross-process resume 的 supported scope。

### 前置条件

- supported runtime input/output contract 100% 达到 E3。
- structural/dynamic/protocol/full restart contract 达到 E4。
- failure contract 达到 F1。
- CPU/GPU rank-1/rank-2 required scenario 均有同 run-id evidence。
- 所有 comparator/oracle mutation 被拒绝。
- source tree clean，source commit 完全一致，无 retry。

### 执行步骤

1. 准备 promotion state、受限范围声明和所有静态门禁，运行不依赖 commit hash 的
   pre-commit 检查，然后创建 PR 74 的唯一 candidate commit。
2. 在该 immutable candidate commit 上运行完整 production contract、CPU、CPU-MPI、
   GPU、GPU-MPI 和 comparator jobs；所有 run 必须精确指向该 commit。
3. 合并 artifact 并由工具生成外部/CI history row；不得为记录 evidence 再修改或提交
   仓库文件，也不得手工编辑布尔值或性能数据。
4. 连续完成三次 clean、retry-free production run；任一失败或 retry 都重新计数。
5. 检查 runtime ratio、finalize fraction、output bytes ratio 和 artifact quota。只有
   evaluator 对 candidate commit 返回 ready，PR 74 才能标记完成并合并。
6. 如果 candidate commit 未通过，关闭或保留该 PR 但不得合并；修复必须从 PR 73
   基线另开新的单提交 PR，不能 amend PR 74，也不能增加第二个 commit。

### 验收

- 三个 run-id、source commit、artifact SHA-256 和 scope exclusions 可追溯。
- evaluator 对缺失 GPU、旧 registry、低证据等级、非有限 performance、retry、dirty
  source 或 VDS 误宣称均返回 blocker。
- promotion 文案使用第 1 节限定结论，不出现“全功能等价”。
- candidate commit 在三连证据完成前不得合并或标记为 promoted；不得为了完成计划
  伪造 history 或跳过失败 run。

### 唯一 commit

```text
ci(bundled-io): promote scoped A/B evidence
```

## 21. 依赖与执行顺序

```text
PR59
  └─ PR60 → PR61 → PR62

PR63 → PR64 → PR65

PR66 → PR67
PR68

PR59–68 → PR69 → PR70
PR70 → PR71 → PR72
PR70 → PR73
PR72/PR73 → PR73b
PR59–73b → PR74
```

PR 63–68 在逻辑上可以独立设计，但实际提交仍按编号串行进行，以保持每个分支基于
上一 PR 的单一已审计 HEAD。

## 22. 每 PR completion log 模板

每个 PR 在总计划 completion log 追加一行：

```markdown
| PR NN: title | Complete | This commit | <commands and counts> | <behavior evidence, mutations, remaining boundary> |
```

交付信息必须报告：

- commit hash；
- 修改文件清单；
- 真实运行与 skip 数；
- 新增/变更 contract 状态；
- evidence level 和 case IDs；
- 未解决 blocker；
- 工作树是否仍只包含用户原有未跟踪目录。
