# SPONGE H5 Input Bundle Implementation Plan

This document is the implementation plan for the input side of the SPONGE H5
bundle work. It complements:

- `docs/sponge_input_bundle_inventory.md`
- `docs/sponge_output_bundle_design.md`
- `schemas/mdin.schema.json`

The output bundle work already provides H5/H5MD writer infrastructure and
runtime output routing. The input side must now provide the matching resolver
and reader architecture for topology, protocol, restart, and rerun trajectory
inputs while preserving current legacy behavior until each H5 path is
implemented and validated.

## 1. Contract

### 1.1 Input Files

The new input model is centered on four launch inputs:

```text
topology.spgt.h5      # topology and Hamiltonian container
protocol.spgp.h5      # protocol and sampling method container
restart.spgr.h5       # concrete launch state container
run.mdin              # human-editable binding and run-policy deck
```

Rerun and future trajectory-analysis runs may also bind:

```text
trajectory.spg.h5md   # trajectory/history input for rerun-style streaming
```

The canonical ownership rule is:

```text
container owns canonical data/state
mdin owns binding + small editable launch policy
resolver owns compatibility checks and runtime assembly
```

### 1.2 mdin Binding Surface

The authoritative SPONGE schema is distributed under `schemas/`.
Downstream tools such as Mokda should sync from this repository instead of
maintaining an independent mdin schema authority.

H5 input bindings:

```toml
[input.h5.topology]
path = "topologies/protein.topology.spgt.h5"

[input.h5.protocol]
path = "protocols/metadyn.protocol.spgp.h5"

[input.h5.restart]
path = "runs/prod_0007.restart.spgr.h5"
load = "structural"

[input.h5.trajectory]
path = "runs/prod_0007.spg.h5md"
particle_stream = "all"
```

Parser-visible flattened keys:

```text
input_h5_topology_path
input_h5_protocol_path
input_h5_restart_path
input_h5_restart_load
input_h5_trajectory_path
input_h5_trajectory_particle_stream
```

Current parser-supported TOML aliases must be treated as first-class contract:

```toml
[write.interval]
trajectory = 5000
information = 500
mdout = 500
restart = 500000

[barostat.monte_carlo]
update_interval = 100
```

They normalize to existing runtime keys:

```text
write_trajectory_interval
write_information_interval
write_mdout_interval
write_restart_file_interval
monte_carlo_barostat_*
```

### 1.3 Restart Load Contract

`restart.spgr.h5` always provides a concrete launch state. It is not optional
for H5 bundle launches, including initial launches; newly built systems should
first be materialized into an initial restart container.

`[input.h5.restart] load` selects which components are restored:

| Value | Required behavior |
|---|---|
| `structural` | Load coordinates, box, step/time, and velocity when required by launch policy. This is the default and cannot be disabled. |
| `dynamic` | Load `structural` plus engine continuation state such as integrator, thermostat, barostat, and RNG state when compatible. |
| `protocol` | Load `structural` plus protocol-owned continuation state such as CV/restraint references, SITS state, and metadynamics state when compatible. |
| `full` | Load `structural`, `dynamic`, and `protocol`. |
| `custom` | Reserved until explicit component-list keys are implemented. It must hard-error until then. |

`fresh` is intentionally invalid.

### 1.4 Rerun H5MD Contract

For `mode = "rerun"`, H5MD trajectory input is a history stream, not restart
state. A rerun launch still binds topology/protocol/restart context, then
streams frames from `[input.h5.trajectory]`.

Required H5MD datasets for ordinary rerun:

```text
/particles/<stream>/position/value
/particles/<stream>/box/edges/value
```

Optional datasets:

```text
/particles/<stream>/velocity/value
/particles/<stream>/force/value
```

Frame selection uses existing rerun controls:

```text
rerun_start
rerun_strip
rerun_frame_limit
rerun_need_box_update
```

The default stream is `all`.

## 2. Principles

1. **Resolver first, readers second, runtime last.**
   Do not let individual modules open H5 inputs ad hoc. All H5 input bindings
   are resolved into a single immutable plan before runtime assembly.

2. **Canonical H5 paths are never exposed in mdin.**
   mdin binds container files and small policies only. HDF5 internal paths
   remain part of schema and reader code.

3. **Legacy behavior remains the fallback.**
   If no `[input.h5.*]` bindings are present, existing Xponge/legacy input and
   legacy rerun behavior must remain unchanged.

4. **H5 and legacy inputs must not silently mix for the same role.**
   Example: `input_h5_trajectory_path` and legacy rerun `crd`/`box`/`vel`
   should be mutually exclusive unless a future explicit override policy is
   introduced.

5. **Restart state is componentized.**
   Coordinates and box are structural. Thermostat/barostat/RNG are dynamic.
   SITS/metadynamics/CV/restraint references are protocol-owned continuation
   state. The interface must match that functional split.

6. **Compatibility failures are hard errors.**
   Atom count, atom ordering, schema version, topology hash, and requested
   component availability must fail early with actionable messages.

7. **H5MD rerun consumes wrappers, not shards.**
   VDS shard paths are writer internals. The input binding is the wrapper H5MD
   file.

8. **Tests define the migration envelope.**
   Every phase must add contract tests for the resolver/reader behavior before
   broad runtime integration.

## 3. Current Implementation Snapshot

### 3.1 Interface Audit Conclusions

The previous interface audit produces these binding decisions:

- The old `[bind]` table is not part of the new contract. Binding paths belong
  under `[input.h5.*]` and `[output.h5.*]`.
- Top-level launch keys such as `mode`, `step_limit`, `dt`,
  `target_temperature`, and `target_pressure` must appear before TOML tables in
  examples, because the current parser flattens tables into controller command
  keys and does not treat later top-level keys as a separate section reset.
- `restart.spgr.h5` supports multiple recovery levels, so
  `[input.h5.restart] load = "..."` is the interface point. Current runtime
  support covers structural state, NHC dynamic state, SITS
  `nk`/`log_norm`/`log_nk`, and native/transitional metadynamics state;
  unsupported component requests fail instead of silently falling back to
  structural-only launch.
- Rerun uses trajectory history, not restart state. H5 rerun input therefore
  belongs under `[input.h5.trajectory]`, with ordinary rerun selection keys
  remaining top-level launch policy.
- Schema authority must live in SPONGE. Mokda copies are downstream sync
  artifacts, not the source of truth.
- Parser aliases are valid only when schema, docs, and executable parser tests
  are updated together.

### Implemented or partially implemented

- TOML table flattening exists in `SPONGE/third_party/toml/toml.cpp`.
- Parser aliases exist for `[write.interval]` and `[barostat.monte_carlo]`.
- SPONGE now distributes `schemas/mdin.schema.json`, `schemas/cv.schema.json`,
  `schemas/catalog.json`, and `schemas/tombi.toml`.
- H5 output contracts and writers exist under `SPONGE/utils/control/` and
  `SPONGE/utils/h5md/`.
- `restart.spgr.h5` output writer exists for structural state plus selected
  module state.
- `sponge.input.v2` is the normative input schema version. XPONGE and XpongeCPP
  write `/identity/uuid` consistently across a generated topology, protocol,
  and initial restart set; SPONGE runtime checkpoints also write the launch
  UUID.
- Topology, protocol, restart, and trajectory readers validate their canonical
  schema names, v2 schema versions, and UUID syntax before runtime assembly.
  Restart validation also checks topology/atom-order lineage and recomputes the
  canonical `/run/state_hash` before accepting the state.
- H5 output plan resolver exists.
- H5 input contract constants, input resolver, and input assembler exist.
- H5 restart and H5MD trajectory readers exist.
- H5 topology/protocol/restart/trajectory metadata readers and no-HDF5
  compatibility helpers exist.
- H5 input validation opens bound topology/protocol/restart/trajectory files
  before MD allocation, checks metadata compatibility, and verifies that
  requested `dynamic`, `protocol`, and `full` restart loads contain a currently
  supported payload.
- Restart reader typed payloads exist for writer-produced NHC, SITS, dynamic
  module state, canonical MetaD grid/scatter/edge/hill/runtime state, and
  transitional metadynamics/protocol sidecar text state.
- Runtime apply hooks exist for NHC dynamic state, full SITS
  `nk`/`log_norm`/`log_nk` protocol state, native named MetaD state, and
  transitional metadynamics text materialization for pre-schema-v1 restart
  files.
- Rerun can consume H5MD trajectory input, including VDS wrappers and explicit
  `particle_stream` selection.

### Remaining input-side pieces

- The topology core (mass, charge, exclusions, bonded terms, LJ/nb14/GB,
  virtual atoms, Urey-Bradley, CMAP, LJ soft-core, custom pairwise/listed
  forces, `bond_soft`, and Ryckaert-Bellemans) is read into runtime structures
  from typed HDF5. QC type, EAM, SW, EDIP, Tersoff, and all eight ReaxFF
  runtime modules now load from module-owned typed staging without a temporary
  legacy input file.
- Dynamic restart state is implemented for NHC, typed portable Bussi thermostat
  Philox RNG, pressure-based barostat `g`/typed Philox RNG, integrator
  mode/step/time metadata, portable
  counter-based Middle Langevin/Andersen Philox state, and Monte Carlo
  barostat SplitMix64/adaptive state. All five portable stochastic components
  have independent E4 production A/B continuation cases in addition to the
  focused process-level checkpoint/resume closure and are promoted to
  `supported` in the registry.
- Named CVs, positional/CV restraints, static steering, constraints, typed
  reference coordinates, and canonical MetaD state have direct readers/apply
  hooks. SITS method fields and atom selection now load directly from
  `/sits/method`, `/sits/atom_indices`, or `/sits/atom_numbers_policy`; the
  compatibility `/sits/config` table is injected in memory and no longer
  creates `sits.txt` or `sits_atom.txt`. Soft-wall names and potential source
  are also passed directly to `SOFT_WALLS` without creating
  `soft_walls.txt`. SITS checkpoints capture and restore versioned `nk`,
  `log_norm`, and `log_nk` directly from runtime device state.
- Transitional MetaD and protocol-sidecar text materialization remains solely
  for compatibility with older bundles; it is not a canonical native path.
- XPONGE produces declarative `/wall/hard/*` objects and SPONGE consumes their
  typed bounds and NPT policy directly. Native input conflicts with legacy
  `hard_wall_*` commands instead of applying ambiguous precedence; `/wall/soft/*`
  is also native.
- Generic `/parameters/restart/protocol_state/cv/*` and
  `/parameters/restart/protocol_state/restraint/*` are reserved extension
  roots, not generic runtime apply hooks. Native continuation exists for the
  concrete reference, SITS, and MetaD schemas documented in the inventory.
- Dynamic continuation uses backend-independent counter addressing for Middle
  Langevin/Andersen and module-owned SplitMix64 state for the Monte Carlo
  barostat. Legacy checkpoints carrying the old `unsupported:*` scalar text
  markers continue to hard-error rather than being guessed or silently reset.
- A v2 checkpoint emitted from an H5-bound launch derives topology,
  atom-order, and producer-protocol lineage from its input bundles. Raw-input
  launches must provide canonical `topology_hash` and `atom_order_hash` under
  `[output.h5.restart]`; SPONGE now fails during restart-output initialization
  instead of emitting an incomplete v2 checkpoint. The raw-input planner is
  responsible for supplying those identities.
- XPONGE's post-run `legacy-outputs-to-bundle` importer now emits a launchable
  `sponge.input.v2` restart only when a matching topology bundle and optional
  protocol bundle are supplied and verified. Without that lineage it writes
  `xponge.legacy_restart.archive` with `restart_launchable=false`; both forms
  carry the recomputable canonical state hash.
- Cross-process reopen-and-append for single-file or VDS trajectory writers
  remains out of scope. A new process creates a new output UUID and epoch.
- End-to-end startup smoke tests now cover H5 structural restart launch with a
  minimal H2 NVE case and H5MD rerun launch.
- H5 structural NPT launch is covered by a minimal H2 NPT startup smoke.
- Restart H5 output round-trip is covered by a real H5-input launch followed by
  `RestartH5Reader` readback.
- Trajectory/observable output round-trip validation is covered by a real
  H5-input launch, `TrajectoryH5Reader` readback for particles, and HDF5 layout
  inspection for the observable-only file.

This section is a snapshot marker, not a final status report. As phases are
implemented, keep the todo list below authoritative for remaining work.

## 4. Architecture

### 4.1 Layer Diagram

```text
run.mdin
  |
  | TOML flatten + parser aliases
  v
CONTROLLER commands
  |
  | Input contract helpers
  v
H5InputResolver
  |
  | validates paths, modes, policy, schema versions, compatibility metadata
  v
ResolvedInputPlan
  |
  +--> TopologyH5Reader   -> Runtime topology/Hamiltonian data
  +--> ProtocolH5Reader   -> CV/restraint/protocol definitions
  +--> RestartH5Reader    -> structural/dynamic/protocol restart components
  +--> TrajectoryH5Reader -> rerun frame stream
  |
  v
RuntimeInputAssembler
  |
  +--> Xponge::system bridge where legacy data structures are still required
  +--> MD_INFORMATION coordinate/velocity/box/runtime state
  +--> Module-specific state loaders
```

### 4.2 Main Types

New namespace:

```cpp
namespace SpongeH5InputPlan
```

Core structs:

```cpp
enum class RestartLoadPolicy
{
    Structural,
    Dynamic,
    Protocol,
    Full,
    Custom
};

struct H5Binding
{
    bool enabled = false;
    std::string path;
    bool has_recommended_suffix = false;
};

struct RestartInputPlan
{
    H5Binding binding;
    RestartLoadPolicy load_policy = RestartLoadPolicy::Structural;
};

struct TrajectoryInputPlan
{
    H5Binding binding;
    std::string particle_stream = "all";
};

struct ResolvedInputPlan
{
    bool valid = true;
    std::string error_message;
    bool any_h5_input_enabled = false;
    H5Binding topology;
    H5Binding protocol;
    RestartInputPlan restart;
    TrajectoryInputPlan trajectory;
    bool legacy_input_allowed = true;
};
```

### 4.3 Reader Interfaces

Readers should be narrow and testable:

```cpp
class RestartH5Reader
{
  public:
    bool Open(const std::string& path);
    bool Read_Metadata(RestartMetadata* out);
    bool Read_Structural_State(RestartStructuralState* out);
    bool Read_Dynamic_State(RestartDynamicState* out);
    bool Read_Protocol_State(RestartProtocolState* out);
    std::string Last_Error() const;
};

class TrajectoryH5Reader
{
  public:
    bool Open(const std::string& path, const std::string& particle_stream);
    bool Read_Metadata(TrajectoryMetadata* out);
    bool Read_Frame(int64_t frame_index, TrajectoryFrame* out);
    int64_t Frame_Count() const;
    std::string Last_Error() const;
};
```

Initial readers can target the existing HighFive backend directly. If read-side
mockability becomes painful, introduce a read backend analogous to the write
backend after Phase 1, not before.

### 4.4 Runtime Assembly Boundary

`MD_INFORMATION::Read_Coordinate_And_Velocity()` currently handles both rerun
allocation and legacy Xponge coordinate import. It should not become a general
HDF5 resolver. The desired split is:

```text
main/controller startup:
  resolve input plan
  load topology/protocol/restart into runtime staging objects

MD_INFORMATION:
  allocate runtime arrays
  accept structural state already materialized by the assembler
  run legacy fallback when no H5 input plan is active

RERUN_information:
  keep legacy FILE* path
  add H5MD trajectory-backed frame provider
```

## 5. File-by-File Design

### `schemas/mdin.schema.json`

Role:

- Authoritative schema for SPONGE mdin TOML.
- Must include `[input.h5.*]`, `[output.h5.*]`, parser aliases, and legacy keys.

Design:

- Keep `[input.h5.topology] path`.
- Keep `[input.h5.protocol] path`.
- Keep `[input.h5.restart] path/load`.
- Keep `[input.h5.trajectory] path/particle_stream`.
- Keep `[write.interval]` alias group.
- Keep `[barostat.monte_carlo]` alias group.

Validation:

- JSON parse test.
- AJV validation of the canonical NPT and rerun examples.

### `SPONGE/third_party/toml/toml.cpp`

Role:

- TOML flattening and parser-level alias normalization.

Design:

- Keep generic underscore flattening.
- Normalize only explicit, documented aliases.
- Continue to hard-error when alias and canonical key map to the same runtime
  command.

Current aliases:

```text
write_interval_* -> write_*_interval
barostat_monte_carlo_* -> monte_carlo_barostat_*
```

Future aliases:

- Add only when schema and docs are updated in the same phase.

### `SPONGE/utils/control/h5_input_contract.hpp` (new)

Role:

- Central constants for input-side parser keys, suffixes, default values, and
  enum parsing.

Design:

- Mirror the output contract style in `h5_output_contract.hpp`.
- Define key strings:
  - `input_h5_topology_path`
  - `input_h5_protocol_path`
  - `input_h5_restart_path`
  - `input_h5_restart_load`
  - `input_h5_trajectory_path`
  - `input_h5_trajectory_particle_stream`
- Define recommended suffixes:
  - `.spgt.h5`
  - `.spgp.h5`
  - `.spgr.h5`
  - `.spg.h5md`
- Define default restart load policy: `structural`.
- Define default trajectory particle stream: `all`.

Tests:

- Parse legal load policy values.
- Reject `fresh`.
- Reject `custom` until component lists exist.

### `SPONGE/utils/h5md/input_plan.hpp` (new)

Role:

- Resolve controller commands into `ResolvedInputPlan`.

Design:

- Read only `CONTROLLER` commands, not HDF5 files.
- Validate path presence and mutual exclusion with legacy keys.
- Validate restart load policy.
- Validate suffix warnings separately from hard errors.
- Derive `any_h5_input_enabled`.

Tests:

- Empty mdin produces legacy plan.
- Topology/protocol/restart bindings produce H5 plan.
- H5 rerun trajectory conflicts with legacy `crd`/`box`/`vel`.
- Missing restart in H5 input mode hard-errors.
- `custom` restart load hard-errors.

### `SPONGE/utils/h5md/h5_input_metadata.hpp` (new)

Role:

- Shared metadata structs and compatibility-check helpers.

Design:

- `TopologyMetadata`: schema version, atom count, atom ordering hash, topology
  hash, force-field hash.
- `ProtocolMetadata`: schema version, topology hash, protocol hash, method
  inventory.
- `RestartMetadata`: schema version, atom count, atom ordering hash, producer
  topology/protocol hashes, available component flags.
- `TrajectoryMetadata`: schema version, particle stream, atom count, frame
  count, atom ordering hash, has velocity, has force, has VDS manifest.

Tests:

- Compatibility success and failure cases are unit-tested without opening HDF5.

### `SPONGE/utils/h5md/restart_h5_reader.hpp` / `.cpp` (new)

Role:

- Read `restart.spgr.h5` structural and continuation state.

Design:

- Phase 2 implements structural reads only:
  - `/particles/all/position/value`
  - `/particles/all/velocity/value`
  - `/particles/all/box/edges/value`
  - `/particles/all/step`
  - `/particles/all/time`
- Phase 6 adds dynamic/protocol payload readers:
  - NHC
  - SITS state
  - metadynamics text state
  Runtime apply hooks for dynamic/protocol/full are still separate work.

Tests:

- Read a restart written by existing `RestartH5Writer`.
- Missing structural datasets hard-error.
- Optional velocity is handled according to launch policy.

### `SPONGE/utils/h5md/trajectory_h5_reader.hpp` / `.cpp` (new)

Role:

- Read H5MD trajectory frames for rerun.

Design:

- Open a single H5MD file or VDS wrapper.
- Resolve `/particles/<stream>`.
- Validate required datasets.
- Provide `Read_Frame(frame_index, out)`.
- Convert box edges into existing SPONGE box representation.

Tests:

- Read a trajectory written by existing trajectory writer.
- Read VDS wrapper produced by VDS writer.
- Missing position/box hard-errors.
- Optional velocity is detected and read only when present.

### `SPONGE/MD_core/rerun.hpp`

Role:

- Current legacy rerun frame streamer.

Design:

- Introduce a frame-provider abstraction:

```cpp
class RerunFrameProvider
{
  public:
    virtual bool Read_Frame(int64_t frame_index, RerunFrame* out) = 0;
};
```

- Keep legacy `FILE*` provider unchanged behind this interface.
- Add H5MD provider backed by `TrajectoryH5Reader`.
- `RERUN_information::Iteration()` should consume the provider and then perform
  the existing MPI broadcast/device copy.

Tests:

- Legacy provider behavior unchanged.
- H5 provider maps start/strip/frame_limit correctly.
- End-of-trajectory behavior matches legacy incomplete-frame stop semantics.

### `SPONGE/MD_core/MD_core.cpp`

Role:

- Coordinate/velocity allocation and import into `MD_INFORMATION`.

Design:

- Add a structural-state injection path for H5 restart data.
- Keep legacy Xponge import path when no H5 input plan is active.
- For `mode = RERUN`, initialize rerun with the selected frame provider after
  topology/restart context has been assembled.

Tests:

- Existing legacy startup tests or smoke cases must still pass.
- H5 structural restart startup loads coordinate/velocity/box correctly.

### `SPONGE/main.cpp`

Role:

- Top-level initialization sequence.

Design:

- Resolve input plan after `CONTROLLER::Initial()`.
- Initialize topology/protocol/restart readers before MD arrays are finalized.
- Preserve existing initialization order for legacy mode.
- Pass `ResolvedInputPlan` or assembled staging state to MD/runtime modules.

Tests:

- Legacy input path with no H5 keys has identical command usage behavior.
- H5 plan initialization fails before expensive GPU allocations when bindings
  are invalid.

### `SPONGE/utils/h5md/input_assembler.hpp` (new)

Role:

- Bridge resolved H5 input data into existing SPONGE runtime structures.

Design:

- Own transient host-side `RuntimeInputState`.
- Materialize structural state for `MD_INFORMATION`.
- Route dynamic/protocol state to module-specific loaders.
- Keep Xponge bridge localized while topology H5 support is incomplete.

Tests:

- Structural state is copied without changing atom order.
- Velocity absent/present behavior matches policy.

### `SPONGE/utils/h5md/h5_legacy_sidecar.hpp` (new)

Role:

- Compatibility-only bridge for imported legacy bundles. It may materialize
  old topology/protocol H5 containers into existing runtime input commands
  while native HDF5 readers are missing, but it is not part of the canonical
  new input contract.

Design:

- Read `/parameters/sponge/files/legacy_sidecars/key` and
  `/parameters/sponge/files/legacy_sidecars/path`.
- Resolve relative sidecar paths relative to the containing H5 file.
- Inject only whitelisted input command keys into `CONTROLLER`.
- Reject conflicts with user-provided mdin commands.
- Keep coordinates and velocities out of this bridge; structural launch state
  remains owned by `restart.spgr.h5`.
- New native-H5 bundles should not depend on this bridge. If a required
  topology/protocol object exists only as a sidecar, validation should report
  "native reader missing" or "legacy compatibility input" instead of treating
  that sidecar as canonical state.

Tests:

- HDF5 reader test for key/path arrays and relative path resolution.
- Injection test for allowed topology/Hamiltonian sidecar keys.
- Negative tests for unsupported keys and command conflicts.
- Runtime smoke where mdin omits legacy topology/Hamiltonian keys and startup
  succeeds through the `.spgt.h5` manifest. This smoke proves compatibility
  coverage only; it does not close the native H5 materialization contract.

### `SPONGE/utils/h5md/topology_native_h5_reader.hpp` (partial)

Role:

- Read canonical typed topology/Hamiltonian datasets from `topology.spgt.h5`
  without using legacy text files.

Design:

- Implemented readers: `/atoms/mass`, `/atoms/charge`,
  `/topology/exclusions/{offset,list}`, `/forcefield/lj/{type,params}`,
  `/forcefield/bond/{atoms,k,r0}`, and
  `/forcefield/angle/{atoms,k,theta0}`,
  `/forcefield/dihedral/{atoms,pk,pn,ipn,gamc,gams}`,
  `/forcefield/improper/{atoms,pk,phi0}`, and
  `/forcefield/nb14/{atoms,params}`, `/forcefield/gb/params`, and
  `/forcefield/virtual_atom/{type,atom,from_offset,from,parameter_offset,parameter}`,
  `/forcefield/urey_bradley/{atoms,angle_k,angle_theta0,bond_k,bond_r0}`,
  `/forcefield/cmap/{atoms,type,resolution,grid_value}`, and
  `/forcefield/lj_soft_core/{atom_type_A,atom_type_B,atom_type_count_A,atom_type_count_B,pair_AA,pair_AB,pair_BA,pair_BB}`
  with optional `/forcefield/subsys_division`.
- Remaining module-scoped readers: residues and many-body groups.
- Validate `schema_version`, atom count, zero-based atom indices, CSR offsets,
  units, and local module `source_hash` before returning runtime staging data.
- Return typed staging structs, not controller command strings.
- Report missing native reader coverage through a structured capability error.

Tests:

- Implemented HDF5 fixtures for atom mass/charge, exclusions, LJ, bond, angle,
  dihedral, improper, nb14, GB, and virtual atom minimal datasets.
- Implemented negative fixtures for atom-count mismatch, non-positive mass, bad
  CSR offset, LJ parameter row-count mismatch, dihedral runtime parameter length
  mismatch, nb14 parameter row-count mismatch, virtual atom arity mismatch, and
  GB atom-count mismatch. Remaining: bad-index coverage for each force-field
  family and unsupported module schema version.

### Native protocol H5 readers (incremental implementation)

Role:

- Read canonical typed CV, restraint, constraint, wall, SITS, metadynamics, and
  steering protocol objects from `protocol.spgp.h5`.

Design:

- Resolve object references by `/cv/<name>` and protocol object names, not by
  legacy file order.
- Validate protocol topology hash, CV dimensions, restraint atom indices, and
  metadynamics grid/domain compatibility before module initialization.
- Return typed protocol staging structs that module-specific initializers can
  consume directly.

Current status:

- Named scalar CV objects are read directly by
  `protocol_cv_h5.hpp`; CV-harmonic restraints are read directly by
  `protocol_restraint_h5.hpp`. Positional restraints are staged as native
  runtime data by `xponge/load/native/restraint_h5.hpp`; atom indices, weights,
  full-system reference coordinates, and defaults are passed to
  `RESTRAIN_INFORMATION` without temporary files.
- One enabled named metadynamics object is read by
  `protocol_metadynamics_h5.hpp`, validated against its CV references, and
  injected into the existing runtime `meta` section without creating `cv.txt`.
  Its stable object name is retained for diagnostic and restart paths.
- Static steering is read by `protocol_steer_h5.hpp`, validated against enabled
  scalar named CVs, and injected without creating a text sidecar.
- Distance constraints are read directly into the runtime constraint staging
  structure. SITS definitions and atom selection are read directly by
  `SITS/sits_h5_input.hpp`, with explicit mdin keys taking precedence per
  field. Soft walls are read directly by the module-owned
  `NativeSoftWallH5Reader`; hard-wall bounds and explicit NPT policy are read
  directly by `NativeHardWallH5Reader`.

Tests:

- HDF5 fixtures for a minimal CV graph, restraint object, SITS protocol object,
  and metadynamics protocol definition.
- Negative fixtures for missing CV refs, dimension mismatch, and topology hash
  mismatch.

### Native Protocol Restart State Readers (incremental implementation)

Role:

- Replace protocol sidecar text materialization with native HDF5 state readers
  for `restart.spgr.h5`.

Design:

- Read restraint/CV references from `/parameters/restart/references`.
- Read SITS state from `/parameters/restart/bias/sits/<module>`.
- Read metadynamics grid/scatter/hills/history from
  `/parameters/restart/bias/meta/<name>` as typed arrays.
- Read module-owned protocol continuation state from
  `/parameters/restart/protocol_state/<module>/<name>`.
- Apply state through module-level hooks only; H5 readers do not mutate modules.

Named restraint and CV reference coordinates use typed restart datasets and are
validated/applied by their owning modules. Metadynamics protocol objects and
their accumulated grid/scatter/edge/hill/runtime state now use named native
objects. SPONGE checkpoints write state schema v1 directly from `META`; startup
validates identity and layout before copying state to host/device buffers.
XPONGE's earlier typed conversion layout remains a transition read format, and
text materialization remains only for older nodes carrying text snapshots.

Tests:

- Writer/reader/apply round trips for SITS, CV references, restraint reference
  coordinates, and metadynamics typed grid/hills state.
- Runtime closure for canonical typed MetaD scatter+sink continuation covers
  scatter point-order preservation, accumulated scatter potential/force,
  edge normalization/force, hill sink history, finite force/thermo execution,
  and typed checkpoint re-export.
- Negative fixtures for protocol hash mismatch and state dimension mismatch.

### `tests/h5_bundle/`

Role:

- Main H5 input/output contract tests.

Design:

- Add input tests alongside current output tests:
  - input plan resolver
  - restart reader using writer-produced files
  - trajectory reader using writer-produced files
  - VDS trajectory input
  - compatibility helpers

### `tests/control/`

Role:

- Parser and schema-adjacent contract tests that do not require HDF5.

Design:

- Keep TOML alias tests here.
- Add restart-load policy parsing tests if implemented in header-only helpers.

## 6. Todo List

Each phase has a matching detailed plan in Section 7. A phase is considered
complete only when its code, schema/docs, and validation gate are all complete.

### Phase 0: Authoritative mdin schema and parser alias baseline

Status: complete for the schema/alias baseline.

Tasks:

- [x] Distribute `schemas/mdin.schema.json`.
- [x] Distribute `schemas/catalog.json`.
- [x] Distribute `schemas/tombi.toml`.
- [x] Add `[write.interval]` parser alias.
- [x] Add `[barostat.monte_carlo]` parser alias.
- [x] Add parser alias tests.
- [x] Add a schema validation test target independent of Mokda.
- [x] Add a script to compare Mokda schema copies against SPONGE authority.

### Phase 1: Input contract constants and resolver

Tasks:

- [x] Add `h5_input_contract.hpp`.
- [x] Add `input_plan.hpp`.
- [x] Parse and validate `[input.h5.*]`.
- [x] Enforce H5 vs legacy mutual exclusion for restart structural state and
      rerun trajectory input.
- [x] Add resolver unit tests.

### Phase 2: Structural restart reader

Tasks:

- [x] Add header-only `restart_h5_reader.hpp`.
- [x] Read metadata and structural datasets.
- [x] Load coordinates, velocities, box, step, and time into a staging struct.
- [x] Add round-trip test source against files produced by `RestartH5Writer`.
      Verified through CTest in the `dev-cuda13` HighFive/HDF5 environment.

### Phase 3: Runtime structural assembly

Tasks:

- [x] Add `input_assembler.hpp`.
- [x] Integrate structural restart load into `MD_INFORMATION`.
- [x] Preserve legacy Xponge path when no H5 restart binding is present.
- [x] Add startup smoke tests for H5 structural restart.
      Verified with a generated H2 H5 fixture: `input_h5_restart_path` starts
      without `coordinate_in_file`/`velocity_in_file` and reads coordinates,
      velocity, box, step, and time from `restart.spgr.h5`.

### Phase 4: H5MD trajectory reader for rerun

Tasks:

- [x] Add header-only `trajectory_h5_reader.hpp`.
- [x] Add explicit VDS wrapper manifest validation.
- [x] Honor `[input.h5.trajectory] particle_stream` when resolving H5MD paths.
- [x] Add H5 rerun frame provider path.
- [x] Integrate with `RERUN_information`.
- [x] Test H5 frame-cursor start/strip and end-of-trajectory selection helper.
- [x] Test runtime start/strip/frame_limit behavior for H5MD rerun.
      Verified with a generated H2 H5MD fixture using
      `input_h5_trajectory_path`, no legacy `crd`/`box`/`vel`, and
      `rerun_frame_limit = 1`.
- [x] Test runtime end-of-trajectory behavior.
      Verified with the generated two-frame H2 H5MD fixture and no
      `rerun_frame_limit`: runtime consumes both frames, stops at EOF, and exits
      successfully without legacy `crd`/`box`/`vel`.

### Phase 5: Topology and protocol H5 metadata readers

Tasks:

- [x] Add minimal topology metadata reader.
- [x] Add minimal protocol metadata reader.
- [x] Add shared metadata structs and no-HDF5 compatibility helpers.
- [x] Implement HDF5-backed topology/protocol/restart/trajectory metadata
      reader sources.
- [x] Add HDF5-backed input validation tests for topology/protocol/restart and
      trajectory compatibility.
- [x] Implement compatibility topology/Hamiltonian materialization through the
      H5 legacy sidecar manifest bridge into existing native loaders.
- [x] Complete canonical native HDF5 topology/Hamiltonian materialization
      through typed `/atoms`, `/topology`, `/forcefield`, `/manybody`, and
      `/qc` datasets. Partial support exists for mass, charge, exclusions, LJ,
      bond, angle, dihedral, improper, nb14, GB, virtual atom, Urey-Bradley,
      CMAP, LJ soft-core, QC type, EAM, SW, EDIP, Tersoff, ReaxFF, custom
      pairwise/listed forces, and `bond_soft`. QC and every canonical
      many-body/custom-force runtime initialize from typed values directly.
      ReaxFF is read once into a shared typed definition and its EEQ,
      bond-order, bond, VDW, over/under, angle, torsion, and hydrogen-bond
      modules stage from that object; native launches do not create
      `.sponge_h5_native_qc`,
      `.sponge_h5_native_manybody/{eam,eam_atom_type,sw,edip,tersoff,reaxff,reaxff_type}.txt`, or
      `.sponge_h5_native_custom_force`; the custom-force and EDIP text
      materializers and the unused ReaxFF compatibility serializer remain only
      as explicit legacy compatibility/test utilities. SW has no runtime text
      materializer.

### Phase 6: Dynamic and protocol restart state

Tasks:

- [x] Add typed reader payload for writer-produced NHC dynamic state.
- [x] Add typed reader payloads for writer-produced SITS and metadynamics
      protocol state.
- [x] Implement Nose-Hoover chain `dynamic` restart load in the runtime module.
- [x] Implement remaining supported `dynamic` restart load components in
      runtime modules: Bussi thermostat RNG, pressure-based barostat `g`/RNG,
      and integrator mode/step/time metadata.
- [x] Implement portable continuation for Bussi, the pressure-based barostat,
      Middle Langevin, and Andersen with a backend-independent counter Philox
      engine, and for the Monte Carlo barostat with module-owned SplitMix64
      plus adaptive counters.
- [x] Add a focused process-level uninterrupted-versus-checkpoint/resume
      closure for Bussi, the pressure-based barostat, Middle Langevin,
      Andersen, and the Monte Carlo barostat.
- [x] Bind those continuation checks into independent production A/B
      cases/assertions and promote their registry entries from `unsupported`
      to `supported`.
- [x] Implement SITS `nk`, `log_norm`, and `log_nk` protocol restart load and
      versioned checkpoint export in the runtime module.
- [x] Implement native `/sits/method`, `/sits/atom_indices`, and
      `/sits/atom_numbers_policy` loading without temporary text files.
- [x] Implement transitional metadynamics text protocol restart materialization
      for existing legacy sidecar loaders.
- [x] Implement transitional protocol restart sidecar text materialization for
      CV/restraint style file-backed state before module initialization.
- [x] Implement native typed protocol restart components for
      `/parameters/restart/references` and `/parameters/restart/bias/meta` with
      module-owned HDF5 schemas/apply hooks.
- [x] Keep generic `/parameters/restart/protocol_state/cv/*` and
      `/parameters/restart/protocol_state/restraint/*` as reserved extension
      roots. Concrete mutable state must gain a module-owned schema/apply hook;
      unsupported payloads fail validation instead of entering a text fallback.
- [x] Implement `full` as dynamic + protocol runtime load for currently
      supported components.
- [x] Extend `full` to include the newly supported Bussi RNG,
      pressure-based barostat, and integrator metadata dynamic components.
- [ ] Extend `full` when future dynamic/protocol component hooks are added.
- [x] Keep `custom` reserved with a clear error.
- [x] Prevent unsupported `dynamic`, `protocol`, and `full` components from
      silently degrading to `structural` before runtime module apply hooks
      exist.
- [x] Add restart reader payload round-trip tests for NHC/SITS/metadyn state.
- [x] Add module-level apply helper tests for NHC dynamic state.
- [x] Add module-level apply helper tests for legacy SITS `nk` protocol state.
- [x] Add an isolated process-level SITS apply/export closure for `nk`,
      `log_norm`, and `log_nk` device restoration. The fixture deliberately
      supplies `log_nk != log(nk)` and requires exact typed checkpoint
      write-back, proving the apply hook restored the device state instead of
      reconstructing it during initialization.
- [x] Add validation tests that reject requested dynamic/protocol loads when
      the restart only contains currently unsupported payloads.
- [x] Add validation test that accepts `full` only when supported dynamic and
      protocol payloads are both present.
- [x] Add writer/reader and materializer tests for protocol sidecar text state.
- [x] Add module-level state loader/apply tests for Bussi thermostat and
      pressure-based barostat dynamic state.
- [ ] Add module-level state loader/apply tests for future dynamic and
      structured protocol components as their schemas/apply hooks are added.

### Phase 7: End-to-end launch integration

Tasks:

- [x] Wire `ResolvedInputPlan` into `main.cpp`.
- [x] Fail early on invalid H5 bindings and incompatible H5 metadata.
- [x] Run H5 structural NPT launch.
      Verified with the generated H2 NPT fixture: H5 structural restart input
      reaches normal NPT startup and initializes middle Langevin plus MC
      barostat without legacy coordinate input.
- [x] Run H5 structural restart launch smoke.
      Verified with the generated H2 NVE fixture. This validates the input-side
      structural restart path but does not close the NPT-specific gate.
- [x] Run H5MD rerun launch.
      Verified with the generated H2 H5MD fixture. The runtime opens
      `trajectory/prod.spg.h5md` and honors `rerun_frame_limit = 1`.
- [x] Verify restart H5 output round-trip.
      Verified by running the generated H2 H5 restart-output fixture and reading
      `output/restart_out.spgr.h5` back with `RestartH5Reader`.
- [x] Verify trajectory/observable H5 output round-trip from an H5-input launch.
      Verified by running the generated H2 H5 trajectory-output fixture,
      reading `output/traj_out.spg.h5md` back with `TrajectoryH5Reader`, and
      inspecting `output/obs_out.obs.spg.h5md` HDF5 observable layout.

### Phase 8: Migration, docs, and downstream sync

Tasks:

- [x] Update `docs/input-reference/*.md`.
- [x] Update Mokda schema sync workflow.
- [x] Add migration notes for legacy `coordinate_in_file`, `velocity_in_file`,
      `crd`, `box`, and `vel`.
- [x] Add example mdin files under a stable examples directory.

## 7. Detailed Phase Plans

### Phase 0 Detailed Plan: Authoritative mdin schema and parser alias baseline

Objective:

Make SPONGE itself the mdin schema authority and keep parser-visible TOML
aliases executable, tested, and documented.

Primary files:

- `schemas/mdin.schema.json`
- `schemas/catalog.json`
- `schemas/tombi.toml`
- `CMakeLists.txt`
- `SPONGE/third_party/toml/toml.cpp`
- `tests/control/test_toml_command_aliases.cpp`
- `scripts/validate_mdin_schema_examples.py`
- `scripts/check_mdin_schema_sync.py`

Implementation steps:

1. Keep schema files under `schemas/`.
2. Install the schema directory through CMake.
3. Maintain parser aliases in `CanonicalCommandKey()` only.
4. Add unit tests in `tests/control`.
5. Add a lightweight schema validation test that parses canonical TOML examples
   and validates them against `schemas/mdin.schema.json`.
6. Add an optional script:

```text
scripts/check_mdin_schema_sync.py ../Mokda
```

The script should compare SPONGE schema files with known Mokda copies and print
actionable sync commands.

Acceptance criteria:

- TOML aliases normalize to canonical commands.
- Alias/canonical conflicts hard-error.
- Schema examples validate.
- SPONGE install includes schema files.

Validation gate:

```text
python3 -m json.tool schemas/mdin.schema.json
python3 scripts/validate_mdin_schema_examples.py --schema schemas/mdin.schema.json
tests/control/test_toml_command_aliases.cpp through direct compile or CTest
```

Risks:

- Full CMake tests currently depend on HighFive configuration. Keep parser and
  schema tests capable of running without HDF5.

### Phase 1 Detailed Plan: Input contract constants and resolver

Objective:

Introduce an input-side plan resolver without opening HDF5 files or changing
runtime behavior.

Primary files:

- `SPONGE/utils/control/h5_input_contract.hpp`
- `SPONGE/utils/h5md/input_plan.hpp`
- `tests/control/test_h5_input_plan.cpp`
- `schemas/mdin.schema.json`
- `docs/sponge_input_bundle_inventory.md`

Implementation steps:

1. Add `SPONGE/utils/control/h5_input_contract.hpp`.
2. Add `SPONGE/utils/h5md/input_plan.hpp`.
3. Implement:
   - `Parse_Restart_Load_Policy()`
   - `Resolve_Input_Plan(CONTROLLER*)`
   - `Has_H5_Input_Binding(CONTROLLER*)`
4. Validate:
   - topology/protocol/restart path keys
   - trajectory path only for rerun or future analysis modes
   - restart load policy
   - legacy conflict sets
5. Add unit tests not requiring HDF5.

Conflict rules:

| H5 key | Conflicts with |
|---|---|
| `input_h5_restart_path` | `coordinate_in_file`, `velocity_in_file` for launch structural state, unless running legacy compatibility mode |
| `input_h5_trajectory_path` | rerun `crd`, `box`, `vel` |
| `input_h5_topology_path` | native force-field/topology file groups once topology H5 materialization is implemented |
| `input_h5_protocol_path` | protocol-owned legacy CV/restraint/bias files once protocol H5 materialization is implemented |

Initial implementation may warn rather than hard-error for topology/protocol
legacy overlap until topology/protocol readers are available, but trajectory and
restart structural conflicts should hard-error.

Acceptance criteria:

- Resolver returns legacy plan when no H5 input keys are present.
- Resolver returns valid H5 plan for standard NPT example.
- Resolver rejects invalid restart load values.
- Resolver rejects H5 trajectory plus legacy rerun files.

Validation gate:

```text
tests/control/test_h5_input_plan.cpp through direct compile or CTest
python3 scripts/validate_mdin_schema_examples.py --schema schemas/mdin.schema.json
```

Done boundary:

This phase is complete when a caller can ask one question, "what H5 input plan
does this controller imply?", without opening any HDF5 file or touching runtime
MD state.

### Phase 2 Detailed Plan: Structural restart reader

Objective:

Read the structural subset of `restart.spgr.h5` and expose it as a staging
object that can later be injected into `MD_INFORMATION`.

Primary files:

- `SPONGE/utils/h5md/h5_structural_state.hpp`
- `SPONGE/utils/h5md/restart_h5_reader.hpp`
- `SPONGE/utils/h5md/h5_input_metadata.hpp`
- `tests/h5_bundle/test_restart_h5_reader.cpp`
- Existing writer files under `SPONGE/utils/h5md/`

Implementation steps:

1. Add `RestartStructuralState`:

```cpp
struct RestartStructuralState
{
    int64_t step = 0;
    double time = 0.0;
    std::size_t atom_count = 0;
    std::vector<float> position_xyz;
    std::vector<float> velocity_xyz;
    std::array<float, 9> box_edges;
    bool has_velocity = false;
};
```

2. Add `RestartH5Reader`. Done as a header-only reader using HighFive directly.
3. Implement structural reads from `/particles/all`. Done for step/time,
   position, optional velocity, and box edges.
4. Validate shapes. Done with hard errors for missing or malformed structural
   datasets:
   - position: `[1,N,3]`
   - velocity: `[1,N,3]` when present
   - box edges: `[1,3,3]`
   - step/time: `[1]`
5. Add tests that write a restart through `RestartH5Writer`, then read it back.
   Done in `tests/h5_bundle/test_restart_h5_reader.cpp`.

Acceptance criteria:

- Reader rejects missing position or box.
- Reader reports atom count.
- Reader preserves step/time.
- Reader treats velocity as optional at read time and lets policy decide
  whether absence is acceptable.

Validation gate:

```text
test_restart_h5_reader writes restart.spgr.h5 with RestartH5Writer
test_restart_h5_reader reads the same file with RestartH5Reader
```

Done boundary:

This phase does not mutate `MD_INFORMATION`. It only proves that restart H5
structural data can be read into a typed staging object with shape validation.

### Phase 3 Detailed Plan: Runtime structural assembly

Objective:

Use structural restart data as the source of coordinates, velocities, box, step,
and time for H5 launches.

Primary files:

- `SPONGE/utils/h5md/input_assembler.hpp`
- `SPONGE/MD_core/MD_core.cpp`
- `SPONGE/MD_core/MD_core.h`
- `SPONGE/main.cpp`
- `tests/control/test_h5_input_assembler.cpp`

Implementation steps:

1. Add `RuntimeInputState` in `input_assembler.hpp`. Implemented as
   `Apply_Restart_Structural_State`, keeping the staging object in
   `h5_structural_state.hpp`.
2. Convert `RestartStructuralState` into:
   - `MD_INFORMATION::coordinate`
   - `MD_INFORMATION::velocity`
   - `MD_INFORMATION::sys.box_length` / cell representation
   - `MD_INFORMATION::sys.start_time`
   Done for coordinate, optional/zero-filled velocity, box length/angle, and
   start time.
3. Decide step semantics:
   - preserve restart step for continuation metadata
   - keep current simulation loop step initialization unless the engine already
     supports continuing from nonzero step
   Current implementation prints the restart step but does not yet move
   `sys.steps`.
4. Add a narrow branch in `MD_INFORMATION::Read_Coordinate_And_Velocity()`:
   - H5 structural state present: allocate and copy from state
   - no H5 state: legacy behavior
   Done for non-rerun structural restart.
5. Keep GPU allocation and copies in the existing MD layer.
   Done by reusing the same allocation/copy pattern as legacy coordinate
   readers.

Acceptance criteria:

- No H5 keys: legacy path unchanged.
- H5 restart structural state: coordinates and box come from restart H5.
- Missing velocity with NVE/NVT policy is either zero-filled or hard-error by an
  explicit rule documented in the resolver.

Validation gate:

```text
tests/control/test_h5_input_assembler.cpp through direct compile or CTest
one H5 structural NPT startup smoke test
one H5 structural restart startup smoke test without legacy coordinate input
one legacy startup smoke test
```

Done boundary:

This phase is complete for the structural restart assembly path when H5 restart
reaches normal MD startup without requiring legacy coordinate or velocity files.
The NPT-specific smoke remains under Phase 7.

### Phase 4 Detailed Plan: H5MD trajectory reader for rerun

Objective:

Let rerun consume `trajectory.spg.h5md` or a VDS wrapper through
`[input.h5.trajectory]`.

Primary files:

- `SPONGE/utils/h5md/trajectory_h5_reader.hpp`
- `SPONGE/MD_core/rerun.h`
- `SPONGE/MD_core/rerun.hpp`
- `SPONGE/utils/h5md/input_plan.hpp`
- `tests/h5_bundle/test_trajectory_h5_reader.cpp`

Implementation steps:

1. Add `TrajectoryFrame`:

```cpp
struct TrajectoryFrame
{
    int64_t step = 0;
    double time = 0.0;
    std::vector<float> position_xyz;
    std::vector<float> velocity_xyz;
    std::array<float, 9> box_edges;
    bool has_velocity = false;
};
```

Implemented by reusing `RestartStructuralState` as the frame staging object,
because rerun needs the same step/time/position/velocity/box payload shape.

2. Add `TrajectoryH5Reader`. Done as a header-only HighFive reader. It accepts
   the resolved `particle_stream` and builds `/particles/<stream>/...` paths
   instead of always reading `/particles/all`.
3. Add legacy and H5 implementations of a rerun frame provider. Done directly
   in `RERUN_information::Iteration()`; a separate provider abstraction can be
   added later if the branch grows.
4. Modify `RERUN_information::Initial()` to choose provider:
   - H5 trajectory path present: H5 provider
   - otherwise: legacy provider
   Done through the resolved input plan.
5. Modify `RERUN_information::Iteration()` to request a frame from provider,
   then keep existing MPI/device-copy logic.
   Done for rank-0 H5 reads plus existing MPI/device-copy flow.
6. Implement end-of-trajectory behavior matching legacy incomplete reads.
   Done for H5MD runtime EOF: out-of-range frame selection sets
   `sys.step_limit`, and the generated two-frame H2 EOF smoke exits cleanly
   after two printed frames.
7. Normalize `rerun_frame_limit` in `system_information::Initial()` to the main
   loop's inclusive `step_limit` convention. Done by mapping an N-frame limit to
   `step_limit = N - 1`, with a runtime smoke covering `rerun_frame_limit = 1`.

Acceptance criteria:

- Rerun reads H5MD position and box frames.
- Requested `particle_stream` is honored and fails if the stream is absent.
- `rerun_start`, `rerun_strip`, and `rerun_frame_limit` select the expected
  frame indices.
- VDS wrappers are accepted when the VDS maps complete frames.
- Legacy rerun still works without H5 trajectory binding.

Validation gate:

```text
test_trajectory_h5_reader reads direct H5MD writer output
test_trajectory_h5_reader reads VDS-wrapper output
test_h5_input_validation rejects missing requested trajectory stream
rerun smoke test with input_h5_trajectory_path and no crd/box/vel
legacy rerun smoke test with crd/box/vel
rerun_frame_limit runtime smoke for H5MD trajectory input
```

Done boundary:

This phase is complete for H5MD runtime launch and frame-limit handling. Full
completion now includes an explicit H5MD end-of-trajectory runtime smoke for the
direct H5MD reader path. Legacy rerun parity should still be watched by existing
legacy tests when touched.

### Phase 5 Detailed Plan: Topology and protocol H5 metadata readers

Objective:

Add enough topology/protocol H5 reading to enforce compatibility checks before
full topology/protocol runtime materialization. This phase is now explicitly
split into Phase 5a metadata/compatibility and Phase 5b native H5
materialization.

Primary files:

- `SPONGE/utils/h5md/topology_h5_reader.hpp`
- `SPONGE/utils/h5md/protocol_h5_reader.hpp`
- `SPONGE/utils/h5md/h5_input_metadata.hpp`
- `SPONGE/utils/h5md/h5_legacy_sidecar.hpp`
- `SPONGE/utils/h5md/topology_native_h5_reader.hpp`
- `SPONGE/utils/h5md/protocol_native_h5_reader.hpp`
- `SPONGE/utils/h5md/input_validation.hpp`
- `SPONGE/utils/h5md/input_plan.hpp`
- `SPONGE/main.cpp`
- `tests/control/test_h5_input_metadata.cpp`
- `tests/h5_bundle/test_h5_input_validation.cpp`

Implementation steps:

1. Add `TopologyH5MetadataReader`.
2. Add `ProtocolH5MetadataReader`.
3. Read schema version, hashes, atom count, and method inventories.
4. Implement compatibility functions:
   - topology vs restart
   - topology vs trajectory
   - topology vs protocol
   - protocol vs restart protocol components
5. Add `input_validation.hpp` to open bound H5 files and apply metadata
   compatibility checks before MD allocation. Done for topology/protocol,
   restart structural state, and trajectory input.
6. Add the H5 legacy sidecar manifest bridge as a compatibility-only bridge:
   - read `/parameters/sponge/files/legacy_sidecars/{key,path}` from topology
     and protocol containers.
   - resolve relative paths against the containing H5 file.
   - inject whitelisted command keys before `Xponge::system.Load_Inputs()`.
   - reject unsupported keys and conflicts with explicit mdin commands.
7. Add native HDF5 topology/Hamiltonian materializers:
   - done for `/atoms/mass`, `/atoms/charge`,
     `/topology/exclusions/{offset,list}`, `/forcefield/lj/{type,params}`,
     `/forcefield/bond/{atoms,k,r0}`, and
     `/forcefield/angle/{atoms,k,theta0}`,
     `/forcefield/dihedral/{atoms,pk,pn,ipn,gamc,gams}`,
     `/forcefield/improper/{atoms,pk,phi0}`, and
     `/forcefield/nb14/{atoms,params}`, `/forcefield/gb/params`, and
     `/forcefield/virtual_atom/{type,atom,from_offset,from,parameter_offset,parameter}`,
     `/forcefield/urey_bradley/{atoms,angle_k,angle_theta0,bond_k,bond_r0}`,
     `/forcefield/cmap/{atoms,type,resolution,grid_value}`, and
     `/forcefield/lj_soft_core/{atom_type_A,atom_type_B,atom_type_count_A,atom_type_count_B,pair_AA,pair_AB,pair_BA,pair_BB}`
     with optional `/forcefield/subsys_division`.
   - residue atom counts are loaded from `/residues/atom_offset`; atom display
     names/types under XPONGE provenance metadata are not runtime Hamiltonian
     inputs; all contract-declared many-body groups have module-specific typed
     readers.
   - reject launch when a required module only has legacy sidecar text and
     native mode is requested.
8. Add native HDF5 protocol materializers:
   - read `/cv`, `/restraint`, `/constraint`, `/sits`, `/meta`, `/steer`, and
     `/wall` typed protocol objects.
   - validate CV references, dimensions, atom indices, and topology hash before
     module initialization.
   - feed module-specific staging structs directly to module initializers.

Acceptance criteria:

- Mismatched atom count hard-errors.
- Mismatched atom ordering hash hard-errors when both sides provide hashes.
- Protocol-owned state is checked only when `load = protocol` or `full`.
- Phase 5a compatibility: a launch can omit mdin `mass_in_file`,
  `charge_in_file`, and `qc_type_in_file` when `topology.spgt.h5` supplies
  those paths through the legacy sidecar manifest.
- Phase 5b native: the same launch can omit those keys because
  `topology.spgt.h5` supplies native `/atoms/mass`, `/atoms/charge`, and
  `/qc/type` datasets, with no legacy text sidecar dependency.

Validation gate:

```text
tests/control/test_h5_input_metadata.cpp through direct compile or CTest
negative fixtures for mismatched atom_count/topology_hash/protocol_hash
test_h5_input_validation verifies HDF5-backed metadata compatibility failures
startup fails before MD allocation when metadata is incompatible
test_h5_legacy_sidecar verifies HDF5 sidecar manifest read/inject behavior
runtime smoke for H5 topology sidecar materialization without mdin native keys
native topology fixture stages mass/charge/qc_type without legacy sidecars or temporary text
native protocol fixture materializes CV/restraint/metadynamics definitions without legacy sidecars
```

Done boundary:

Phase 5a stops at metadata validation plus compatibility sidecar-backed runtime
materialization. Phase 5b is complete only when the H5 input path can launch
from native topology/protocol datasets without relying on embedded or external
legacy text for required topology, Hamiltonian, CV, restraint, or bias
definitions.

### Phase 6 Detailed Plan: Dynamic and protocol restart state

Objective:

Implement restart component policies beyond structural state.

Primary files:

- `SPONGE/utils/h5md/restart_h5_reader.hpp`
- `SPONGE/utils/h5md/input_assembler.hpp`
- `SPONGE/MD_core/MD_core.cpp`
- `SPONGE/SITS/SITS.cpp`
- `SPONGE/SITS/SITS.h`
- Barostat/thermostat module files that own continuation state
- Metadynamics/CV module files that own protocol continuation state

Implementation steps:

1. Add typed restart payload structs. Done for NHC dynamic state, generic
   dynamic RNG/integrator/thermostat/barostat state maps, and writer-produced
   SITS/metadynamics/protocol-sidecar protocol state.
2. For `dynamic` reader payloads:
   - read NHC state. Done.
   - read Bussi thermostat lambda plus typed portable Philox RNG state. Done;
     the old standard-library RNG scalar remains read-only compatibility.
   - read pressure-based barostat `g` plus typed portable Philox RNG state.
     Done; the old standard-library RNG scalar remains read-only
     compatibility.
   - read integrator mode/step/time metadata. Done.
   - read typed portable RNG groups for Bussi/pressure-based barostat/Middle
     Langevin/Andersen counter Philox state and Monte Carlo barostat SplitMix64
     state. Legacy unsupported text markers remain rejected; the older Bussi
     and pressure standard-library scalar states are compatibility-only.
3. For `protocol` reader payloads:
   - read SITS state. Done for float state datasets under the writer path.
   - read metadynamics text state. Done as a transitional compatibility path.
   - read CV/restraint and other file-backed protocol references through
     `/parameters/restart/protocol_sidecars/<command_key>` text datasets.
     Done for whitelisted protocol command keys as a transitional
     compatibility path.
   - read native restraint/CV reference state from
     `/parameters/restart/references/*`. Done for positional-restraint and RMSD
     CV reference coordinates, including named checkpoint write-back.
   - read native metadynamics grid/scatter/edge/hills/runtime state from
     `/parameters/restart/bias/meta/<name>`. Done for canonical state schema v1,
     with read compatibility for XPONGE's earlier typed conversion layout.
   - read module-owned protocol continuation state from
     `/parameters/restart/protocol_state/*`. Planned.
4. For runtime `dynamic`:
   - add module-level apply hooks for engine continuation state. Done for
     Nose-Hoover chain, Bussi thermostat, pressure-based barostat, Middle
     Langevin, Andersen, and Monte Carlo barostat state.
   - check method compatibility before applying.
   - hard-error for legacy Middle Langevin/Andersen/Monte Carlo unsupported
     text markers and for malformed or incompatible typed RNG schemas.
5. For runtime `protocol`:
   - add module-level apply hooks for protocol-owned continuation state. Done
     for full SITS `nk`/`log_norm`/`log_nk`, positional-restraint references,
     RMSD CV references, and native named MetaD state.
   - check protocol compatibility before applying.
   - materialize metadynamics text state to the legacy sidecar files consumed by
     `META::Initial()`. Done for hills/history/edge/potential/direct text
     components as a transitional bridge.
   - materialize generic protocol sidecar text state before CV/protocol module
     initialization and inject the resulting files as controller commands.
     Done for whitelisted keys such as `cv_in_file`, `restrain_atom_id`,
     `restrain_coordinate_in_file`, and `restrain_weight_in_file` as a
     transitional bridge.
   - replace metadynamics/CV/restraint text materialization with native typed
     module apply hooks. Done for CV/restraint definitions and named references,
     and for named MetaD grid/scatter/edge/hill/runtime state. Legacy MetaD text
     materialization remains a compatibility fallback only.
   - hard-error when metadynamics text state is requested but the `meta` module
     is not initialized.
6. For runtime `full`:
   - compose dynamic + protocol for currently supported components. Done for
     NHC, all five portable stochastic module states, integrator metadata,
     full SITS `nk`/`log_norm`/`log_nk`, native/transitional MetaD state,
     references, and compatibility protocol-sidecar text payloads.
   - extend the composition when future dynamic/protocol component apply hooks
     are introduced.
7. For `custom`:
   - hard-error until component list keys exist.
8. Add module-level apply hooks rather than letting readers mutate modules
   directly.
9. Until those hooks exist, validation/runtime must hard-error for unsupported
   `dynamic`, `protocol`, and `full` payloads instead of silently loading only
   structural state. Validation now rejects dynamic requests with malformed or
   legacy unsupported RNG payloads or no supported payload, and protocol
   requests with no supported SITS `nk`, typed/legacy-compatible MetaD,
   reference, or protocol sidecar payload.

Acceptance criteria:

- `load = structural` ignores protocol/dynamic state even if present.
- `load = protocol` loads protocol continuation state and checks protocol
  compatibility. Currently supported natively: full SITS
  `nk`/`log_norm`/`log_nk`, positional-restraint references, RMSD CV
  references, and named MetaD
  grid/scatter/edge/hill/runtime state. Transitional MetaD text and generic
  protocol sidecar materialization remain for older file-backed state. New
  native support must use typed
  `/parameters/restart/references`, `/parameters/restart/bias`, and
  `/parameters/restart/protocol_state` datasets.
- `load = dynamic` loads supported engine continuation state and checks method
  compatibility. Runtime support includes Nose-Hoover chain state, Bussi
  thermostat RNG/lambda, pressure-based barostat `g`/RNG, integrator
  mode/step/time metadata, Middle Langevin/Andersen counter Philox state, and
  Monte Carlo barostat SplitMix64/adaptive state. All five portable stochastic
  components have independent E4 production continuation gates and focused
  process-level checkpoint/resume closure, and are registered as `supported`.
- `load = full` composes supported dynamic and protocol components. Future
  serialized module states must add their own reader/apply hooks and then join
  the same validation path.

Validation gate:

```text
restart reader writer/reader round-trip for NHC/SITS/metadyn payloads
restart reader writer/reader round-trip for protocol sidecar text payloads
restart reader writer/reader round-trip for native CV/restraint/meta state payloads
module-level writer/reader/apply round-trip tests for each supported dynamic component
negative tests for requested component absent from restart.spgr.h5
negative tests for requested component present but unsupported by current runtime
runtime or validation test that dynamic/protocol/full never silently degrade to structural
```

Done boundary:

The current reader/apply subset is complete when readers return typed component
payloads for writer-produced supported state, validation rejects unsupported
payload-only requests, and module-specific apply hooks consume the supported
payloads. Direct module mutation from the H5 reader is not acceptable. Future
Structured future CV/restraint extension state remains a separate module phase.
Portable stochastic state is serialized and applied by the runtime, and its
independent production A/B continuation evidence is bound to the contract
registry.

### Phase 7 Detailed Plan: End-to-end launch integration

Objective:

Wire H5 input into normal SPONGE launch paths.

Primary files:

- `SPONGE/main.cpp`
- `SPONGE/MD_core/MD_core.cpp`
- `SPONGE/MD_core/rerun.hpp`
- `SPONGE/utils/h5md/input_plan.hpp`
- `SPONGE/utils/h5md/input_assembler.hpp`

Implementation steps:

1. Resolve input plan in `main.cpp` after controller initialization.
2. Read metadata and fail early before GPU allocation. Done through
   `input_validation.hpp`.
3. Assemble structural state before MD coordinate allocation.
   Done for H5 restart structural launch; Xponge native coordinate loading now
   skips text coordinate requirements when a valid H5 restart or H5 rerun
   trajectory binding owns the structural source.
4. Initialize rerun provider if `mode = rerun`.
   Done for H5MD trajectory input, including runtime launch without legacy
   `crd`/`box`/`vel`.
5. Ensure output H5 writers can consume state from H5 input path.
6. Add small end-to-end cases:
   - initial restart H5 -> NPT -> output H5 restart
   - H5MD trajectory -> rerun -> observable H5 output

Acceptance criteria:

- Legacy examples still run.
- H5 structural NPT starts without Xponge coordinate input.
- H5MD rerun consumes trajectory input without `crd`/`box`/`vel`.
- Errors mention the binding key and HDF5 path that failed validation.

Validation gate:

```text
legacy startup smoke
H5 structural NPT smoke
H5 structural restart smoke without legacy coordinate input
H5MD rerun smoke
restart output from H5 launch can be read back by RestartH5Reader
trajectory output from H5 launch can be read back by TrajectoryH5Reader
```

Done boundary:

This phase is complete when the canonical test fixtures are executable launch
decks, not only schema examples.

### Phase 8 Detailed Plan: Migration, docs, and downstream sync

Objective:

Make the interface usable by users and downstream tools.

Primary files:

- `docs/input-reference/h5-input-migration.md`
- `docs/input-reference/io.md`
- `docs/input-reference/barostat.md`
- `docs/sponge_input_bundle_inventory.md`
- `scripts/check_mdin_schema_sync.py`

Implementation steps:

1. Update input reference docs with:
   - `[input.h5.*]`
   - restart load policies
   - H5MD rerun input
   - parser aliases
2. Keep end-to-end H5 launch coverage in `tests/h5_bundle/`.
3. Add a Mokda sync note:
   - source: `SPONGE/schemas/`
   - destination: Mokda schema locations
   - validation command
4. Add migration table:
   - `coordinate_in_file` -> initial `restart.spgr.h5`
   - `velocity_in_file` -> initial `restart.spgr.h5`
   - rerun `crd`/`box`/`vel` -> `trajectory.spg.h5md`
   - protocol continuation text files -> `restart.spgr.h5` protocol state

Acceptance criteria:

- A user can write a new H5 NPT mdin from docs alone.
- A user can write a new H5 rerun mdin from docs alone.
- Mokda can sync schemas without editing schema structure by hand.

Validation gate:

```text
python3 scripts/validate_mdin_schema_examples.py --schema schemas/mdin.schema.json
python3 scripts/check_mdin_schema_sync.py ../Mokda
manual TOML parse check for every fenced toml block in touched docs
```

Done boundary:

This phase is complete when the docs, schema, examples, and downstream sync
script describe the same interface names and defaults.

## 8. Validation Matrix

| Area | Test type | Required by phase |
|---|---|---|
| TOML alias normalization | Direct parser unit | Phase 0 |
| mdin schema examples | JSON Schema validation | Phase 0 |
| input resolver | Unit, no HDF5 | Phase 1 |
| input metadata compatibility helpers | Unit, no HDF5 | Phase 1/5 |
| restart reader | HDF5 writer/reader round-trip | Phase 2 |
| structural runtime assembly | CPU/GPU smoke | Phase 3 |
| H5MD trajectory reader | HDF5/VDS round-trip | Phase 4 |
| rerun H5 input | End-to-end rerun smoke | Phase 4/7 |
| topology/protocol compatibility | Unit metadata fixtures | Phase 5 |
| dynamic/protocol restart state | Module-level round-trip | Phase 6 |
| legacy compatibility | Existing examples/smoke | Every phase |

## 9. Open Decisions

1. Whether missing velocity in structural restart should zero-fill or hard-error
   for each mode.
2. Whether initial H5 topology/protocol support should materialize full
   topology/Hamiltonian data or only validate metadata while legacy/Xponge data
   still supplies parameters.
3. How nonzero restart step should interact with current `sys.steps` and output
   cadence.
4. Whether protocol state paths should be standardized now or introduced module
   by module as readers are implemented.
5. Whether a generic HDF5 read backend is needed before more readers are added.
