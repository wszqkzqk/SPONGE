# Bundled I/O A/B Behavior Alignment Test Plan

## 1. Goal and Boundary

The production A/B gate must prove that every declared supported bundled I/O
contract has all of the following:

1. A legacy reference surface.
2. A bundled input or output surface.
3. Runtime evidence that the contract was exercised.
4. A deterministic or statistical equivalence assertion.
5. A traceable result in the gate evidence report.

A passing gate proves equivalence only for contracts explicitly registered as
supported. Cross-process VDS reopen-and-append remains unsupported until the
writer implements an append/resume open mode. Unsupported and deferred
contracts must be recorded explicitly and must not contribute to the supported
coverage percentage.

## 2. Mandatory PR and Commit Rule

Implementation is divided into the independent PR scopes in this document.
The following rule is mandatory for all work under this plan:

> Complete exactly one independent PR scope, verify its acceptance criteria,
> and then create exactly one corresponding commit before starting the next PR.

Additional constraints:

- Do not mix changes belonging to different PR scopes in one commit.
- Do not mark a PR complete or commit it before all acceptance checks for that
  PR pass.
- The commit must contain only files owned by that PR and must preserve
  unrelated worktree changes.
- If a PR must be split because its scope is no longer independently
  reviewable, update this plan first and give each new scope its own commit.
- Follow-up fixes discovered after a completed PR belong to a new, explicitly
  named PR scope and a new commit; do not amend an already verified boundary
  unless the user explicitly requests it.
- Record `This commit` and the verification commands in the PR completion log
  within that PR's commit. Report the resulting hash after the commit is
  created; a Git commit cannot contain its own hash.

Recommended commit subjects:

```text
test(bundled-io): add executable A/B contract registry
test(bundled-io): complete contract manifests and fixtures
test(bundled-io): harden statistical equivalence checks
test(bundled-io): close bundled output behavior coverage
test(bundled-io): close bundled input behavior coverage
test(bundled-io): cover rerun and failure semantics
ci(bundled-io): add executable shadow promotion guards
test(bundled-io): complete and promote the A/B gate matrix
```

## 3. Evidence Model

Create `benchmarks/bundled_io/contracts/ab_contracts.json` as the single source
of truth. Every contract entry must contain:

```json
{
  "contract_id": "trajectory.output.velocity",
  "direction": "output",
  "component": "trajectory",
  "status": "supported",
  "minimum_evidence": "E3",
  "legacy_surface": "vel",
  "bundled_surface": "/particles/all/velocity/value",
  "case_ids": ["normal_particle_payload_vds_off"],
  "assertion_ids": ["trajectory_velocity_statistical_equivalence"]
}
```

Evidence levels are normative:

| Level | Required evidence |
|---|---|
| E0 | Schema, key, or dataset path exists. |
| E1 | Converted payload is semantically equal to its legacy source. |
| E2 | Runtime proves that the owning module consumed or produced the payload. |
| E3 | Legacy and bundled runtime behavior is deterministically or statistically equivalent. |
| E4 | Output is reloaded and continued execution remains equivalent. |
| F1 | Invalid input, conflicts, exit status, and stable diagnostics are equivalent. |

Every runtime assertion must emit an `EvidenceRecord`. Case declarations do
not count as evidence. The gate must write `ab_evidence.json` containing case,
contract, assertion, profile, seed, backend, thread/rank configuration,
comparison policy, result, and artifact paths.

## 4. PR 1: Executable Contract Registry

### Scope

- Add the canonical JSON registry and a typed Python loader.
- Replace free-form `AbCase.covers` labels with registered `contract_ids` and
  `assertion_ids`.
- Make assertions return executed evidence records.
- Replace source-string searches in the production-gate manifest with registry
  and case-model validation.
- Generate a supported/deferred/unsupported coverage report.

### Required guards

- Contract IDs, case IDs, and assertion IDs are unique.
- Every supported contract has at least one case and assertion.
- Every referenced contract exists.
- Every executed assertion reports evidence for a declared contract.
- Evidence below `minimum_evidence` fails the gate.
- E0/E1 path checks cannot satisfy an E3 behavior contract.

### Acceptance

- Removing a case, assertion, or evidence report makes the manifest test fail.
- Adding an unmapped supported contract makes the manifest test fail.
- Existing A/B cases still run and produce the new evidence report.

## 5. PR 2: Contract Enumeration and Fixtures

### Enumeration to add

- Topology: improper, nb14, and the `virtual_atoms_in_file` alias.
- H5 output controls: trajectory chunk size and repair policy.
- Legacy output: mdinfo, force trajectory, restart, and QC SCF output.
- Trajectory input particle-stream selection.
- Rerun start, strip, frame limit, and box-update behavior.
- Restart load policies: structural, dynamic, protocol, and full.
- Output families: trajectory-only, observable-only, restart-only, and all.
- Sidecar default suppression, explicit coexistence, provenance, same-path
  idempotence, and different-path conflict.

### Fixture design

- Split the broad synthetic fixture into targeted fixtures whose owning module
  can be proven active.
- Give every targeted module a non-zero, distinguishable sentinel payload.
- Require an owning mdout column, module stream, state transition, or explicit
  runtime activation record for E2.
- Retain the broad fixture for schema and conversion coverage only; do not use
  path presence from that fixture as E3 evidence.

### Acceptance

- Every registered E1 input has a legacy-to-H5 semantic comparison.
- Every supported module fixture reaches E2 and has a non-empty, non-trivial
  runtime result.
- Pure bundled fixtures contain no undeclared legacy dependencies.

## 6. PR 3: Statistical and Deterministic Comparators

### Deterministic policy

- Compare all rerun and short deterministic NVE rows and datasets.
- Require exact step, time, frame count, schema, and non-finite pattern.
- Use contract-specific numeric tolerances; do not use one global tolerance.

### Statistical policy

- Preserve independent replica and block-mean analysis for stochastic MD.
- Use equivalence confidence intervals for means and declared bounds for
  variance ratios.
- Apply Holm correction across observable families.
- Define margins per physical quantity and unit.
- Reject mismatched NaN, positive infinity, and negative infinity patterns.

### Structured trajectory statistics

- Position: center of mass, per-atom MSD, bond/distance distributions, and RDF
  where applicable.
- Velocity: component means/variances, kinetic temperature, and norm
  distribution.
- Force: net force, RMS, quantiles, and per-atom norm distribution.
- Box: matrix components, volume, lengths, and angles.

### Comparator mutation tests

Inject an energy offset, variance increase, atom permutation, local force
offset, missing non-finite value, and frame-schedule change. Every injected
defect must be rejected reliably.

## 7. PR 4: Output Behavior Closure

Add these independent case groups:

| Case group | Required behavior |
|---|---|
| `normal_particle_payload_vds_off/on` | Position, velocity, force, box, step, and time. |
| `normal_output_routing` | All eight legacy output routes, suppression, explicit coexistence, and provenance. |
| `normal_restart_roundtrip_structural` | Legacy and H5 restart reload and E4 continuation. |
| `normal_restart_roundtrip_dynamic` | RNG, thermostat, barostat, and NHC state continuation. |
| `module_output_nhc/sits/meta/qc/reaxff` | Generic observable and dedicated module streams agree. |
| `output_family_modes` | Trajectory-only, observable-only, restart-only, and all enabled. |
| `vds_chunk_boundaries` | Frame counts at chunk-1, chunk, chunk+1, and 2*chunk+1. |
| `vds_repair` | Strict, complete-prefix, failure status, and repair metadata. |

Compare mdinfo as structured key/value content. Normalize line endings before
exactly comparing QC SCF and metadynamics diagnostic text.

## 8. PR 5: Input Behavior Closure

Create focused A/B coverage for:

- Basic topology: mass, charge, residue, exclusions, LJ, and nb14.
- Bonded force fields: bond, angle, dihedral, improper, Urey-Bradley, and CMAP.
- Specialized topology: GB, softcore, subsystem division, and all virtual-atom
  types including periodic-boundary cases.
- Many-body force fields: EAM, SW, EDIP, Tersoff, and ReaxFF.
- Custom pairwise and listed forces.
- Protocol behavior: constraints, CV, restraint, steering, and soft wall.
- Bias/state behavior: SITS, metadynamics, and NHC.
- QC type, energy, spin-square, and SCF text.

Each case must prove that the legacy and bundled branches activate the same
module and compare a module-owned result. Initialization messages or H5 path
existence alone are insufficient.

## 9. PR 6: Rerun and Failure Semantics

### Rerun matrix

- Start/strip: 0/0, 1/0, 0/1, and 1/1.
- Frame limit: 1, exact EOF, beyond EOF, and default unlimited.
- Box update: false and true.
- Input trajectory: velocity present and absent.
- Particle stream: default `all` and one non-default stream.
- Restart: absent, structural, dynamic, protocol, and full.
- Output: VDS off and on.

Use a risk-driven pairwise matrix plus explicit boundary cases instead of the
complete Cartesian product.

### F1 failure matrix

- Atom count, shape, dtype, and schema-version mismatch.
- Missing topology, protocol, restart, or trajectory binding.
- Mixed legacy/H5 restart or trajectory inputs.
- Unsupported sidecar key, key/path length mismatch, and path conflict.
- Invalid chunk size, VDS value, repair policy, and restart policy.
- Dynamic/protocol/full state requested without its owning runtime module.

Compare exit code, error category, and stable diagnostic tokens. Full stderr
byte equality is not required.

## 10. PR 7: Execution Matrix and CI Promotion

PR 7 is split into three independently reviewable commits because the promotion
mechanism can be verified while the runtime matrix is still red. PR 7a must not
claim behavior coverage or enable a release gate; PR 7b owns portable runtime
cases and CPU/MPI evidence; PR 7c owns GPU evidence and eventual promotion.

### PR 7a: Executable Shadow Promotion Guard

- Add a machine-readable matrix containing every required axis and risk-driven
  combination below.
- Validate the matrix structurally and reject missing axes, combinations,
  scenario evidence, or non-executable case references.
- Add a promotion evaluator for supported-contract evidence, scenario status,
  consecutive retry-free runs, and performance budgets.
- Run contract, medium, and production tiers in an explicitly non-release,
  shadow CI workflow. Known runtime failures must remain visible.

Acceptance:

- Mutation tests reject a removed axis, removed risk combination, unresolved
  scenario, missing evidence, fewer than three retry-free runs, and each
  performance-budget violation.
- The current matrix evaluates as not ready while deferred scenarios or
  behavior gaps remain.
- Release workflows do not consume the shadow result.

### PR 7b: Portable Runtime Matrix and CPU/MPI Evidence

- Add pinned same-semantic TIP3P legacy/bundled fixtures for orthogonal and
  nonorthogonal boxes without requiring an online converter during the gate.
- Define executable cases for every declared CPU/GPU, OMP, and MPI scenario.
- Run and compare all CPU rank-1 and rank-2 scenarios, including a single
  rank-0 H5 finalize report and output ownership under MPI.
- Make deterministic position comparisons and statistical position features
  periodic-box aware. Bound pair-distribution work with deterministic sampling.
- Keep GPU scenarios executable but unproven when no GPU runner is available.

Acceptance:

- CPU rank-1 covers deterministic NVE; middle, Bussi, and NHC NVT; orthogonal
  and nonorthogonal NPT; all three constraint modes; and OMP 1/4.
- CPU rank-2 covers deterministic NVE and nonorthogonal statistical NPT, with
  rank-0 output ownership asserted.
- Periodic-image-only position changes pass, while material force/position
  mutations remain rejected.
- Contract, comparator, matrix manifest, and workflow validation tests pass.

### PR 7c: GPU Evidence and Final Promotion

Add a risk-driven execution matrix:

- Deterministic NVE.
- Statistical middle-Langevin, Bussi, and NHC NVT.
- NPT/barostat with orthogonal and non-orthogonal boxes.
- SHAKE, SETTLE, and unconstrained systems.
- CPU and GPU, comparing legacy versus bundled within the same backend.
- `OMP_NUM_THREADS=1` and `4`.
- MPI rank counts 1 and 2, including rank-0 output ownership.

Gate tiers:

| Tier | Contents | Target |
|---|---|---|
| contract | Registry, schema, conversion, and negative contract tests. | Under 5 minutes. |
| medium | Deterministic matrix and shortened statistical A/B. | 30-45 minutes. |
| production | Full replicas, continuation, module, and environment matrix. | Nightly and release gate. |
| soak | Repeated production runs for flake and resource analysis. | Scheduled. |

Production promotion requires:

- 100% of supported contracts at their minimum evidence level.
- All runtime input/output contracts at E3.
- Restart and dynamic state contracts at E4.
- All declared invalid combinations at F1.
- Three consecutive production runs without retry.
- All comparator mutation tests reject their injected defects.
- Every PR 7a matrix scenario is backed by executable case IDs whose evidence
  metadata proves the declared backend, thread count, MPI rank count, and
  rank-0 output ownership; a declaration-only scenario cannot be promoted.

## 11. PR 8: VDS Chunk-Boundary Behavior Closure

This follow-up scope closes the deferred runtime behavior contract for
`output_h5_trajectory_chunk_size` without weakening the unsupported
cross-process append boundary.

- Run same-semantic legacy-input and bundled-input deterministic NVE branches
  with a fixed chunk size.
- Cover frame counts at chunk-1, chunk, chunk+1, and 2*chunk+1.
- Compare mdout plus every trajectory/observable H5 dataset across branches.
- Assert the VDS wrapper frame count and physical shard count at each boundary.
- Emit E3 evidence for the chunk-size contract from the executed cases.

Acceptance:

- The four boundary cases produce exactly 3, 4, 5, and 9 frames for chunk size
  4 and exactly 1, 1, 2, and 3 shards respectively.
- A wrong frame count, shard count, or cross-branch numeric payload fails the
  gate.
- Registry, manifest, comparator, and full contract tests pass.

## 12. PR 9: Focused Non-Zero EDIP Input Behavior Gate

This follow-up replaces the broad fixture's all-zero EDIP sentinel with a
minimal deterministic case whose owning module has distinguishable behavior.

- Generate a two-atom legacy input with `EDIP_in_file`, finite separation, and
  no competing force-field module.
- Convert the same input to native `/manybody/edip` H5 datasets and remove the
  converter's unused legacy-sidecar directory before execution.
- Require non-zero, cross-branch equivalent EDIP energy and force values.
- Compare complete deterministic mdout rows and the non-zero legacy force route
  across the legacy-input and bundled-input branches.
- Keep H5 force/output schedule behavior outside this input-only evidence scope;
  that independent contract must be repaired and gated separately.

Acceptance:

- The bundled mdin has no `EDIP_in_file` and the runtime directory has no
  `legacy_sidecars`; all required EDIP pair/triple/type datasets are present.
- Both branches emit finite `EDIP` values above the declared sentinel floor and
  non-zero force payloads.
- Zero energy, zero force, force mismatch, or mdout mismatch fails the gate.
- Registry, manifest, comparator, real focused A/B, Ruff, and diff checks pass.

## 13. PR 10: H5 Sidecar-Table Failure Semantics

Close the deferred F1 contract for malformed bundled sidecar tables with real
process execution rather than reader-only checks.

- Mutate the topology H5 sidecar table with an unsupported command key.
- Create unequal key/path dataset lengths.
- Duplicate an allowed key with a different path to trigger command conflict.
- Require non-zero exit, `spongeErrorValueErrorCommand`, and stable diagnostic
  tokens for every mutation.
- Keep these cases bundled-only because legacy text inputs have no equivalent
  H5 key/path table representation.

Acceptance:

- All three mutations reach the production sidecar materialization path and are
  rejected before MD execution.
- Unsupported-key evidence names the topology binding and invalid key;
  length-mismatch evidence names the key/path cardinality error; conflict
  evidence names the key and both existing/H5 paths.
- Registry, manifest, all three real CPU failure cases, Ruff, and diff checks
  pass.

## 14. PR 11: Restart Owner-State Failure Semantics

Close the deferred F1 owner-state contract by exercising each non-structural
restart load policy through the production process.

- Request `dynamic` state containing Nose-Hoover-chain payload while the
  current run has no Nose-Hoover-chain thermostat owner.
- Request `protocol` state containing metadynamics payload while the current
  run has no initialized metadynamics owner.
- Request `full` state with the same missing dynamic owner and require it to
  reject before partially applying later protocol state.
- Require non-zero exit, `spongeErrorConflictingCommand`, and owner-specific
  stable diagnostic tokens for all three cases.

Acceptance:

- The three cases use the same bundled restart payload and differ only in the
  declared `input_h5_restart_load` policy.
- Dynamic/full diagnostics name the Nose-Hoover-chain state and missing
  thermostat; protocol diagnostics name metadynamics state and its missing
  module.
- Registry, manifest, all three real CPU failure cases, Ruff, and diff checks
  pass.

## 15. PR 12: Focused Non-Zero Custom-Pair Input Gate

Close the custom-pair input contract with a minimal deterministic fixture whose
module-owned energy and force cannot pass trivially.

- Generate a two-atom legacy input with a named custom pair potential and a
  non-zero parameter payload.
- Convert it to `/forcefield/custom_force/pairwise`, then remove both legacy
  custom-force mdin keys, the H5 sidecar table, and the sidecar directory from
  the bundled branch.
- Require non-zero, cross-branch equivalent `custom_pair` energy and force.
- Require both native materialized inputs to be present and non-empty.

Acceptance:

- The bundled branch has no `pairwise_force_in_file` or
  `custom_pair_in_file`, yet materializes both inputs from native H5 data.
- Complete deterministic mdout rows match, `custom_pair` is non-zero in both
  branches, and all legacy force-route values match within the force tolerance.
- Registry, manifest, real focused CPU A/B, Ruff, and diff checks pass.

## 16. PR 13: Focused Native Exclusion Semantics Gate

Close the exclusions contract with a three-atom, non-periodic Coulomb fixture
whose excluded-pair behavior has an independent analytic oracle.

- Exclude the close `(0, 1)` pair while retaining two non-zero interactions.
- Convert to `/topology/exclusions/{offset,list}` and remove all bundled
  sidecar routes.
- Require cross-branch deterministic mdout and force equivalence.
- Independently require energy `1/4 - 1/3 = -1/12` and the corresponding nine
  Cartesian force values; the unexcluded energy would instead be `-13/12`.

Acceptance:

- Bundled offsets/list are exactly `[0, 1, 1, 1]` and `[1]`, with no
  `exclude_in_file` or legacy sidecar table.
- Both branches emit non-zero Coulomb behavior, match each other, and satisfy
  the analytic energy/force oracle.
- Energy and force mutations representing ignored or corrupted exclusions are
  rejected by unit tests.
- Registry, manifest, real focused CPU A/B, Ruff, and diff checks pass.

## 17. PR 14: Focused Native LJ Soft-Core Behavior Gate

Close the LJ soft-core input contract with a two-atom periodic fixture whose
A-state interaction is absent and B-state interaction produces distinguishable
soft-core energy and force at `lambda_lj = 0.5`.

- Convert the legacy `LJ_soft_core_in_file` payload to native
  `/forcefield/lj_soft_core`, then remove the legacy key, sidecar table, and
  sidecar directory from the bundled branch.
- Preserve source-semantic pair values in H5 and apply the legacy loader's
  `12/6` A/B normalization exactly once when pure native H5 owns the payload.
- Require non-zero, cross-branch equivalent `LJ_soft` energy and force values,
  plus complete deterministic mdout equivalence.
- Keep `/forcefield/subsys_division` absent so this case proves only the
  soft-core parameter contract. Subsystem division remains deferred because
  the current force kernel copies but does not consume its mask.

Acceptance:

- The bundled branch has no `LJ_soft_core_in_file`, subsystem input, or legacy
  sidecars, and contains the exact native atom-type and pair datasets.
- Both branches emit the same non-zero `LJ_soft` result and six-value force
  payload; zero or mismatched energy/force mutations are rejected.
- Legacy text, H5 fixture-equivalence, and native-reader tests retain their
  source-semantic parameter values while runtime behavior matches.
- Registry, manifest, real focused CPU A/B, Ruff, and diff checks pass.

## 18. PR 15: Focused Native Virtual-Atom Behavior Gates

Close the virtual-atom reconstruction contracts with one all-types fixture and
one boundary-crossing fixture. Both bundled branches must be pure native H5,
with no legacy key or sidecar fallback.

- Enumerate types 0, 1, 2, and 3 in one ragged native payload and require exact
  offsets, parent indices, parameters, reconstructed coordinates, non-zero PM
  behavior, real-atom force signal, and full legacy/bundled force equivalence.
- Make type 3 depend on a type 2 site so the runtime case also exercises ordered
  multi-layer reconstruction.
- Put the type 1 parents on opposite sides of a periodic box and require the
  minimum-image coordinate `9.75`; a non-periodic interpolation would produce
  `7.25` and must fail the oracle.
- Refresh virtual coordinates from global records before domain decomposition,
  then switch to local records after `Get_Local`; copy type 3 device records
  from their populated host source.

Acceptance:

- The bundled topology contains the exact `/forcefield/virtual_atom` ragged
  datasets and retains neither virtual-atom legacy key nor sidecar metadata.
- Both branches match the per-type and PBC coordinate oracles, emit finite
  non-zero PM/force behavior, and match complete force and mdout payloads.
- Coordinate, PBC, and trivial-real-force mutations are rejected by unit tests.
- Registry, manifest, real focused CPU A/B, related H5 CTest, Ruff, and diff
  checks pass.

## 19. PR 16: Virtual-Atom Plural Alias Behavior Gate

Close the `virtual_atoms_in_file` compatibility contract independently from the
canonical singular key.

- Reuse the boundary-crossing type 1 fixture, but expose the source only through
  `virtual_atoms_in_file` in the legacy branch.
- Let the legacy native loader resolve the plural alias when the canonical key
  is absent, without changing custom virtual-atom module key resolution.
- Require conversion to the exact `/forcefield/virtual_atom` ragged payload and
  remove both spellings, the sidecar table, and the sidecar directory from the
  bundled branch.
- Apply the same `9.75` PBC coordinate, non-zero PM, real-atom force, complete
  force, and mdout equivalence assertions as the canonical PBC contract.

Acceptance:

- The legacy branch contains only `virtual_atoms_in_file`; the bundled branch
  contains neither spelling and has no fallback sidecar route.
- Both branches satisfy the same runtime coordinate/energy/force semantics.
- Registry, required-semantic inventory, manifest, real focused CPU A/B,
  related H5 CTest, Ruff, clang-format, and diff checks pass.

## 20. Restart-Absent Same-Bootstrap Rerun Gate

Goal: prove that omitting restart input has the same behavior in legacy and
bundled rerun paths, rather than merely proving that both processes start.

Fixture and route requirements:

- Use byte-identical legacy coordinate and velocity files as bootstrap state
  for both branches.
- Keep legacy text trajectory input and bundled H5 trajectory input as the
  only branch-specific state route.
- Remove `input_h5_restart_path`, `input_h5_restart_load`, and the bundled
  `restart.spgr.h5` file so an accidental structural fallback cannot satisfy
  the test.
- Run two selected frames with box update disabled to isolate restart absence
  from the known legacy box-update invalid-free defect.

Behavior requirements:

- Compare every mdout row and column deterministically.
- Check the selected rerun frame indices against an independent oracle.
- Compare bundled trajectory and observable output against the legacy branch,
  including position, box, velocity, force, and scalar observables.
- When no restart exists, derive the expected fixed box from the shared
  orthogonal bootstrap coordinate file; existing structural cases continue to
  use the restart H5 box.

Acceptance:

- The route assertions prove identical bootstrap bytes and complete absence of
  restart bindings and payloads.
- All three E3 assertions pass in a real CPU run.
- Registry, manifest, restart/input CTest, Ruff, smoke, and diff checks pass.

## 21. PR 18: NHC Dynamic Restart Continuation Gate

Close the dynamic restart input, NHC-owned state, and dynamic restart output
contracts with one deterministic producer checkpoint rather than two
independently evolved branches.

- Run one NVT Nose-Hoover-chain producer with chain length three and write the
  particle structural state plus NHC state through both legacy text and H5
  restart routes at the same restart event.
- Require the producer NHC state to be finite and non-zero, and compare the
  legacy and H5 structural/NHC payloads before either continuation starts.
- Fork that single checkpoint into a legacy continuation that binds coordinate,
  velocity, and NHC text restart files and a bundled continuation that binds
  H5 topology/protocol plus `input_h5_restart_load = "dynamic"`.
- Continue both deterministic NHC runs for two steps with identical integrator,
  thermostat, topology, and output settings.
- Compare every mdout column, particle position/velocity/force/box streams, NHC
  coordinate/velocity trajectories, observable-only output, final structural
  restart, final NHC restart, and integrator mode/step/time state.
- Compare legacy NHC trajectory/restart text against the H5 module datasets in
  each branch; account only for the legacy `%f` representation boundary and do
  not relax the existing position or velocity tolerances.

Acceptance:

- The producer state and continued NHC state are non-trivial.
- The bundled mdin contains no legacy coordinate, velocity, or NHC restart
  binding; the legacy mdin contains no H5 binding.
- Any NHC component swap, particle-state mismatch, mdout mismatch, schedule
  mismatch, missing output family, or output-to-input continuation mismatch
  fails the gate.
- The three contracts emit E4 evidence from the same executed case, and
  registry, manifest, real CPU A/B, related CTest, Ruff, and diff checks pass.

## 22. PR 19: Protocol/Full Metadynamics Continuation Gate

Close the protocol and full restart-load policies together with the
metadynamics input contract, using one deterministic checkpoint and three
same-semantic continuations.

- Generate a two-atom periodic NVT producer with a distance collective
  variable, one-dimensional well-tempered metadynamics, and a three-element
  Nose-Hoover chain. Require multiple finite, non-zero hills and a finite
  potential snapshot in the H5 restart.
- Treat the H5 restart event as the checkpoint boundary. Project its stored
  hill and potential text into the legacy branch instead of copying the
  producer's terminal sidecar, because the terminal sidecar contains the next
  hill written after that restart event.
- Fork the checkpoint into: a pure legacy structural/NHC/metadynamics text
  continuation; a `protocol` continuation that loads H5 structural and
  metadynamics state while retaining only the legacy NHC dynamic state; and a
  `full` continuation that loads structural, NHC dynamic, and metadynamics
  state entirely from H5.
- Compare every mdout column; particle position, velocity, force, and box;
  NHC coordinate/velocity trajectories and final state; metadynamics
  `meta/rbias/rct` streams; final hill histories; observable-only output; and
  final restart/integrator state.
- Apply existing deterministic numeric tolerances only at legacy text
  precision boundaries. Reject missing/non-finite/trivial bias, a different
  hill count, a route that mixes owners, or any schedule/shape mismatch.

Acceptance:

- All three continuations use the same producer checkpoint and yield finite,
  non-zero metadynamics work with equivalent particle, thermostat, bias, and
  restart behavior.
- The `protocol` branch proves that only protocol state is added beyond the H5
  structural baseline; the `full` branch contains no legacy restart binding.
- `input.restart_load.protocol`, `input.restart_load.full`, and
  `input.bias.metadynamics` emit E4 evidence from the executed case.
- Registry, manifest, real CPU A/B, related CTest, Ruff, and diff checks pass,
  followed by exactly one PR-scoped commit.

## 23. PR 20: Unrestricted QC Observable and SCF Text Gate

Close the non-trivial QC spin-square and non-empty SCF text contracts with a
deterministic unrestricted-QC rerun under both VDS modes.

- Start from one two-atom `H/N` full-contract fixture and change both branches
  from singlet/restricted to triplet/unrestricted QC. Keep atom indices,
  symbols, charge, coordinates, and every non-QC input unchanged.
- Keep the bundled sidecar QC type payload and its typed `/qc/type` metadata
  synchronized at multiplicity three. This case proves output behavior only;
  it does not promote `input.qc.type`, because the pure bundled runtime still
  does not materialize `/qc/type` into the QC input consumer.
- Require a present, finite, non-zero `QC_S_sq` value in every branch and exact
  deterministic legacy/bundled mdout equivalence. Bridge the same values to
  `/observables/all/qc/spin_square/{step,time,value}`.
- Enable `qc_scf_print_iter`, require non-empty legacy and bundled
  `qc_scf_output`, normalize line endings, and compare the complete text.
  Require both bundled trajectory and observable files to archive exactly the
  same text at `/parameters/sponge/qc/scf_output`.
- Run the focused case with VDS disabled and enabled so both trajectory writer
  implementations carry the QC parameter text while observable behavior stays
  invariant.
- Remove `qc_scf_exact_equivalence` evidence from the old full-contract cases,
  whose SCF files are empty, so successful process startup can no longer
  satisfy this contract.

Acceptance:

- Both VDS modes produce the same two finite, non-zero spin-square frames in
  legacy mdout and bundled observable output.
- SCF output is non-empty and exactly equal across legacy text, bundled text,
  trajectory metadata, and observable metadata.
- Zero spin-square, missing observable mapping, empty SCF text, text mismatch,
  or missing H5 archive fails the gate.
- `input.qc.spin_square`, `input.qc.scf_text`, and
  `output.legacy.qc_scf_output` emit E3 evidence only from the focused cases;
  `input.qc.type` remains explicitly deferred.
- Registry, manifest, both real CPU A/B cases, Ruff, and diff checks pass,
  followed by exactly one PR-scoped commit.

## 24. PR 21: Constraint Sidecar Projection Gate

Close the executable constraint sidecar contract with a deterministic,
non-trivial two-atom projection while keeping native typed `/constraint`
consumption explicitly separate.

- Generate two equal-mass atoms at the target distance `1.5` with opposing
  radial velocities `-1.0` and `1.0`, so merely parsing or initializing the
  constraint module cannot satisfy the gate.
- Convert the same legacy input into a bundle. Require the converter's typed
  `/constraint/default/pairs/{atoms,r0}` payload and the exact
  `constrain_in_file` sidecar binding to coexist, while the bundled mdin does
  not retain a direct legacy key.
- Run four deterministic NVE steps with SHAKE enabled. For every emitted
  frame, require the pair-distance residual to stay below `1e-5` and the
  relative radial-velocity residual to stay below `1e-4` in both branches.
- Compare the complete legacy/bundled position and velocity trajectories with
  the existing deterministic tolerances, in addition to full mdout equality.
- Record this evidence for `input.protocol.constraint.sidecar` only. Keep
  `input.protocol.constraint` deferred until the runtime consumes typed
  `/constraint` without the sidecar table.
- Mutation tests must reject both a wrong constrained distance and an
  unprojected radial velocity, preventing module-startup evidence from passing.

Acceptance:

- Both branches reduce the initial relative radial speed from `2.0` to at most
  `1e-4` and preserve the target distance in every one of four frames.
- The sidecar payload is byte-identical to the legacy constraint file and is
  bound only through the protocol H5 sidecar table in the bundled branch.
- Removing constraint execution, changing the target geometry, retaining
  radial motion, changing the schedule, or diverging either trajectory fails
  the gate.
- `input.protocol.constraint.sidecar` emits E3 evidence from the focused case;
  native typed `input.protocol.constraint` remains explicitly deferred.
- Registry, manifest, real CPU A/B, Ruff, and diff checks pass, followed by
  exactly one PR-scoped commit.

## 25. PR 22: Native Improper Runtime Gate

At the PR 22 boundary, separate SPONGE's executable native improper behavior
from the then-broken legacy-to-bundle integration, and close only the former
with a non-trivial four-atom case. PR 43 below closes the retained conversion
gap.

- Generate four non-coplanar equal-mass atoms and one harmonic improper with
  `pk=10.0` and `phi0=0.2`. Require the legacy canonical
  `improper_dihedral_in_file` route to produce finite, non-zero
  `improper_dihedral` energy and force.
- Convert the input, then explicitly normalize the converter's legacy-schema
  `/forcefield/improper/k` dataset to the runtime reader's
  `/forcefield/improper/pk`. Remove the complete H5 sidecar table and sidecar
  directory so the bundled result can only come from typed native state.
- Require the bundled mdin to contain neither `improper_dihedral_in_file` nor
  the `improper_in_file` alias. Validate the exact one-row atom, `pk`, and
  `phi0` payload before execution.
- Compare the module-owned energy row, complete mdout row, and complete binary
  force payload deterministically. Reject trivial energy, trivial force, any
  shape change, or any cross-branch value mismatch.
- Track the passing behavior as `input.topology.improper.native_runtime`. At
  this boundary, keep `input.topology.improper` deferred for end-to-end
  conversion because the converter writes `/k` while the runtime requires
  `/pk`.

Acceptance:

- Legacy canonical text and pure native typed H5 produce the same finite,
  non-zero improper energy and force without bundled sidecars.
- Removing or ignoring the typed payload makes the semantic gate fail even if
  the process and improper module initialization complete.
- The registry and documentation do not use the native runtime evidence to
  claim that legacy-to-bundle improper conversion is closed.
- Registry, manifest, real CPU A/B, Ruff, and diff checks pass, followed by
  exactly one PR-scoped commit.

## 26. PR 23: GB Native-State/Sidecar-Activation Gate

Close the executable GB hybrid route without using its evidence to claim that
pure native `/forcefield/gb` initialization is supported.

- Generate two non-periodic atoms with charges `+1/-1`, separation `2.0`, and
  identical GB radius/scale parameters `1.5/0.8`. Require finite, non-zero
  module-owned `gb` energy.
- Convert the same input and retain typed `/forcefield/gb/params`. Reduce the
  H5 sidecar table to the single `gb_in_file` binding so mass and charge are
  necessarily consumed from typed native state.
- Treat the remaining GB sidecar as an activation binding, not GB payload
  evidence: native topology state is loaded before sidecar injection, and the
  later binding causes `main` to call `gb.Initial` against that native state.
- Require the bundled mdin to contain no direct `gb_in_file`, while the bound
  sidecar remains byte-identical to the legacy input and the typed parameter
  matrix is validated independently.
- Compare complete mdout and force payloads. In addition to non-zero GB energy,
  require the total force oracle `+/-0.10313021`; reject the Coulomb-only
  `+/-0.25` result so process startup or an ignored GB force cannot pass.
- Track this route as `input.topology.gb.hybrid_activation`. Keep pure native
  `input.topology.gb` deferred until initialization no longer depends on the
  legacy-key sidecar.

Acceptance:

- Legacy text and the typed-state/sidecar-activation bundle produce equivalent
  non-zero GB energy and GB-influenced force.
- Missing or trivial GB energy, Coulomb-only force, sidecar-table contamination,
  typed parameter drift, or any cross-branch mdout/force mismatch fails.
- The registry explicitly distinguishes the executable hybrid route from the
  deferred pure native contract.
- Registry, manifest, real CPU A/B, Ruff, and diff checks pass, followed by
  exactly one PR-scoped commit.

## 27. PR 29: VDS Complete-Prefix Repair Behavior Gate

Close the VDS repair-policy contract with separate evidence for successful
production finalization and terminal-shard failure. A normal no-op is required
behavior, but is not by itself evidence that repair works.

- Run a deterministic five-frame, chunk-size-four SPONGE case from both the
  legacy-input and bundled-input routes with
  `output_h5_trajectory_repair_policy=complete_prefix`.
- Compare every mdout column and every numeric/string H5 dataset across the two
  routes. Require two complete shards, five materialized frames,
  `repair_status=not_applied`, and `repaired_shard_count=0` in both wrappers.
- In the same production gate, execute the real HighFive terminal/resume smoke
  helper with a backend that fails finalization of the terminal shard. Require
  `Finalize_With_Repair` to retain only the valid prefix, rewrite completion
  metadata to one frame, report `repair_status=applied`, and count one repaired
  shard.
- Mutation-test the production no-op oracle against strict policy metadata,
  false applied status, non-zero repaired count, wrong completion frame count,
  and an open manifest entry.
- Keep `output.vds.cross_process_append_resume` unsupported. Complete-prefix
  repair finalizes the current process output; it does not reopen an existing
  wrapper for append.

Acceptance:

- Legacy-input and bundled-input production runs retain identical complete
  five-frame semantics and do not report a repair when every shard finalizes.
- The real backend failure helper proves that a failed terminal shard is
  removed and the wrapper exposes only the valid complete prefix.
- Either half of the evidence being absent, any metadata mutation, an
  incomplete retained manifest, or any cross-branch H5 difference fails the
  gate.
- Registry, manifest, real CPU A/B, real-backend VDS smoke, Ruff, and diff
  checks pass, followed by exactly one PR-scoped commit.

## 28. PR 31: Residue Membership/COM-Virial Follow-up Gate

Strengthen PR 30 from residue-count and whole-molecule evidence to a consumer
whose result changes with residue membership while keeping the declared residue
count fixed.

- Retain PR 30's `normal_residue_sidecar_pbc_mapping` case and its residue-count,
  runtime-partition, force, and whole-molecule PBC mapping mutations. Add the
  COM-virial case as a second consumer of the same sidecar contract instead of
  replacing the existing gate.
- Use two disconnected two-atom molecules with the same `[2,2]` residue
  partition in the legacy and bundled routes. Isolate `residue_in_file` in the
  topology H5 sidecar table and retain only the constraint and restraint
  support bindings in protocol H5.
- Enable `restrain_refcoord_scaling="com_res"`, restraint virial calculation,
  and pressure output. Compare complete mdout and force output, and require
  finite non-trivial `bond`, `restrain`, `pressure`, and `Pxx` behavior.
- Execute a bundled control with `[1,3]`: the declared residue count remains
  two, but the disconnected second residue is split by the runtime. Require
  unchanged bond/restraint energy and force together with a distinct,
  deterministic pressure/Pxx fingerprint. This prevents sidecar parsing or
  module initialization alone from satisfying the gate.
- Fix CPU ownership of atom-to-group maps used by `com_ug`, `com_res`, and
  `com_mol`: allocate independent device storage before copying and then free
  the temporary host map. This avoids a dangling CPU alias that is overwritten
  by later COM scratch allocation.
- Keep typed `/atoms/residue_index` deferred; this PR closes only the isolated
  topology-sidecar behavior route.

Acceptance:

- Legacy and bundled `[2,2]` routes produce equivalent finite mdout,
  trajectory, and non-trivial force behavior at the COM-residue consumer.
- The `[1,3]` control keeps two declared residues and the same energy/force but
  changes both pressure and `Pxx`, proving membership-sensitive virial use.
- Mutation tests reject wrong residue count, pressure, `Pxx`, and trivial force
  without searching production-test source text for implementation tokens.
- Registry, manifest, real CPU A/B, Ruff, clang-format, and diff checks pass,
  followed by exactly one PR-scoped commit.

## 29. PR 32: Pure Native GB Behavior Gate

Close the remaining pure typed generalized-Born input contract without using a
legacy activation sidecar.

- Retain PR 23's native-state plus `gb_in_file` sidecar-activation case as a
  separate compatibility route.
- Add a second deterministic two-atom case whose bundled branch contains
  `/forcefield/gb/params` but no `gb_in_file` command, H5 sidecar table, or
  `legacy_sidecars` directory. The legacy branch continues to consume the same
  radii and scale factors through `gb_in_file`.
- Initialize GB in non-periodic mode when either the legacy command exists or
  `Xponge::system.generalized_born` already contains native parameters.
- Require complete mdout equality, `Coulomb=-0.50`, `gb=-0.25`,
  `potential=-0.75`, and the complete six-value force oracle with maximum
  magnitude `0.10313020646572113`. Reject zero GB output and the Coulomb-only
  force fingerprint.

Acceptance:

- Both the pure typed and hybrid GB cases pass real CPU legacy/bundled A/B.
- The pure typed route cannot satisfy the gate through a retained text binding
  or sidecar payload.
- Registry, manifest, related native-reader CTest, full smoke, Ruff,
  clang-format, and diff checks pass, followed by exactly one PR-scoped commit.

## 30. PR 33: Core Topology Payload-Sensitivity Gate

Replace the existing shared statistical evidence for mass, charge, and
Lennard-Jones input with an additional deterministic, payload-specific behavior
gate.

- Use one two-atom NVE fixture with `dt=0`, no PBC, explicit non-zero velocity,
  charge, and LJ interactions. The legacy branch reads `mass_in_file`,
  `charge_in_file`, and `LJ_in_file`; the bundled branch contains only typed
  `/atoms/mass`, `/atoms/charge`, and `/forcefield/lj` datasets, with no topology
  sidecar table or `legacy_sidecars` directory.
- Require complete one-frame mdout and force equivalence between the two routes.
- Run three independent bundled controls which mutate only one typed owner:
  doubling both masses changes temperature from `26.21` to `52.42` while force
  remains unchanged; zeroing charge changes Coulomb energy from `-0.50` to zero
  and changes force; zeroing both LJ pair arrays changes LJ energy from `-0.02`
  to zero and changes force.
- Make non-periodic LJ initialization consume the same materialized
  `Xponge::system.classical_force_field.lj` state as periodic LJ, while retaining
  named-module legacy loading.

Acceptance:

- The real CPU A/B and all three payload controls pass with finite, non-trivial
  module-owned responses.
- Manifest mutation tests reject trivial correct values, unchanged observables,
  missing required force changes, unexpected mass-only force changes, and
  non-finite force payloads.
- The original statistical TIP3P evidence remains active as a broader regression.
- Registry, full manifest and smoke suites, SPONGE rebuild, related native-reader
  CTests, Ruff, clang-format, and diff checks pass, followed by exactly one
  PR-scoped commit.

## 31. PR 34: Complete H5 Metadata Failure Semantics

Close the complete topology metadata F1 contract by adding the missing schema
version rejection to the existing atom-count, shape, and dtype process gates.

- Preserve the three schema versions emitted by current supported producers:
  legacy unit fixtures (`0`), native H5 fixtures (`1`), and XPONGE conversion
  (`xponge.legacy_to_bundle.v1`).
- Add a bundled-only mutation that replaces `/schema/version` with
  `unsupported.topology.v999` and require a non-zero process exit,
  `spongeErrorValueErrorCommand`, and route-specific stable diagnostics.
- Run validation from the common `Xponge::System::Load_Inputs` entry before any
  native text input loader. Legacy-only inputs remain unaffected.
- Merge the partial three-mutation `failure.h5_metadata.runtime_rejections`
  registry record into the complete `failure.h5_metadata` contract so one owner
  enumerates all four failure semantics.

Acceptance:

- All four real bundled failure cases reject their mutation with the exact
  expected error category and stable diagnostic tokens.
- Known supported schema versions continue to pass existing real A/B and H5
  input-validation tests.
- Registry, full manifest and smoke suites, SPONGE rebuild, related CTests,
  Ruff, clang-format, and diff checks pass, followed by exactly one PR-scoped
  commit.

## 32. PR 35: Subsystem-Division Energy-Partition Behavior Gate

Close the typed subsystem-division contract with a partition-sensitive runtime
observable while preserving the established soft-core dynamics.

- Treat `subsys_division` as an ownership label for short-range soft-core LJ
  energy accounting. A pair whose atom masks are equal contributes to
  `LJ_soft_intra`; a pair whose masks differ contributes to `LJ_soft_inter`.
- Keep the existing pair energy, total energy, force, virial, and lambda
  derivative formulas unchanged. Require the explicit conservation invariant
  `LJ_soft_inter + LJ_soft_intra = LJ_soft_short`.
- Use a deterministic two-atom soft-core fixture. The legacy branch consumes
  `subsys_division_in_file`; the bundled branch consumes only typed
  `/forcefield/subsys_division`, with no retained subsystem command, topology
  sidecar table, or `legacy_sidecars` directory.
- Run the same-semantic `[0,1]` A/B pair and require non-zero inter-system
  energy. Run a second same-semantic `[0,0]` control A/B pair and require the
  same energy to move entirely to the intra-system observable.
- Require complete mdout equivalence for both semantic inputs, unchanged total
  soft-core energy and force across the mask mutation, and exact legacy/bundled
  partition values. This proves runtime consumption without assigning an
  unsupported force-filtering meaning to the historical mask.

Acceptance:

- The `[0,1]` routes match at `LJ_soft_inter=-0.06`,
  `LJ_soft_intra=0.00`, and `LJ_soft_short=-0.06`; the `[0,0]` routes match at
  `LJ_soft_inter=0.00`, `LJ_soft_intra=-0.06`, and the same total.
- Both controls retain the complete non-zero six-value force fingerprint, and
  the bundled branch has only the native subsystem dataset.
- Mutation tests reject non-finite, trivial, mixed, non-conserving, and
  total-energy-changing partition evidence.
- Registry, manifest, real CPU A/B, full smoke, related native-reader CTests,
  SPONGE rebuild, Ruff, clang-format, and diff checks pass, followed by exactly
  one PR-scoped commit.

## 33. PR 36: Typed Residue Runtime-Behavior Gate

Close the native typed residue contract by materializing
`/atoms/residue_index` into `Xponge::system.residues` before MD-core residue
consumers initialize.

- Apply the typed route only to native H5 input. Reject simultaneous
  `residue_in_file` ownership instead of silently selecting either source.
- Require a one-dimensional atom-sized residue-index payload whose labels start
  at zero and remain contiguous and nondecreasing. Convert those labels to the
  per-residue atom counts consumed by the native runtime. When
  `/residues/atom_offset` is also present, require it to start at zero, end at
  atom count, remain nondecreasing, and describe exactly the same mapping.
- Run two deterministic same-semantic A/B pairs. The legacy branches consume
  `residue_in_file`; the bundled branches contain typed
  `/atoms/residue_index=[0,0,1,1]` and
  `/residues/atom_offset=[0,2,4]`, with no topology sidecar table or residue
  sidecar payload.
- Retain only non-target protocol sidecars required to drive the behavior
  oracles. The PBC case proves whole-molecule coordinate mapping; the COM case
  proves residue-membership-sensitive restraint virial and pressure.
- Reuse the `[1,3]` wrong-partition control by mutating the typed payload to
  `[0,1,1,1]`. Require unchanged bond/restraint energy and complete force, but
  a large pressure/Pxx change, proving that the typed membership reaches its
  runtime consumer.

Acceptance:

- Typed and legacy PBC routes both load two residues, produce non-zero
  `bond=2.0`, retain the complete 12-value force, and emit the mapped position
  fingerprint `(19,21,25,28)`.
- Typed and legacy COM routes both produce `bond=2.0`, `restrain=2.0`,
  `pressure=0.04`, and `Pxx=0.11`. The `[1,3]` control preserves energy and
  force while changing pressure/Pxx to `-11.53/-34.60`.
- Both existing residue-sidecar behavior gates continue to pass, demonstrating
  that typed materialization does not override the sidecar route.
- A residue-index/offset disagreement exits with `spongeErrorBadFileFormat`;
  simultaneous typed and `residue_in_file` ownership exits with
  `spongeErrorConflictingCommand`, retaining route-specific diagnostics.
- Registry, manifest, real CPU A/B, full smoke, related native-reader CTests,
  SPONGE rebuild, Ruff, clang-format, and diff checks pass, followed by exactly
  one PR-scoped commit.

## 34. PR 37: Typed SW Pair/Three-Body Behavior Gate

Close `input.manybody.sw` with a same-semantic typed/legacy A/B gate rather
than treating successful SW initialization as behavioral evidence.

- The legacy branch consumes the canonical `SW_in_file`; the bundled branch
  contains only `/manybody/sw` atom-type, pair, and triple datasets, with no SW
  command, H5 sidecar table, or sidecar directory.
- Materialize the typed payload into the existing SW parser so both routes
  share the established validation and force kernel. Validate dataset ranks,
  type ranges, duplicate/missing parameter rows, and complete pair/triple
  parameter coverage before runtime initialization.
- Compare the complete one-step mdout and nine-value force payload. Require
  both branches to produce `SW=194.50` with maximum absolute force
  `404.27862548828125`.
- Mutate only the typed triple lambda from `32.5` to `0`. Require the control
  to produce `SW=158.79`, maximum absolute force `343.52105712890625`, and a
  force delta greater than `1`, proving that the typed three-body parameters
  reach the runtime consumer.

Acceptance:

- The typed bundled route materializes
  `.sponge_h5_native_manybody/sw.txt`; no legacy SW route survives in its
  command file or H5 sidecar metadata.
- The typed and legacy branches match the complete deterministic mdout and
  force behavior, while the lambda-zero counterfactual changes both energy and
  force as expected.
- The existing SW sidecar pair/three-body gate still passes, demonstrating
  that native typed materialization does not override the sidecar route.
- Registry, manifest, real CPU A/B, full smoke, related native-reader CTests,
  SPONGE rebuild, Ruff, clang-format, and diff checks pass, followed by exactly
  one PR-scoped commit.

## 35. PR 38: Typed Constraint Projection Behavior Gate

Close `input.protocol.constraint` by materializing typed protocol pairs into
`Xponge::system.classical_force_field.constraints` before the SHAKE/SETTLE
runtime builds its constraint list.

- Preserve migration compatibility: when `constrain_in_file` is bound by the
  protocol sidecar table, the sidecar remains the owner even if the converted
  H5 file also retains typed constraint datasets.
- For a pure typed route, require both
  `/constraint/default/pairs/atoms` with shape `[n,2]` and
  `/constraint/default/pairs/r0` with shape `[n]`. Reject empty payloads,
  mismatched lengths, out-of-range or self pairs, and non-positive/non-finite
  distances.
- Reuse the deterministic two-atom projection fixture. The legacy branch
  consumes `constrain_in_file`; the bundled branch contains only typed pair
  `[0,1]` with `r0=1.5`, with no constraint command, sidecar table, or sidecar
  directory.
- Compare all four position and velocity frames. Require zero distance error
  and reduction of the initial relative radial speed from `2.0` to at most
  `1e-4` through both routes.
- Change only the typed target and matching bootstrap distance to `2.0` while
  retaining non-zero relative motion. Require all four frames to converge to
  the new target; an implementation that ignores the typed payload must fail.

Acceptance:

- Legacy and pure typed routes have identical complete mdout, position, and
  velocity behavior at `r0=1.5`; both report four projected frames and maximum
  distance residual `0`.
- The `r0=2.0` control reports four projected frames, zero distance residual,
  and radial velocity residual below `1e-4`.
- An out-of-range typed pair `[0,2]` exits with
  `spongeErrorBadFileFormat` and stable dataset/row diagnostics.
- The existing constraint sidecar projection gate remains green, proving that
  typed materialization does not override the migration route.
- Registry, manifest, real CPU A/B, full smoke, related input CTests, SPONGE
  rebuild, Ruff, clang-format, and diff checks pass, followed by exactly one
  PR-scoped commit.

## 36. PR 39: Typed Tersoff Angular Behavior Gate

Close `input.manybody.tersoff` with a pure typed same-semantic behavior gate.

- The legacy branch consumes canonical `TERSOFF_in_file`; the bundled branch
  contains `/manybody/tersoff` atom types, names, map, entry types, raw
  14-parameter rows, and runtime 18-parameter rows, with no Tersoff command,
  sidecar table, or sidecar directory.
- Materialize the typed payload through the established Tersoff text parser.
  Validate atom/type counts, names, entry shapes and indices, duplicate rows,
  map consistency, finite raw values, and agreement between raw and converted
  runtime parameters before force initialization.
- Compare the complete deterministic mdout and nine-value force payload.
  Require both routes to produce `potential=-173.23`,
  `eff_pot=-173.23468`, and maximum absolute force
  `135.94906616210938`.
- Change gamma from `1` to `0` in both typed parameter representations and
  rerun the bundled branch. Require `potential=-196.06`,
  `eff_pot=-196.05984`, maximum absolute force `144.7831268310547`, and
  maximum force delta greater than `25`.

Acceptance:

- Pure typed and legacy routes match the complete energy and force behavior,
  while the gamma-zero counterfactual proves the angular bond-order parameter
  reaches the runtime kernel.
- The bundled route materializes
  `.sponge_h5_native_manybody/tersoff.txt`; no legacy ownership remains.
- A map/entry conflict and a raw/runtime-parameter conflict both exit through
  `spongeErrorBadFileFormat` with dataset-specific diagnostics.
- The existing Tersoff sidecar angular gate remains green, proving typed
  materialization does not override migration-sidecar precedence.
- Registry, manifest, real CPU A/B, full smoke, native input CTests, SPONGE
  rebuild, Ruff, clang-format, and diff checks pass, followed by exactly one
  PR-scoped commit.

## 37. PR 40: Typed QC Type-Sensitive Behavior Gate

Close `input.qc.type` with an unrestricted, type-sensitive rerun gate rather
than an initialization check.

- The legacy branch consumes `qc_type_in_file`; the bundled branch removes
  only that binding and payload and consumes `/qc/type` with two atoms,
  symbols `H/N`, charge `0`, and multiplicity `3`. The fixture also removes
  the duplicate residue sidecar owner while retaining its typed residue data,
  so QC is the only newly varied behavior.
- Materialize the typed QC payload through the existing `qc_type_in_file`
  parser. Validate count, atom-index and symbol shapes, unique in-range atom
  indices, non-empty symbols, and positive multiplicity. Preserve explicit
  `qc_type_in_file` precedence for migration bundles.
- Compare both rerun frames across the complete mdout column set, including
  non-trivial `QC` and `QC_S_sq`, and require normalized-exact SCF text
  equivalence.
- Change only typed multiplicity from triplet `3` to singlet `1` and rerun the
  bundled branch. Require a non-trivial response in QC energy or spin square,
  proving typed electronic state reaches the SCF runtime.

Acceptance:

- Legacy and typed bundled routes match both complete deterministic mdout
  frames and exact SCF output. The typed route materializes
  `.sponge_h5_native_qc/qc_type.txt` and contains no QC type sidecar owner.
- The singlet control changes maximum QC energy by about `89.27` and maximum
  spin square by about `2.0065` relative to the triplet baseline; both runs
  remain finite and successful.
- Existing unrestricted VDS-off/on sidecar QC gates remain green.
- Registry, manifest, real CPU A/B, full smoke, related native input CTests,
  SPONGE rebuild, Ruff, changed-line clang-format, and diff checks pass,
  followed by exactly one PR-scoped commit.

## 38. PR 41: Typed Steering CV Behavior Gate

Close `input.protocol.steering` with a pure typed same-semantic behavior gate,
using the protocol surface actually consumed by the CV controller.

- Correct the contract inventory: working legacy steering is a `steer` block
  plus its referenced CV definitions in `cv_in_file`; the equivalent bundled
  surface is `/cv/config`. Do not promote the unconsumed
  `steer_cv_in_file`/`/steer` conversion surface.
- Materialize `/cv/config` through the established CV text parser. Require a
  positive section count, one-dimensional section/key/value vectors,
  consistent monotonic key offsets, and config-safe section names, keys, and
  values. Preserve an existing `cv_in_file` owner so converted sidecar bundles
  keep migration-route precedence.
- Reuse the isolated two-atom distance fixture. The legacy branch consumes a
  `steer` section with weight `2.0` and a distance CV over atoms `0 1`; the
  bundled branch contains the same definitions only in `/cv/config`, with no
  `cv_in_file` command, protocol sidecar table, or sidecar directory.
- Compare the complete deterministic mdout and six-value force payload. Then
  change only bundled `/cv/config/value` weight from `2.0` to `0.0` and run a
  process-level control, proving the typed parameter reaches the steering
  energy and force consumers.

Acceptance:

- Legacy and pure typed routes both produce `steer_cv=3.0` for distance `1.5`
  and weight `2.0`, with analytic force `(2,0,0,-2,0,0)` and identical complete
  mdout behavior.
- The typed route materializes `.sponge_h5_native_protocol/cv.txt`; no legacy
  CV ownership remains in its command file, H5 sidecar metadata, or files.
- The weight-zero control produces steering energy `0`, maximum absolute force
  `0`, and maximum force delta `2.0` from baseline.
- The existing steering `cv_in_file` protocol-sidecar gate remains green,
  proving that typed materialization does not override the migration route.
- Registry, manifest, real CPU A/B, full smoke, related native input CTests,
  SPONGE rebuild, Ruff, clang-format, and diff checks pass, followed by exactly
  one PR-scoped commit.

## 39. PR 42: Typed SITS Configuration/Selection Behavior Gate

Close `input.protocol.sits` with pure typed configuration, atom-selection, and
restart-state behavior rather than treating SITS initialization as evidence.

- Use the deterministic two-atom selective LJ/Coulomb fixture. The legacy
  branch consumes inline `SITS_*` configuration, `SITS_atom_numbers=2`, and
  `SITS_nk_in_file`; the bundled branch removes every inline `SITS_*` owner and
  consumes `/sits/config`, `/sits/atom_indices=[0,1]`, and typed Nk protocol
  restart state.
- Materialize `/sits/config` through the established controller command
  interface and `/sits/atom_indices` through the established atom-list parser.
  Require exactly the active SITS section, consistent section offsets and
  key/value lengths, unique config keys, config-safe tokens, a mode entry, and
  unique in-range atom indices. Preserve explicit inline configuration and atom
  routes for migration compatibility.
- Compare the complete deterministic mdout row and six-value force payload.
  Require non-trivial `SITS_AA_kAB`, `SITS_bias`, and `SITS_fb` behavior and the
  existing asymmetric `Nk=(1,4)` oracle.
- Change only typed `pe_a` from `1.0` to `0.5`; require bias `-1.4591`, force
  factor `0.6471`, and maximum force delta above `0.01`. Separately change the
  typed atom set from `[0,1]` to `[0]`; require `SITS_AA_kAB=-0.61`, bias
  `-0.7295`, force factor `0.6471`, and maximum force delta above `0.02`.
- Supply a complete inline SITS configuration with `pe_a=0.5` while leaving
  typed config at `1.0`. Require the inline fingerprint, no materialized typed
  config, and a materialized typed atom list, proving config and atom ownership
  precedence independently.
- Run process-level duplicate-index `[0,0]` and duplicate-config-key controls.
  Require `spongeErrorBadFileFormat`, `Materialize_H5_SITS_Input`, and
  payload-specific diagnostics for both.

Acceptance:

- Legacy and pure typed routes both produce `SITS_AA_kAB=-1.22`,
  `SITS_bias=-0.5317`, `SITS_fb=0.7049`, `eff_pot=-1.2829471`, and the same
  non-zero six-value force with maximum magnitude `0.2489061951637268`.
- The typed branch materializes `.sponge_h5_native_protocol/sits.txt` and
  `sits_atom.txt`; it contains no inline SITS configuration or protocol
  sidecar owner.
- Both typed counterfactuals match their deterministic behavior fingerprints,
  so ignored config or ignored atom selection cannot pass. The inline control
  proves explicit config precedence without suppressing typed atom selection;
  duplicate atom and config-key inputs are rejected before SITS initialization.
- The existing typed Nk restart case remains green, proving the new
  configuration path does not weaken restart-state behavior.
- Registry, manifest, real CPU A/B, full smoke, related native input CTests,
  SPONGE rebuild, Ruff, clang-format, and diff checks pass, followed by exactly
  one PR-scoped commit.

## 40. PR 43: End-to-End Improper Conversion Behavior Gate

Close `input.topology.improper` with unmodified converter output for both
canonical `improper_dihedral_in_file` and alias `improper_in_file` inputs.

- Change XPONGE's typed improper schema from `/forcefield/improper/k` to the
  SPONGE-native `/forcefield/improper/pk`; retain reverse export fallback for
  older bundles that still contain `/k`.
- Route legacy `improper_in_file` through the same native improper loader as
  `improper_dihedral_in_file`, so the alias branch consumes the same force-field
  payload instead of merely starting the process without an improper term.
- Generate two independent four-atom cases with identical `pk=10.0` and
  `phi0=0.2` semantics, differing only in the canonical versus alias legacy
  key. Run each legacy input against the exact `legacy-to-bundle` output
  without renaming datasets, deleting metadata, or rewriting the bundled mdin.
- Require both conversion manifests to report `typed_converted` at
  `/forcefield/improper`, require neither improper legacy key nor improper
  sidecar ownership in the bundle, and validate exact atoms, count, `pk`, and
  `phi0` values before execution.
- Compare the complete deterministic mdout row and 12-value force payload.
  Change only typed `pk` from `10` to `5` and require both energy and every
  force component to scale by one half, proving the typed dataset reaches the
  force kernel rather than merely initializing the module.

Acceptance:

- Canonical, alias, and bundled routes each produce
  `improper_dihedral=31.36` and maximum force
  `35.415924072265625`, with complete cross-branch equality.
- The `pk=5` controls produce `improper_dihedral=15.68`, maximum force
  `17.707962036132812`, and exact half-scale force within deterministic
  float32 tolerances.
- XPONGE parser/exporter focused tests, both converted real CPU A/B cases, the
  existing pure-native regression, full smoke, related native input CTests,
  Ruff, and diff checks pass. XPONGE and SPONGE are committed separately, once
  per independent PR scope.

## 41. PR 44: Execution-Matrix Residue Ownership Fixture Repair

Repair the pinned TIP3P matrix fixture after the typed residue runtime gate made
legacy and typed residue inputs mutually exclusive.

- Make XPONGE's canonical `residue_in_file` contract `typed_required`, so
  `legacy-to-bundle` emits `/atoms/residue_index` and
  `/residues/atom_offset` without retaining an active residue sidecar owner.
- Prove reverse conversion consumes typed H5 state by changing one residue into
  two residues before export and checking the exported legacy partition.
- Regenerate the pinned topology and matching protocol with the repaired
  converter, remove the stale residue sidecar, and pin converter provenance and
  new SHA-256 values in the fixture manifest.
- Parse the pinned H5 in the manifest gate and require both typed residue
  datasets, no `residue_in_file` sidecar key, and no external residue sidecar.
- Re-run same-semantic CPU/GPU deterministic matrix cases and at least one GPU
  statistical thermostat/constraint case; initialization alone is not
  acceptance evidence.

Acceptance:

- XPONGE bundle regression tests pass, the real TIP3P conversion has one
  residue owner, and XPONGE is committed separately from SPONGE.
- The pinned-fixture hash/ownership gate passes.
- CPU and GPU deterministic cases pass complete mdout and H5 comparisons, and
  a four-replica GPU statistical case passes the statistical mdout/H5 gates.
- The former GPU startup failure naming
  `Xponge::Read_H5_Residue_Atom_Numbers` is absent.

## 42. PR 45: CUDA+MPI NCCL Discovery Prerequisite

Make the existing CUDA+MPI build path capable of resolving NCCL through the
project's CMake module search path.

- Add a standard `FindNCCL.cmake` that accepts `NCCL_ROOT`, searches both
  `lib` and `lib64`, and exports the plural include/library variables already
  consumed by `cmake/utils/mpi.cmake`.
- Keep package installation and runner provisioning outside the finder; the
  finder only resolves an available NCCL installation.
- Configure and build the real `SPONGE` target with CUDA, MPI, OpenMPI, and
  NCCL enabled before accepting the prerequisite.

Acceptance:

- CMake reports both MPI and NCCL found under `PARALLEL=cuda` and `MPI=ON`.
- The complete CUDA+MPI `SPONGE` executable links successfully with `libmpi`
  and `libnccl` and has no unresolved runtime libraries.
- CMake formatting and diff checks pass, followed by one PR-scoped commit.

## 43. PR 46: Portable GPU MPI Device Mapping

Make the GPU rank-2 matrix executable on both one-GPU and multi-GPU runners
without changing the simulated system or branch semantics.

- For GPU cases with more than one MPI rank, default both ranks to device `0`,
  which permits local single-GPU validation.
- Allow `SPONGE_BUNDLED_IO_AB_GPU_DEVICES` to provide an explicit per-rank map
  such as `0 1` on a multi-GPU runner.
- Inject the same device map into legacy and bundled mdin files and remove any
  stale fixture device key before insertion.
- Keep GPU rank-1 behavior unchanged.

Acceptance:

- Static tests prove the shared-device default, explicit two-device override,
  and absence of an injected rank-1 device command.
- A real two-rank CUDA+MPI nonorthogonal NPT+SETTLE case runs four replicas and
  passes statistical mdout, trajectory, observable, and restart comparisons.
- Evidence reports two MPI ranks and rank-0-only output ownership.
- Smoke, Ruff, and diff checks pass, followed by one PR-scoped commit.

## 44. PR 47: Production-Fixture Residue Ownership Repair

Bring the production full-contract and focused residue fixtures in line with
the converter's typed-only residue contract without weakening the dedicated
legacy-sidecar behavior coverage.

- Keep `/atoms/residue_index` and `/residues/atom_offset` as the sole residue
  owner in both full-contract bundles. Remove the stale `residue_in_file`
  binding, materialized payload, and manifest sidecar metadata from the
  compatibility bundle while retaining all unrelated compatibility sidecars.
- Construct the focused residue-sidecar cases explicitly from current typed
  converter output: copy the legacy residue payload into a bundle-local
  sidecar, replace the topology binding table with that single owner, and
  remove both typed residue datasets. The typed focused cases then promote the
  same semantic partition back to one typed owner.
- Reject partial typed residue state, no owner, and dual ownership during
  full-contract preparation. Add a static fixture gate that parses both H5
  topologies, verifies typed ownership, checks the manifest, and rejects an
  external residue sidecar directory.
- Remove the typed-QC fixture's obsolete assumption that converter output
  contains a residue sidecar. QC mutation must not change residue ownership.

Acceptance:

- Sidecar and typed residue PBC mapping plus COM/virial behavior cases pass
  against legacy with the complete semantic oracles.
- Full-contract compatibility bundles pass both VDS modes, all QC cases, all
  failure-semantics cases, and the CMAP potential VDS invariant without the
  `Read_H5_Residue_Atom_Numbers` ownership conflict.
- The production sweep may expose later independent behavior gaps only after
  residue initialization succeeds; typed SITS configuration, no-velocity
  rerun temperature, strip/frame numbering, and exact-EOF legacy cleanup are
  retained as failures for follow-up PRs rather than hidden by the fixture.
- Contract smoke, fixture manifest/equivalence CTests, Ruff, formatting, and
  diff checks pass, followed by one PR-scoped commit.

## 45. PR 48: Inactive Typed SITS Configuration Semantics

Preserve the legacy meaning of a syntactically valid SITS section that omits
`mode`: the configuration exists, but the SITS module remains inactive.

- Allow `/sits/config` without a `mode` key to be materialized. Keep all
  structural, token, duplicate-key, and atom-ownership validation unchanged;
  the existing SITS initializer remains the authority that requires `mode`
  before activating the module.
- Add a focused two-atom A/B case whose legacy `SITS_in_file` and pure typed
  `/sits/config` both contain only `CV=distance` and `nk=2`. Remove the SITS
  sidecar owner and embedded restart copy from the bundled branch.
- Require exact typed payload values, no `mode`, no active SITS sidecar, exact
  materialized text, no SITS mdout columns, and complete deterministic mdout
  equivalence.
- Keep the active typed SITS gate mapped to the same contract. It must still
  produce non-zero bias and force behavior with a complete production-mode
  configuration and typed Nk restart state.

Acceptance:

- Active and inactive typed SITS A/B cases both pass, proving both sides of
  the mode-presence contract rather than accepting initialization alone.
- The previous full-contract pure-native bad-format failure at
  `Materialize_H5_SITS_Input` is absent; later module-consumption differences
  remain visible for independent follow-up.
- SPONGE rebuild, 141-test contract smoke, registry/manifest mutation checks,
  Ruff, formatting, clang-format, and diff checks pass, followed by one
  PR-scoped commit.

## 46. PR 49: Canonical and Extra NB14 Behavior Contracts

Separate the two legacy NB14 surfaces and consume both corresponding typed H5
representations with their distinct parameter semantics.

- Keep `nb14_in_file` mapped to `/forcefield/nb14`: two-column parameters are
  LJ and Coulomb scale factors and require typed LJ atom/pair data.
- Enumerate `nb14_extra_in_file` independently at
  `/forcefield/nb14_extra`: three-column parameters are raw A, B, and Coulomb
  scale values, so runtime materialization applies the legacy `12`/`6`
  derivative multipliers before NB14 initialization.
- Reject simultaneous canonical and extra typed roots, partial atoms/params,
  and two-column `nb14_extra` payloads. Preserve the existing three-column
  canonical materialized-A/B representation for compatibility.
- Keep the pre-existing untracked native topology reader outside this PR.
  Consume `nb14_extra` through a small independently owned H5 helper and one
  tracked `Xponge::System::Load_Inputs` integration point, so the commit does
  not absorb unrelated H5 foundation work.
- Add independent two-atom pure typed A/B cases. Require exact non-zero NB14
  LJ/electrostatic energy and complete six-value force fingerprints, then zero
  only the typed LJ parameter(s) and require zero NB14 LJ energy, unchanged
  electrostatic energy, and a non-trivial force response.

Acceptance:

- Canonical scaled input matches at `nb14_LJ=-0.75`, `nb14_EE=-0.25`, and
  maximum force `4.875`; its zero-LJ control changes force by `1.5`.
- Extra input matches at `nb14_LJ=-0.03`, `nb14_EE=-0.25`, and maximum force
  `0.46728515625`; its zero-LJ control changes force by `0.09228515625`.
- The pure full-contract branch now emits matching `nb14_LJ`/`nb14_EE`
  columns and advances to independent EAM/restraint/soft-wall gaps; the
  sidecar full-contract VDS-off/on regressions remain green. Both fixture
  routes assert exactly one NB14 owner.
- Real process failure gates reject dual canonical/extra roots and two-column
  extra parameters with stable bad-file-format categories and diagnostics;
  the complete 21-case failure sweep remains green.
- Registry coverage advances to 83 supported contracts with no deferred
  contracts and one explicitly unsupported contract.
- SPONGE rebuild, native topology reader CTest, 142-test smoke, Ruff,
  formatting, clang-format, and diff checks pass, followed by one PR-scoped
  commit.

## 47. PR 50: Typed EAM Funcfl/Setfl Behavior Contracts

Consume typed `/manybody/eam` input through the same legacy EAM runtime parser
without weakening payload validation or retaining two active owners.

- Materialize both declared EAM formats. Funcfl reconstructs the raw embedding,
  Z, and density tables; setfl reconstructs element metadata, embedding,
  density, and lower-triangular raw pair tables from native pair potentials.
- Preserve optional `/manybody/eam/atom_type` as the exact
  `EAM_atom_type_in_file` sequence. Reject typed/legacy ownership conflicts.
- Validate scalar controls, table dimensions, finite values, supported format,
  and funcfl/setfl-specific metadata before `EAM::Initial` consumes the
  materialized files.
- Add isolated two-atom funcfl and two-type setfl A/B cases. Compare all mdout
  columns and complete forces, then zero only the typed pair source and require
  an exact energy/force response.
- Make full-contract pure and sidecar fixtures select exactly one EAM owner;
  sidecar inventory no longer requires dormant typed EAM data.

Acceptance:

- Focused funcfl matches at `EAM=267.33` and force
  `(-165.90681,0,0,165.90681,0,0)`; its zero-Z control has `EAM=46.12` and
  zero force.
- Focused setfl matches at `EAM=61.49` and force
  `(-11.530274,0,0,11.530274,0,0)`; its zero-pair control has `EAM=46.12` and
  zero force.
- Pure full-contract VDS-off/on restores the legacy `EAM=23.06` column and
  advances to only restraint/soft-wall gaps. Sidecar VDS-off/on remains green.
- Unknown EAM format and invalid embedding-table shape exit through stable
  bad-file-format diagnostics; the complete failure sweep remains green.
- SPONGE rebuild, fixture-equivalence CTest, contract smoke, Ruff,
  clang-format, and diff checks pass, followed by one PR-scoped commit.

## 48. PR 51: Typed Positional-Restraint Behavior Contract

Treat positional restraint as one composite runtime contract rather than
conflating it with CV-restraint configuration. Its typed owner spans protocol
atom IDs and anisotropic weights plus restart reference coordinates.

- Materialize `/restraint/default/atom_indices` and
  `/restraint/default/weight` together with
  `/parameters/restart/references/restraint/default/coordinate` into the
  established positional-restraint parser inputs.
- Validate atom-index bounds, exact vector/matrix shapes, atom-count agreement,
  and finite values before the legacy parser consumes any materialized text.
- Reject typed/legacy dual ownership. Make pure and sidecar full-contract
  fixtures select exactly one complete owner across the protocol and restart
  sidecar tables.
- Add an isolated two-atom A/B case with exact energy and complete-force
  oracles. Independently swap atom IDs, zero weights, and move references to
  the current positions; each control must produce its predicted response.
- Correct the registry surface from the inactive `restrain_in_file` to
  `/restraint` shorthand to the three-part positional contract. Separate typed
  and sidecar positional-restraint evidence, and track typed CV-restraint and
  typed soft-wall config as deferred until each has a non-zero runtime behavior
  gate.

Acceptance:

- The focused case matches at `restrain=42.5` with force
  `(-20,0,0,0,-30,0)`. Swapped-ID, zero-weight, and matching-reference
  controls each have zero restraint energy and force.
- Pure full-contract VDS-off/on restores the legacy `restrain=22.50/11.25`
  frames and advances to only the typed soft-wall gap. Sidecar VDS-off/on and
  both residue COM/virial restraint regressions remain green.
- Dual positional-restraint owners and a two-column weight matrix exit through
  stable conflict/bad-file-format diagnostics; all 25 failure cases pass.
- Registry coverage records 84 supported contracts, two explicitly deferred
  typed protocol contracts, and one explicitly unsupported VDS contract.
- SPONGE rebuild, fixture-equivalence CTest, contract smoke, Ruff, formatting,
  clang-format, and diff checks pass, followed by one PR-scoped commit.

## 49. PR 52: Typed Soft-Wall Behavior Contract

Close the final pure full-contract runtime gap by materializing typed
`/wall/soft/{count,name,potential}` through the established `SOFT_WALLS`
configuration parser.

- Require a positive scalar count and exact `[count]` name/potential vectors.
  Reject empty, duplicate, overlong, or configuration-delimiter-bearing module
  names, empty potentials, and potential lines that can inject a legacy
  configuration section.
- Materialize each typed entry as one `[[[ name ]]]` section with exactly one
  `[[ potential ]]` value, then inject `soft_walls_in_file` only after rejecting
  a typed/legacy dual owner.
- Add an isolated two-atom pure typed A/B case. Require the exact non-zero
  `z_wall` energy and complete six-value force fingerprint, then change only
  the typed potential to a zero expression and require zero energy and force.
- Make full-contract pure and sidecar fixtures select exactly one soft-wall
  owner. Pure retains `/wall/soft`; sidecar removes the dormant typed group and
  retains `soft_walls_in_file` in the protocol sidecar table.
- Add process-level failures for dual ownership and invalid count, name, and
  potential dataset shapes. Promote typed soft wall only after the focused and
  full-contract runtime gates pass.

Acceptance:

- The focused case matches at `z_wall=10.0` with force
  `(0,0,-2,0,0,-6)`. Its typed zero-potential control has zero soft-wall energy
  and force.
- Pure full-contract VDS-off/on now matches every legacy mdout column,
  including `z_wall`, and passes complete rerun trajectory and observable
  equivalence. Sidecar VDS-off/on remains green.
- Dual owners and all three typed shape mutations exit through stable
  conflict/bad-file-format diagnostics; the complete 29-case failure sweep
  remains green.
- Registry coverage advances to 85 supported contracts, one explicitly
  deferred typed CV-restraint contract, and one explicitly unsupported VDS
  contract.
- SPONGE rebuild, fixture-equivalence CTest, contract smoke, Ruff, formatting,
  clang-format, and diff checks pass, followed by one PR-scoped commit.

## 50. PR 53: Typed and Sidecar CV-Restraint Behavior Contracts

Close the final deferred input contract by making the historical
`restrain_in_file` and `restrain_cv_in_file` surfaces reach the same
`RESTRAIN_CV` consumer as typed protocol configuration.

- Extend the CV controller to merge `cv_in_file`, `restrain_in_file`, and
  `restrain_cv_in_file`. Identical repeated section/key/value definitions are
  accepted; conflicting definitions fail instead of silently selecting one
  owner.
- Read `/cv/config`, `/restraint/config`, and `/restraint/cv/config` through one
  structured validator. Treat `/restraint/config` and
  `/restraint/cv/config` as independent historical input surfaces, validate
  each present root's scalar count/vector offsets/key/value structure, reject
  duplicates, merge identical shared CV definitions, and materialize one
  canonical `cv_in_file` for the established runtime parser.
- Add independent two-atom typed and H5-sidecar A/B cases. Both must exercise
  the real legacy restraint keys, emit an exact non-zero `restrain_cv` energy
  and complete force fingerprint, then zero only their own weight payload and
  produce zero energy and force.
- Make pure and sidecar full-contract fixtures select exactly one composite
  CV-restraint owner. Include the non-zero module column in the formal rerun
  input-semantic evidence under both VDS modes.
- Add process-level failures for typed/legacy dual ownership, inconsistent key
  offsets, and conflicting shared CV definitions. A single valid typed
  restraint root remains a supported input, matching the corresponding legacy
  surface.

Acceptance:

- Both focused routes match at `restrain_cv=1.0` with force
  `(4,0,0,-4,0,0)`. Their independent weight-zero controls have zero
  CV-restraint energy and force. A reordered duplicate typed CV definition
  preserves the exact energy and force, proving key order is not semantic.
- Pure and sidecar full-contract VDS-off/on match every mdout column, complete
  rerun trajectory, and observable output while exercising
  `restrain_cv=12.06` in both frames.
- All three ownership/structure/definition mutations exit through stable
  conflict/bad-file-format diagnostics, and the complete 33-case failure sweep
  remains green.
- Registry coverage reaches 87 supported contracts, zero deferred contracts,
  and one explicitly unsupported cross-process VDS contract.
- SPONGE rebuild, fixture-equivalence CTest, contract smoke, Ruff, formatting,
  clang-format, and diff checks pass, followed by one PR-scoped commit.

## 51. PR 54: Attested Production Promotion Evidence

Close the evidence-integrity gap between passing tests and the PR 7c promotion
evaluator. A production-history row must be derived from the actual contract,
execution-matrix, and comparator-mutation artifacts rather than authored as a
set of trusted booleans and performance numbers.

- Give contract, matrix, and comparator reports one shared production run ID.
- Record all eight frozen comparator-mutation test families through a pytest
  hook. Parameterized cases count as passed only when every selected node
  passes; missing, skipped, or failed nodes and a non-zero pytest session exit
  status cannot attest the run.
- Merge independently produced CPU rank-1, CPU rank-2, GPU rank-1, and GPU
  rank-2 matrix reports only when their run IDs match and duplicate scenario
  payloads are identical.
- Give the scheduled production CPU rank-1 and rank-2 jobs the same GitHub run
  ID and upload the real comparator report with the main production artifact.
  GPU evidence remains an explicit separate-runner input, not a fabricated CI
  placeholder.
- Revalidate 100% supported-contract coverage at the declared minimum evidence
  level and every matrix scenario's backend, OMP, MPI, and rank-0 ownership
  metadata before deriving a production run.
- Derive conservative performance values from the maximum observed scenario
  runtime, finalize, and output-size ratios. Reject missing, negative, NaN, or
  infinite performance values.
- Append history atomically and reject duplicate run IDs. Every persisted run
  records its source commit plus SHA-256 provenance for the contract report,
  canonical merged matrix report, and comparator report. The history loader
  rejects unattested handwritten rows.

Acceptance:

- A complete split matrix, complete supported-contract report, and real
  mutation report derive one loadable production-history row with exact hashes.
- Missing mutation evidence, below-minimum contract evidence, false environment
  metadata, stale/conflicting run IDs, non-finite performance values, duplicate
  history rows, and provenance-free handwritten rows are rejected.
- The real statistics suite writes all eight mutation families under one run
  ID, including all parameterizations of non-finite kind/sign checks.
- This PR does not change `promotion_state=shadow`: final PR 7c promotion still
  requires three complete, consecutive, retry-free production runs and the
  cross-process VDS contract remains explicitly unsupported.
- Contract smoke, the real mutation hook, Ruff, formatting, and diff checks
  pass, followed by one PR-scoped commit.

## 52. PR 55: Clean Source Provenance Gate

Prevent production evidence from being attributed to a commit when the tested
source tree contains tracked modifications or untracked implementation files.
This is required before the three-run PR 7c streak can be meaningful: a GPU
result from a dirty worktree is useful debugging evidence but not release
evidence for the named commit.

- Before reading any behavior artifact, resolve the repository `HEAD` and
  require it to equal the requested source commit.
- Run `git status --porcelain=v1 --untracked-files=all` and reject any tracked
  modification, staged change, deletion, rename, or untracked path.
- Persist `source_tree_state=clean` with the source commit and artifact hashes.
  The production-history loader rejects rows without this attestation.
- Keep generated production outputs under ignored `outputs/` paths so evidence
  collection itself does not dirty an otherwise clean checkout.

Acceptance:

- A temporary clean Git repository passes exact-HEAD attestation.
- Commit mismatch, modified tracked content, and untracked content each fail
  before contract, matrix, or comparator artifacts are read.
- The real current worktree is rejected as dirty, proving that the existing
  local GPU-capable binaries cannot yet establish a release evidence streak.
- Contract smoke, focused manifest tests, Ruff, formatting, and diff checks
  pass, followed by one PR-scoped commit.

## 53. PR 56: Track Bundled I/O Implementation Baseline

Put the implementation, schemas, fixtures, tests, examples, and documentation
that the A/B gates exercise under version control so clean-source provenance
can name an actual bundled I/O commit.

- Track the audited baseline scope: 36 modified tracked files and 94 previously
  untracked implementation/test assets. Exclude `.codegraph`, build outputs,
  and all `.tmp-*` evidence directories.
- Make aggregate CMake targets build every executable selected by their CTest
  labels, and order staged VDS execution after both mock and real backend
  dependencies are available.
- Keep legacy `rst` classified as an output, distinguish H5 path bindings from
  the rerun particle-stream selector, and align static matrix expectations with
  the generated full-contract fixtures.
- Preserve independent typed CV restraint roots and let injected H5 legacy
  sidecars retain ownership after the pre-injection conflict validator has
  accepted the input plan. Repair restart-smoke `mdinfo` rewriting so it does
  not create duplicate TOML keys.
- Record opt-in runtime results separately from compile/file-level coverage.
  Passing registry and CTest labels does not imply production promotion or full
  bundled/legacy behavioral equivalence.

Acceptance:

- A fresh CPU build completes the `SPONGE` target and the staged contract,
  coverage, matrix, mock-backend, VDS, backend-I/O, rerun, and default smoke
  labels without failures; opt-in runtime skips are reported explicitly.
- The opt-in normal input/output A/B matrix and VDS terminal/resume smoke pass.
  The EDIP non-triviality gate and dynamic NHC restart-load closure remain
  explicit failing evidence to be fixed in later behavior-alignment PRs.
- The 152-test bundled I/O contract smoke, schema examples, Ruff, lockfile,
  formatting, and diff checks pass.
- Stage only the audited baseline paths and create exactly one PR-scoped commit
  authored by `Xiaoxuan Yu <xiaoxuan_yu@pku.edu.cn>`.

## 54. PR 57: Non-Trivial REAXFF/EDIP Runtime Parity

Repair the combined many-body runtime gate so its EDIP assertion exercises an
interacting pair rather than the full-contract fixture's excluded pair outside
the EDIP cutoff.

- Rewrite only temporary test copies of the two rerun frames, using identical
  float32 values in the legacy binary trajectory and bundled H5 trajectory.
- Remove the same pair exclusion from legacy, native H5, and legacy-sidecar
  input routes before execution, while retaining REAXFF and the complete
  sidecar/pure-H5 route matrix.
- Assert the prepared binary/H5 coordinates and exclusion ownership before
  launching SPONGE; retain the existing finite, non-zero EDIP and REAXFF
  checks, full mdout-column parity, H5 observable checks, and VDS off/on runs.
- Keep the focused production `normal_edip_nonzero` energy/force gate as an
  independent pure-H5 regression.

Acceptance:

- `test_h5_reaxff_edip_runtime_parity` passes with runtime smoke enabled and
  reports finite non-zero EDIP and REAXFF values for every input/output route.
- The focused production EDIP A/B case passes without sidecars and preserves
  its non-zero energy and force comparison.
- Static staged H5 bundle labels, formatting, and diff checks pass, followed
  by exactly one PR-scoped commit authored by
  `Xiaoxuan Yu <xiaoxuan_yu@pku.edu.cn>`.

## 55. PR 58: Dynamic NHC Restart-Load Runtime Closure

Close the dynamic Nose-Hoover-chain restart-load runtime gate with a fixture
that exercises normal NVT continuation without unrelated full-contract rerun
modules.

- Preserve the real bundled topology, protocol binding, and restart state, but
  remove force-field, many-body, QC, enhanced-sampling, restraint, steering,
  constraint, and wall roots only from the temporary dynamic-NHC test copy.
- Keep the other protocol/full success and unsupported-policy failure cases on
  their existing broad fixtures so their contract coverage does not shrink.
- Run one NVT step at `dt=0`, require NHC initialization, reject non-finite
  `mdout` and restart state, require a non-zero input state, and require exact
  input/output NHC-state equality.
- Record the observed crash correctly: the former fixture reached the QC Fock
  builder during normal dynamics and failed before it could establish NHC
  restart continuation behavior.

Acceptance:

- The focused restart-load closure passes with runtime smoke enabled and the
  output restart contains a finite, non-empty, non-zero, exactly preserved NHC
  state.
- The complete opt-in runtime group passes the normal input/output matrix,
  non-trivial REAXFF/EDIP parity, restart-load closure, and VDS terminal/resume
  smoke together.
- Static staged H5 bundle labels, the 152-test contract smoke, formatting, and
  diff checks pass, followed by exactly one PR-scoped commit authored by
  `Xiaoxuan Yu <xiaoxuan_yu@pku.edu.cn>`.

## 56. Artifacts and Temporary-Space Policy

The detailed next-round strengthening scopes for PR 59 through PR 74 are
maintained in
[`bundled_io_ab_next_round_strengthening_plan.md`](bundled_io_ab_next_round_strengthening_plan.md).
That plan explicitly keeps cross-process VDS reopen-and-append/resume
unsupported for this round and requires every promotion statement to exclude
that capability.

- Use `SPONGE_BUNDLED_IO_AB_RUN_ROOT` for all heavy runs.
- Put production artifacts on the workspace filesystem rather than `/tmp`.
- Delete successful per-replica payloads after evidence extraction.
- Retain failed case inputs, stdout/stderr, mdout, H5 metadata, evidence, and a
  minimal reproduction command.
- Record artifact byte counts and enforce a configurable per-case quota.
- Cleanup must only remove directories created by the current gate run.

## 57. PR Completion Log

Append one row immediately after completing and committing each PR.

| PR | Status | Commit | Verification | Notes |
|---|---|---|---|---|
| PR 1: Executable contract registry | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; sidecar rerun VDS-off A/B | 21 contract tests and one real A/B case pass; pure-native rerun exposes a strict column-set gap retained for input behavior closure. |
| PR 2: Contract enumeration and fixtures | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; H5 manifest/native-reader/fixture-equivalence CTest; sidecar rerun VDS-off A/B | Runtime key inventory has exact contract ownership; full-contract mdin/manifest now enumerate output controls, all legacy routes, particle stream, and rerun controls; missing topology runtime cases are explicitly deferred to PR 5. |
| PR 3: Comparators | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real sidecar rerun VDS-off A/B (expected QC failure) | 41 contract/comparator tests pass; TOST plus Holm correction, quantity-specific margins, structured position/velocity/force/box statistics, exact non-finite patterns, and all six comparator mutations are enforced. The real gate still rejects the known QC first-frame behavior gap (`-34183.44` vs `-34139.92`) without weakening tolerances. |
| PR 4: Output behavior closure | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; five output/VDS CTest targets; two-step NVE restart continuation; real TIP3P A/B (expected force-schedule failure) | 45 contract tests and all output writer/backend targets pass. Normal VDS-off/on cases now require position, velocity, force, box, structured mdinfo, explicit route provenance/content, and E4 structural restart continuation. The real gate exposes legacy `frc` 51-frame vs H5 50-frame schedule skew; QC SCF non-empty/exact comparison and dynamic restart continuation remain explicit behavior gaps rather than relaxed assertions. |
| PR 5: Input behavior closure | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; native-topology reader and fixture-equivalence CTest; saved TIP3P and sidecar-rerun production artifacts | 52 contract tests require every supported input contract to produce a present, non-trivial module-owned result and then pass deterministic or statistical A/B comparison. TIP3P mass/charge/LJ pass; the rerun gate rejects the known QC first-frame mismatch (`-34183.44` vs `-34139.92`). Reader/path checks remain E1/E2 only. Residue/exclusions, improper, GB/softcore/subsystem/virtual atoms, SW/EDIP/Tersoff, custom pairwise, constraint/steering, SITS/meta/NHC, and QC type/spin/SCF stay explicitly deferred until dedicated non-trivial fixtures exist. |
| PR 6: Rerun and failure semantics | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; nine process-level F1 cases; five real rerun boundary cases | 60 contract/manifest tests and all nine textual/binding F1 cases pass. The harness now keeps `crd/box/vel` as legacy rerun inputs and inserts overrides before TOML module tables; this removes the earlier QC first-frame skew (`-33917.02` now matches). The 0/0 limit-1 boundary passes, while the gate exposes optional-velocity potential (`1951.16` vs `1947.41`), strip frame counter (`2` vs `1`), and box-update legacy invalid-free failures. H5 atom/shape/dtype/schema and sidecar-table process mutations, restart absent/dynamic/protocol/full continuation, and owner-state F1 cases remain explicitly deferred. |
| PR 7a: Executable shadow promotion guard | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; execution-matrix CLI; workflow YAML parse; `git diff --check` | 70 contract/manifest tests pass. The typed matrix and mutation-tested promotion evaluator reject missing axes/combinations/evidence, unproven environment metadata, retries, comparator-mutation gaps, and runtime/finalize/output-size budget overruns. Medium and production jobs are explicit `continue-on-error` shadow jobs and are absent from the release workflow; runtime evidence remained intentionally deferred at this boundary. |
| PR 7b: Portable runtime matrix and CPU/MPI evidence | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; six CPU rank-1 fast A/B cases; two CPU rank-2 fast A/B cases; Ruff; workflow YAML parse; `git diff --check` | 74 contract/manifest tests and all eight CPU/MPI runtime cases pass. The matrix pins same-semantic fixtures, records backend/OMP/MPI/rank-0 ownership evidence, and compares complete H5 trajectory/observable/restart behavior with periodic-aware deterministic and statistical position checks. Four GPU cases are executable but lack real-device evidence, so promotion remains shadow-only for PR 7c. |
| PR 7c: GPU evidence and final promotion | Pending | | | |
| PR 8: VDS chunk-boundary behavior closure | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; four real CPU VDS chunk-boundary A/B cases; Ruff; `git diff --check` | 76 contract/manifest/comparator tests pass. With chunk size 4, same-semantic deterministic NVE branches produce and compare 3/4/5/9 frames with 1/1/2/3 shards. Frame-count and shard-count mutations are rejected. The cases use `dt=0.0001` to isolate chunk/finalize behavior without relaxing deterministic numeric tolerances. |
| PR 9: Focused non-zero EDIP input behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_edip_nonzero` CPU A/B; Ruff; `git diff --check` | 79 contract/manifest/comparator tests and the focused real A/B pass. The pure bundled branch removes the sidecar table and directory, materializes `/manybody/edip`, and matches legacy at non-zero `EDIP=52.98159` and maximum force `400.30548`. The independently exposed H5 zeroth-frame/force-payload schedule gap remains outside this input-only contract. |
| PR 10: H5 sidecar-table failure semantics | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; three real bundled CPU failure cases; Ruff; `git diff --check` | 79 contract/manifest/comparator tests and all three F1 process cases pass. Unsupported key, key/path length mismatch, and same-key different-path conflict each exit with code 238 and `spongeErrorValueErrorCommand`, while retaining route-specific stable diagnostic tokens. |
| PR 11: Restart owner-state failure semantics | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; three real bundled CPU failure cases; Ruff; `git diff --check` | 79 contract/manifest/comparator tests and all three F1 process cases pass. Dynamic, protocol, and full policies each exit with code 235 and `spongeErrorConflictingCommand`; diagnostics distinguish missing Nose-Hoover-chain and metadynamics owners. |
| PR 12: Focused non-zero custom-pair input gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_custom_pair_nonzero` CPU A/B; Ruff; `git diff --check` | 80 contract/manifest/comparator tests and the focused real A/B pass. The pure bundled branch has no legacy custom-force keys or sidecars, materializes both native inputs, and matches legacy at `custom_pair=31.57` with maximum force `252.554443359375`. |
| PR 13: Focused native exclusion semantics gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_exclusions_coulomb_oracle` CPU A/B; Ruff; `git diff --check` | 83 contract/manifest/comparator tests and the focused real A/B pass. The pure bundled branch uses offsets/list `[0,1,1,1]`/`[1]`, and both branches match the `-1/12` energy and analytic nine-value force oracle with maximum force `0.1111111119389534`. |
| PR 14: Focused native LJ soft-core behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; `ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R 'test_topology_native_h5_reader|test_h5_input_fixture_equivalence'`; real `normal_lj_soft_core_nonzero` CPU A/B; Ruff; `git diff --check` | 85 contract/manifest/comparator tests, both H5 input contract tests, and the focused real A/B pass. The gate first exposed legacy `potential=-0.06` versus bundled `-0.01`; pure native runtime normalization now matches legacy at `LJ_soft=-0.06`, `eff_pot=-0.056708537`, maximum force `0.1920633763074875`, and zero full-mdout error. Subsystem division remains deferred because its mask is not consumed by the force kernel. |
| PR 15: Focused native virtual-atom behavior gates | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; `ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R 'test_topology_native_h5_reader|test_h5_input_fixture_equivalence'`; real `normal_virtual_atoms_all_types` and `normal_virtual_atoms_pbc_boundary` CPU A/B; Ruff; clang-format; `git diff --check` | 88 contract/manifest/comparator tests, both H5 input contract tests, and both focused real A/B cases pass. The all-types gate covers exact type 0/1/2/3 ragged payloads, layered type 2-to-3 reconstruction, 24-value coordinate oracles, non-zero PM, and complete force parity; the PBC gate distinguishes `9.75` from the incorrect non-periodic `7.25`. It exposed and fixed the type 3 host/device copy source and the missing global first-frame coordinate refresh. Registry coverage advances to 52 supported, 22 deferred, and 1 unsupported contract. |
| PR 16: Virtual-atom plural alias behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; `ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R 'test_topology_native_h5_reader|test_h5_input_fixture_equivalence'`; real `normal_virtual_atoms_plural_alias` CPU A/B; Ruff; clang-format; `git diff --check` | 89 contract/manifest/comparator tests, both H5 input contract tests, and the focused real A/B pass. The legacy branch resolves only `virtual_atoms_in_file`; the bundled branch contains the exact native payload with neither key nor sidecars. Both paths match the boundary-crossing coordinate `9.75`, non-zero PM, complete force, and mdout behavior. Registry coverage advances to 53 supported, 21 deferred, and 1 unsupported contract. |
| PR 17: Restart-absent same-bootstrap rerun gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; `ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R 'test_(restart_h5_reader\|h5_input_validation\|h5_restart_load_runtime_closure)$'`; real `rerun_restart_absent_same_bootstrap_vds_off` CPU A/B; Ruff; `git diff --check` | 90 contract/manifest/comparator tests, two related CTest passes (runtime closure skipped by its build guard), and the real two-frame A/B pass. Both branches use byte-identical coordinate/velocity bootstrap files, the bundled branch has no restart binding or payload, and all mdout columns, frame selection, trajectory, and observable semantics pass E3 comparison. Registry coverage advances to 54 supported, 20 deferred, and 1 unsupported contract. |
| PR 18: NHC dynamic restart continuation gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; `ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R 'test_(restart_h5_reader\|h5_restart_load_runtime_closure\|module_h5_mappings_with_mock_backend)$'`; real `normal_nhc_dynamic_restart_continuation` CPU A/B; Ruff; `git diff --check` | 91 contract/manifest/comparator tests, two related CTest passes (runtime closure skipped by its build guard), and the real E4 A/B pass. One producer emits non-zero three-chain state through legacy and H5 restart routes, then a text-restart continuation and `dynamic` H5 continuation compare all mdout columns, particle position/velocity/force/box, NHC trajectories, observable output, and final restart state. Registry coverage advances to 57 supported, 17 deferred, and 1 unsupported contract. |
| PR 19: Protocol/full metadynamics continuation gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; `ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R 'test_(restart_h5_reader\|h5_restart_load_runtime_closure\|module_h5_mappings_with_mock_backend)$'`; real `normal_meta_protocol_full_restart_continuation` CPU A/B; Ruff; `git diff --check` | 92 contract/manifest/comparator tests, two related CTest passes (runtime closure skipped by its build guard), and the real E4 A/B pass. One producer checkpoint with three non-zero hills forks into pure legacy, H5 `protocol` plus legacy NHC, and H5 `full` routes. All mdout columns, particle/NHC/metadynamics streams, final hill/potential state, observable output, and restart state agree within existing text-precision tolerances. Registry coverage advances to 60 supported, 14 deferred, and 1 unsupported contract. |
| PR 20: Unrestricted QC observable and SCF text gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `rerun_qc_unrestricted_sidecar_vds_off/on` CPU A/B; Ruff; `git diff --check` | Triplet unrestricted QC produces two finite non-zero spin-square frames in both branches and archives a non-empty exact SCF trace in trajectory and observable output under both VDS modes. Empty-output evidence was removed from the old full-contract cases. `input.qc.type` remains deferred because pure `/qc/type` is not consumed. Registry coverage advances to 62 supported, 12 deferred, and 1 unsupported contract. |
| PR 21: Constraint sidecar projection gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_constraint_sidecar_projection` CPU A/B; Ruff; `git diff --check` | 97 contract/manifest/comparator tests and the real four-frame A/B pass. Both branches preserve distance `1.5` with zero measured residual and reduce relative radial speed from `2.0` to `3.411315e-6`; complete position and velocity trajectories are identical. The supported sidecar contract is separated from deferred native typed `/constraint` consumption. Registry coverage advances to 63 supported, 12 deferred, and 1 unsupported contract. |
| PR 22: Native improper runtime gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; `ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R 'test_topology_native_h5_reader\|test_h5_input_fixture_equivalence'`; real `normal_improper_native_nonzero` CPU A/B; Ruff; `git diff --check` | 99 contract/manifest/comparator tests, both related H5 input CTests, and the focused real A/B pass. Pure typed `/forcefield/improper/{atoms,pk,phi0}` with no sidecars matches canonical legacy input at `improper_dihedral=31.36`, maximum force `35.415924072265625`, and zero full-mdout error. End-to-end conversion remains deferred for the documented `/k` versus `/pk` and sidecar-key gaps. Registry coverage advances to 64 supported, 12 deferred, and 1 unsupported contract. |
| PR 23: GB native-state/sidecar-activation gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_gb_hybrid_nonzero` CPU A/B; Ruff; `git diff --check` | 101 contract/manifest/comparator tests and the focused real A/B pass. Typed GB state plus the isolated H5 activation binding matches legacy at `gb=-0.25`, total potential `-0.75`, and maximum force `0.10313020646572113`; the gate rejects the Coulomb-only `0.25` force. Pure native GB without the activation sidecar remains deferred. Registry coverage advances to 65 supported, 12 deferred, and 1 unsupported contract. |
| PR 24: SW sidecar pair/three-body behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_sw_sidecar_pair_three_body` CPU A/B; Ruff; `git diff --check` | The isolated H5 `SW_in_file` sidecar route matches legacy at non-zero `SW=194.50` and the complete nine-value force payload. A `lambda=0` control changes the energy to `158.79` and the force fingerprint, so the gate proves both pair and three-body behavior rather than initialization alone. Typed `/manybody/sw` consumption remains deferred. Registry coverage advances to 66 supported, 12 deferred, and 1 unsupported contract. |
| PR 25: Tersoff sidecar angular behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_tersoff_sidecar_angular` CPU A/B; Ruff; `git diff --check` | The isolated H5 `TERSOFF_in_file` sidecar route matches legacy at `potential=-173.23`, `eff_pot=-173.23468`, and the complete nine-value force payload. A `gamma=0` control changes the potential to `-196.06` and the force fingerprint, proving angular bond-order behavior rather than initialization alone. Typed `/manybody/tersoff` consumption remains deferred. Registry coverage advances to 67 supported, 12 deferred, and 1 unsupported contract. |
| PR 26: H5 topology metadata runtime-rejection gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; three real bundled CPU F1 cases; Ruff; `git diff --check` | Production process gates reject topology atom-count mismatch with exit 238/`spongeErrorValueErrorCommand`, and reject non-vector or string `/atoms/mass` payloads with exit 234/`spongeErrorBadFileFormat`, retaining route-specific stable diagnostic tokens. An incompatible `/schema/version` value is currently accepted, so the complete four-mutation metadata contract remains deferred. Registry coverage advances to 68 supported, 12 deferred, and 1 unsupported contract. |
| PR 27: Steering CV protocol-sidecar behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_steering_cv_sidecar_nonzero` CPU A/B; Ruff; `git diff --check` | An isolated distance CV at `1.5` with weight `2.0` produces `steer_cv=3.0` and the analytic six-value force `(2,0,0,-2,0,0)` through both legacy `cv_in_file` and the H5 protocol sidecar route. The gate removes typed `/cv`, rejects zero-weight energy and force mutations, and keeps the unconsumed `steer_cv_in_file`/typed `/steering` contract deferred. Registry coverage advances to 69 supported, 12 deferred, and 1 unsupported contract. |
| PR 28: SITS typed Nk restart behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_sits_nk_typed_restart_nonzero` CPU A/B; Ruff; `git diff --check` | A two-temperature selective LJ/Coulomb SITS fixture consumes legacy `SITS_nk_in_file` with `Nk=(1,4)` versus typed `/parameters/restart/bias/sits/SITS/nk` under protocol restart load. Both routes produce `SITS_AA_kAB=-1.22`, `SITS_bias=-0.5317`, `SITS_fb=0.7049`, `eff_pot=-1.2829471`, and the same non-zero six-value force with maximum magnitude `0.2489061951637268`; the gate removes embedded/external bundled Nk text and rejects initialization-only output, unscaled force, and the symmetric `Nk=(1,1)` control. Typed `/sits` configuration materialization remains deferred. Registry coverage advances to 70 supported, 12 deferred, and 1 unsupported contract. |
| PR 29: VDS complete-prefix repair behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_vds_complete_prefix_noop` CPU A/B; `SPONGE_H5_ENABLE_RUNTIME_SMOKE=1 ctest --test-dir build-dev-cpu-h5-2 --output-on-failure -R '^test_h5_vds_terminal_resume_smoke$'`; Ruff; `git diff --check` | Legacy-input and bundled-input production routes each retain five deterministic frames in two complete shards with `repair_status=not_applied` and zero repaired shards. The same gate runs a real HighFive helper that injects terminal-shard finalize failure and verifies one-shard complete-prefix repair with rewritten completion metadata. Five metadata mutations are rejected. Cross-process reopen/append remains unsupported. Registry coverage advances to 71 supported, 11 deferred, and 1 unsupported contract. |
| PR 30: Residue topology-sidecar runtime behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_residue_sidecar_pbc_mapping` CPU A/B; Ruff; `git diff --check` | An isolated H5 `residue_in_file` topology-sidecar route matches the legacy `[2,2]` residue partition through runtime domain state, non-zero `bond=2.0`, the complete 12-value force payload with maximum magnitude `4.0`, and the cross-PBC whole-molecule position fingerprint `(19,21,25,28)`. Partition and mapping mutations are rejected. Typed `/atoms/residue_index` remains deferred because it is not materialized into `Xponge::system.residues`. Registry coverage advances to 72 supported, 11 deferred, and 1 unsupported contract. |
| PR 31: Residue membership/COM-virial follow-up gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (119 tests); 79-test production manifest; real `normal_residue_sidecar_pbc_mapping` and `normal_residue_sidecar_com_res_virial` CPU A/B; SPONGE rebuild; Ruff; clang-format; `git diff --check` | PR 30's residue-count, runtime-partition, force, and whole-molecule PBC mapping gate remains active. The second isolated `[2,2]` residue-sidecar case drives `com_res` restraint virial behavior with equivalent `bond=2.0`, `restrain=2.0`, `pressure=0.04`, `Pxx=0.11`, and 12-value force output. A same-count `[1,3]` control preserves energy and force but changes pressure/Pxx to `-11.53/-34.60`, proving membership-sensitive consumption. CPU atom-to-group maps now own their storage instead of aliasing freed host buffers. Typed residue input remains deferred. Registry coverage stays at 72 supported, 11 deferred, and 1 unsupported contract. |
| PR 32: Pure native GB behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (120 tests); 80-test production manifest; real `normal_gb_native_nonzero` and `normal_gb_hybrid_nonzero` CPU A/B; two native topology reader CTests; SPONGE rebuild; Ruff; changed-line clang-format; `git diff --check` | The pure bundled branch consumes `/forcefield/gb/params` with no `gb_in_file`, sidecar table, or sidecar files. Both routes match `Coulomb=-0.50`, `gb=-0.25`, `potential=-0.75`, the complete six-value force oracle, and all mdout columns. Non-periodic GB initialization now recognizes preloaded native state while preserving the hybrid activation route. Registry coverage advances to 73 supported, 10 deferred, and 1 unsupported contract. |
| PR 33: Core topology payload-sensitivity gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (128 tests); 88-test production manifest; real `normal_core_topology_payload_sensitivity` CPU A/B plus three typed-payload controls; two native topology reader CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Pure typed mass, charge, and LJ routes contain no topology sidecars and match legacy across the complete one-frame mdout and force payload. Independent controls produce temperature `26.21 -> 52.42` with unchanged force, Coulomb `-0.50 -> 0` with maximum force delta `0.25`, and LJ `-0.02 -> 0` with maximum force delta `0.04541015625`, proving each payload reaches its consumer. The broader `normal_core_h5_output` case remains active and still reports its pre-existing legacy/H5 force schedule mismatch (`51` versus `50` frames), which is not weakened by this input gate. Registry coverage remains 73 supported, 10 deferred, and 1 unsupported contract. |
| PR 34: Complete H5 metadata failure semantics | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (130 tests); 90-test production manifest; four real `failure_h5_topology_*` CPU process gates; three supported-schema real controls; three metadata/input CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Atom-count, mass-shape, mass-dtype, and unknown-schema mutations now share one complete F1 contract. Schema versions `0`, `1`, and `xponge.legacy_to_bundle.v1` each start successfully, while `unsupported.topology.v999` exits 238 with `spongeErrorValueErrorCommand` and route-specific diagnostics. The partial metadata registry record is removed; coverage remains 73 supported, advances to 9 deferred, and retains 1 unsupported contract. |
| PR 35: Subsystem-division energy-partition behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (130 tests); 90-test production manifest; real `normal_lj_soft_core_nonzero` and `normal_subsystem_division_partition` CPU A/B plus all-intra legacy/bundled controls; two native topology reader CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | The legacy `[0,1]` sidecar and pure typed bundled dataset both produce `LJ_soft_inter=-0.06`, `LJ_soft_intra=0.00`, and the same complete force. Same-semantic `[0,0]` controls move the full `-0.06` to the intra observable while preserving total soft-core energy and force, proving mask consumption without changing dynamics. Registry coverage advances to 74 supported, 8 deferred, and 1 unsupported contract. |
| PR 36: Typed residue runtime-behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (131 tests); 91-test production manifest; real typed PBC and COM-virial CPU A/B plus both sidecar residue regressions; two native topology reader CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Pure typed `/atoms/residue_index=[0,0,1,1]` and `/residues/atom_offset=[0,2,4]` materialize into the native `[2,2]` runtime partition without a topology sidecar. They match legacy at `bond=2.0`, the whole-molecule PBC fingerprint `(19,21,25,28)`, `restrain=2.0`, `pressure=0.04`, and `Pxx=0.11`. A typed `[1,3]` control preserves the complete force and energy while changing pressure/Pxx to `-11.53/-34.60`; inconsistent typed datasets and legacy/typed ownership conflicts are rejected. Registry coverage advances to 75 supported, 7 deferred, and 1 unsupported contract. |
| PR 37: Typed SW pair/three-body behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract`; real `normal_sw_typed_pair_three_body` and `normal_sw_sidecar_pair_three_body` CPU A/B cases; related native topology CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Pure typed `/manybody/sw` with no SW command or sidecars matches the legacy pair/three-body route at `SW=194.50` and maximum force `404.27862548828125`. A typed lambda-zero control changes energy to `158.79` and maximum force to `343.52105712890625`, proving that the typed three-body payload reaches the force kernel. Registry coverage advances to 76 supported, 6 deferred, and 1 unsupported contract. |
| PR 38: Typed constraint projection behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (133 tests); 93-test production manifest; real `normal_constraint_typed_projection` and `normal_constraint_sidecar_projection` CPU A/B cases; input-validation and fixture-equivalence CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Pure typed pair `[0,1]` with `r0=1.5` and no constraint sidecar matches the legacy route across four complete position/velocity frames with zero distance residual and radial velocity residual `3.411315e-6`. A typed `r0=2.0` control follows the new target with residual `0`; out-of-range pair `[0,2]` is rejected as `spongeErrorBadFileFormat`. Sidecar precedence remains unchanged. Registry coverage advances to 77 supported, 5 deferred, and 1 unsupported contract. |
| PR 39: Typed Tersoff angular behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (134 tests); 94-test production manifest; real `normal_tersoff_typed_angular` and `normal_tersoff_sidecar_angular` CPU A/B cases plus two typed failure controls; two native input CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Pure typed `/manybody/tersoff` with no Tersoff command or sidecars matches legacy at `potential=-173.23`, `eff_pot=-173.23468`, and maximum force `135.94906616210938`. A typed gamma-zero control produces `potential=-196.06`, `eff_pot=-196.05984`, and maximum force `144.7831268310547`, proving angular parameter consumption. Conflicting map/entry and raw/runtime parameter payloads are rejected as bad-file-format. Registry coverage advances to 78 supported, 4 deferred, and 1 unsupported contract. |
| PR 40: Typed QC type-sensitive behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (135 tests); 95-test production manifest; real `rerun_qc_type_typed_unrestricted_vds_off` CPU A/B plus singlet control; two native input CTests; SPONGE rebuild; Ruff; changed-line clang-format; `git diff --check` | Typed `/qc/type` with no QC sidecar matches legacy across two complete mdout frames, non-trivial `QC`/`QC_S_sq`, and exact SCF text. Changing only multiplicity `3 -> 1` changes QC by up to `89.27` and spin square by `2.0065`, proving type-sensitive SCF consumption. Registry coverage advances to 79 supported, 3 deferred, and 1 unsupported contract. |
| PR 41: Typed steering CV behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (136 tests); 96-test production manifest; real `normal_steering_cv_typed_nonzero` and `normal_steering_cv_sidecar_nonzero` CPU A/B cases plus typed weight-zero control; two native input CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Pure typed `/cv/config` with no CV command or sidecars matches the legacy `cv_in_file` route at `steer_cv=3.0` and analytic force `(2,0,0,-2,0,0)`. Changing only weight `2.0 -> 0.0` zeros energy and force with maximum force delta `2.0`, proving typed config consumption. The registry now reflects the real `/cv/config` surface instead of the unconsumed `/steer` conversion. Coverage advances to 80 supported, 2 deferred, and 1 unsupported contract. |
| PR 42: Typed SITS configuration/selection behavior gate | Complete | This commit | `pixi run -e dev-cpu smoke-bundled-io-contract` (138 tests); 98-test production manifest; real `normal_sits_typed_configuration_nonzero` and `normal_sits_nk_typed_restart_nonzero` CPU A/B cases plus typed config/selection/precedence controls; two native input CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | Pure typed `/sits/config`, `/sits/atom_indices`, and Nk restart state with no inline SITS owner match legacy at `SITS_AA_kAB=-1.22`, `SITS_bias=-0.5317`, `SITS_fb=0.7049`, and maximum force `0.2489061951637268`. Typed `pe_a=0.5` and atom set `[0]` produce distinct deterministic bias, factor, energy, and force fingerprints; complete inline config wins over typed config without suppressing typed atoms. Duplicate atom and duplicate config-key payloads exit as bad-file-format. Registry coverage advances to 81 supported, 1 deferred, and 1 unsupported contract. |
| PR 43: End-to-end improper conversion behavior gate | Complete | XPONGE `45e52be`, `b890fbf`; SPONGE this commit | XPONGE 93-test bundle regression; `pixi run -e dev-cpu smoke-bundled-io-contract` (138 tests); 103-test registry/manifest suite; real canonical, alias, and pure-native improper CPU A/B cases plus typed `pk=5` controls; two native input CTests; SPONGE rebuild; Ruff; clang-format; `git diff --check` | XPONGE emits native `/forcefield/improper/pk`, retains `/k` reverse-export fallback, and omits improper sidecars for typed conversion. SPONGE consumes the declared `improper_in_file` alias when canonical input is absent. Unmodified canonical and alias conversions match legacy at `improper_dihedral=31.36` and maximum force `35.415924072265625`; independent `pk=5` controls halve energy and every force component. Coverage reaches 82 supported, 0 deferred, and 1 explicitly unsupported contract. |
| PR 44: Execution-matrix residue ownership fixture repair | Complete | XPONGE `2c0dc3e`; SPONGE this commit | XPONGE 101-test bundle regression; pinned-fixture hash/ownership gate; CPU and GPU deterministic matrix A/B; four-replica GPU Middle+SETTLE statistical A/B; Ruff; `git diff --check` | XPONGE now emits typed residue state without an active residue sidecar and reverse-exports mutated typed offsets. The regenerated TIP3P fixture has exactly one residue owner. Complete CPU/GPU deterministic and GPU statistical behavior comparisons pass; the former conflicting-owner startup failure is closed. |
| PR 45: CUDA+MPI NCCL discovery prerequisite | Complete | This commit | CUDA 13 + OpenMPI + NCCL configure; 107-step CUDA `SPONGE` target build; `ldd`; cmake-format; `git diff --check` | The project CMake module path can now resolve an available NCCL installation through `NCCL_ROOT`. The real CUDA+MPI executable links HDF5, CUDA 13, OpenMPI, and NCCL without unresolved runtime libraries, enabling the GPU rank-2 behavior gate. |
| PR 46: Portable GPU MPI device mapping | Complete | This commit | shared/explicit/rank-1 device-map mutation test; real four-replica CUDA+MPI rank-2 NPT+SETTLE statistical A/B; 139-test smoke; Ruff; `git diff --check` | GPU rank-2 defaults to a portable single-GPU shared map and accepts an explicit per-rank multi-GPU map. The real case passes complete mdout and all three H5 families with `mpi_rank_count=2` and rank-0-only output ownership. |
| PR 47: Production-fixture residue ownership repair | Complete | This commit | four focused residue CPU A/B cases; full-contract sidecar VDS-off/on; 24-case sidecar QC/rerun/failure sweep; CMAP VDS invariant; 140-test smoke; H5 manifest and fixture-equivalence CTests; Ruff; format; `git diff --check` | Full-contract bundles now have exactly one typed residue owner, while focused sidecar cases explicitly replace typed state with one bundle-local sidecar owner. The stale startup conflict is closed across compatibility paths. Independent typed-SITS and rerun boundary behavior failures remain visible for subsequent alignment PRs. |
| PR 48: Inactive typed SITS configuration semantics | Complete | This commit | focused active and inactive typed SITS CPU A/B; full-contract pure-native progression check; SPONGE rebuild; 141-test smoke; Ruff; format; clang-format; `git diff --check` | Typed SITS config without `mode` now preserves the valid legacy inactive state: exact text is materialized, no SITS columns appear, and complete mdout is equal. The existing production-mode case remains non-trivial. Full-contract pure-native execution advances to independent NB14/EAM/restraint/soft-wall consumption gaps. |
| PR 49: Canonical and extra NB14 behavior contracts | Complete | This commit | canonical/extra focused CPU A/B plus zero-LJ controls; pure full-contract VDS-off/on progression; sidecar full-contract VDS-off/on; dual-root and extra-shape process failures; 21-case failure sweep; native topology reader and fixture-equivalence CTests; SPONGE rebuild; 142-test smoke; Ruff; format; clang-format; `git diff --check` | The registry now distinguishes `nb14_in_file` scaling from raw `nb14_extra_in_file` A/B semantics. A small independently owned H5 helper preserves the untracked topology-reader boundary. Both pure typed routes match exact non-zero energies and complete force fingerprints, parameter-zero controls change force, and full-contract fixtures enforce exactly one NB14 owner. Pure full-contract NB14 columns are restored; EAM/restraint/soft-wall remain visible. |
| PR 50: Typed EAM funcfl/setfl behavior contracts | Complete | This commit | focused funcfl/setfl CPU A/B plus zero-pair controls; pure full-contract VDS-off/on progression; sidecar full-contract VDS-off/on; unknown-format and table-shape process failures; complete failure sweep; fixture-equivalence CTest; SPONGE rebuild; contract smoke; Ruff; format; clang-format; `git diff --check` | Typed EAM input now materializes both legacy formats and optional atom types through the established EAM runtime parser. Exact energy and full-force oracles plus zero-pair controls prove table consumption, full-contract fixtures enforce exactly one EAM owner, and pure rerun now advances to only restraint/soft-wall gaps. |
| PR 51: Typed positional-restraint behavior contract | Complete | This commit | focused typed positional-restraint CPU A/B plus three independent payload controls; pure full-contract VDS-off/on progression; sidecar VDS-off/on; two residue COM/virial regressions; dual-owner and weight-shape process failures; 25-case failure sweep; fixture-equivalence CTest; SPONGE rebuild; 142-test smoke; Ruff; format; clang-format; `git diff --check` | Protocol atom IDs/weights and restart reference coordinates now form one validated typed owner and materialize through the established positional-restraint parser. Exact energy/full-force oracles prove all three payload fields are consumed. The registry separates typed and sidecar evidence, retains CV restraint and typed soft wall as deferred, and pure rerun advances to only the soft-wall gap. |
| PR 52: Typed soft-wall behavior contract | Complete | This commit | focused typed soft-wall CPU A/B plus zero-potential control; pure and sidecar full-contract VDS-off/on; dual-owner and three shape process failures; 29-case failure sweep; fixture-equivalence CTest; SPONGE rebuild; 142-test smoke; Ruff; format; clang-format; `git diff --check` | Typed `/wall/soft/{count,name,potential}` now validates and materializes through the established soft-wall parser. The focused gate matches `z_wall=10.0` and force `(0,0,-2,0,0,-6)` while the typed zero-potential control zeros both. Full-contract fixtures enforce exactly one owner, and pure VDS-off/on now pass the complete rerun behavior gate. Registry coverage reaches 85 supported, one deferred, and one unsupported contract. |
| PR 53: Typed and sidecar CV-restraint behavior contracts | Complete | This commit | focused typed/sidecar CV-restraint CPU A/B plus weight-zero and reordered-definition controls; pure and sidecar full-contract VDS-off/on E3 gates; dual-owner, partial-owner, offset, and definition-conflict process failures; 33-case failure sweep; fixture-equivalence CTest; SPONGE rebuild; 142-test smoke; Ruff; format; clang-format; `git diff --check` | The CV controller now merges all three legacy configuration surfaces, while typed protocol roots pass one structured, conflict-aware materializer. Focused routes match `restrain_cv=1.0` and force `(4,0,0,-4,0,0)`; full-contract routes exercise `restrain_cv=12.06`. Registry coverage reaches 87 supported, zero deferred, and one unsupported contract. |
| PR 54: Attested production promotion evidence | Complete | This commit | eight derived-evidence positive/negative guards; real 24-test comparator mutation hook; 150-test contract smoke; workflow YAML parse; execution-matrix CLI; Ruff; format; `git diff --check` | Production history is now derived from same-run contract, split CPU/GPU/MPI matrix, and comparator artifacts. Scheduled CPU jobs share one run ID and upload the comparator report. History records conservative scenario performance, source commit, and three SHA-256 provenance values; stale/conflicting IDs, incomplete evidence, failed mutation sessions, non-finite metrics, duplicates, and unattested handwritten rows fail. Promotion remains shadow pending three complete retry-free runs and cross-process VDS closure. |
| PR 55: Clean source provenance gate | Complete | This commit | clean/mismatched/modified/untracked temporary-Git guards; real dirty-tree CLI rejection; 152-test contract smoke; Ruff; format; `git diff --check` | Production derivation now requires exact `HEAD`, an empty porcelain status including untracked files, and persisted `source_tree_state=clean`. The current worktree is correctly rejected before artifact loading, so its available RTX 4090 cannot be used to claim a commit-qualified PR 7c streak until the bundled I/O source baseline is tracked. |
| PR 56: Track bundled I/O implementation baseline | Complete | This commit | fresh CPU `SPONGE` build; staged H5 bundle CTest labels; opt-in runtime smoke audit; 152-test contract smoke; schema examples; Ruff; `pixi lock --check`; `git diff --check` | The audited bundled I/O implementation, schemas, fixtures, examples, and tests are now tracked. Static/file-level gates pass; opt-in normal A/B and VDS terminal/resume pass. EDIP remains all-trivial in the manybody gate and dynamic NHC restart-load exits after initialization, so full behavioral equivalence and PR 7c promotion remain unproven. |
| PR 57: Non-trivial REAXFF/EDIP runtime parity | Complete | This commit | focused production `normal_edip_nonzero` CPU A/B; runtime `test_h5_reaxff_edip_runtime_parity`; staged H5 bundle labels; formatting; `git diff --check` | The combined gate now gives legacy binary, bundled sidecar, and pure H5 routes the same interacting two-atom frames and empty pair exclusions. EDIP and REAXFF remain subject to finite non-zero assertions and full mdout/H5 observable parity across legacy output and bundled VDS off/on output. The independent dynamic NHC restart-load failure remains open. |
| PR 58: Dynamic NHC restart-load runtime closure | Complete | This commit | focused runtime `test_h5_restart_load_runtime_closure`; complete four-test opt-in runtime group; staged H5 bundle labels; 152-test contract smoke; formatting; `git diff --check` | The dynamic NHC case now uses an isolated normal-NVT fixture instead of entering unrelated QC and rerun modules. It executes one finite `dt=0` step, confirms NHC initialization, rejects non-finite mdout/restart values, and requires exact non-zero NHC-state roundtrip. Existing protocol/full success and unsupported-policy failure branches remain on their broad fixtures. |
| PR 59: Structural restart producer-source closure | Complete | This commit | 154-test contract smoke; real CUDA `normal_structural_restart_continuation` E4 A/B; execution-matrix CLI; Ruff; `git diff --check` | Independent legacy-input and bundled-input producers each run one non-zero step and emit their own restart. The two-step non-zero-dt NVE continuations bind only the branch-owned text or H5 source, record source directory/payload/SHA-256, and compare all mdout columns plus trajectory position/velocity/force/box and final structural restart state. A source-swap mutation is rejected before launch. `input.restart_load.structural` now requires E4 from this focused case; cross-process VDS reopen-and-append/resume remains unsupported. |
| PR 60: Dynamic restart contract enumeration | Complete | This commit | 155-test contract smoke; 13 staged H5 contract tests plus staged smoke/backend groups; four-test runtime smoke; `test_h5_input_validation`; `test_h5_dynamic_state_modules`; real CUDA NHC E4 A/B; execution-matrix CLI; Ruff; clang-format; `git diff --check` | Dynamic restart evidence is split into integrator, NHC, Bussi, pressure-barostat, Middle, Andersen, and MC-barostat contracts. Integrator and NHC are supported/E4; Bussi and pressure barostat remain deferred/E4 until PR 61/62; three non-restorable RNG states are unsupported/E4 with module-specific reader diagnostics. Coverage reports supported/deferred/unsupported layers separately, and an exact inventory mutation rejects a missing module contract. Registry status is 89 supported, 2 deferred, and 4 unsupported. |
| PR 61: Bussi dynamic restart continuation closure | Complete | This commit | 156-test contract smoke; staged H5 contract/runtime groups; `test_h5_dynamic_state_modules`; real CUDA `normal_bussi_dynamic_restart_continuation` E4 A/B plus three mutation controls; execution-matrix CLI; Ruff; clang-format; `git diff --check` | A continuous legacy-input Bussi NVT reference is compared with an independent bundled-input producer prefix and a two-step continuation loaded directly from the bundled producer H5 checkpoint. Checkpoint/final RNG text and lambda, continuous time, relative step progression, all mdout columns, and trajectory/restart position, velocity, force, and box are checked. Corrupt RNG and missing lambda fail with stable diagnostics; a valid alternate RNG state runs but changes final RNG and velocity. The CUDA-observed velocity replay envelope is recorded locally as `rtol=1e-5, atol=1.5e-4`; other deterministic tolerances are unchanged. Registry status is 90 supported, 1 deferred, and 4 unsupported. |
| PR 62: Pressure-based barostat dynamic restart continuation closure | Complete | This commit | 157-test contract smoke; staged H5 contract/runtime groups; `test_h5_dynamic_state_modules`; real CUDA `normal_pressure_barostat_dynamic_restart_continuation` E4 A/B plus five mutation controls; execution-matrix CLI; Ruff; clang-format; `git diff --check` | An Andersen pressure barostat consumes its RNG and six-value `g` state from an independent bundled-input checkpoint, then runs a two-step non-zero-dt NPT continuation against a continuous reference. The gate checks prefix parity, checkpoint/final RNG and `g`, continuous time, relative step progression, pressure/Pxx/Pyy/Pzz, particle state, box evolution, and final restart. Instantaneous pressure at the restart zeroth frame is explicitly excluded because virial is not restart state; both post-restart pressure rows use recorded local tolerances. Missing/corrupt/owner-mismatched state fails with stable diagnostics, while swapped `g` and alternate RNG change box and final `g`. NPT pressure output registration is de-duplicated so mdout/H5 observable columns remain parseable. Registry status is 91 supported, 0 deferred, and 4 unsupported. |
| PR 63: Bonded and custom-listed payload-sensitivity oracles | Complete | This commit | 158-test contract smoke; real CUDA `normal_{bond,angle,urey_bradley,dihedral,listed}_payload_sensitivity` E3 A/B cases plus five zero-stiffness controls; execution-matrix CLI; Ruff; `git diff --check` | Five isolated fixtures now give bond, angle, Urey-Bradley, dihedral, and custom-listed contracts independent analytic energy and complete force oracles. Each bundled branch owns the corresponding typed H5 payload without its legacy route or sidecars; custom listed definitions and parameters are materialized from H5 at runtime. Mutating only the target stiffness to zero must zero that module's energy and complete force, excluding default-value, wrong-field, and cross-module false positives. Full-contract reruns remain composition regressions rather than the sole E3 evidence. Registry status remains 91 supported, 0 deferred, and 4 unsupported. |
| PR 64: Remaining payload-sensitivity oracles | Complete | This commit | 159-test contract smoke; real CUDA positional-restraint sidecar, soft-wall sidecar, typed CV, typed QC, ReaxFF, and EDIP E3 A/B cases plus independent payload controls; `test_h5_input_matrix_contract`; execution-matrix CLI; Ruff; clang-format; `git diff --check` | Six remaining functional contracts now have focused evidence instead of relying only on broad full-contract reruns. Structural restart loading restores only the restraint reference sidecar needed by positional restraint and cannot override unrelated typed protocol owners. ReaxFF uses a periodic interacting O-H pair; CV mutates only its definition; QC mutates only the typed element symbol and compares the nondeterministic GPU SCF trace by identical step/iteration structure plus bounded finite values. Registry status remains 91 supported, 0 deferred, and 4 unsupported. |
| PR 65: Nonzero-dt input evolution gates | Complete | This commit | 161-test contract smoke; real CUDA SITS, EDIP, typed-residue, and sidecar-residue three-frame nonzero-dt A/B evolution; frozen-coordinate and schedule mutations; execution-matrix CLI; Ruff; `git diff --check` | The registry/case audit now derives input-only, single-step, and zero-dt sets instead of maintaining counts by hand. Of 23 current input-only supported contracts, five are exercised through real nonzero-dt evolution spanning topology/PBC, manybody/custom, and protocol/stateful cohorts; the other 18 have reviewable registry reasons tied to independent frame-local payload mutations. Every evolution compares three mdout rows and position/velocity/force/box schedules, requires finite outputs and coordinate change, and rejects a frozen trajectory or altered schedule. Registry status remains 91 supported, 0 deferred, and 4 unsupported. |
| PR 66: Process output-family combinations | Complete | This commit | 162-test contract smoke; all 14 real CUDA output-family A/B cases; execution-matrix CLI; Ruff; `git diff --check` | The process harness now exercises every non-empty combination of trajectory, observable, and restart writers under both default legacy suppression and explicit legacy coexistence. Each branch requires exactly the enabled H5 files, no disabled-writer files or temporary artifacts, readable finalized outputs, and matching sidecar provenance; restart-only additionally reloads the bundled-produced checkpoint in a structural continuation. A manifest mutation guard rejects any missing combination or legacy mode. VDS remains off in this matrix and cross-process VDS append/resume remains explicitly unsupported. Registry status advances to 92 supported, 0 deferred, and 4 unsupported. |
| PR 67: Writer finalize failure isolation | Complete | This commit | 163-test contract smoke; six real CUDA and CPU append/finalize fault processes; all 14 CUDA and CPU normal output-family regressions; `test_h5_output_plan`, `test_h5md_writers_with_mock_backend`, and `test_highfive_backend_io`; execution-matrix CLI; Ruff; clang-format; `git diff --check` | Test builds expose a compile-gated `SPONGE_H5_TEST_FAULT` hook for trajectory, observable, and restart append/finalize points; production builds do not compile the hook. A failed writer persists `failed` plus a stable reason, closes without leaving a temporary artifact, and records a deferred process failure. `Main_Clear` finalizes every still-healthy writer before emitting one nonzero `spongeErrorValueErrorCommand`, so the two unaffected families remain readable, `finalized`, and completion-metadata consistent. The six-point manifest and process oracles reject swallowed failures, missing status/error metadata, or skipped peer finalization. Registry status advances to 93 supported, 0 deferred, and 4 unsupported. |
| PR 68: Rerun velocity/force/optional-field direct closure | Complete | This commit | 165-test contract smoke; real CPU and CUDA rerun direct-payload A/B cases covering VDS off/on, velocity/force present/absent, non-default particle stream, and exact-EOF box update; payload/schedule/presence mutation guards; `test_h5_output_plan` and `test_trajectory_h5_reader`; execution-matrix CLI; Ruff; clang-format; `git diff --check` | Rerun inputs and legacy outputs now use distinct paths. The direct legacy-sidecar-to-H5 bridge compares position, box, velocity, and force presence, shape, frame count, step/time schedule, and values, and records all four payload decisions in assertion evidence. H5 rerun now zeroes absent velocity, preserves selected-frame numbering, writes the loaded box, gathers distributed force before H5 append, and keeps force scheduling aligned with trajectory scheduling. Two broad CUDA full-contract regressions remain rejected before this bridge by their existing QC mdout differences; no comparator tolerance was relaxed. Registry status remains 93 supported, 0 deferred, and 4 unsupported. |
| PR 69: Normal-success finite-output gate | Complete | This commit | 170-test contract smoke; real CPU focused-EDIP, full-contract rerun, and trajectory-only output-family cases; real CUDA focused-EDIP, trajectory-only output-family, and repaired full-contract rerun; matching NaN/+Inf/-Inf mutation guards; execution-matrix CLI; Ruff; `git diff --check` | Every successful SPONGE process now scans all mdout values, numeric trajectory/observable/restart H5 datasets (including VDS shards), and legacy coordinate/velocity/force/box outputs before equivalence or oracle evaluation. Each successful evidence record carries the `normal_success` classification plus run-directory, file, dataset, and value counts; F1 cases are separately classified, and non-finite propagation cannot attach supported-functionality evidence. The gate exposed the full-contract rerun soft-wall singularity, so temporary legacy/H5 second-frame inputs are now synchronized away from that exact singular coordinate rather than whitelisting matching infinities. CUDA structural and Bussi attempts remained rejected later by their existing finite force-tolerance differences; no tolerance was relaxed. Registry status remains 93 supported, 0 deferred, and 4 unsupported. |
| PR 70: Oracle-aware evidence and layered coverage | Complete | This commit | 175-test contract smoke; real CPU and elevated CUDA `normal_core_topology_payload_sensitivity` E3 A/B; forged assertion, empty details, wrong-contract oracle, handwritten pass, and E0-as-E3 mutation guards; 12 execution-matrix tests collected/skipped outside the opt-in runtime; Ruff; `git diff --check` | Passing records are now revalidated against contract-owned assertions, minimum evidence levels, non-empty details, and contract-scoped input oracles before they enter coverage or promotion. A shared input assertion emits one oracle/result per contract with the contract ID, control mutation, expected/actual delta, and compared payload. Reports persist an explicit registry summary plus independent supported denominators for inventory (1), conversion (0), runtime behavior (77), continuation (10), and failure semantics (5); E0 inventory can no longer inflate E3/E4/F1 coverage. Promotion rejects malformed or handwritten records. Registry status remains 93 supported, 0 deferred, and 4 unsupported; cross-process VDS append/resume remains excluded. |
| PR 71: Feature × CPU/MPI execution matrix | Complete | This commit | 176-test contract smoke; 26-scenario matrix validation; real CPU rank-1 matrix covering EDIP, ReaxFF, typed QC, SITS with VDS, metadynamics E4 continuation, typed positional/CV restraint, soft wall, virtual atoms, and typed constraint; real CPU-MPI rank-2 EDIP, three-step nonzero-dt SITS, virtual-atom PBC, and constraint-sidecar cases; feature/VDS metadata and required-combination mutations; execution-matrix CLI; Ruff; `git diff --check` | The matrix now carries explicit `feature_family` and `vds` axes instead of proving only TIP3P environment presence. Every required CPU rank-1 feature and the manybody/SITS/virtual-atom/constraint rank-2 risk set has a concrete source case, rank-0 ownership check, finite performance metadata, and module-owned mdout/H5 evidence. Rank-2 SITS uses deterministic module-owned mdout/observable comparison over three nonzero-dt frames; stochastic particle-velocity comparison remains covered by the core environment matrix rather than being silently tolerance-widened. The metadynamics run exposed an ASCII NHC `.crd` sidecar, so the finite scanner now distinguishes numeric text from binary float32 by content. GPU feature crossings remain for PR 72, and cross-process VDS append/resume remains excluded. |
| PR 72: Formal GPU shadow evidence workflow | Complete | This commit | 179-test contract smoke; 39-scenario matrix validation; elevated RTX 4090 GPU rank-1 feature matrix (10 cases) and GPU-MPI rank-2 risk matrix (3 cases); CUDA+MPI Pixi task dry-run; lock check; workflow manifest/provenance mutation guards; Ruff; `git diff --check` | Dedicated self-hosted CUDA rank-1 and CUDA-MPI rank-2 jobs now share the production run ID and source attestation with CPU, contract, and comparator artifacts. GPU evidence records device map, model, driver/runtime, ranks, OMP, and rank-0 ownership; merge rejects absent GPU artifacts, dirty or mixed commits, missing hardware metadata, skipped scenarios, and conflicting duplicate records. The feature matrix crosses all rank-1 feature families plus rank-2 manybody, virtual-atom, and constraint risks. Promotion remains shadow-only pending three clean retry-free production runs; cross-process VDS append/resume remains explicitly excluded. |
| PR 72a: Pinned shadow converter toolchain | Complete | This commit | pinned-base patch apply/check and byte-identical runtime-source hashes; converter CLI smoke under CUDA Python; 179-test contract smoke; workflow YAML parse; Ruff; `git diff --check` | Every shadow job checks out remote XPONGE commit `d607b89` and applies one reviewed runtime patch whose SHA-256 is locked by the workflow manifest test. This removes dependence on an unpushed sibling worktree and makes `legacy-to-bundle` available under CPU, MPI, GPU, and GPU-MPI Pixi Python. The immutable PR 72 commit passes all 10 GPU rank-1 feature crossings and all three GPU-MPI rank-2 risk crossings with one run ID; the full GPU selector separately exposed two existing core statistical blockers, which remain failures rather than tolerance changes. |
| PR 73: VDS resume scope boundary | Complete | This commit | 189-test contract smoke; scope-exclusion/report/history positive and mutation guards; same-process-repair evidence forgery rejection; execution-matrix CLI; Ruff; `git diff --check` | Contract evidence and production summaries now carry the exact unsupported-contract inventory plus one fixed scoped-equivalence statement. Registry validation keeps cross-process VDS resume unsupported/E4 with no case or assertion, promotion rejects missing/stale exclusions, unbounded claims, silently supported VDS resume, and run-history rows without the boundary. Existing same-process shard, tail-repair, and complete-prefix no-op evidence cannot satisfy this contract. No writer append/open mode or cross-process runtime case was added. |
| PR 73b: Independent-replica statistical matrix | Complete | This commit | 199-test contract smoke; 48-replica RTX 4090 Middle and NHC production GPU A/B plus retained-artifact comparator replay; correlated-block, sample-plan, H5-family, feature-boundary, and std-threshold mutations; execution-matrix CLI; Ruff; `git diff --check` | Production statistical evidence now estimates SEM across at least 48 independent replica means (384 raw paired blocks) instead of treating within-run blocks or flat H5 components as independent. Structured trajectory comparison uses 4096-pair PBC and global physical distributions; per-atom mutation sensitivity remains in the lower-level gate. Single-frame restart is excluded from statistical H5 families and remains covered by deterministic/E4 gates. Both previously blocked GPU core scenarios pass mdout plus trajectory/observable H5 comparison with the existing 3-sigma margins. Promotion remains shadow-only pending immutable clean-source matrices and three retry-free production runs. |
| PR 74: Scoped promotion candidate | Candidate | This commit | Static contract, matrix, workflow, comparator, provenance, scope-boundary, and declaration-only blocker gates before the immutable commit; production evidence is external and must be generated afterward. | The matrix enters `candidate`, but remains outside the release workflow and cannot become ready from declarations alone. This commit makes no production-pass claim: it may be merged only after three consecutive clean, retry-free runs on this exact commit derive attested history rows and the final evaluator returns ready. Cross-process VDS reopen-and-append/resume remains excluded. |
