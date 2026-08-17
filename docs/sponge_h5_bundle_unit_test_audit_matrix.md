# H5 Bundle Unit Test Audit Matrix

This audit maps the staged plan in
`tests/h5_bundle/INPUT_OUTPUT_MATRIX_TEST_PLAN.md` to concrete tests. It is a
coverage index, not a replacement for the tests themselves.

## Phase Coverage

| Plan item | Evidence | Default status |
|---|---|---|
| H5 bundle CTest target and build-only aggregate inventories stay synchronized with the human-readable target manifest. | `test_h5_test_targets_manifest.py` parses `tests/h5_bundle/CMakeLists.txt` and `tests/h5_bundle/TEST_TARGETS.md`, then requires every CTest target name, source file, label set, and `h5_bundle` label to match exactly. It also checks the `sponge_h5_bundle_*` build-only aggregate targets against the manifest's Build-only targets table, verifies explicit aggregate dependencies are real compiled CTest targets, and rejects Python-only CTest scripts as aggregate build dependencies. | Runs under `h5_bundle;contract`. |
| The legal matrix table in `INPUT_OUTPUT_MATRIX_TEST_PLAN.md` stays synchronized with both the static matrix test and runtime smoke inventory. | `test_h5_matrix_plan_manifest.py` parses the plan's `Minimum legal case set` table, `test_h5_io_matrix_spec.cpp::Legal_Matrix_Cases`, `test_h5_input_output_smoke_matrix.cpp::expected_runtime_case_names`, and the runtime case vectors themselves, then requires the case names and mode/input/output/VDS dimensions to match. It also checks that runtime execution still calls `Run_Legacy_Sidecar_Smoke_Cases` from `Run_Normal_Mode_Matrix`, calls `Run_Rerun_Frame_Selection_Smoke_Cases` from `Run_Rerun_Mode_Matrix`, binds those execution paths to `Sidecar_Smoke_Cases` and `Rerun_Selection_Cases`, and keeps mdout/H5 output assertions attached to every successful smoke family so process success alone is insufficient. | Runs under `h5_bundle;contract;matrix`. |
| CI labels, staged script commands, and the runtime smoke gate stay synchronized with the plan. | `test_h5_ci_plan_manifest.py` parses the plan's `CI Labels` block, CMake CTest labels, `run_h5_bundle_tests.sh` command dispatch and staged label order, `TEST_TARGETS.md` expected validation sequence, the smoke CTest properties, `test_h5_input_output_smoke_matrix.cpp`, `test_h5_matrix_plan_manifest.py`, and `test_h5_io_contract_manifest.py`; it requires every planned label to be registered, staged labels to have script commands, script `run_staged()` to match the documented validation sequence, `test-smoke-runtime` to set `SPONGE_H5_ENABLE_RUNTIME_SMOKE=1`, the runtime smoke CTest to use skip code 77, runtime smoke `main()` to validate preparation before checking the runtime gate and before running normal/rerun matrices, and the manifest sources to retain the runtime-dimension, branch-assertion, and semantic-equivalence guards. | Runs under `h5_bundle;contract;matrix`. |
| The plan's Implementation File Map remains a tested phase-to-file ownership index instead of passive documentation. | `test_h5_audit_matrix_manifest.py` parses the `Implementation File Map` table in `INPUT_OUTPUT_MATRIX_TEST_PLAN.md`, checks the exact Phase 0-7 order, requires each phase's primary files to stay listed, and requires the verification boundary keywords for fixture relocation, input sidecars, normal smoke, sidecar smoke, rerun/VDS, legal matrix, full-contract coverage, and native parity. It also self-checks the parser and guard names so this phase-to-file ownership test cannot be removed silently. | Runs under `h5_bundle;contract;coverage`. |
| Input fixture generation emits stable `legacy_input`, sidecar-stripped `bundled_input`, sidecar-preserving `bundled_input_with_legacy_sidecar`, manifest schema/version, required H5 payload files, sidecar materialization roots, and the core override fixtures. | `test_h5_io_contract_manifest.py::require_fixture_group_structure` checks both `core_structural` and `full_contract_rerun` fixture groups; it verifies required directories, mdin files, topology/protocol/restart H5 files, rerun trajectory H5 files, pure bundled absence of external `legacy_sidecars`, sidecar bundle presence of external `legacy_sidecars`, manifest presence, exact schema/schema_version, and core-only `mdin.override_same_path.spg.toml` plus `mdin.override_conflict.spg.toml`. The same manifest test parses TOML mdin input references, requires every legacy input file reference to exist, requires pure bundled mdin to retain no legacy input file references, requires sidecar-bundled mdin references to stay under materialized `legacy_sidecars/`, verifies the sidecar-preserving bundled mdin does not explicitly spell H5 sidecar-table keys that must be injected from the H5 key/path tables, and requires every materialized sidecar payload plus every H5 sidecar-table entry to have explicit manifest `sidecar_key` and `sidecar_path` fields matching `legacy_sidecars/<key>/<basename>`. | Runs under `h5_bundle;contract;coverage`. |
| Phase 0 fixture helper resolves fixture roots, copies fixture cases into temp workspaces, describes mdin, bundle input, and expected normal/rerun output paths, and reports missing fixture paths with explicit path-bearing errors. | `tests/h5_bundle/h5_input_matrix_fixture.hpp`, exercised by `test_h5_io_matrix_spec::Test_Fixture_Helper_Copies_Case_To_Temp`, `Test_Fixture_Helper_Describes_Case_Paths`, and `Test_Fixture_Helper_Missing_Path_Has_Explicit_Message`. | Runs under `h5_bundle;contract;matrix`. |
| Phase 1 legacy input leaves H5 input disabled and allows legacy input. | `test_h5_input_matrix_contract::Test_Legacy_Input_Allows_Legacy_Path`; `Test_Empty_H5_Input_Paths_Do_Not_Disable_Legacy_Fallback` locks that empty `input_h5_*_path` commands do not enable H5 input, do not disable legacy fallback, and keep `Has_H5_Input_Binding` false until a non-empty H5 path appears. | Runs under `h5_bundle;input;contract`. |
| Phase 1 pure bundled input resolves topology, protocol, restart, restart-load policy, and injects no sidecars. | `test_h5_input_matrix_contract::Test_Pure_Bundled_Input_Has_No_Legacy_Sidecars`; `Test_H5_Restart_Load_Policy_Validation` locks the default `structural` restart-load policy, accepted `structural`/`dynamic`/`protocol`/`full` values, invalid-value rejection, and reserved `custom` rejection; `Test_H5_Input_Path_Suffix_Flags_Are_Non_Fatal` locks that non-recommended H5 input suffixes only clear `has_recommended_suffix` flags and do not invalidate the resolved input plan; pure fixture tree checks in `test_h5_io_matrix_spec::Test_Pure_Bundled_Fixtures_Are_Sidecar_Free`, including topology/protocol/restart H5 files and rerun trajectory H5 files; `test_h5_input_fixture_equivalence.py` auto-discovers each pure/sidecar H5 file set and compares every matched file after excluding legacy-sidecar tables. | Runs under `h5_bundle;input;contract;matrix`. |
| Phase 1 any H5 input requires the required topology, protocol, and normal-mode restart bindings. | `test_h5_input_matrix_contract::Test_H5_Input_Requires_Topology_And_Protocol_Bindings` locks missing `input_h5_topology_path` and `input_h5_protocol_path` diagnostics after any H5 input is enabled; `Test_Normal_H5_Input_Requires_Restart_Binding` locks the normal-mode `input_h5_restart_path` requirement. | Runs under `h5_bundle;input;contract`. |
| Phase 1 bundled input with legacy sidecars stores key/path tables, resolves sidecar paths, and injects allowed keys. | `test_h5_input_matrix_contract::Test_Bundled_Input_With_Sidecar_Injects_Allowed_Keys`; `Test_Full_Contract_Sidecar_Key_Sets_Are_Injectable` locks full-contract topology/protocol/restart sidecar key sets and verifies every key is accepted by the production allowlists; `Test_Full_Contract_Sidecars_Match_Legacy_Source_Files` verifies sidecar paths point to `legacy_sidecars/<key>/<file>` and byte-match the original legacy inputs; sidecar table inventory in `test_h5_io_matrix_spec::Test_Bundled_With_Sidecar_Fixtures_Carry_Sidecar_Tables` now verifies topology/protocol/restart sidecar paths resolve to absolute, existing files under `legacy_sidecars`, locks both core structural and full-contract sidecar key distribution by H5 container, verifies the materialized `legacy_sidecars` file set equals the union of H5 sidecar-table references and mdin `legacy_sidecars/...` references, and verifies rerun trajectory H5 remains sidecar-free; `test_h5_input_fixture_equivalence.py` now requires each sidecar-preserving fixture manifest to contain embedded sidecar entries and byte-compares every embedded sidecar H5 scalar against its legacy source text; `test_h5_io_contract_manifest.py` now checks the core structural manifest and verifies every materialized `legacy_sidecars/<key>/<file>` payload byte-matches the relocated legacy source. | Runs under `h5_bundle;input;contract;matrix` and `h5_bundle;contract;coverage`. |
| Phase 1 legacy override behavior covers absent key, same-path idempotence, and different-path conflict. | `test_h5_legacy_sidecar::Test_Accepts_Relative_Existing_Path_For_Same_Sidecar` locks relative/absolute same-file idempotence in the sidecar injector; `test_h5_input_matrix_contract::Test_Legacy_Sidecar_Override_Conflict`; core fixture checks for `mdin.override_same_path.spg.toml` and `mdin.override_conflict.spg.toml`; runtime inventory and override fixture preparation in `test_h5_input_output_smoke_matrix::Validate_Runtime_Smoke_Preparation`. | Contract runs by default; runtime smoke is opt-in. |
| Phase 1 H5 restart input rejects legacy restart inputs. | `test_h5_input_matrix_contract::Test_H5_Restart_Rejects_Legacy_Restart_Inputs`. | Runs under `h5_bundle;input;contract`. |
| Phase 1 H5 input validation covers both resolved-plan validation and the controller-level parser path used by `SPONGE/main.cpp`. | `test_h5_input_validation::Test_Validates_Controller_Input_Bindings` builds topology/protocol/restart/trajectory H5 files, sets `input_h5_*` commands on `CONTROLLER`, requests `input_h5_restart_load = "full"`, and calls `SpongeH5InputValidation::Validate_Input_Bindings(&controller)` directly; `Test_Controller_Input_Bindings_Propagate_Resolver_Errors` verifies the same controller entry point propagates missing required H5 binding and reserved `input_h5_restart_load = custom` resolver errors before file-level validation. | Runs under `h5_bundle;input;backend-io;contract`. |
| Phase 1 v2 artifact identity and restart integrity are enforced rather than treated as descriptive metadata. | `test_h5_input_validation` rejects an invalid topology schema name, invalid protocol UUID, restart topology-lineage mismatch, and state-data mutation without a matching `/run/state_hash`. `RestartH5Reader::Compute_State_Hash` recomputes the canonical digest over `/particles/**` and `/parameters/restart/**`; the runtime fixture set is regenerated by XPONGE with `sponge.input.v2`, `/identity/uuid`, `/run/topology_hash`, `/run/atom_order_hash`, `/run/producer_protocol_hash`, and `/run/state_hash`. The enabled input/output, ReaxFF/EDIP, and restart-load smoke tests therefore exercise the same cross-language integrity path used in production. | Runs under `h5_bundle;input;backend-io;contract`; cross-language execution is runtime-smoke gated. |
| Phase 2 minimal normal-mode smoke covers legacy, bundled, and bundled-with-sidecar inputs against legacy and bundled outputs. | `test_h5_output_plan::Test_Resolve_Legacy_Output_Plan_Matrix` and `Test_Explicit_Legacy_Sidecar_Collection` lock the exact legacy output sidecar key set (`mdout`, `mdinfo`, `crd`, `box`, `vel`, `frc`, `rst`, `qc_scf_output`), default-enable behavior, H5-output suppression of defaults, and explicit sidecar collection order; `test_h5_input_output_smoke_matrix::Run_Normal_Mode_Matrix`; `test_h5_input_output_smoke_matrix::Require_Normal_Prepared_Mdin`; `test_h5_matrix_plan_manifest.py` locks the normal baseline run, the bundled-input core-mdout comparison branch, the strict text comparison branch for legacy/sidecar-compatible rows, `Require_Normal_Legacy_Restart_Output` for normal legacy-output runtime rows (`restart_coordinate.txt` and `restart_velocity.txt`), the bundled-output H5 restart/trajectory/observable readback branch, bundled-output preparation that keeps explicit `mdout`/`mdinfo` comparison files while rejecting default legacy trajectory/restart sidecar keys such as `crd`, `box`, `vel`, `frc`, `rst`, and `qc_scf_output`, and `Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs`, which checks that bundled-output runtime cases do not materialize default legacy `mdcrd.dat`, `mdbox.txt`, or `restart` files. | Registered by default but skips unless `SPONGE_H5_ENABLE_RUNTIME_SMOKE=1`. |
| Phase 2 bundled restart output checks position, velocity, box, step, and time. | `test_h5_input_output_smoke_matrix::Require_H5_Restart_Matches_Core_State`. | Runtime smoke only. |
| Phase 2 bundled trajectory output checks frame count, step/time, first-frame position, and box. | `test_h5_input_output_smoke_matrix::Require_H5_Trajectory_Has_Frames` and `Require_H5_Trajectory_First_Frame_Matches_Core_State`. | Runtime smoke only. |
| Phase 2 bundled observable output checks observable frame metadata, mdout column metadata, and per-column values against the same case's `mdout.txt`; the same observable-stream check is applied to both trajectory H5 and observable-only H5 outputs. | `test_h5_input_output_smoke_matrix::Require_H5_Observable_Stream_Matches_Mdout`. | Runtime smoke only. |
| Phase 3 sidecar smoke covers injected sidecars, same-path override, different-path conflict, and pure bundled input without sidecar files. | `test_h5_input_output_smoke_matrix::Run_Legacy_Sidecar_Smoke_Cases`; `Validate_Runtime_Smoke_Preparation` verifies the absent-key injection fixture does not explicitly spell representative H5 sidecar-table keys in mdin while those keys resolve from topology/protocol H5 sidecar tables, verifies the different-path fixture has `mass_in_file`, `override_mass.txt`, an H5 `mass_in_file` sidecar binding, and native `/atoms/mass`, so the failure remains tied to native H5 topology mass ownership, and also locks the same-path `qc_type_in_file` binding and pure-bundled sidecar-directory removal. For sidecar runtime fixtures, `Require_Materialized_Sidecars_Are_Exactly_H5_Referenced` requires the materialized `legacy_sidecars` file set to equal the union of H5 sidecar-table references, so runtime preparation cannot silently carry extra or missing sidecar payloads. Successful sidecar legacy-output smoke branches also call `Require_Normal_Legacy_Restart_Output`. `test_h5_matrix_plan_manifest.py` verifies that `Run_Normal_Mode_Matrix` still executes this function, that it remains bound to `Sidecar_Smoke_Cases`, that each switch branch retains its own run/output/failure assertions, and that the Phase 3 preparation tokens stay present. | Runtime smoke only. |
| Phase 4 rerun smoke covers bundled and bundled-with-sidecar rerun input to legacy output. | `test_h5_input_output_smoke_matrix::Run_Rerun_Mode_Matrix`; `test_h5_matrix_plan_manifest.py` locks the rerun legacy baseline run, the `Rerun_Legacy_Output_Cases` execution path, pure-bundled core mdout comparison, and sidecar/full rerun mdout comparison. | Runtime smoke only. |
| Phase 4 rerun smoke covers legacy, bundled, and bundled-with-sidecar rerun input to bundled output. | `test_h5_input_output_smoke_matrix::Rerun_Bundled_Output_Cases`; `test_h5_matrix_plan_manifest.py` locks the bundled-output branch for VDS on/off, pure-bundled versus sidecar/legacy comparison selection, trajectory frame readback, observable readback from both trajectory H5 and observable-only H5 outputs, and the same `Require_No_Default_Legacy_Trajectory_Or_Restart_Outputs` guard against default legacy `mdcrd.dat`, `mdbox.txt`, or `restart` leakage. | Runtime smoke only. |
| Phase 4 full-contract rerun H5 input plan resolves topology, protocol, restart, and trajectory bindings, rejects H5 trajectory input outside `mode = rerun`, allows rerun H5 trajectory input without a restart binding, defaults the H5 trajectory particle stream, treats rerun mode case-insensitively, and rejects legacy rerun `crd`/`box`/`vel` overrides. | `test_h5_input_matrix_contract::Test_Full_Contract_Rerun_H5_Trajectory_Input_Resolves`, `Test_H5_Trajectory_Is_Rerun_Only`, `Test_Rerun_H5_Trajectory_Input_Does_Not_Require_Restart_Binding`, `Test_Rerun_H5_Trajectory_Default_Stream_And_Mode_Case`, and `Test_H5_Trajectory_Rejects_Legacy_Rerun_Inputs`. | Runs under `h5_bundle;input;contract`. |
| Phase 4 full-contract bundled mdin preserves rerun mode, H5 input bindings, H5 output bindings, and does not mix legacy restart/rerun keys. | `test_h5_io_matrix_spec::Test_Full_Contract_Bundled_Mdin_Locks_H5_Input_Output_Contract`. | Runs under `h5_bundle;contract;matrix;rerun;vds`. |
| Phase 4 rerun frame selection obeys `rerun_start`, `rerun_strip`, and `rerun_frame_limit`. | Static fixture reader checks in `test_h5_io_matrix_spec::Test_Rerun_Frame_Selection_Uses_Fixture_Trajectory`; runtime row-count and stable-column mdout comparison in `test_h5_input_output_smoke_matrix::Run_Rerun_Frame_Selection_Smoke_Cases`; `test_h5_matrix_plan_manifest.py` treats the legacy and sidecar second-frame runtime cases as explicit extra smoke inventory outside the legal support matrix and verifies that `Run_Rerun_Mode_Matrix` still executes the frame-selection smoke through `Rerun_Selection_Cases`. The selection smoke intentionally excludes currently unstable aggregate columns while the runtime scrub forces H5 restart loading to `structural`. | Static test runs by default; runtime smoke is opt-in. |
| Phase 4 rerun trajectory input step/time, frame count, position, and box are readable; runtime bundled-output readback checks the stable rerun output position/box state for every bundled-output case while H5 rerun input rows keep distinct second-frame time expectations. | `test_h5_io_matrix_spec::Test_Rerun_Bundled_Fixtures_Contain_Trajectory_Input`; `test_trajectory_h5_reader`; runtime bundled-output readback in `test_h5_input_output_smoke_matrix::Require_H5_Trajectory_Frame_Matches_Rerun_Runtime_State`; `test_h5_matrix_plan_manifest.py` locks that assertion into the rerun bundled-output branch. | Static/backend tests run by default; runtime smoke is opt-in. |
| Phase 5 legal matrix is explicit instead of blind Cartesian product, with the full plan-table case-name set and mode/input/output/VDS combination counts locked. | `test_h5_io_matrix_spec::Legal_Matrix_Cases` and `Test_Legal_Matrix_Cases_Are_Explicit_And_Fixture_Backed` assert the 15 names from `INPUT_OUTPUT_MATRIX_TEST_PLAN.md`, including `normal_sidecar_in_*` and `rerun_sidecar_in_*`, plus exact mode/input/output/VDS counts; `test_h5_matrix_plan_manifest.py` also derives the runtime matrix dimensions from `Normal_Smoke_Cases`, `Rerun_Legacy_Output_Cases`, and `Rerun_Bundled_Output_Cases` and requires them to equal the plan table. | Runs under `h5_bundle;contract;matrix;rerun;vds`. |
| Phase 5 fixture-backed bundled rerun mdin aligns with the legal matrix dimensions: H5 input keys are present, H5 output keys are present, VDS is enabled for the broad fixture, and legacy rerun keys are absent. | `test_h5_io_matrix_spec::Test_Full_Contract_Bundled_Mdin_Locks_H5_Input_Output_Contract`. | Runs under `h5_bundle;contract;matrix;rerun;vds`. |
| Phase 5 VDS on/off is covered for rerun bundled trajectory output. | Static legal matrix in `test_h5_io_matrix_spec`; runtime cases in `test_h5_input_output_smoke_matrix::Rerun_Bundled_Output_Cases`; `test_h5_input_output_smoke_matrix::Require_VDS_Shards_Are_Complete` checks VDS-on runtime wrapper finalization/status metadata, repair metadata, wrapper completion datasets, shard manifest ranges/status, materialized shard files, shard step/time arrays, and shard position/box payload slices against the wrapper VDS datasets (`VDS shard position` and `VDS shard box`); `test_h5_vds_terminal_resume_smoke.cpp` checks complete-prefix terminal tail-shard repair and the no-op resume-policy path where all terminal shards finalize; VDS writer contracts in `test_vds_trajectory_writer_with_mock_backend` lock wrapper-relative VDS source paths for particle, observable, and module virtual datasets across both shard segments, and `test_highfive_backend_io` validates real HDF5 VDS readback; `test_h5_matrix_plan_manifest.py` verifies the runtime rerun matrix still calls the VDS shard-completion check from the `spec.vds` branch, keeps the VDS trajectory-observable readback branch, keeps the non-VDS trajectory-observable mdout comparison branch, and keeps the wrapper status assertions attached. | Static/mock/backend tests run by default; runtime smoke is opt-in. |
| Phase 6 broad fixture covers manifest-level buckets, sidecar key/path tables, and bundled output plan paths. | `test_h5_io_contract_coverage.cpp` checks sidecar tables with exact topology/protocol/restart key sets and relocation-safe `legacy_sidecars/<key>/...` paths, representative manifest buckets, manifest entry field completeness via `Require_Manifest_Entry_Fields`, legacy output sidecar preservation via `Test_Full_Contract_Rerun_Legacy_Output_Sidecar_Plan_Is_Preserved`, explicit bucket labels via `Test_Full_Contract_Rerun_Manifest_Buckets_Are_Explicit`, plan-bucket evidence via `Test_Phase6_Plan_Buckets_Are_Represented`, and required HDF5 bundle paths; `test_h5_io_contract_manifest.py` checks required manifest fields for both `core_structural` and `full_contract_rerun`, exact manifest schema/schema_version, top-level `bundled_mdin`/`case_root` relocation, exact core and full-contract protocol-sidecar sets, full-contract output-sidecar sets, component/direction/payload metadata, status-specific direction/payload-kind/override-policy semantics, source relocation, materialized sidecar byte provenance, bundle-file existence, and a semantic-evidence map that binds every required converted/typed input contract to `test_h5_input_fixture_equivalence.py`; `test_h5_input_matrix_contract::Test_Full_Contract_Sidecar_Key_Sets_Are_Injectable` verifies the broad fixture's sidecar key sets are injectable; `Test_Full_Contract_Sidecars_Match_Legacy_Source_Files` verifies sidecar file provenance against the legacy source files; `test_h5_io_matrix_spec::Test_Full_Contract_Bundled_Mdin_Locks_H5_Input_Output_Contract` locks the broad fixture's bundled mdin input/output keys. | Runs under `h5_bundle;contract;coverage`. |
| Phase 6 broad fixture covers topology typed datasets with semantic equivalence to legacy source files. | `test_h5_input_fixture_equivalence.py` compares H5 topology payloads against legacy inputs for mass, charge, QC type, bonds, angles, LJ, dihedrals, GB, Urey-Bradley, NB14, exclusions, virtual atom, LJ soft core/subsystem, and CMAP; `test_h5_io_contract_coverage.cpp` locks required topology, forcefield, QC, and `/manybody` HDF5 paths. | Runs under `h5_bundle;input;contract;matrix` and `h5_bundle;contract;coverage`. |
| Phase 6 broad fixture covers restart structural/dynamic/protocol payloads. | `test_h5_input_fixture_equivalence.py` compares restart coordinate, velocity, box, Nose-Hoover chain state, hills, metadynamics scatter/potential, every manifest-declared embedded sidecar text payload, and typed restraint reference coordinates; `test_h5_io_contract_manifest.py` requires `restart.restrain_coordinate`; `test_h5_io_contract_coverage.cpp` requires `/parameters/restart/references/restraint/default/coordinate`. | Runs under `h5_bundle;input;contract;matrix` and `h5_bundle;contract;coverage`. |
| Phase 6 broad fixture covers protocol typed datasets and protocol sidecars. | `test_h5_input_fixture_equivalence.py` compares protocol configs for CV, SITS, restraint, restraint-CV, steering-CV, metadynamics grid, soft walls, SITS atom indices, SITS `nk`, constraints, restraint atom ids, and restraint weights; manifest and coverage tests require `protocol.restrain_atom_id`, `protocol.restrain_weight`, and protocol-sidecar buckets. | Runs under `h5_bundle;input;contract;matrix` and `h5_bundle;contract;coverage`. |
| Native hard-wall protocol bounds are validated without a legacy text bridge. | `test_h5_input_validation::Test_Protocol_Reader_Loads_Native_Hard_Wall` accepts typed `/wall/hard/bounds_low`, `/wall/hard/bounds_high`, and matching `allow_npt` policy, then rejects invalid intervals. `NativeHardWallH5Reader` also rejects NaN, all-unbounded definitions, non-boolean policy values, and dataset/attribute disagreement; `HARD_WALL::Initial` rejects simultaneous native and legacy ownership. | Runs under `h5_bundle;input;backend-io;contract`; XPONGE producer validation runs in `test_1_io_bundle_saver.py`. |
| Phase 6 broad fixture covers custom pairwise/listed force payloads and many-body and specialized force fields. | `test_h5_input_fixture_equivalence.py` compares custom pair/listed force payloads against `pairwise_force.txt`, `custom_pair.txt`, `listed_forces.txt`, and `custom_bond.txt`; it also checks EAM, SW, EDIP, TERSOFF, and ReaxFF legacy semantics against `/manybody/*` HDF5 payloads. | Runs under `h5_bundle;input;contract;matrix`. |
| Phase 7 full-contract native parity expansion is explicit about Current runtime native-parity exceptions, the REAXFF/EDIP pure bundled closure, and restart-load runtime closure. | `INPUT_OUTPUT_MATRIX_TEST_PLAN.md` lists only `input_h5_restart_load` as the current broad-matrix runtime native-parity exception, and lists the current pure bundled rerun native core mdout columns: `temperature`, `LJ_short`, `LJ_long`, `LJ`, `LJ_soft`, `LJ_soft_short`, `LJ_soft_long`, `PM`, `bond`, `angle`, `urey_bradley`, and `dihedral`; `test_h5_reaxff_edip_runtime_parity.cpp` runs legacy, sidecar-compatible, and pure-bundled input against legacy and VDS on/off bundled output, compares `EDIP` plus all `REAXFF_*` mdout columns, verifies bundled observable output, and requires that pure bundled runs create none of `.sponge_h5_native_manybody/{edip,reaxff,reaxff_type}.txt`; `TopologyManybodyH5Materializer::Read_ReaxFF` validates fixed-width typed sections once and `native_init.cpp` stages all eight ReaxFF modules from that shared object; custom pairwise/listed forces and QC likewise initialize directly; `test_h5_restart_load_runtime_closure.cpp` requires native QC/EAM/SW/Tersoff/ReaxFF diagnostics and absence of `.sponge_h5_native_qc`, `.sponge_h5_native_custom_force`, and all former many-body text bridge files; the same closure runs supported `input_h5_restart_load = "protocol"`, `"dynamic"`, and `"full"` cases, proves SITS `nk`/`log_norm`/`log_nk` device restoration by exact typed write-back without a protocol-sidecar directory, and locks expected rejection diagnostics for unsupported continuation state; manifest tests bind these runtime cases to the plan and required coverage buckets. | Runs under `h5_bundle;input;smoke;matrix;rerun;vds;manybody`, `h5_bundle;input;smoke;rerun;restart`, `h5_bundle;contract;matrix`, and `h5_bundle;contract;coverage`. |

The Current runtime native-parity exceptions remain explicitly bound to
`Scrub_Runtime_Unstable_Rerun_Features`,
`rerun_runtime_scrub_prepare_check`, and
`Require_Runtime_Unstable_Rerun_Features_Scrubbed`. The unsimplified evidence
cases remain `manybody_legacy_in_legacy_out_no_manybody_scrub`,
`manybody_sidecar_in_legacy_out_no_manybody_scrub`,
`manybody_sidecar_in_bundled_out_vds_*_no_manybody_scrub`,
`manybody_pure_bundled_in_legacy_out_no_sidecar`, and
`manybody_pure_bundled_in_bundled_out_vds_*_no_sidecar`. The restart closure
now covers both the protocol SITS sidecar compatibility bridge and isolated native
SITS `nk`/`log_norm`/`log_nk` device-state restoration, NHC dynamic state, initialized
metadynamics restart loading, and pure-bundled native custom-force direct
initialization. Compatibility coverage remains anchored by
`topology_custom_force_h5_materializer.hpp`; terminal VDS coverage remains in
`test_h5_vds_terminal_resume_smoke.cpp`.
`test_h5_matrix_plan_manifest.py` requires the exception table to remain
aligned with the executable smoke matrix, and runtime preparation builds
`rerun_runtime_scrub_prepare_check` before applying those scoped exceptions.

## Acceptance Evidence

The `INPUT_OUTPUT_MATRIX_TEST_PLAN.md` acceptance criteria are backed by these
test families:

- Internal plan and validation tests: `test_h5_output_plan`,
  `test_h5_input_matrix_contract`, and `test_h5_input_validation`.
- Sidecar injection/conflict handling: `Test_Bundled_Input_With_Sidecar_Injects_Allowed_Keys`
  and `Test_Legacy_Sidecar_Override_Conflict`.
- Runtime input/output smoke: `Run_Normal_Mode_Matrix`,
  `Run_Rerun_Mode_Matrix`, `Require_Core_Mdout_Equivalent`, and
  `Require_H5_Observable_Stream_Matches_Mdout`.
- Rerun bundled trajectory VDS on/off: `Rerun_Bundled_Output_Cases` and the
  Phase 5 VDS on/off rows above.
- Broad full-contract fixture checks: `test_h5_io_contract_manifest.py` and
  `test_h5_io_contract_coverage.cpp`.
- Converter sidecar provenance: byte-match checks and relocation-safe sidecar
  path checks.
- Legacy override behavior: absent key, same-path, and different-path conflict
  cases.

## Phase 6 Buckets

`test_h5_audit_matrix_manifest.py` keeps the Phase 6 Required coverage buckets
in the plan bound to `test_h5_io_contract_coverage.cpp`. The required buckets
are: Topology typed datasets; Restart structural state; Restart dynamic state;
Protocol sidecars; SITS state and sidecars; Metadynamics state and sidecars;
Custom pairwise/listed force payloads; QC/ReaxFF sidecars; Rerun trajectory
input; Legacy sidecar key/path tables; Bundled output
trajectory/restart/observable paths.

## Runtime Gate

`test_h5_input_output_smoke_matrix` always validates its inventory, generated
mdin preparation, input-family key sets, and helper-derived normal/rerun output
paths before checking the runtime gate. Its preparation phase now also
constructs the full 15-case legal runtime inventory from the normal baseline,
normal matrix, rerun baseline, rerun legacy-output cases, and rerun
bundled-output VDS on/off cases, then compares that set to the plan table. The
same preparation checks also lock that normal bundled-output cases request H5
restart/trajectory/observable output, while rerun bundled-output cases request
H5 trajectory/observable output with the selected VDS flag and do not request
H5 restart output. It returns CTest skip code 77 unless
`SPONGE_H5_ENABLE_RUNTIME_SMOKE=1` is set. A machine that can initialize the
selected SPONGE backend must run:

```bash
SPONGE_H5_ENABLE_RUNTIME_SMOKE=1 \
  ctest --test-dir build-h5-tests -R 'test_h5_input_output_smoke_matrix|test_h5_reaxff_edip_runtime_parity|test_h5_restart_load_runtime_closure|test_h5_vds_terminal_resume_smoke' \
  --output-on-failure
```

The current default suite proves parser, contract, fixture, reader, writer,
static matrix, sidecar provenance, and manifest behavior. End-to-end legacy
alignment remains proven only when the runtime smoke gate succeeds in a
runnable SPONGE environment. The equivalent helper command is:

```bash
tests/h5_bundle/run_h5_bundle_tests.sh test-smoke-runtime
```

## Production Gate

`pixi run -e dev-cuda13 smoke-bundled-io-production` is the bundled I/O
production-readiness gate. It depends on:

- `smoke-bundled-io-contract`: XPONGE-generated bundle contract smoke plus the
  production-gate manifest test.
- `smoke-bundled-io-runtime`: gated runtime matrix, REAXFF/EDIP parity,
  restart-load closure, and VDS terminal/resume smoke.
- `smoke-bundled-io-staged`: full staged H5 bundle contract, mock,
  backend-io, matrix, rerun, VDS, and smoke-label validation.
- `smoke-gpu-medium`: CUDA runtime preflight plus selected real
  validation/performance scenarios on a GPU environment, covering
  thermostat/CV/softcore, ReaxFF throughput, and nonorthogonal NVT long-run RDF.

The VDS terminal/resume smoke covers complete-prefix terminal tail-shard repair
and the no-op complete-prefix policy path. Cross-process reopen-and-append
resume remains outside the current production contract until the writer
implements an explicit append/resume open mode.

Legacy/bundled behavior A/B is a separate opt-in release gate, not part of the
default smoke task:

- `pixi run -e dev-cuda13 ab-bundled-io-medium`: ten independently seeded
  paired NVT runs for each normal-MD input, with a denser output cadence for
  block statistics, plus deterministic full-contract rerun checks.
- `pixi run -e dev-cuda13 ab-bundled-io-production`: longer five-replica NVT
  runs with the same statistical contract and the same deterministic rerun
  matrix.

Normal MD comparison deliberately does not require frame-by-frame physical
agreement. For every emitted mdout column, it removes configured warmup frames,
forms non-overlapping block means for each replica, and requires both a
confidence-bounded practical mean equivalence and a bounded block-standard-
deviation ratio. Step/frame/time schedules remain exact because they are not
stochastic. Non-finite values are never silently skipped: their complete
per-frame pattern must match exactly. The shorter medium profile uses a 5%
relative practical margin; the longer production profile tightens that margin
to 1% with five longer replicas. The medium profile uses ten replicas to keep
the shorter run statistically powered. Their absolute practical
margins are respectively 0.15 and 0.05 in native observable units, so
near-zero energy terms are judged against a meaningful noise floor rather than
an undefined relative percentage.

Each branch also validates its own H5 output contract before the cross-branch
comparison. Trajectory, restart, and observable files must have complete frame
metadata and shapes; the observable H5 original-name/HDF5-name mapping must
cover every mdout column in order, and each H5 value stream must match that
branch's mdout stream. The cross-branch gate then compares every H5 object path
and dataset schema. Numeric trajectory/restart/observable payloads use the
same statistical policy (including frame mean and RMS summaries when shaped as
frames); schedules and string metadata remain exact.

The A/B metrics include wall runtime, H5 trajectory/observable/restart file
sizes, total H5 bytes, VDS shard counts, block-statistics summaries, output
schema counts, and runtime-reported H5 finalize timing for trajectory,
observable, restart, and total finalize cost. This is flush-through-finalize
evidence, not a per-frame or per-writer flush breakdown.

The rerun A/B cases use both the pure bundled fixture and the
bundled-with-legacy-sidecar fixture, each with VDS off and on. They retain
`input_h5_restart_load = "structural"` because a rerun does not initialize
dynamic NHC/metadynamics modules. They verify the broad full-contract input H5
path inventory before execution, and compare every active mdout column and
every emitted H5 dataset value/schema exactly after execution. The separately
required `test_h5_restart_load_runtime_closure` smoke exercises the supported
`dynamic`, `protocol`, and `full` restart-load combinations with the modules
that own those states initialized; it also locks the expected rejection when a
full restart state is incompatible with the requested runtime mode. This
closes the old gap where rerun output was merely checked for path existence.
Legacy-text rerun currently reaches its frame loop before H5 output writers are
initialized, so it has no peer H5 artifact. The gate records that limitation
explicitly and bridges it semantically: legacy and bundled mdout rows must be
identical, bundled observable H5 must match bundled mdout for every emitted
column, and bundled trajectory H5 payloads must match the selected bundled H5
input frames. Removing this limitation requires a runtime/output implementation
change, not a weaker A/B assertion.

The A/B gate also includes the existing SITS `ala2_sits` static system as a
normal-run ff19SB CMAP case. That system is a solvated ACE-ALA-NME peptide and
contains one normal ff19SB 24x24 CMAP entry (`4 6 8 14 16 0`), so it covers a
real force-field CMAP input path rather than only the synthetic full-contract
fixture. Its complete finite `cmap` stream is evaluated with the same paired
block-statistical equivalence policy as every other MD observable. The native H5 reader is
aligned with legacy `nb14_in_file` semantics for this path: two-column nb14
params are interpreted as LJ and Coulomb scale factors and materialized through
the already-loaded LJ pair parameters, while three-column nb14 params remain
direct A/B/Coulomb-scale materialized values.

Those A/B gates are required evidence before promoting bundled I/O from an
opt-in/shadow path to a legacy-behavior replacement candidate. Together with
the contract, reader, matrix, and restart-closure suites, a passing
`ab-bundled-io-production` proves equivalence for the explicitly enumerated
supported input/output contract and runtime scenarios. It does not claim to
prove behavior outside that contract, including cross-process VDS
reopen-and-append resume, which remains unsupported.

Equivalence is limited to explicitly enumerated supported contracts and runtime
scenarios and excludes cross-process VDS reopen-and-append/resume.
