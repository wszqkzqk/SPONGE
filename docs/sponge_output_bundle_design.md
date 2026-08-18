# SPONGE Output Bundle Design

## 0. Goal

SPONGE output should support two HDF5-native modes:

```text
single-file H5MD                 # default local/simple path
chunked H5MD shards + HDF5 VDS   # cloud/segment-sync path
```

The design goal is not to make restart files carry trajectories. Full trajectory
history can be large, especially for cloud production runs. Restart remains a
lightweight single-state file, while the output bundle owns observation history.

Core rule:

```text
restart.spgr.h5 owns one launchable state.
output bundle owns trajectory and observable history.
H5MD is the logical output schema.
HDF5 VDS is an optional assembly layer for chunked output.
```

The normative schema versions are:

```text
sponge.input.v2   # topology, protocol, and restart/checkpoint artifacts
sponge.output.v2  # trajectory, observable, VDS wrapper, and VDS shards
```

Every process launch generates a new UUID. The same value is written to
`/identity/uuid` in every output artifact created by that launch. Publication
epochs start at zero for the new UUID; they are not continued from files left
by an earlier process.

## 1. Mode Selection

Use a run-level output flag to select the physical output mode:

```toml
[output.h5.trajectory]
path = "prod_0007.spg.h5md"
vds = false
chunk_size = 20
repair_policy = "strict"
```

With the current TOML flattening parser, this becomes:

```text
output_h5_trajectory_path = "prod_0007.spg.h5md"
output_h5_trajectory_vds = false
output_h5_trajectory_chunk_size = 20
output_h5_trajectory_repair_policy = "strict"
```

Recommended modes:

| `output_h5_trajectory_vds` | Runtime write target | Final artifact | Use case |
|---:|---|---|---|
| `false` | one H5MD file | same H5MD file | local runs, small/medium outputs, simple archive |
| `true` | segmented H5MD shard files | H5MD VDS wrapper + shard directory | cloud runs, partial sync, large trajectories |

`chunk_size` is the VDS file-level shard size in trajectory frames. It is not an
HDF5 dataset internal chunk shape, and it does not affect striping, compression,
or `write_trajectory_interval`. The default is 20 trajectory frames.

`repair_policy` controls VDS wrapper finalization only. The default `strict`
mode fails on incomplete or non-contiguous shard manifests. The explicit
`complete_prefix` mode may finalize from the complete contiguous shard prefix;
it does not delete orphan shard files and is invalid when VDS output is off.

## 2. Recommended File Names

Single-file mode:

```text
prod_0007.spg.h5md
```

VDS mode:

```text
prod_0007.spg.h5md                 # small VDS wrapper / logical H5MD entrypoint
prod_0007.spg.shards/              # writer-derived internal shard directory
  segment_000000.spg.h5md
  segment_000001.spg.h5md
  segment_000002.spg.h5md
```

The VDS wrapper is the user-facing H5MD file. It exposes the same logical layout
as the single-file output, but its large datasets are HDF5 virtual datasets that
map to shard files. The shard directory is derived by the writer from the
trajectory H5MD path and is not configurable from `run.mdin`.

Observable-only mode:

```text
prod_0007.obs.spg.h5md
```

The observable-only file is a lightweight H5MD artifact for analysis and sync. It
contains `/h5md`, `/observables`, and `/parameters`, but no trajectory particle
fields under `/particles`.

## 3. Ownership Boundary

The output bundle owns time-dependent observations:

```text
- positions over time
- velocities over time, if requested
- forces over time, if requested
- box over time
- energies
- temperature
- pressure
- volume
- CV values
- bias values
- wall/custom force energy outputs
- mdout-like scalar time series
- output metadata and labels
```

It does not own:

```text
- canonical topology
- canonical force-field definition
- canonical protocol definition
- launchable restart state
- resolver binding policy
```

It may record hashes, UUIDs, and resolved launch metadata for provenance, but
those records are derived/cache data.

If any canonical H5 output is enabled, legacy output files are disabled by
default. Legacy files are written only when their legacy path keys are
explicitly set.

## 4. Logical H5MD Layout

Both output modes must expose the same logical H5MD layout.

Recommended top level:

```text
/
  h5md/
  particles/
  observables/
  parameters/
```

### 4.1 `/h5md`

```text
/h5md
  attrs:
    version = [1, 1]

/h5md/creator
  attrs:
    name = "SPONGE"
    version = "<engine version>"

/h5md/author
  attrs:
    name = "<optional>"
    email = "<optional>"
```

### 4.2 `/particles`

For ordinary MD:

```text
/particles/all/
```

For replica/PIMD-like output:

```text
/particles/replica0/
/particles/replica1/
...
```

or, if compatibility with MindSPONGE naming is desired:

```text
/particles/trajectory0/
/particles/trajectory1/
...
```

The group name should be stable within one output file and recorded in
`/parameters/sponge/output/particle_streams`.

For each particle stream:

```text
/particles/all/step
/particles/all/time
/particles/all/position/value
/particles/all/velocity/value
/particles/all/force/value
/particles/all/box/edges/value
/particles/all/species
```

Recommended shapes:

```text
/particles/all/step                 int64[n_frame]
/particles/all/time                 float64[n_frame]
/particles/all/position/value       float32[n_frame, N, 3]
/particles/all/velocity/value       float32[n_frame, N, 3]
/particles/all/force/value          float32[n_frame, N, 3]
/particles/all/box/edges/value      float32[n_frame, 3, 3]
/particles/all/species              int32[N]
```

Units should be stored as attributes on each value array:

```text
position/value.attrs["unit"] = "Angstrom" or "nm"
velocity/value.attrs["unit"] = "Angstrom ps-1" or "nm ps-1"
force/value.attrs["unit"] = engine force unit
box/edges/value.attrs["unit"] = same as position
time.attrs["unit"] = "ps"
```

Use one shared `step` and `time` dataset per particle stream. Element-local
`step` and `time` paths should be HDF5 hard links:

```text
/particles/all/position/step   -> hard link to /particles/all/step
/particles/all/position/time   -> hard link to /particles/all/time
/particles/all/velocity/step   -> hard link to /particles/all/step
/particles/all/velocity/time   -> hard link to /particles/all/time
/particles/all/box/edges/step  -> hard link to /particles/all/step
/particles/all/box/edges/time  -> hard link to /particles/all/time
```

In VDS mode, `value` datasets are virtual datasets in the wrapper file. `step`
and `time` may also be VDS datasets if they span multiple shards.

### 4.3 `/observables`

For ordinary MD:

```text
/observables/all/
```

For replica/PIMD-like output, mirror the particle stream naming:

```text
/observables/replica0/
/observables/replica1/
```

or:

```text
/observables/trajectory0/
/observables/trajectory1/
```

Recommended scalar layout:

```text
/observables/all/step
/observables/all/time
/observables/all/potential_energy/value
/observables/all/kinetic_energy/value
/observables/all/total_energy/value
/observables/all/temperature/value
/observables/all/pressure/value
/observables/all/volume/value
```

Recommended enhanced-sampling layout:

```text
/observables/all/cv/value
/observables/all/cv/labels
/observables/all/biases/value
/observables/all/biases/labels
/observables/all/bias_potential/value
```

Recommended shapes:

```text
/observables/all/step                         int64[n_frame]
/observables/all/time                         float64[n_frame]
/observables/all/potential_energy/value       float32[n_frame]
/observables/all/temperature/value            float32[n_frame]
/observables/all/cv/value                     float32[n_frame, n_cv]
/observables/all/cv/labels                    string[n_cv]
/observables/all/biases/value                 float32[n_frame, n_bias]
/observables/all/biases/labels                string[n_bias]
```

Use one shared `step` and `time` dataset per observable stream. Element-local
`step` and `time` paths should be HDF5 hard links:

```text
/observables/all/potential_energy/step   -> hard link to /observables/all/step
/observables/all/potential_energy/time   -> hard link to /observables/all/time
/observables/all/temperature/step        -> hard link to /observables/all/step
/observables/all/temperature/time        -> hard link to /observables/all/time
/observables/all/cv/step                 -> hard link to /observables/all/step
/observables/all/cv/time                 -> hard link to /observables/all/time
```

### 4.4 `/parameters`

`/parameters` is for application-specific custom data and metadata that is not
standard particle or observable time-series data.

Recommended:

```text
/parameters/sponge/
/parameters/output_labels/
/parameters/finalization/
```

Do not put ordinary trajectory positions, velocities, forces, or box histories
under `/parameters`.

### 4.5 `/parameters/sponge`

`/parameters/sponge` stores SPONGE-specific runtime bookkeeping that is not part of H5MD.

```text
/parameters/sponge/schema/name = "sponge.output.h5md"
/parameters/sponge/schema/version
/parameters/sponge/output/mode = "single" | "vds"
/parameters/sponge/output/status
/parameters/sponge/output/frame_count
/parameters/sponge/output/publication_epoch
/parameters/sponge/output/streams/<stream>/committed_count
/parameters/sponge/output/particle_streams
/parameters/sponge/output/observable_streams
/parameters/sponge/output/chunk_policy
/parameters/sponge/output/compression
/parameters/sponge/output/trajectory_chunk_size
/parameters/sponge/output/shard_manifest
/parameters/sponge/provenance/topology_hash
/parameters/sponge/provenance/protocol_hash
/parameters/sponge/provenance/restart_hash
/parameters/sponge/provenance/mdin_hash
/parameters/sponge/provenance/engine_version
/parameters/sponge/provenance/command_line
```

`status` should be one of:

```text
open
closing
finalized
failed
```

## 5. Single-File H5MD Mode

When `output_h5_trajectory_vds = false`, SPONGE writes one ordinary H5MD file
directly:

```text
prod_0007.spg.h5md
```

Properties:

```text
- simplest reader compatibility
- no shard manifest needed
- no VDS dependency
- larger single file transfer for cloud workflows
- best default for local runs and modest outputs
```

Large datasets should still be chunked and compressed internally using HDF5
chunking filters.

Recommended dataset chunking:

```text
position/value       [frame_chunk, atom_chunk, 3]
velocity/value       [frame_chunk, atom_chunk, 3]
force/value          [frame_chunk, atom_chunk, 3]
box/edges/value      [frame_chunk, 3, 3]
scalar/value         [frame_chunk]
cv/value             [frame_chunk, n_cv]
biases/value         [frame_chunk, n_bias]
```

## 6. Chunked H5MD + VDS Mode

When `output_h5_trajectory_vds = true`, SPONGE writes completed frame segments
into shard files and maintains a small VDS wrapper:

```text
prod_0007.spg.h5md
prod_0007.spg.shards/segment_000000.spg.h5md
prod_0007.spg.shards/segment_000001.spg.h5md
...
```

Each shard contains at most `output_h5_trajectory_chunk_size` trajectory frames.
The default is 20. The shard index is determined by trajectory frame index:

```text
shard_index = trajectory_frame_index / output_h5_trajectory_chunk_size
```

Each shard is itself an HDF5/H5MD-like segment containing a contiguous frame
range:

```text
/parameters/sponge/shard/index
/parameters/sponge/shard/frame_start
/parameters/sponge/shard/frame_count
/particles/all/step
/particles/all/time
/particles/all/position/value
/observables/all/step
/observables/all/time
/observables/all/<name>/value
```

The wrapper exposes the full logical H5MD layout with VDS datasets:

```text
/particles/all/position/value   VDS over shard position/value datasets
/particles/all/velocity/value   VDS over shard velocity/value datasets
/particles/all/box/edges/value  VDS over shard box/edges/value datasets
/observables/all/<name>/value   VDS over shard observable datasets
```

Recommended manifest:

```text
/parameters/sponge/output/shard_manifest/table
```

Manifest fields:

```text
index                  int64
path                   string
frame_start            int64
frame_count            int64, trajectory frame count
step_start             int64
step_end               int64
time_start             float64
time_end               float64
content_hash           string, optional
status                 enum: open/complete/failed
```

The wrapper also serializes each completed shard's independent stream counts:

```text
/parameters/sponge/output/shard_manifest/stream_counts/particles
/parameters/sponge/output/shard_manifest/stream_counts/observables
/parameters/sponge/output/shard_manifest/stream_counts/nose_hoover_chain
/parameters/sponge/output/shard_manifest/stream_counts/sits
/parameters/sponge/output/shard_manifest/stream_counts/metadynamics
/parameters/sponge/output/shard_manifest/stream_counts/qc
/parameters/sponge/output/shard_manifest/stream_counts/reaxff
```

Every array is aligned with `shard_manifest/index`. Counts may differ between
streams because trajectory, mdout, and module cadences are independent.

Cloud sync should transfer only completed shard files plus the updated VDS wrapper
and manifest. A shard is syncable only after its status is `complete`.
Observable streams in a shard may have a different frame count from trajectory
streams because `write_mdout_interval` and `write_trajectory_interval` are
separate cadence domains.

On Windows, a reader that keeps the published wrapper open may prevent replacing
that pathname. Intermediate publication therefore keeps the last complete wrapper
visible and retries the newest complete prefix on the next publication. Initial
and final publication remain strict and report an error if the pathname is still
busy.

`shard_manifest/path` is relative to the wrapper directory. Each completed
entry also records `shard_manifest/byte_size`; a downloader must verify that
size before atomically promoting a downloaded shard. A future minor revision
may add a cryptographic content hash without changing these v2 semantics.

## 7. VDS Finalization

VDS mode does not require copying all trajectory data into one monolithic HDF5
file. Finalization means:

```text
1. Close the current shard.
2. Mark all complete shards immutable.
3. Validate the shard manifest.
4. Update the private same-directory wrapper working file.
5. Create HDF5 hard links for element-local step/time paths in the wrapper.
6. Mark /parameters/sponge/output/status = finalized.
7. Flush the working file, copy it to a temporary snapshot, and atomically
   replace the public wrapper.
```

The public wrapper is never modified in place. Intermediate completed-prefix
publications use the same temporary-snapshot and atomic-replace sequence, so a
polling synchronizer observes either the previous complete wrapper or the next
complete wrapper.

The finalized user-facing artifact is:

```text
prod_0007.spg.h5md
```

but it depends on:

```text
prod_0007.spg.shards/
```

If a fully portable single-file archive is required later, an explicit repack
operation may materialize the VDS into a standalone H5MD file. That is an archive
export step, not required for normal VDS finalization.

## 8. Relation to `restart.spgr.h5`

Output finalization must not require transferring the full trajectory to start a
new run.

Production continuation should use:

```text
restart.spgr.h5
```

not:

```text
*.spg.h5md:selected frame
```

The restart path is a replaceable checkpoint slot rather than an append-only
history. Each checkpoint is first written and finalized at a same-directory
temporary path, then closed and atomically renamed over the public path. A
failed write, finalization, close, or rename leaves the previously published
checkpoint untouched. Stale temporary files are never treated as restart
artifacts and are removed before the next attempt.

Each published restart contains:

```text
/parameters/sponge/output/restart_generation
/parameters/sponge/output/last_complete_step
/parameters/sponge/output/last_complete_time
/identity/uuid
/schema/name = "sponge.restart.h5"
/schema/version = "sponge.input.v2"
```

`restart_generation` is positive and increases once per successful atomic
replacement. The producer advances its in-memory generation only after the
rename succeeds, so readers can combine generation, step, time, and file
identity to detect a newly available checkpoint without a push channel.

A polling consumer recognizes a checkpoint produced by the current launch only
when all of the following hold:

```text
/identity/uuid equals the expected launch UUID
/parameters/sponge/schema/name equals "sponge.restart.h5"
/parameters/sponge/schema/version equals "sponge.input.v2"
/parameters/sponge/output/mode equals "checkpoint"
/parameters/sponge/output/status equals "finalized"
restart_generation is positive and greater than the accepted generation
frame_count is 1 and the required structural datasets contain exactly one frame
```

An existing restart loaded as input belongs to a parent launch. The new process
does not continue its generation: it creates a new UUID and publishes generation
1 for its first checkpoint. File name, mtime and inode are not restart identity.

The last trajectory frame may be used to generate a restart file, but that is a
conversion operation. The next launch should bind the produced `restart.spgr.h5`
directly.

## 9. Failure and Resume Semantics

Frame-prefix safety and cross-process writer resume are separate capabilities.
Version 2 implements safe published prefixes, but does **not** implement
cross-process continuation of an existing single-file or VDS output.

Required behavior:

```text
- Every appended frame must be either complete or ignored.
- A frame is complete only when all required arrays for that stream are written.
- frame_count must only advance after completing required writes.
- VDS finalization refuses inconsistent shard manifests unless repair is explicitly requested.
```

After a process failure, a new SPONGE process creates a new `/identity/uuid`,
starts publication epoch at zero, and overwrites or creates its configured
outputs. It does not reopen an old wrapper, scan an old shard directory, clean
orphan shards, or continue an old shard index. Cross-process VDS resume remains
explicitly out of scope for this revision.

Recommended markers:

```text
/parameters/sponge/output/status
/parameters/sponge/output/frame_count
/parameters/sponge/output/last_complete_step
/parameters/sponge/output/last_complete_time
```

Incremental readers should prefer the append-only publication contract:

```text
/parameters/sponge/output/publication_epoch
/parameters/sponge/output/streams/<stream>/committed_count
```

The writer flushes frame data first, appends all dirty stream watermarks and one
publication epoch, and then flushes the publication metadata. One publication
therefore does not add a separate flush for every module. The legacy
`output/frame_count` remains available for compatibility, but it represents the
file's primary stream rather than every embedded cadence domain.

In VDS mode, shard-level markers are also required:

```text
/parameters/sponge/shard/status
/parameters/sponge/shard/frame_start
/parameters/sponge/shard/frame_count
/parameters/sponge/shard/last_complete_step
/parameters/sponge/shard/last_complete_time
```

## 10. Minimal Single-File H5MD Example

```text
prod_0007.spg.h5md
  /h5md
  /particles/all/step
  /particles/all/time
  /particles/all/position/value
  /particles/all/position/step   -> hard link
  /particles/all/position/time   -> hard link
  /particles/all/velocity/value
  /particles/all/velocity/step   -> hard link
  /particles/all/velocity/time   -> hard link
  /particles/all/box/edges/value
  /particles/all/box/edges/step  -> hard link
  /particles/all/box/edges/time  -> hard link
  /observables/all/...
  /parameters/...
  /parameters/sponge/...
```

## 11. Minimal VDS H5MD Example

```text
prod_0007.spg.h5md
  /h5md
  /particles/all/step                 VDS over shards
  /particles/all/time                 VDS over shards
  /particles/all/position/value       VDS over shards
  /particles/all/position/step        -> hard link
  /particles/all/position/time        -> hard link
  /observables/all/temperature/value  VDS over shards
  /parameters/sponge/output/shard_manifest

prod_0007.spg.shards/
  segment_000000.spg.h5md
  segment_000001.spg.h5md
  segment_000002.spg.h5md
```

## 12. Observable-only H5MD

An optional observable-only H5MD artifact may be requested with:

```toml
[output.h5.observable]
path = "prod_0007.obs.spg.h5md"
```

The parser-visible key is:

```text
output_h5_observable_path
```

The observable-only file must contain:

```text
/h5md
/observables
/parameters
/parameters/sponge
```

It must not contain trajectory particle fields:

```text
/particles/all/position
/particles/all/velocity
/particles/all/force
/particles/all/box
```

Its cadence follows observable streams such as `write_mdout_interval`; it is not
controlled by trajectory frame cadence. It may be produced together with the
trajectory H5MD file, or as the only H5 analysis output.

## 13. Current SPONGE Output Inventory and Mapping

This section records the currently observed SPONGE output contract and how each
output should be represented in the finalized `*.spg.h5md` bundle or in
`restart.spgr.h5`.

The main rule is:

```text
trajectory-like time series -> *.spg.h5md
analysis/scalar time series -> *.spg.h5md
restart/checkpoint state    -> restart.spgr.h5
SPONGE-only logs/debug data -> *.spg.h5md:/parameters/sponge
```

### 13.1 Output Cadence Domains

SPONGE currently uses separate output intervals. The H5MD layout must preserve
that separation instead of forcing all datasets onto one shared frame axis.

| Existing control | Meaning | Bundle axis |
|---|---|---|
| `write_information_interval` | Default base interval for information output. | Used only as fallback when a more specific interval is absent. |
| `write_trajectory_interval` | Coordinate, box, velocity, and force trajectory cadence. | `/particles/<stream>/*/step` and `/particles/<stream>/*/time`. |
| `write_mdout_interval` | Scalar `mdout` / `Step_Print` cadence. | `/observables/<stream>/*/step` and `/observables/<stream>/*/time`. |
| `write_restart_file_interval` | Legacy restart export cadence. | `restart.spgr.h5`, not the output trajectory bundle. |

TOML decks may use `[write.interval] information`, `trajectory`, `mdout`, and
`restart`/`restart_file` aliases. The parser normalizes those aliases to the
four existing `write_*_interval` runtime commands before output planning.

`/particles` and `/observables` may therefore have different frame counts.

### 13.2 Legacy File Outputs

| Current output | Existing control/key | Current format | New contract |
|---|---|---|---|
| `mdinfo.txt` or `<prefix>.info` | `mdinfo` | Text log containing command and mdin echo. | `/parameters/sponge/log/mdinfo_text` as UTF-8 text, or recorded as an external sidecar under `/parameters/sponge/files`. |
| `mdout.txt` or `<prefix>.out` | `mdout` | Text table generated by `Step_Print`. | Canonical data goes to `/observables/<stream>/<column>/value`; original labels and ordering go to `/parameters/sponge/mdout/columns`. |
| `mdcrd.dat` or `<prefix>.dat` | `crd` | Binary `VECTOR[N]` frames. | `/particles/all/position/value` with shape `float32[n_frame,N,3]`. |
| `mdbox.txt` or `<prefix>.box` | `box` | Text box length and angle per frame. | `/particles/all/box/edges/value` as `float32[n_frame,3,3]`; original length/angle may be mirrored under `/parameters/sponge/box_legacy` if needed. |
| Velocity trajectory | `vel` | Optional binary `VECTOR[N]` frames. | `/particles/all/velocity/value` with shape `float32[n_frame,N,3]`. |
| Force trajectory | `frc` | Optional binary `VECTOR[N]` frames. | `/particles/all/force/value` with shape `float32[n_frame,N,3]`. |
| Legacy restart export | `rst`, `max_restart_export_count` | Amber `*.rst7` or native coordinate/velocity text files. | Belongs to `restart.spgr.h5`; the output bundle may only record legacy export metadata under `/parameters/sponge/restart_exports`. |

The finalized H5MD file should not preserve legacy filenames as canonical data
owners. Legacy names are provenance only.

If canonical H5 outputs are enabled, these legacy files are disabled by default
unless their legacy path keys are explicitly set.

### 13.3 H5MD Particle Data

| H5MD path | Type/shape | Source |
|---|---|---|
| `/particles/all/position/step` | `int64[n_frame]` | Trajectory write step. |
| `/particles/all/position/time` | `float64[n_frame]` | Trajectory write time. |
| `/particles/all/position/value` | `float32[n_frame,N,3]` | `crd`. |
| `/particles/all/velocity/step` | `int64[n_frame]` | Velocity trajectory write step. |
| `/particles/all/velocity/time` | `float64[n_frame]` | Velocity trajectory write time. |
| `/particles/all/velocity/value` | `float32[n_frame,N,3]` | `vel`. |
| `/particles/all/force/step` | `int64[n_frame]` | Force trajectory write step. |
| `/particles/all/force/time` | `float64[n_frame]` | Force trajectory write time. |
| `/particles/all/force/value` | `float32[n_frame,N,3]` | `frc`. |
| `/particles/all/box/edges/step` | `int64[n_frame]` | Box trajectory write step. |
| `/particles/all/box/edges/time` | `float64[n_frame]` | Box trajectory write time. |
| `/particles/all/box/edges/value` | `float32[n_frame,3,3]` | `box`, converted from length/angle. |

The H5MD `particles` group is reserved for physical particle time series. It
must not contain SPONGE-specific logs, mdin text, restart policy, or force-field
tables.

### 13.4 `mdout` / `Step_Print` Scalar Observables

Every active `Step_Print` column should be treated as an observable. Column
names may be dynamic, so the file must store both the original SPONGE label and
the sanitized HDF5 dataset name.

Recommended schema:

| HDF5 path | Type/shape | Meaning |
|---|---|---|
| `/observables/all/<name>/step` | `int64[n_mdout_frame]` | Step index for this observable stream. |
| `/observables/all/<name>/time` | `float64[n_mdout_frame]` | Time for this observable stream. |
| `/observables/all/<name>/value` | `float64[n_mdout_frame]` or `float64[n_mdout_frame,k]` | Scalar or small-vector observable values. |
| `/parameters/sponge/mdout/columns/original_name` | string array | Original `Step_Print` labels. |
| `/parameters/sponge/mdout/columns/hdf5_name` | string array | Sanitized dataset names. |
| `/parameters/sponge/mdout/columns/category` | string array | Optional category such as `core`, `bonded`, `bias`, `qc`, or `reaxff`. |

Currently observed scalar categories include:

| Category | Observed columns |
|---|---|
| Core run state | `step`, `time`, `temperature`, `frame`, `density`, `potential`, `eff_pot`. |
| Pressure | `pressure`, `Pxx`, `Pyy`, `Pzz`, `Pxy`, `Pxz`, `Pyz`. |
| Nonbonded | `LJ_short`, `LJ_long`, `LJ`, `LJ_soft`, `LJ_soft_short`, `LJ_soft_long`, `PM`, `PM_direct`, `PM_reciprocal`, `PM_self`, `PM_correction`, `Coulomb`, `GB`. |
| Bonded | `bond`, `angle`, `urey_bradley`, `dihedral`, `improper_dihedral`, `cmap`, `nb14_LJ`, `nb14_EE`, plus module-specific bonded aliases. |
| Many-body and custom forces | `SW`, `EAM`, `EDIP`, `TERSOFF`, custom force section names, pairwise force names, and soft wall force names. |
| Restraint and CV | `restrain`, dynamic CV print names, `steer_cv`, `restrain_cv`. |
| Metadynamics | `meta`, `rbias`, `rct`. |
| SITS | Dynamic SITS enhancing energy, bias, and factor names. |
| Quantum chemistry | `QC`, `QC_S_sq`. |
| ReaxFF | `REAXFF_BOND`, `REAXFF_VDW`, `REAXFF_EEQ`, `REAXFF_ELP`, `REAXFF_OVUN`, `REAXFF_ANG`, `REAXFF_PEN`, `REAXFF_COA`, `REAXFF_TOR`, `REAXFF_CONJ`, `REAXFF_HB`, `REAXFF`. |

`step` and `time` may be represented as standard H5MD step/time arrays rather
than duplicated as ordinary observable values. If the original `mdout` table is
also embedded for byte-level provenance, place it below
`/parameters/sponge/mdout/raw_text`; it must not become the canonical numeric
owner.

### 13.5 Thermostat and Barostat Outputs

| Current output | Existing control/key | Current format | New contract |
|---|---|---|---|
| Nose-Hoover chain coordinate trajectory | `nose_hoover_chain_crd` | Text rows of chain coordinates. | `/observables/all/thermostat/nose_hoover_chain/coordinate/value` as `float32[n_frame,n_chain]`. |
| Nose-Hoover chain velocity trajectory | `nose_hoover_chain_vel` | Text rows of chain velocities. | `/observables/all/thermostat/nose_hoover_chain/velocity/value` as `float32[n_frame,n_chain]`. |
| Nose-Hoover chain restart | `nose_hoover_chain_restart_output` | Text restart state. | `restart.spgr.h5:/parameters/restart/thermostat/nose_hoover_chain`; only provenance belongs in output. |
| Barostat internal state | Future or implementation-specific. | Runtime state. | `restart.spgr.h5`; observable pressure and box history remain in `*.spg.h5md`. |

Thermostat/barostat internal continuation state is restart state. Only
time-series diagnostics belong in the output bundle.

### 13.6 Enhanced Sampling and Bias Outputs

Enhanced sampling has two different kinds of output: observables and state.
Observables belong in `*.spg.h5md`; state needed to continue a run belongs in
`restart.spgr.h5`.

| Current output | Existing control/key | Current format | New contract |
|---|---|---|---|
| CV print values | CV `print` settings and dynamic names. | `mdout` columns. | `/observables/all/cv/<name>/value` or `/observables/all/<safe_column_name>/value`. |
| Restraint energy | `restrain` `Step_Print`. | `mdout` column. | `/observables/all/restrain/value`. |
| Steered CV diagnostics | `steer_cv`. | `mdout` column. | `/observables/all/steer_cv/value` or per-CV subgroup if dimensionality requires it. |
| Restrained CV diagnostics | `restrain_cv`. | `mdout` column. | `/observables/all/restrain_cv/value` or per-CV subgroup if dimensionality requires it. |
| Metadynamics scalar energy | `meta`. | `mdout` column. | `/observables/all/metadynamics/meta/value`. |
| Metadynamics reweighting bias | `rbias`. | `mdout` column. | `/observables/all/metadynamics/rbias/value`. |
| Metadynamics `rct` | `rct`. | `mdout` column. | `/observables/all/metadynamics/rct/value`. |
| Metadynamics potential export | `meta_potential_out_file`, default `Meta_Potential.txt`. | Text grid or scatter potential. | `/parameters/sponge/metadynamics/<name>/potential_export`. Continuation state belongs in `restart.spgr.h5`. |
| Metadynamics direct export | default `Meta_directly.txt`. | Text direct potential/force export. | `/parameters/sponge/metadynamics/<name>/direct_export`. |
| Metadynamics hill log | hardcoded `myhill.log`. | Text append log of hill centers and hill metadata. | `/parameters/sponge/metadynamics/<name>/hills`. |
| Metadynamics hill/history log | hardcoded or derived `history.log`. | Text log. | `/parameters/sponge/metadynamics/<name>/history`. |
| Metadynamics edge/sumhill log | `meta_edge_in_file` or default `sumhill.log`. | Text edge-effect/sumhill output. | `/parameters/sponge/metadynamics/<name>/edge`. |
| SITS nk trajectory | `SITS_nk_traj_file`. | Binary float records. | `/observables/all/sits/<module>/nk/value` as `float32[n_frame,k_numbers]`. |
| SITS nk restart | `SITS_nk_rest_file`. | Text restart state. | `restart.spgr.h5:/parameters/restart/bias/sits/<module>/nk`. |
| SITS scalar diagnostics | Dynamic SITS print names. | `mdout` columns. | `/observables/all/sits/<module>/<quantity>/value`. |

Metadynamics text exports under `/parameters/sponge/metadynamics` are
diagnostics/provenance. They are not the canonical continuation state for a new
H5 bundle. Restartable metadynamics state must be written to typed restart
datasets:

| Restart HDF5 path | Type/shape | Meaning |
|---|---|---|
| `/parameters/restart/bias/meta/<name>/grid/min` | `float32[ndim]` | Bias grid lower bound. |
| `/parameters/restart/bias/meta/<name>/grid/max` | `float32[ndim]` | Bias grid upper bound. |
| `/parameters/restart/bias/meta/<name>/grid/count` | `int64[ndim]` | Bias grid shape. |
| `/parameters/restart/bias/meta/<name>/potential/value` | `float32[grid_count...]` | Accumulated grid potential. |
| `/parameters/restart/bias/meta/<name>/scatter/position` | `float32[n,ndim]` | Scatter/history positions. |
| `/parameters/restart/bias/meta/<name>/scatter/weight` | `float32[n]` | Scatter/history weights. |
| `/parameters/restart/bias/meta/<name>/hills/center` | `float32[n_hill,ndim]` | Hill centers. |
| `/parameters/restart/bias/meta/<name>/hills/height` | `float32[n_hill]` | Hill heights. |
| `/parameters/restart/bias/meta/<name>/hills/sigma` | `float32[n_hill,ndim]` | Hill widths. |
| `/parameters/restart/bias/meta/<name>/history/<field>` | typed datasets | Module-specific typed history fields. |

Wall bias/protocol outputs should be represented as protocol observables when
they are time series. Wall definitions remain in `protocol.spgp.h5`; wall
diagnostic values belong in `*.spg.h5md`.

### 13.7 Quantum Chemistry Outputs

| Current output | Existing control/key | Current format | New contract |
|---|---|---|---|
| Quantum chemistry energy | `QC`. | `mdout` column. | `/observables/all/qc/energy/value`. |
| Spin square | `QC_S_sq`. | `mdout` column when enabled. | `/observables/all/qc/spin_square/value`. |
| SCF detailed output | `qc_scf_output`. | Text SCF iteration log. | `/parameters/sponge/qc/scf_output`. |

The SCF text log is diagnostic/provenance data. It is not a standard H5MD
observable unless it is converted into explicit numeric iteration datasets in a
future schema.

### 13.8 ReaxFF Outputs

| Current output | Existing control/key | Current format | New contract |
|---|---|---|---|
| ReaxFF energy terms | ReaxFF `Step_Print`. | `mdout` columns. | `/observables/all/reaxff/<term>/value`. |
| EEQ charge snapshot | hardcoded `eeq_charges.txt`. | Text file overwritten with latest charges. | `/parameters/sponge/reaxff/eeq_charges/value` as a fixed-shape `float32[N]` latest debug snapshot with unit `e`. |

If future output policy records EEQ charge history, use
`/observables/all/reaxff/eeq_charges/value` with shape `float32[n_frame,N]`.
The current native implementation preserves latest-snapshot semantics: it
pre-creates the fixed-shape dataset before SWMR publication and overwrites it
at the mdout cadence. A bundle-enabled run writes the snapshot into every
enabled trajectory/observable H5 sink (the VDS wrapper for VDS mode) and does
not create `eeq_charges.txt`; a legacy-only run retains the text file.

### 13.9 SPONGE Parameters Extension

The following groups are reserved for SPONGE-specific non-H5MD-standard data:

| HDF5 path | Type/shape | Meaning |
|---|---|---|
| `/parameters/sponge/schema_version` | string scalar | SPONGE output extension schema version. |
| `/parameters/sponge/files` | group/table | Original sidecar filenames, if produced or imported. |
| `/parameters/sponge/log/mdinfo_text` | string scalar or string array | Embedded `mdinfo` text. |
| `/parameters/sponge/mdout/columns` | group/table | `mdout` column metadata and name mapping. |
| `/parameters/sponge/mdout/raw_text` | optional string scalar | Raw `mdout` text for provenance only. |
| `/parameters/sponge/restart_exports` | group/table | Legacy restart export paths and intervals. |
| `/parameters/sponge/metadynamics` | group | Potential, direct, hills, edge, and history text/table exports. |
| `/parameters/sponge/qc` | group | QC diagnostic logs. |
| `/parameters/sponge/reaxff` | group | ReaxFF diagnostic snapshots. |
| `/parameters/sponge/thermostat` | group | Non-continuation thermostat diagnostics, if not represented as observables. |

No SPONGE extension data should be placed at H5MD top level. Do not introduce
`/sponge` or `/vmd_structure`.
