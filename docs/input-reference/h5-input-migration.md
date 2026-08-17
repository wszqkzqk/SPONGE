# H5 Input Migration

SPONGE's H5 input bundle model keeps large data and restart state in HDF5
containers while keeping `mdin.spg.toml` as the editable binding and launch
policy file.

## Canonical Bindings

Use these TOML tables for H5 bundle inputs:

```toml
[input.h5.topology]
path = "topologies/protein.topology.spgt.h5"

[input.h5.protocol]
path = "protocols/metadyn.protocol.spgp.h5"

[input.h5.restart]
path = "runs/prod_0007.restart.spgr.h5"
load = "structural"
```

For rerun:

```toml
mode = "rerun"
rerun_frame_limit = 1000
rerun_start = 0
rerun_strip = 0

[input.h5.topology]
path = "topologies/protein.topology.spgt.h5"

[input.h5.protocol]
path = "protocols/analysis.protocol.spgp.h5"

[input.h5.trajectory]
path = "runs/prod.spg.h5md"
particle_stream = "all"
```

## Legacy Input Mapping

| Legacy key | H5 bundle replacement |
|---|---|
| `coordinate_in_file` | `input.h5.restart.path`, structural position state |
| `velocity_in_file` | `input.h5.restart.path`, structural velocity state |
| `amber_rst7` / `rst7` | `input.h5.restart.path` |
| `crd` | `input.h5.trajectory.path`, position frames |
| `box` | `input.h5.trajectory.path`, box edge frames |
| `vel` | `input.h5.trajectory.path`, optional velocity frames |
| `frame_limit` | keep as top-level flat key, or use `rerun_frame_limit` |
| `rerun_start` | unchanged top-level flat key |
| `rerun_strip` | unchanged top-level flat key |
| `rerun_need_box_update` | unchanged top-level flat key |

For structural restart input, `/particles/all/velocity/value` is optional. If
it is absent, SPONGE initializes every atomic velocity component to zero before
copying the launch state to the device. Thermostats may introduce stochastic
motion during later integration steps, but the input assembly stage does not
sample a Maxwell-Boltzmann distribution.

H5 and legacy inputs for the same role are mutually exclusive. For example,
`input_h5_restart_path` must not be combined with `coordinate_in_file`,
`velocity_in_file`, or `rst7`; `input_h5_trajectory_path` must not be combined
with rerun `crd`, `box`, or `vel`.

## Restart Load Policy

`input.h5.restart.load = "structural"` restores coordinates, velocity, box,
step, and time. `dynamic` additionally restores compatible integrator,
thermostat, barostat, and RNG state; `protocol` restores compatible SITS,
MetaD, restraint-reference, and CV-reference continuation state; `full`
composes both. Each enabled module validates ownership and schema compatibility
before applying state. `custom` remains reserved until explicit component-list
keys are defined.
