# SinkMeta/meta Input Parameters

This list is derived from `META::Initial()` in `SPONGE/bias/sinkmeta.cpp` and defaults in `SPONGE/bias/sinkmeta.h`.

## Parameters (module: `meta`)

| Parameter | Type / Shape | Required | Default (if any) | Meaning / Notes |
| --- | --- | --- | --- | --- |
| `CV` | string list | Yes | None | CV module names used by Meta; determines dimensionality. If missing, META is not initialized. |
| `dip` | float | No | `0.0` | Extra dip term for submarine/sink behavior (adds to bias shift via `kB*T`). |
| `welltemp_factor` | float | No | `1e9` | Well-tempered bias factor; an explicit value must be finite, normal, and > 1. |
| `Ndim` | int | No | `CV` count | Explicit dimensionality; must match `CV` list length. |
| `subhill` | flag | No | `false` | Enable sub-hill (Gaussian) behavior. Presence-only (no value read). |
| `kde` | int | No | `0` | Nonzero enables KDE mode and also sets `subhill=true`. Also switches sigma scaling to `1.414/sigma`. |
| `mask` | int | No | `0` | Enable mask mode (n-dim area exit label). |
| `max_force` | float | No | `0.1` | Edge-force criterion for the exit label; only read when `mask` is set and must be a finite positive normal float. |
| `sink` | int | No | `0` | Nonzero enables negative-hill (sink/submarine) behavior. |
| `sumhill_freq` | int | No | `0` | History frequency for `sumhill` accumulation (affects Rbias/RCT). |
| `convmeta` | int | No | `0` | ConvolutionMeta flag; also sets `do_negative=true`. |
| `grw` | int | No | `0` | GRW flag; also sets `do_negative=true`. |
| `CV_period` | float array (ndim) | Yes | None in META | Periodic box length per CV. Always requested in `META::Initial()`. |
| `CV_sigma` | float array (ndim) | Yes | None in META | Gaussian width per CV (must be > 0). Stored internally as inverse sigma. |
| `cutoff` | float array (ndim) | No | `3 * CV_sigma` | Neighbor cutoff for lookup and border wall; if present enables `do_cutoff`. |
| `potential_in_file` | string | No | None | Read potential from file; if set, `Read_Potential()` is called (grid/scatter settings below are bypassed). |
| `scatter_in_file` | string | No | None | Read scatter potential from file; sets `use_scatter=true`, `usegrid=false`, and calls `Read_Potential()`. |
| `edge_in_file` | string | No | `sumhill.log` | Explicit V1 edge-effect cache. A missing, unversioned, malformed, or configuration-mismatched explicit file is rejected. The implicit default cache is regenerated when missing, unversioned, or stale for the active configuration. |
| `scatter` | int | No | `0` | Number of scatter points; if > 0 uses scatter points instead of grid. |
| `CV_minimal` | float array (ndim) | Conditionally | None in META | Grid minimum per CV; required when not using `potential_in_file` or `scatter_in_file`. |
| `CV_maximum` | float array (ndim) | Conditionally | None in META | Grid maximum per CV; must be > `CV_minimal`. Required when not using `potential_in_file` or `scatter_in_file`. |
| `CV_grid` | int array (ndim) | Conditionally | None in META | Grid points per CV; must be > 1. Required when not using `potential_in_file` or `scatter_in_file`. |
| `height` | float | No | `1.0` | Initial hill height (`height_0`). |
| `wall_height` | float | No | None | Enables border wall and sets `border_potential_height`. |
| `potential_out_file` | string | No | `Meta_Potential.txt` | Output file name for writing potential. |
| `potential_update_interval` | int | No | `write_information_interval` or `1000` | Hill-deposition interval; it must be positive. |

## Parameters used from other modules / global

| Parameter | Scope | Type / Shape | Required | Default (if any) | Meaning / Notes |
| --- | --- | --- | --- | --- | --- |
| `CV_point` | Each CV module | float array (scatter size) | Conditionally | None in META | Scatter point coordinates per CV; read only when `scatter > 0`. |
| `write_information_interval` | Controller | int | No | `1000` | Global controller interval; used as the default for `potential_update_interval` and must be positive when Meta is enabled. |

## Notes

- `potential_in_file` and `scatter_in_file` are mutually exclusive; specifying both is an error.
- `CV_sigma` is inverted internally; `kde` mode uses `1.414 / sigma`, otherwise `1.0 / sigma`.
- Default file names set in code: `read_potential_file_name` and `write_potential_file_name` both start as `Meta_Potential.txt`.
- Temperature-dependent bias, history loading, and reweighting use the current top-level `target_temperature` (including the step-zero schedule value), not the instantaneous kinetic temperature. The temperature and all derived inverse-thermal/well-tempered factors must be finite positive normal floats; failures are reported instead of leaving stale reweighting state.
- Meta file paths use dynamically sized strings; whitespace, nested paths longer than 512 bytes, and `#` characters inside quoted TOML paths are preserved.

## Persistent formats

- Direct potential V1 starts with `SPONGE_META_POTENTIAL_V1 <ndim> <grid|scatter> <subhill|d_force>`. The reader checks all four fields, grid shape, record count, payload width, numeric tokens, and grid coordinate order before committing any state. Legacy potential files are accepted only when their dimension, payload, and representation can be determined without ambiguity; an all-grid-point legacy scatter file must be regenerated with a V1 header.
- Edge cache V1 starts with `SPONGE_META_EDGE_V1 <ndim> <grid-size> <scatter-size> <normal_lse|normal_lse_force> <fingerprint>`. The fingerprint covers every grid/scatter value that affects the cache. Records store `normal_lse` directly in log space, followed by log-normalization derivatives when required; they never round-trip through `exp`.
- `lnbias.dat` starts with `SPONGE_META_SCATTER_LOG_V1 <scatter-size>` and stores the scatter index plus `normal_lse`, also without exponentiating it.
- Persisted floats use enough decimal digits to round-trip a `float`. Hills and potential inputs require complete finite numeric tokens; trailing characters, undeclared records, invalid indices, and nonrepresentable derived values are errors.
- `history.log`, edge, scatter-log, and potential outputs are written to completed temporary files and then atomically replace their destinations, so validation failures do not truncate an existing result.
