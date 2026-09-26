# Vulkan Backend Porting Guide

How to port a SPONGE module to the Vulkan compute backend. The bond module
(`SPONGE/bond/bond.cpp` + `SPONGE/vulkan_backend/glsl/bond_force.comp`)
is the reference implementation — read it first.

## Architecture

- `SPONGE/third_party/device_backend/vulkan_api.h` — macro surface included from
  `common.h` when `USE_VULKAN` is defined. Provides the device API
  (deviceMalloc/deviceMemcpy/deviceMemset/deviceFree/streams), `VK_LAUNCH`,
  and the `Launch_Device_Kernel` fallback (runs the C++ kernel on CPU after a
  `sponge_vk::HostBarrier()`).
- `SPONGE/vulkan_backend/vulkan_runtime.cpp` — runtime: device
  init, host-visible coherent buffers (device pointer == host mapped pointer),
  stream-ordered command recording, glslang-based runtime GLSL compilation.
- `SPONGE/vulkan_backend/glsl/*.comp` — GLSL compute kernels.
  Every `.comp` file is one kernel named by its file name (e.g.
  `bond_force.comp` → kernel `bond_force`). At build time
  `cmake/utils/embed_vulkan_shaders.cmake` embeds all sources; at runtime each
  kernel is compiled with `#version 450` + device-feature defines + the
  contents of `glsl/sponge_common.glsl` prepended.
- Unported kernels keep working: `Launch_Device_Kernel` falls back to the CPU
  (OpenMP) version of the same C++ kernel. Ported kernels get an explicit
  `#ifdef USE_VULKAN ... #else ... #endif` dispatch at the launch site.

## Kernel dispatch pattern (per launch site)

```cpp
#ifdef USE_VULKAN
    MyVkParams params{...};                 // plain scalars + Boundary only
    const void* buffers[] = {d_a, d_b, d_c, vk_or_dummy(maybe_null, d_a)};
    VK_LAUNCH(my_kernel, grid_x, 1, block_x, 1, buffers, &params, NULL);
#else
    Launch_Device_Kernel(My_Kernel, grid, block, 0, NULL, ...);
#endif
```

- `VK_LAUNCH(name, grid_x, grid_y, block_x, block_y, buffers, params, stream)`:
  `buffers` is a C array of device pointers bound as storage buffers 0..N-1 in
  declaration order; `params` is a pointer to a POD struct uploaded as push
  constants. Grid Z is always 1; block is (block_x, block_y, 1).
- Every buffer entry must be a valid allocation — substitute a dummy
  (`vk_or_dummy`, copy it from bond.cpp) for nullable pointers. The GLSL side
  must not read buffers disabled by flag params.
- Push constant struct: only `int`/`float`/`Boundary`/fixed arrays. NO
  `VECTOR`/`float3` members (std430 vec3 has alignment 16, C++ VECTOR is 12).
  Use three floats or six-float `LTMatrix3`-style members instead.
- Match the std430 layout exactly: `Boundary` is `uint` + 12 floats = 52
  bytes and matches C++ `Boundary` (policy@0, cell@4, rcell@28). Add a
  `static_assert(sizeof(MyVkParams) == N, ...)` next to the struct.

## GLSL conventions

Common helpers available in every kernel (from `sponge_common.glsl`):

- Types: `Vec3` (float x,y,z — 12 bytes, matches `VECTOR`), `LTMat`
  (a11,a21,a22,a31,a32,a33 — matches `LTMatrix3`), `Boundary`.
- Vector ops: `v3_add/sub/scale/mul/neg/dot/len/floor`, `v3_mul_lt` (row-vector
  times lower-triangular matrix, same as C++ `operator*(VECTOR, LTMatrix3)`).
- PBC: `get_displacement(a, b, bd)`, `wrap_coordinate(c, bd)`.
- Virial: `virial_from_force_dis(f, dr)`.
- Atomics (MACROS — GLSL forbids atomics through function parameters):
  `afadd(slot, v)` (float atomicAdd), `aadd_v3(slot, v)`, `aadd_lt(slot, m)`.
  For ints use `atomicAdd` directly.
- RNG: `philox_normal4(uvec2 offset, uvec2 stream, uvec2 seed)`,
  `philox_uint32(...)` — bit-compatible with `SPONGE_PHILOX4X32_10`;
  64-bit values are passed as `uvec2(lo, hi)`, helpers `u64_add`, `u64_shr2`.
- Reductions: `block_sum(float)`, `block_sum_int(int)` (subgroup + shared),
  `Kahan`/`kahan_add` for compensated sums.
- `layout(local_size_x_id = 0, local_size_y_id = 1, local_size_z_id = 2) in;`
  is already declared in the common header — do not redeclare.

Buffer declarations:

```glsl
layout(set = 0, binding = 0, std430) buffer Crd { Vec3 crd[]; };
layout(set = 0, binding = 1, std430) buffer K { float bond_k[]; };
layout(push_constant) uniform Params { int n; Boundary boundary; } pc;
```

- `char` arrays (e.g. labels) are read as `uint label[]` with
  `(label[i >> 2] >> ((i & 3) * 8)) & 0xFF` (see `bond_get_local.comp`).
- `vec3`/`vec4` must NOT appear in buffer/push-constant layouts (alignment).
  `Vec3` arrays in std430 have stride 12, matching C++ `VECTOR`.

## GLSL gotchas (all bitten during bring-up)

- `do { ... } while (0)` is illegal: GLSL requires a boolean in `while`.
  Use `do { ... } while (false)`.
- Atomics only on direct buffer/shared l-values: wrap them in macros, never
  in helper functions.
- `if` conditions must be `bool` — no `if (int_flag)`.
- No implicit `int`→`uint` comparisons: cast with `int(i)`/`uint(i)`.
- Integer literals: use `1u`, `0xFFu` suffixes where needed.
- `erfcf`, `sincosf` etc. do not exist — add approximations to
  `sponge_common.glsl` when a module needs them (keep them precise).

## Offline shader check (fast, do this before building)

```bash
ENV=.pixi/envs/dev-vulkan
{ echo "#version 450"; echo "#extension GL_EXT_shader_atomic_float : enable";
  echo "#extension GL_KHR_shader_subgroup_basic : enable";
  echo "#extension GL_KHR_shader_subgroup_arithmetic : enable";
  echo "#extension GL_KHR_shader_subgroup_ballot : enable";
  echo "#extension GL_KHR_shader_subgroup_shuffle : enable";
  echo "#define SPONGE_VK_ATOMIC_FLOAT 1"; echo "#define SPONGE_VK_FP64 0";
  echo "#define SPONGE_VK_SUBGROUP_SIZE 32";
  cat SPONGE/vulkan_backend/glsl/sponge_common.glsl;
  cat SPONGE/vulkan_backend/glsl/my_kernel.comp; } > /tmp/t.comp
$ENV/bin/glslangValidator --target-env vulkan1.1 -V /tmp/t.comp -o /dev/null
```

## Build and end-to-end validation

```bash
LDFLAGS="-B/tmp/apple-ld-shim" pixi run -e dev-vulkan compile 8
```

Reference CPU binary (same branch): `.pixi/envs/dev-cpu/bin/SPONGE`.
Vulkan binary: `.pixi/envs/dev-vulkan/bin/SPONGE`.

Test cases (ready-made):

- `/tmp/sponge_vk_case/mdin.spg.toml` — TIP3P water NVE (exercises
  neighbor_list, LJ, PME, bond, integrator). `mdin50.toml` = 50 steps.
- `/tmp/sponge_ala_case/sponge.mdin` — alanine dipeptide step-0 energies
  (bond, urey_bradley, dihedral, improper_dihedral, nb14_LJ/EE, LJ, PM,
  pressure tensor).

Run both binaries in the case dir (mdout.txt is overwritten each run — save a
copy per backend) and diff the mdout tables. fp32 noise up to ~1e-4 relative
is acceptable; systematic shifts are not.
