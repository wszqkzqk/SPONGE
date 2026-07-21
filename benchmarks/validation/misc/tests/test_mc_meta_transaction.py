import ctypes
import json
import math
import os
import shlex
import shutil
import stat
import struct
import subprocess
import sys
import time
from pathlib import Path

import pytest

BOX_LENGTH = 10.0
PAIR_DISTANCE = 2.0
RESTRAINT_WEIGHT = 2.0
RESTRAINT_REFERENCE = 1.5
MC_RATIO = 0.05
REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
META_IO_PROBE_SOURCE = Path(__file__).with_name("meta_io_transaction_probe.cpp")
META_ATOMIC_PROBE_SOURCE = Path(__file__).with_name(
    "meta_atomic_output_probe.cpp"
)


def _meta_probe_compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = shutil.which("c++")
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the meta I/O probes")
    return [compiler]


def _build_meta_probe(
    tmp_path_factory, name, source, *, extra_flags=(), needs_dependencies=False
):
    build_dir = tmp_path_factory.mktemp(name)
    executable = build_dir / name
    command = [
        *_meta_probe_compiler_command(),
        "-std=c++17",
        "-pthread",
        *extra_flags,
        "-DUSE_CPU",
        f"-I{REPOSITORY_ROOT / 'SPONGE'}",
    ]
    if needs_dependencies:
        command.append(f"-I{REPOSITORY_ROOT / '.pixi/envs/dev-cpu/include'}")
    result = subprocess.run(
        [*command, str(source), "-o", str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    return executable, result


@pytest.fixture(scope="module")
def meta_io_transaction_probe(tmp_path_factory):
    executable, result = _build_meta_probe(
        tmp_path_factory, "meta_io_transaction_probe", META_IO_PROBE_SOURCE
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return executable


@pytest.fixture(scope="module")
def meta_atomic_output_probe(tmp_path_factory):
    executable, result = _build_meta_probe(
        tmp_path_factory,
        "meta_atomic_output_probe",
        META_ATOMIC_PROBE_SOURCE,
        needs_dependencies=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    return executable


@pytest.fixture(scope="module")
def meta_atomic_output_tsan_probe(tmp_path_factory):
    if sys.platform == "win32":
        pytest.skip("ThreadSanitizer probe requires POSIX threads")
    executable, result = _build_meta_probe(
        tmp_path_factory,
        "meta_atomic_output_tsan_probe",
        META_ATOMIC_PROBE_SOURCE,
        extra_flags=("-fsanitize=thread", "-g", "-O1"),
        needs_dependencies=True,
    )
    if result.returncode != 0:
        pytest.skip(f"ThreadSanitizer is unavailable: {result.stderr}")
    return executable


def _write_counted_values(path, values):
    path.write_text(
        str(len(values))
        + "\n"
        + "\n".join(f"{value:.12g}" for value in values)
        + "\n",
        encoding="utf-8",
    )


def _proposal_and_draw(initial_ratio):
    libc = ctypes.CDLL(None)
    libc.srand(0)
    proposal = 2.0 * libc.rand() / 2_147_483_647 - 1.0
    acceptance_draw = libc.rand() / 2_147_483_647
    return 1.0 + proposal * initial_ratio, acceptance_draw


def _proposal_scale(initial_ratio):
    return _proposal_and_draw(initial_ratio)[0]


def _float32(value):
    return struct.unpack("=f", struct.pack("=f", value))[0]


def _float32_ulp(value):
    value = _float32(value)
    bits = struct.unpack("=I", struct.pack("=f", value))[0]
    adjacent = struct.unpack("=f", struct.pack("=I", bits + 1))[0]
    return abs(adjacent - value)


def _forcing_pressure(scale, *, accept):
    delta_volume = BOX_LENGTH**3 * (scale**3 - 1.0)
    magnitude = 1.0e8
    return (
        -math.copysign(magnitude, delta_volume)
        if accept
        else math.copysign(magnitude, delta_volume)
    )


def _write_case(
    case_dir,
    *,
    target_pressure,
    initial_ratio,
    positional_weight=10.0,
    pm_mpi_size=1,
    target_temperature=300.0,
    meta_extra_lines=(),
    scatter_cv_point=None,
    schedule_lines=(),
):
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", (12.0, 12.0))
    _write_counted_values(case_dir / "charge.txt", (0.0, 0.0))
    (case_dir / "coordinate.txt").write_text(
        "\n".join(
            (
                "2",
                "4 5 5",
                "6 5 5",
                "10 10 10",
                "90 90 90",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    (case_dir / "velocity.txt").write_text(
        "2\n0 0 0\n0 0 0\n", encoding="utf-8"
    )
    (case_dir / "restrain_atoms.txt").write_text("0\n1\n", encoding="utf-8")
    distance_block = [
        "distance",
        "{",
        "    CV_type = distance",
        "    atom = 0 1",
    ]
    if scatter_cv_point is not None:
        distance_block.append(f"    CV_point = {scatter_cv_point}")
    distance_block.append("}")
    (case_dir / "cv.txt").write_text(
        "\n".join(
            (
                *distance_block,
                "print",
                "{",
                "    CV = distance",
                "}",
                "restrain",
                "{",
                "    CV = distance",
                f"    weight = {RESTRAINT_WEIGHT}",
                f"    reference = {RESTRAINT_REFERENCE}",
                "}",
                "meta",
                "{",
                "    Ndim = 1",
                "    CV = distance",
                "    CV_minimal = 0",
                "    CV_maximum = 5",
                "    CV_period = 0",
                "    CV_grid = 500",
                "    CV_sigma = 0.2",
                "    height = 100",
                "    potential_update_interval = 2",
                *meta_extra_lines,
                "}",
            )
        )
        + "\n",
        encoding="utf-8",
    )

    settings = {
        "md_name": case_dir.name,
        "mode": "npt",
        "step_limit": 1,
        # Keep the post-force integration displacement below one float ulp so
        # the coordinate trajectory is also an exact restore check.
        "dt": 1.0e-8,
        "cutoff": 4.0,
        "skin": 0.01,
        "PM.MPI_size": pm_mpi_size,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "velocity_in_file": "velocity.txt",
        "cv_in_file": "cv.txt",
        "restrain_atom_id": "restrain_atoms.txt",
        "restrain_refcoord_scaling": "all",
        "restrain_single_weight": positional_weight,
        "thermostat": "berendsen_thermostat",
        "thermostat_tau": 1.0,
        "barostat": "monte_carlo_barostat",
        "target_temperature": target_temperature,
        "target_pressure": target_pressure,
        "monte_carlo_barostat_initial_ratio": initial_ratio,
        "monte_carlo_barostat_update_interval": 2,
        "monte_carlo_barostat_check_interval": 100,
        "monte_carlo_barostat_couple_dimension": "XYZ",
        "mdout": "mdout.txt",
        "crd": "mdcrd.dat",
        "frc": "frc.dat",
        "box": "mdbox.txt",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
        "dont_check_input": 1,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + ("\n" + "\n".join(schedule_lines) if schedule_lines else "")
        + "\n",
        encoding="utf-8",
    )


def _run_case(
    case_dir,
    *,
    target_pressure,
    initial_ratio,
    positional_weight=10.0,
    mpi_np=None,
    target_temperature=300.0,
    meta_extra_lines=(),
    scatter_cv_point=None,
    schedule_lines=(),
):
    _write_case(
        case_dir,
        target_pressure=target_pressure,
        initial_ratio=initial_ratio,
        positional_weight=positional_weight,
        # Keep one reciprocal-space process so MPI runs execute the CV/Meta
        # branch.  A two-rank run has one PP and one PM process and therefore
        # does not claim to exercise multi-PP DD remeshing.
        pm_mpi_size=1,
        target_temperature=target_temperature,
        meta_extra_lines=meta_extra_lines,
        scatter_cv_point=scatter_cv_point,
        schedule_lines=schedule_lines,
    )
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]
    if mpi_np is not None:
        command = [
            "mpirun",
            "--oversubscribe",
            "-np",
            str(mpi_np),
            *command,
        ]
    result = subprocess.run(
        command,
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"SPONGE failed with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    mdout_lines = (
        (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    )
    mdout = {
        name: float(value)
        for name, value in zip(mdout_lines[0].split(), mdout_lines[1].split())
    }
    box = tuple(
        float(value)
        for value in (case_dir / "mdbox.txt")
        .read_text(encoding="utf-8")
        .split()[:3]
    )
    coordinate_bytes = (case_dir / "mdcrd.dat").read_bytes()
    coordinates = struct.unpack(
        f"={len(coordinate_bytes) // 4}f", coordinate_bytes
    )[:6]
    force_bytes = (case_dir / "frc.dat").read_bytes()
    forces = struct.unpack(f"={len(force_bytes) // 4}f", force_bytes)[:6]
    hills = [
        line.split()
        for line in (case_dir / "myhill.log")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]
    return mdout, box, coordinates, forces, hills


def _write_high_dimensional_cv(path, *, ndim, use_scatter, grid_points=2):
    names = [f"distance_{index}" for index in range(ndim)]
    blocks = []
    for name in names:
        entries = ["    CV_type = distance", "    atom = 0 1"]
        if use_scatter:
            entries.append("    CV_point = 2.5")
        blocks.extend((name, "{", *entries, "}"))
    repeated = lambda value: " ".join([str(value)] * ndim)
    blocks.extend(
        (
            "meta",
            "{",
            f"    Ndim = {ndim}",
            f"    CV = {' '.join(names)}",
            f"    CV_minimal = {repeated(0)}",
            f"    CV_maximum = {repeated(5)}",
            f"    CV_period = {repeated(0)}",
            f"    CV_grid = {repeated(grid_points)}",
            f"    CV_sigma = {repeated(1.0)}",
            "    height = 1",
            "    potential_update_interval = 1",
        )
    )
    if use_scatter:
        blocks.append("    scatter = 1")
    blocks.append("}")
    path.write_text("\n".join(blocks) + "\n", encoding="utf-8")


def _run_sponge_process(case_dir):
    return subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


def _build_fsync_path_substitution_interposer(build_dir):
    source = build_dir / "replace_temp_after_fsync.c"
    library = build_dir / (
        "replace_temp_after_fsync.dylib"
        if sys.platform == "darwin"
        else "replace_temp_after_fsync.so"
    )
    source.write_text(
        r"""
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <fcntl.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef int (*fsync_function)(int);

static int injected = 0;

static int descriptor_path(int descriptor, char *path, size_t capacity)
{
#ifdef __APPLE__
    (void)capacity;
    return fcntl(descriptor, F_GETPATH, path);
#else
    char link_name[64];
    int written = snprintf(link_name, sizeof(link_name), "/proc/self/fd/%d",
                           descriptor);
    if (written < 0 || (size_t)written >= sizeof(link_name)) return -1;
    ssize_t path_size = readlink(link_name, path, capacity - 1);
    if (path_size < 0) return -1;
    path[path_size] = '\0';
    return 0;
#endif
}

static int injected_fsync(int descriptor)
{
#ifdef __APPLE__
    int status = fsync(descriptor);
#else
    static fsync_function real_fsync = NULL;
    if (real_fsync == NULL)
    {
        real_fsync = (fsync_function)dlsym(RTLD_NEXT, "fsync");
        if (real_fsync == NULL)
        {
            errno = EIO;
            return -1;
        }
    }
    int status = real_fsync(descriptor);
#endif
    if (status != 0 || injected) return status;
    const char *sentinel = getenv("SPONGE_TEST_TEMP_SENTINEL");
    if (sentinel == NULL) return status;
    char path[PATH_MAX];
    if (descriptor_path(descriptor, path, sizeof(path)) != 0 ||
        strstr(path, "/.sponge-tmp.") == NULL)
    {
        return status;
    }
    if (unlink(path) == 0 && symlink(sentinel, path) == 0) injected = 1;
    return status;
}

#ifdef __APPLE__
__attribute__((used)) static struct
{
    const void *replacement;
    const void *replacee;
} fsync_interpose __attribute__((section("__DATA,__interpose"))) = {
    (const void *)injected_fsync, (const void *)fsync};
#else
int fsync(int descriptor) { return injected_fsync(descriptor); }
#endif
""",
        encoding="utf-8",
    )
    command = ["cc"]
    if sys.platform == "darwin":
        command.extend(
            (
                "-dynamiclib",
                "-undefined",
                "dynamic_lookup",
                str(source),
                "-o",
                str(library),
            )
        )
    else:
        command.extend(
            ("-shared", "-fPIC", str(source), "-ldl", "-o", str(library))
        )
    result = subprocess.run(
        command, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        pytest.skip(f"cannot build fsync fault interposer: {result.stderr}")
    return library


def _build_rename_failure_path_substitution_interposer(build_dir):
    source = build_dir / "replace_temp_during_rename.c"
    library = build_dir / (
        "replace_temp_during_rename.dylib"
        if sys.platform == "darwin"
        else "replace_temp_during_rename.so"
    )
    source.write_text(
        r"""
#define _GNU_SOURCE
#include <dlfcn.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

typedef int (*renameat_function)(int, const char *, int, const char *);

static int injected = 0;

static int injected_renameat(int old_directory, const char *old_path,
                             int new_directory, const char *new_path)
{
    static renameat_function real_renameat = NULL;
    if (real_renameat == NULL)
    {
        real_renameat = (renameat_function)dlsym(RTLD_NEXT, "renameat");
        if (real_renameat == NULL)
        {
            errno = EIO;
            return -1;
        }
    }
    const char *sentinel = getenv("SPONGE_TEST_TEMP_SENTINEL");
    if (!injected && sentinel != NULL && old_path != NULL &&
        strncmp(old_path, ".sponge-tmp.", 12) == 0 &&
        unlinkat(old_directory, old_path, 0) == 0 &&
        symlinkat(sentinel, old_directory, old_path) == 0)
    {
        injected = 1;
        errno = EIO;
        return -1;
    }
    return real_renameat(old_directory, old_path, new_directory, new_path);
}

#ifdef __APPLE__
__attribute__((used)) static struct
{
    const void *replacement;
    const void *replacee;
} renameat_interpose __attribute__((section("__DATA,__interpose"))) = {
    (const void *)injected_renameat, (const void *)renameat};
#else
int renameat(int old_directory, const char *old_path, int new_directory,
             const char *new_path)
{
    return injected_renameat(old_directory, old_path, new_directory, new_path);
}
#endif
""",
        encoding="utf-8",
    )
    command = ["cc"]
    if sys.platform == "darwin":
        command.extend(
            (
                "-dynamiclib",
                "-undefined",
                "dynamic_lookup",
                str(source),
                "-o",
                str(library),
            )
        )
    else:
        command.extend(
            ("-shared", "-fPIC", str(source), "-ldl", "-o", str(library))
        )
    result = subprocess.run(
        command, capture_output=True, text=True, check=False
    )
    if result.returncode != 0:
        pytest.skip(f"cannot build renameat fault interposer: {result.stderr}")
    return library


def _write_two_point_potential(path, representation, *, malformed=False):
    final_token = "1junk" if malformed else "1"
    path.write_text(
        "\n".join(
            (
                f"SPONGE_META_POTENTIAL_V1 1 {representation} d_force",
                "0 2 1",
                "2 2",
                "0.5 0.5 0 0.5",
                f"1.5 1 0 {final_token}",
            )
        )
        + "\n",
        encoding="utf-8",
    )


@pytest.mark.parametrize("use_scatter", [False, True])
def test_meta_supports_more_than_eight_dimensions(tmp_path, use_scatter):
    case_dir = tmp_path / ("scatter_9d" if use_scatter else "grid_9d")
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    _write_high_dimensional_cv(
        case_dir / "cv.txt", ndim=9, use_scatter=use_scatter
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    hills = [
        line.split()
        for line in (case_dir / "myhill.log")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]
    assert len(hills) == 2
    assert all(len(hill) >= 9 for hill in hills)
    mdout_lines = (
        (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    )
    names = mdout_lines[0].split()
    final_values = {
        name: float(value)
        for name, value in zip(names, mdout_lines[-1].split())
    }
    offset = 0.5 if use_scatter else 0.75
    expected_potential = math.exp(-0.5 * (offset / 1.0) ** 2 * 9)
    expected_potential_float = _float32(expected_potential)
    assert final_values["meta"] == float(f"{expected_potential_float:.6f}")

    force_values = struct.unpack(
        f"={(case_dir / 'frc.dat').stat().st_size // 4}f",
        (case_dir / "frc.dat").read_bytes(),
    )
    final_forces = force_values[-6:]
    grid_coordinate = 2.5 if use_scatter else 1.25
    expected_atom_0_x = _float32(
        -9.0 * (grid_coordinate - PAIR_DISTANCE) * expected_potential_float
    )
    force_tolerance = 4.0 * _float32_ulp(expected_atom_0_x)
    assert final_forces == pytest.approx(
        (expected_atom_0_x, 0.0, 0.0, -expected_atom_0_x, 0.0, 0.0),
        abs=force_tolerance,
    )


def test_meta_rejects_overflowing_grid_product(tmp_path):
    case_dir = tmp_path / "overflowing_meta_grid"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    _write_high_dimensional_cv(
        case_dir / "cv.txt",
        ndim=2,
        use_scatter=False,
        grid_points=50_000,
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode != 0
    assert "grid shape overflows its int index space" in (
        result.stdout + result.stderr
    )


@pytest.mark.parametrize(
    ("representation", "input_keyword"),
    [("grid", "potential_in_file"), ("scatter", "scatter_in_file")],
)
def test_meta_reweighting_uses_only_real_partition_terms(
    tmp_path, representation, input_keyword
):
    case_dir = tmp_path / f"two_point_{representation}"
    potential_name = "two_point.meta"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=(f"    {input_keyword} = {potential_name}",),
    )
    _write_two_point_potential(
        case_dir / potential_name, representation=representation
    )
    cv_path = case_dir / "cv.txt"
    cv_path.write_text(
        cv_path.read_text(encoding="utf-8").replace(
            "    height = 100", "    height = 0"
        ),
        encoding="utf-8",
    )

    result = _run_sponge_process(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    mdout_lines = (
        (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    )
    names = mdout_lines[0].split()
    final = dict(zip(names, map(float, mdout_lines[-1].split())))

    thermal_energy = 0.00198716 * 300.0
    bias_factor = _float32(1.0e9)
    beta_f_plus_v = 1.0 / ((bias_factor - 1.0) * thermal_energy)
    beta_f = bias_factor * beta_f_plus_v

    def logsumexp(values):
        maximum = max(values)
        return maximum + math.log(
            sum(math.exp(value - maximum) for value in values)
        )

    potentials = (0.5, 1.0)
    expected_rct = thermal_energy * (
        logsumexp([beta_f * value for value in potentials])
        - logsumexp([beta_f_plus_v * value for value in potentials])
    )
    assert final["rct"] == float(f"{_float32(expected_rct):.6f}")


def test_potential_input_is_strict_and_transactional(tmp_path):
    case_dir = tmp_path / "malformed_potential"
    potential_name = "malformed.meta"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=(f"    scatter_in_file = {potential_name}",),
    )
    _write_two_point_potential(
        case_dir / potential_name, representation="scatter", malformed=True
    )

    result = _run_sponge_process(case_dir)
    assert result.returncode != 0
    assert "invalid potential/force value token '1junk'" in (
        result.stdout + result.stderr
    )
    assert not (case_dir / "sumhill.log").exists()
    assert not (case_dir / "lnbias.dat").exists()


def test_legacy_direct_potential_accepts_writer_rounding(tmp_path):
    case_dir = tmp_path / "legacy_direct_rounding"
    potential_name = "legacy_three_point.meta"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=(f"    potential_in_file = {potential_name}",),
    )
    fixture = (
        os.path.dirname(__file__) + "/../statics/legacy_meta_three_point.txt"
    )
    shutil.copyfile(fixture, case_dir / potential_name)

    result = _run_sponge_process(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr


def test_v1_potential_rejects_inconsistent_grid_metadata(tmp_path):
    case_dir = tmp_path / "inconsistent_potential_grid"
    potential_name = "inconsistent.meta"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=(f"    scatter_in_file = {potential_name}",),
    )
    _write_two_point_potential(
        case_dir / potential_name, representation="scatter"
    )
    potential_path = case_dir / potential_name
    contents = potential_path.read_text(encoding="utf-8")
    potential_path.write_text(
        contents.replace("0 2 1\n", "0 3 1\n"), encoding="utf-8"
    )

    result = _run_sponge_process(case_dir)
    assert result.returncode != 0
    assert "bounds and 2-point extent require 1.5" in (
        result.stdout + result.stderr
    )
    assert not (case_dir / "lnbias.dat").exists()


def test_edge_cache_round_trips_without_recomputation(tmp_path):
    case_dir = tmp_path / "edge_restart"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        scatter_cv_point=2.5,
        meta_extra_lines=("    scatter = 1",),
    )

    first = _run_sponge_process(case_dir)
    assert first.returncode == 0, first.stdout + first.stderr
    edge_path = case_dir / "sumhill.log"
    first_edge = edge_path.read_bytes()
    assert first_edge.startswith(b"SPONGE_META_EDGE_V1 ")
    assert (
        (case_dir / "lnbias.dat")
        .read_text(encoding="utf-8")
        .startswith("SPONGE_META_SCATTER_LOG_V1 ")
    )

    second = _run_sponge_process(case_dir)
    assert second.returncode == 0, second.stdout + second.stderr
    assert edge_path.read_bytes() == first_edge


def test_implicit_edge_cache_recomputes_when_configuration_changes(tmp_path):
    case_dir = tmp_path / "stale_implicit_edge"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        scatter_cv_point=2.5,
        meta_extra_lines=("    scatter = 1",),
    )

    first = _run_sponge_process(case_dir)
    assert first.returncode == 0, first.stdout + first.stderr
    edge_path = case_dir / "sumhill.log"
    first_edge = edge_path.read_bytes()

    cv_path = case_dir / "cv.txt"
    contents = cv_path.read_text(encoding="utf-8")
    assert "    CV_point = 2.5" in contents
    cv_path.write_text(
        contents.replace("    CV_point = 2.5", "    CV_point = 1.5"),
        encoding="utf-8",
    )

    second = _run_sponge_process(case_dir)
    assert second.returncode == 0, second.stdout + second.stderr
    assert edge_path.read_bytes() != first_edge


def test_explicit_edge_cache_rejects_configuration_mismatch(tmp_path):
    case_dir = tmp_path / "stale_explicit_edge"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        scatter_cv_point=2.5,
        meta_extra_lines=("    scatter = 1",),
    )

    first = _run_sponge_process(case_dir)
    assert first.returncode == 0, first.stdout + first.stderr
    edge_path = case_dir / "sumhill.log"
    first_edge = edge_path.read_bytes()

    cv_path = case_dir / "cv.txt"
    contents = cv_path.read_text(encoding="utf-8")
    contents = contents.replace("    CV_point = 2.5", "    CV_point = 1.5")
    contents = contents.replace(
        "    scatter = 1", "    scatter = 1\n    edge_in_file = sumhill.log"
    )
    cv_path.write_text(contents, encoding="utf-8")

    second = _run_sponge_process(case_dir)
    assert second.returncode != 0
    diagnostic = second.stdout + second.stderr
    assert "active configuration requires" in diagnostic
    assert "fingerprint=" in diagnostic
    assert edge_path.read_bytes() == first_edge


def test_edge_logsumexp_has_no_finite_sentinel(tmp_path):
    case_dir = tmp_path / "edge_below_old_sentinel"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        scatter_cv_point=0.0,
        meta_extra_lines=("    scatter = 1",),
    )
    cv_path = case_dir / "cv.txt"
    contents = cv_path.read_text(encoding="utf-8")
    contents = contents.replace("    CV_maximum = 5", "    CV_maximum = 5e15")
    cv_path.write_text(contents, encoding="utf-8")

    result = _run_sponge_process(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    last_record = (
        (case_dir / "sumhill.log").read_text(encoding="utf-8").splitlines()[-1]
    )
    coordinate, stored_lse = map(float, last_record.split())
    inverse_sigma = 1.0 / 0.2
    expected_lse = _float32(-0.5 * (coordinate * inverse_sigma) ** 2)
    assert coordinate > 4.9e15
    assert stored_lse == pytest.approx(
        expected_lse, abs=4.0 * _float32_ulp(expected_lse)
    )
    assert stored_lse < -1.0e30


def test_unrepresentable_edge_logsumexp_has_no_output_side_effects(tmp_path):
    case_dir = tmp_path / "edge_log_overflow"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        scatter_cv_point=0.0,
        meta_extra_lines=("    scatter = 1",),
    )
    cv_path = case_dir / "cv.txt"
    contents = cv_path.read_text(encoding="utf-8")
    contents = contents.replace("    CV_maximum = 5", "    CV_maximum = 1e20")
    cv_path.write_text(contents, encoding="utf-8")

    result = _run_sponge_process(case_dir)
    assert result.returncode != 0
    assert "stable double log-space" in result.stdout + result.stderr
    assert not (case_dir / "sumhill.log").exists()
    assert not (case_dir / "lnbias.dat").exists()


def test_hills_parser_rejects_trailing_garbage_before_history_output(tmp_path):
    case_dir = tmp_path / "malformed_hills"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=("    sumhill_freq = 1",),
    )
    (case_dir / "myhill.log").write_text("2 1junk\n", encoding="utf-8")

    result = _run_sponge_process(case_dir)
    assert result.returncode != 0
    assert "invalid hill height token '1junk'" in result.stdout + result.stderr
    assert not (case_dir / "history.log").exists()


def test_incomplete_hills_tail_is_ignored_and_repaired(tmp_path):
    case_dir = tmp_path / "incomplete_hills_tail"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    first = _run_sponge_process(case_dir)
    assert first.returncode == 0, first.stdout + first.stderr

    hills_path = case_dir / "myhill.log"
    complete_prefix = hills_path.read_bytes()
    assert complete_prefix.endswith(b"\n")
    with hills_path.open("ab") as output:
        output.write(b"2.0 100 partial-record")

    second = _run_sponge_process(case_dir)
    assert second.returncode == 0, second.stdout + second.stderr
    repaired = hills_path.read_bytes()
    assert repaired.startswith(complete_prefix)
    assert repaired.endswith(b"\n")
    assert b"partial-record" not in repaired


def test_grid_sink_history_replay_is_supported(tmp_path):
    case_dir = tmp_path / "grid_sink_restart"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=("    sink = 1", "    sumhill_freq = 1"),
    )

    first = _run_sponge_process(case_dir)
    assert first.returncode == 0, first.stdout + first.stderr
    second = _run_sponge_process(case_dir)
    assert second.returncode == 0, second.stdout + second.stderr
    history_rows = [
        line.split()
        for line in (case_dir / "history.log")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.strip()
    ]
    assert history_rows
    assert all(len(row) == 4 for row in history_rows)


def test_meta_output_path_is_unbounded_and_preserves_hash(tmp_path):
    case_dir = tmp_path / "long_meta_path"
    segments = [f"segment {index}_" + "x" * 80 for index in range(6)]
    relative_prefix = "/".join((*segments, "potential#roundtrip"))
    assert len(relative_prefix) > 512
    assert " " in relative_prefix
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        schedule_lines=(
            f"default_out_file_prefix = {json.dumps(relative_prefix)}",
        ),
    )
    case_dir.joinpath(*segments).mkdir(parents=True)

    result = _run_sponge_process(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (case_dir / f"{relative_prefix}_Meta_Potential.txt").is_file()


def test_meta_atomic_output_supports_near_name_max_target(tmp_path):
    case_dir = tmp_path / "near_name_max"
    output_suffix = "_Meta_Potential.txt"
    output_prefix = "p" * 230
    assert len(output_prefix + output_suffix) == 249
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        schedule_lines=(f'default_out_file_prefix = "{output_prefix}"',),
    )

    result = _run_sponge_process(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    assert (case_dir / f"{output_prefix}{output_suffix}").is_file()


@pytest.mark.skipif(
    not hasattr(os, "symlink") or sys.platform == "win32",
    reason="requires POSIX exec PID and symbolic-link semantics",
)
def test_meta_atomic_output_does_not_reuse_predictable_temporary_symlink(
    tmp_path,
):
    case_dir = tmp_path / "atomic_output_symlink"
    output_prefix = "atomic_meta"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        schedule_lines=(f'default_out_file_prefix = "{output_prefix}"',),
    )
    output_path = case_dir / f"{output_prefix}_Meta_Potential.txt"
    output_path.write_text("previous complete output\n", encoding="utf-8")
    sentinel = case_dir / "temporary-symlink-target.txt"
    sentinel.write_text("must remain unchanged\n", encoding="utf-8")

    pid_path = case_dir / "launcher.pid"
    gate_path = case_dir / "launcher.gate"
    launcher = (
        "import os, pathlib, sys, time; "
        "pathlib.Path(sys.argv[1]).write_text(str(os.getpid()), "
        "encoding='ascii'); "
        "gate = pathlib.Path(sys.argv[2]); "
        "deadline = time.monotonic() + 10.0; "
        "\nwhile not gate.exists():\n"
        "    if time.monotonic() >= deadline: raise SystemExit(124)\n"
        "    time.sleep(0.01)\n"
        "os.execvp(sys.argv[3], [sys.argv[3], '-mdin', 'mdin.spg.toml'])"
    )
    sponge_binary = os.environ.get("SPONGE_BIN", "SPONGE")
    process = subprocess.Popen(
        [
            sys.executable,
            "-c",
            launcher,
            str(pid_path),
            str(gate_path),
            sponge_binary,
        ],
        cwd=case_dir,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    deadline = time.monotonic() + 5.0
    while not pid_path.exists() and time.monotonic() < deadline:
        time.sleep(0.01)
    assert pid_path.exists(), "launcher did not publish its pid"
    child_pid = int(pid_path.read_text(encoding="ascii"))
    reserved_temporary = case_dir / f".sponge-tmp.{child_pid}.0"
    try:
        reserved_temporary.symlink_to(sentinel.name)
    except OSError as error:
        gate_path.touch()
        process.communicate(timeout=120)
        pytest.skip(f"cannot create a test symlink: {error}")
    gate_path.touch()
    stdout, stderr = process.communicate(timeout=120)

    assert process.returncode == 0, stdout + stderr
    assert (
        output_path.read_text(encoding="utf-8") != "previous complete output\n"
    )
    assert sentinel.read_text(encoding="utf-8") == "must remain unchanged\n"
    assert reserved_temporary.is_symlink()
    assert set(case_dir.glob(".sponge-tmp.*")) == {reserved_temporary}


@pytest.mark.skipif(sys.platform == "win32", reason="checks POSIX mode bits")
def test_meta_atomic_replace_preserves_existing_target_mode(tmp_path):
    case_dir = tmp_path / "atomic_output_mode"
    output_prefix = "private_meta"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        schedule_lines=(f'default_out_file_prefix = "{output_prefix}"',),
    )
    output_path = case_dir / f"{output_prefix}_Meta_Potential.txt"
    output_path.write_text("previous complete output\n", encoding="utf-8")
    output_path.chmod(0o600)

    result = _run_sponge_process(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    assert stat.S_IMODE(output_path.stat().st_mode) == 0o600


@pytest.mark.skipif(
    sys.platform == "win32", reason="uses POSIX fsync interposition"
)
def test_meta_atomic_output_rejects_post_open_path_substitution(tmp_path):
    case_dir = tmp_path / "atomic_output_identity"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    baseline = _run_sponge_process(case_dir)
    assert baseline.returncode == 0, baseline.stdout + baseline.stderr

    sentinel = case_dir / "attacker-owned-target.txt"
    sentinel.write_text("must remain unchanged\n", encoding="utf-8")
    interposer_dir = tmp_path / "interposer"
    interposer_dir.mkdir()
    interposer = _build_fsync_path_substitution_interposer(interposer_dir)
    environment = os.environ.copy()
    environment["SPONGE_TEST_TEMP_SENTINEL"] = str(sentinel.resolve())
    if sys.platform == "darwin":
        environment["DYLD_INSERT_LIBRARIES"] = str(interposer)
        environment["DYLD_FORCE_FLAT_NAMESPACE"] = "1"
    else:
        environment["LD_PRELOAD"] = str(interposer)

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "verify the identity of" in diagnostic
    substituted_paths = list(case_dir.glob(".sponge-tmp.*"))
    assert len(substituted_paths) == 1
    assert substituted_paths[0].is_symlink()
    assert substituted_paths[0].resolve() == sentinel.resolve()
    assert sentinel.read_text(encoding="utf-8") == "must remain unchanged\n"


@pytest.mark.skipif(
    sys.platform == "win32", reason="uses POSIX renameat interposition"
)
def test_meta_atomic_cleanup_rechecks_identity_after_failed_rename(tmp_path):
    case_dir = tmp_path / "atomic_output_failed_rename_identity"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    baseline = _run_sponge_process(case_dir)
    assert baseline.returncode == 0, baseline.stdout + baseline.stderr
    output_path = case_dir / "Meta_Potential.txt"
    previous_output = output_path.read_bytes()

    sentinel = case_dir / "attacker-owned-rename-target.txt"
    sentinel.write_text("must remain unchanged\n", encoding="utf-8")
    interposer_dir = tmp_path / "rename_interposer"
    interposer_dir.mkdir()
    interposer = _build_rename_failure_path_substitution_interposer(
        interposer_dir
    )
    environment = os.environ.copy()
    environment["SPONGE_TEST_TEMP_SENTINEL"] = str(sentinel.resolve())
    if sys.platform == "darwin":
        environment["DYLD_INSERT_LIBRARIES"] = str(interposer)
        environment["DYLD_FORCE_FLAT_NAMESPACE"] = "1"
    else:
        environment["LD_PRELOAD"] = str(interposer)

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "Unable to atomically replace" in diagnostic
    assert output_path.read_bytes() == previous_output
    substituted_paths = list(case_dir.glob(".sponge-tmp.*"))
    assert len(substituted_paths) == 1
    assert substituted_paths[0].is_symlink()
    assert substituted_paths[0].resolve() == sentinel.resolve()
    assert sentinel.read_text(encoding="utf-8") == "must remain unchanged\n"


@pytest.mark.skipif(
    sys.platform == "win32", reason="uses POSIX descriptor cleanup"
)
@pytest.mark.parametrize(
    "mode",
    (
        "commit_fatal",
        "close_fatal",
        "close_close",
        "claimed_close_close",
    ),
)
def test_meta_atomic_cleanup_record_interleavings(
    meta_atomic_output_probe, tmp_path, mode
):
    result = subprocess.run(
        [str(meta_atomic_output_probe), mode, str(tmp_path / f"{mode}.out")],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.skipif(sys.platform == "win32", reason="uses ThreadSanitizer")
@pytest.mark.parametrize(
    "mode",
    (
        "commit_fatal",
        "close_fatal",
        "close_close",
        "claimed_close_close",
    ),
)
def test_meta_atomic_cleanup_record_interleavings_are_tsan_clean(
    meta_atomic_output_tsan_probe, tmp_path, mode
):
    environment = os.environ.copy()
    environment["TSAN_OPTIONS"] = "halt_on_error=1"
    result = subprocess.run(
        [
            str(meta_atomic_output_tsan_probe),
            mode,
            str(tmp_path / f"tsan-{mode}.out"),
        ],
        env=environment,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.skipif(sys.platform == "win32", reason="uses POSIX RLIMIT_FSIZE")
def test_hills_partial_payload_write_rolls_back_complete_prefix(tmp_path):
    import resource
    import signal

    case_dir = tmp_path / "hills_short_write"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    file_limit = 1 << 20
    complete_prefix = (
        b"existing uninterpreted prefix"
        + b"x" * (file_limit - 4 - len("existing uninterpreted prefix"))
        + b"\n"
    )
    assert len(complete_prefix) == file_limit - 3
    hills_path = case_dir / "myhill.log"
    hills_path.write_bytes(complete_prefix)

    def constrain_output_size():
        signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
        resource.setrlimit(resource.RLIMIT_FSIZE, (file_limit, file_limit))

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
        preexec_fn=constrain_output_size,
    )
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "record payload" in diagnostic
    assert "previous complete prefix was restored" in diagnostic
    assert hills_path.read_bytes() == complete_prefix


def test_hills_loader_skips_large_unterminated_tail_without_line_allocation(
    tmp_path,
):
    case_dir = tmp_path / "large_incomplete_hills_tail"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=("    sumhill_freq = 1",),
    )
    hills_path = case_dir / "myhill.log"
    hills_path.write_bytes(b"2 100\n" + b"x" * (8 << 20))

    result = _run_sponge_process(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr
    assert "Ignoring incomplete final hills record" in result.stdout
    repaired = hills_path.read_bytes()
    assert repaired.startswith(b"2 100\n")
    assert repaired.endswith(b"\n")
    assert b"x" * 1024 not in repaired


def test_hills_loader_rejects_overlong_committed_record(tmp_path):
    case_dir = tmp_path / "overlong_committed_hills_record"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=("    sumhill_freq = 1",),
    )
    (case_dir / "myhill.log").write_bytes(b"1 " + b"0" * 5000 + b"\n")

    result = _run_sponge_process(case_dir)
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "configuration-derived limit" in diagnostic
    assert not (case_dir / "history.log").exists()


@pytest.mark.skipif(sys.platform == "win32", reason="uses a POSIX FIFO")
def test_hills_loader_rejects_fifo_without_blocking(tmp_path):
    case_dir = tmp_path / "fifo_hills_input"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=("    sumhill_freq = 1",),
    )
    os.mkfifo(case_dir / "myhill.log")

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "is not a regular file" in diagnostic


@pytest.mark.skipif(sys.platform == "win32", reason="uses a POSIX FIFO")
def test_hills_append_rejects_fifo_without_blocking(tmp_path):
    case_dir = tmp_path / "fifo_hills_output"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    os.mkfifo(case_dir / "myhill.log")

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=10,
    )
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "is not a stable regular file" in diagnostic


@pytest.mark.skipif(sys.platform == "win32", reason="uses POSIX file types")
@pytest.mark.parametrize("operation", ("append", "load"))
def test_hills_io_rejects_directory(tmp_path, operation):
    case_dir = tmp_path / f"directory_hills_{operation}"
    extra_lines = ("    sumhill_freq = 1",) if operation == "load" else ()
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=extra_lines,
    )
    (case_dir / "myhill.log").mkdir()

    result = _run_sponge_process(case_dir)
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    expected = (
        "is not a stable regular file"
        if operation == "append"
        else "is not a regular file"
    )
    assert expected in diagnostic


@pytest.mark.skipif(
    sys.platform == "win32" or not hasattr(os, "symlink"),
    reason="uses POSIX no-follow semantics",
)
def test_hills_append_rejects_symlink_without_touching_target(tmp_path):
    case_dir = tmp_path / "symlink_hills_output"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    sentinel = case_dir / "sentinel.txt"
    sentinel.write_text("must remain unchanged\n", encoding="utf-8")
    hills_path = case_dir / "myhill.log"
    hills_path.symlink_to(sentinel.name)

    result = _run_sponge_process(case_dir)
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "is not a stable regular file" in diagnostic
    assert hills_path.is_symlink()
    assert sentinel.read_text(encoding="utf-8") == "must remain unchanged\n"


@pytest.mark.skipif(sys.platform == "win32", reason="uses POSIX file locks")
def test_hills_publication_survives_second_reserve_failure(
    meta_io_transaction_probe,
):
    result = subprocess.run(
        [str(meta_io_transaction_probe), "publish_failure"],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode == 0, result.stdout + result.stderr


@pytest.mark.skipif(sys.platform == "win32", reason="uses POSIX file locks")
def test_same_process_threads_serialize_hills_append_and_load(
    meta_io_transaction_probe, tmp_path
):
    result = subprocess.run(
        [
            str(meta_io_transaction_probe),
            "thread_io",
            str(tmp_path / "hills.log"),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_concurrent_hill_appends_leave_only_complete_records(tmp_path):
    case_dir = tmp_path / "concurrent_hill_append"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]
    processes = [
        subprocess.Popen(
            command,
            cwd=case_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        for _ in range(4)
    ]
    results = [process.communicate(timeout=120) for process in processes]
    for process, (stdout, stderr) in zip(processes, results):
        assert process.returncode == 0, stdout + stderr
    hills_bytes = (case_dir / "myhill.log").read_bytes()
    assert hills_bytes.endswith(b"\n")
    records = [
        line.split() for line in hills_bytes.splitlines() if line.strip()
    ]
    assert len(records) == 4
    assert all(len(record) == 2 for record in records)


@pytest.mark.skipif(sys.platform == "win32", reason="uses POSIX advisory locks")
def test_hills_loader_waits_for_writer_lock_before_snapshot(tmp_path):
    import fcntl

    case_dir = tmp_path / "locked_hills_snapshot"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        meta_extra_lines=("    sumhill_freq = 1",),
    )
    hills_path = case_dir / "myhill.log"
    hills_path.write_text("2 100\n", encoding="ascii")
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]

    with hills_path.open("r+b", buffering=0) as locked_hills:
        fcntl.lockf(locked_hills, fcntl.LOCK_EX)
        process = subprocess.Popen(
            command,
            cwd=case_dir,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        time.sleep(0.5)
        assert process.poll() is None
        assert not (case_dir / "history.log").exists()
        fcntl.lockf(locked_hills, fcntl.LOCK_UN)

    stdout, stderr = process.communicate(timeout=120)
    assert process.returncode == 0, stdout + stderr
    assert (case_dir / "history.log").is_file()


def test_meta_atomic_replace_failure_preserves_existing_target(tmp_path):
    case_dir = tmp_path / "atomic_output_replace_failure"
    output_prefix = "blocked_meta"
    _write_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        schedule_lines=(f'default_out_file_prefix = "{output_prefix}"',),
    )
    output_path = case_dir / f"{output_prefix}_Meta_Potential.txt"
    output_path.mkdir()
    marker = output_path / "previous-output-marker"
    marker.write_text("preserve me\n", encoding="utf-8")

    result = _run_sponge_process(case_dir)
    diagnostic = result.stdout + result.stderr
    assert result.returncode != 0, diagnostic
    assert "Unable to atomically replace" in diagnostic
    assert marker.read_text(encoding="utf-8") == "preserve me\n"
    assert not list(case_dir.glob(".sponge-tmp.*"))


@pytest.mark.parametrize(
    ("run_tag", "mutate", "diagnostic"),
    [
        (
            "write_interval_zero",
            "write_interval",
            "write_information_interval must be positive",
        ),
        (
            "potential_interval_zero",
            "potential_interval",
            "potential_update_interval must be positive",
        ),
        (
            "runtime_inverse_temperature_overflow",
            "temperature",
            "inverse thermal energy is not representable",
        ),
    ],
)
def test_invalid_meta_runtime_contract_is_rejected(
    tmp_path, run_tag, mutate, diagnostic
):
    case_dir = tmp_path / run_tag
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    if mutate == "potential_interval":
        path = case_dir / "cv.txt"
        contents = path.read_text(encoding="utf-8")
        assert "potential_update_interval = 2" in contents
        path.write_text(
            contents.replace(
                "potential_update_interval = 2",
                "potential_update_interval = 0",
            ),
            encoding="utf-8",
        )
    else:
        path = case_dir / "mdin.spg.toml"
        contents = path.read_text(encoding="utf-8")
        if mutate == "write_interval":
            original = "write_information_interval = 1"
            replacement = "write_information_interval = 0"
        else:
            original = "target_temperature = 300.0"
            replacement = "target_temperature = 1e-37"
        assert original in contents
        path.write_text(
            contents.replace(original, replacement), encoding="utf-8"
        )

    if mutate in {"write_interval", "potential_interval"}:
        cv_path = case_dir / "cv.txt"
        cv_contents = cv_path.read_text(encoding="utf-8")
        assert "    atom = 0 1" in cv_contents
        assert "    CV_sigma = 0.2" in cv_contents
        cv_path.write_text(
            cv_contents.replace(
                "    atom = 0 1",
                "    atom = 0 1\n    CV_point = 2.0",
                1,
            ).replace(
                "    CV_sigma = 0.2",
                "    CV_sigma = 0.2\n    scatter = 1",
            ),
            encoding="utf-8",
        )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode != 0
    assert diagnostic in result.stdout + result.stderr
    if mutate in {"write_interval", "potential_interval"}:
        assert not (case_dir / "sumhill.log").exists()
        assert not (case_dir / "lnbias.dat").exists()


@pytest.mark.parametrize("height", [1.0e30, -1.0e30])
def test_well_tempered_height_overflow_or_underflow_is_rejected(
    tmp_path, height
):
    case_dir = tmp_path / (
        "height_underflow" if height > 0 else "height_overflow"
    )
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    cv_path = case_dir / "cv.txt"
    contents = cv_path.read_text(encoding="utf-8")
    contents = contents.replace("    height = 100", f"    height = {height}")
    contents = contents.replace(
        "    potential_update_interval = 2",
        "    potential_update_interval = 1\n    welltemp_factor = 1.0000001",
    )
    cv_path.write_text(contents, encoding="utf-8")

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode != 0
    assert "unrepresentable well-tempered hill height" in (
        result.stdout + result.stderr
    )


def test_extreme_well_tempered_factor_rejects_subnormal_derived_beta(tmp_path):
    case_dir = tmp_path / "subnormal_derived_beta"
    _write_case(case_dir, target_pressure=1.0, initial_ratio=0.0)
    cv_path = case_dir / "cv.txt"
    contents = cv_path.read_text(encoding="utf-8")
    contents = contents.replace(
        "    potential_update_interval = 2",
        "    potential_update_interval = 2\n    welltemp_factor = 3e38",
    )
    cv_path.write_text(contents, encoding="utf-8")

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode != 0
    assert "representable as positive normal floats" in (
        result.stdout + result.stderr
    )


def test_step_zero_temperature_schedule_initializes_meta(tmp_path):
    case_dir = tmp_path / "meta_step_zero_schedule"
    _, _, _, _, hills = _run_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        target_temperature=1.0e-37,
        schedule_lines=(
            'target_temperature_schedule_mode = "step"',
            "target_temperature_schedule_steps = [{step = 0, value = 300.0}]",
        ),
    )
    assert len(hills) == 1


def test_mc_trial_evaluations_do_not_commit_mask_history(tmp_path, mpi_np):
    case_dir = tmp_path / "mask_history_transaction"
    _, _, _, _, hills = _run_case(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        mpi_np=mpi_np,
        scatter_cv_point=2.5,
        meta_extra_lines=(
            "    scatter = 1",
            "    mask = 1",
            "    max_force = 1e-30",
            "    sink = 1",
        ),
    )
    assert len(hills) == 1
    # center, height, sink maximum/index, scatter index, committed exit_tag
    assert len(hills[0]) == 6
    assert float(hills[0][-1]) == pytest.approx(1.0, abs=1.0e-7)


def _assert_final_diagnostics(mdout, forces, distance):
    expected_energy = RESTRAINT_WEIGHT * (distance - RESTRAINT_REFERENCE) ** 2
    expected_force = 2.0 * RESTRAINT_WEIGHT * (distance - RESTRAINT_REFERENCE)
    assert mdout["distance"] == pytest.approx(distance, abs=5.1e-5)
    assert mdout["meta"] == pytest.approx(0.0, abs=1.0e-7)
    # Positional restraint references participate in the same box
    # transaction: accepted trials keep their once-scaled references and
    # rejected trials restore the exact snapshot.
    assert mdout["restrain"] == pytest.approx(0.0, abs=1.0e-7)
    assert mdout["restrain_cv"] == pytest.approx(expected_energy, abs=5.1e-3)
    assert mdout["eff_pot"] == pytest.approx(expected_energy, abs=3.0e-6)
    assert forces == pytest.approx(
        (expected_force, 0.0, 0.0, -expected_force, 0.0, 0.0),
        abs=2.0e-5,
    )


def test_zero_volume_mc_and_meta_commit_one_history_sample(tmp_path, mpi_np):
    case_dir = tmp_path / "zero_volume"
    mdout, box, coordinates, forces, hills = _run_case(
        case_dir, target_pressure=1.0, initial_ratio=0.0, mpi_np=mpi_np
    )

    # A zero proposal has identical old/trial Hamiltonians.  In particular,
    # the hill deposited by the final committed force cannot leak backward
    # into its acceptance test or be deposited by either trial evaluation.
    assert len(hills) == 1
    assert float(hills[0][0]) == pytest.approx(PAIR_DISTANCE, abs=2.0e-6)
    assert box == pytest.approx((BOX_LENGTH,) * 3, abs=1.0e-7)
    assert coordinates == pytest.approx(
        (4.0, 5.0, 5.0, 6.0, 5.0, 5.0), abs=1.0e-7
    )
    _assert_final_diagnostics(mdout, forces, PAIR_DISTANCE)


def test_accepted_mc_trial_commits_meta_at_trial_cv(tmp_path, mpi_np):
    scale = _proposal_scale(MC_RATIO)
    expected_distance = PAIR_DISTANCE * scale
    case_dir = tmp_path / "accepted"
    mdout, box, coordinates, forces, hills = _run_case(
        case_dir,
        target_pressure=_forcing_pressure(scale, accept=True),
        initial_ratio=MC_RATIO,
        mpi_np=mpi_np,
    )

    assert len(hills) == 1
    # This checks the force-evaluation cache boundary: old and trial share a
    # physical step, but the committed hill and printed CV must use trial crd.
    assert float(hills[0][0]) == pytest.approx(expected_distance, abs=2.0e-6)
    assert box == pytest.approx((BOX_LENGTH * scale,) * 3, abs=1.1e-6)
    assert coordinates == pytest.approx(
        (
            4.0 * scale,
            5.0 * scale,
            5.0 * scale,
            6.0 * scale,
            5.0 * scale,
            5.0 * scale,
        ),
        abs=2.0e-6,
    )
    _assert_final_diagnostics(mdout, forces, expected_distance)


def test_rejected_mc_trial_restores_exact_box_coordinates_and_cv(
    tmp_path, mpi_np
):
    scale = _proposal_scale(MC_RATIO)
    case_dir = tmp_path / "rejected"
    mdout, box, coordinates, forces, hills = _run_case(
        case_dir,
        target_pressure=_forcing_pressure(scale, accept=False),
        initial_ratio=MC_RATIO,
        mpi_np=mpi_np,
    )

    assert len(hills) == 1
    # The old `g = -g` path left L*(1+a)*(1-a).  The rejected transaction now
    # restores the authoritative cell and global-order coordinates exactly,
    # even though this skin would make the proposal eligible for remeshing.
    assert box == pytest.approx((BOX_LENGTH,) * 3, abs=1.0e-7)
    assert coordinates == pytest.approx(
        (4.0, 5.0, 5.0, 6.0, 5.0, 5.0), abs=1.0e-7
    )
    assert float(hills[0][0]) == pytest.approx(PAIR_DISTANCE, abs=2.0e-6)
    _assert_final_diagnostics(mdout, forces, PAIR_DISTANCE)
    assert all(math.isfinite(value) for value in (*box, *coordinates, *forces))


def test_mc_acceptance_evaluates_scaled_positional_references(tmp_path, mpi_np):
    scale, acceptance_draw = _proposal_and_draw(MC_RATIO)
    thermal_energy = 0.00198716 * 300.0
    old_cv_energy = (
        RESTRAINT_WEIGHT * (PAIR_DISTANCE - RESTRAINT_REFERENCE) ** 2
    )
    new_cv_energy = (
        RESTRAINT_WEIGHT * (PAIR_DISTANCE * scale - RESTRAINT_REFERENCE) ** 2
    )
    jacobian = -2.0 * thermal_energy * math.log(scale**3)
    correct_delta = new_cv_energy - old_cv_energy + jacobian
    unscaled_reference_delta = (
        correct_delta + 100.0 * 152.0 * (scale - 1.0) ** 2
    )
    decision_threshold = -thermal_energy * math.log(acceptance_draw)
    assert correct_delta < decision_threshold < unscaled_reference_delta
    case_dir = tmp_path / "scaled_reference_acceptance"
    _, box, _, _, hills = _run_case(
        case_dir,
        target_pressure=0.0,
        initial_ratio=MC_RATIO,
        positional_weight=100.0,
        mpi_np=mpi_np,
    )

    # Coordinates and `refcoord_scaling=all` references undergo the same trial
    # transform, so the positional restraint contributes zero to Delta H and
    # this deterministic proposal is accepted.  The inequality above proves
    # that evaluating trial crd against the unscaled reference rejects the
    # same platform-specific libc draw, whether the proposal expands or
    # contracts the box.
    assert box == pytest.approx((BOX_LENGTH * scale,) * 3, abs=1.1e-6)
    assert len(hills) == 1
    assert float(hills[0][0]) == pytest.approx(
        PAIR_DISTANCE * scale, abs=2.0e-6
    )


def _run_pme_rerun(case_dir, update_interval, mpi_np=None):
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", (12.0, 12.0))
    _write_counted_values(case_dir / "charge.txt", (0.3, -0.3))
    (case_dir / "coordinate.txt").write_text(
        "2\n4 5 5\n6 5 5\n10 10 10\n90 90 90\n", encoding="utf-8"
    )
    trajectory = []
    for distance in (2.0, 2.5, 3.0):
        trajectory.extend(
            (5.0 - 0.5 * distance, 5.0, 5.0, 5.0 + 0.5 * distance, 5.0, 5.0)
        )
    (case_dir / "trajectory.dat").write_bytes(
        struct.pack(f"={len(trajectory)}f", *trajectory)
    )
    (case_dir / "box.txt").write_text(
        "10 10 10 90 90 90\n" * 3, encoding="utf-8"
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "rerun",
        "cutoff": 4.0,
        "PM.MPI_size": 1,
        "PME.update_interval": update_interval,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "crd": "trajectory.dat",
        "box": "box.txt",
        "mdout": "mdout.txt",
        "write_information_interval": 1,
        "write_mdout_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
        "dont_check_input": 1,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]
    if mpi_np is not None:
        command = [
            "mpirun",
            "--oversubscribe",
            "-np",
            str(mpi_np),
            *command,
        ]
    result = subprocess.run(
        command,
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"SPONGE failed with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    names = lines[0].split()
    return [
        {name: float(value) for name, value in zip(names, line.split())}
        for line in lines[1:]
    ]


def test_pme_energy_is_current_between_mts_force_updates(tmp_path, mpi_np):
    every_step = _run_pme_rerun(tmp_path / "pme_every_step", 1, mpi_np)
    mts = _run_pme_rerun(tmp_path / "pme_mts", 3, mpi_np)
    assert len(every_step) == len(mts) == 3
    # PME.update_interval controls reciprocal force impulses, not which
    # coordinates define an explicitly requested Hamiltonian.  Step 1 is not
    # an MTS force-update step and used to report the mesh from step 0.
    assert [row["eff_pot"] for row in mts] == pytest.approx(
        [row["eff_pot"] for row in every_step], abs=3.0e-6
    )
    assert [row["PM"] for row in mts] == pytest.approx(
        [row["PM"] for row in every_step], abs=5.1e-3
    )
    assert abs(mts[1]["eff_pot"] - mts[0]["eff_pot"]) > 1.0e-3
