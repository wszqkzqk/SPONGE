import os
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("fatal_cleanup_probe.cpp")
MPI_INIT_PROBE_SOURCE = Path(__file__).with_name("fatal_mpi_init_probe.cpp")
FAKE_MPI_INCLUDE = Path(__file__).with_name("fake_mpi")
SPONGE_MALLOC_ERROR = 1008
SPONGE_SIMULATION_BREAKDOWN = 1007
SPONGE_VALUE_ERROR = 1006


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = (
        shutil.which("c++")
        or shutil.which("clang++")
        or shutil.which("g++")
        or shutil.which("clang-cl")
        or shutil.which("cl")
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the fatal-cleanup probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
        candidates.append(
            Path(os.environ["CONDA_PREFIX"]) / "Library" / "include"
        )
    candidates.append(
        REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include"
    )
    candidates.append(
        REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "Library" / "include"
    )
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (
            candidate / "fftw3.h"
        ).is_file():
            return candidate
    pytest.skip(
        "OpenMP and FFTW headers are required for the fatal-cleanup probe"
    )


def _thread_flags():
    return [] if sys.platform == "win32" else ["-pthread"]


def _dead_strip_flags():
    if sys.platform == "darwin":
        return ["-Wl,-dead_strip"]
    if sys.platform == "win32":
        return []
    return ["-Wl,--gc-sections"]


def _is_msvc_driver(compiler_command):
    driver = Path(compiler_command[0]).name.lower()
    return driver in {"cl", "cl.exe", "clang-cl", "clang-cl.exe"}


def _msvc_flags(flags):
    translated = []
    for flag in flags:
        if flag.startswith("-D"):
            translated.append(f"/D{flag[2:]}")
        elif flag.startswith("-I"):
            translated.append(f"/I{flag[2:]}")
        else:
            translated.append(flag)
    return translated


def _executable_path(build_dir, name):
    suffix = ".exe" if sys.platform == "win32" else ""
    return build_dir / f"{name}{suffix}"


def _build_probe(tmp_path_factory, name, extra_flags=()):
    build_dir = tmp_path_factory.mktemp(name)
    executable = _executable_path(build_dir, name)
    compiler = _compiler_command()
    if _is_msvc_driver(compiler):
        command = [
            *compiler,
            "/nologo",
            "/std:c++17",
            "/EHsc",
            "/utf-8",
            "/DUSE_CPU",
            "/DSPONGE_FATAL_CLEANUP_TESTING",
            "/w",
            *_msvc_flags(extra_flags),
            f"/I{REPOSITORY_ROOT / 'SPONGE'}",
            f"/I{_dependency_include()}",
            str(PROBE_SOURCE),
            f"/Fo{build_dir}{os.sep}",
            f"/Fe{executable}",
        ]
    else:
        command = [
            *compiler,
            "-std=c++17",
            "-DUSE_CPU",
            "-DSPONGE_FATAL_CLEANUP_TESTING",
            "-w",
            *_thread_flags(),
            *extra_flags,
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(PROBE_SOURCE),
            "-o",
            str(executable),
        ]
    compile_result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile fatal-cleanup probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )
    return executable


def _build_mpi_init_probe(tmp_path_factory, name, fake_mpi_flags=()):
    build_dir = tmp_path_factory.mktemp(name)
    executable = _executable_path(build_dir, name)
    compiler = _compiler_command()
    if _is_msvc_driver(compiler):
        command = [
            *compiler,
            "/nologo",
            "/std:c++17",
            "/EHsc",
            "/utf-8",
            "/DUSE_CPU",
            "/DUSE_MPI",
            '/DSPONGE_VERSION_STR="test"',
            "/Gy",
            "/w",
            *_msvc_flags(fake_mpi_flags),
            f"/I{FAKE_MPI_INCLUDE}",
            f"/I{REPOSITORY_ROOT / 'SPONGE'}",
            f"/I{_dependency_include()}",
            str(MPI_INIT_PROBE_SOURCE),
            str(REPOSITORY_ROOT / "SPONGE" / "control.cpp"),
            f"/Fo{build_dir}{os.sep}",
            f"/Fd{build_dir / (name + '.pdb')}",
            f"/Fe{executable}",
            "/link",
            "/OPT:REF",
        ]
    else:
        command = [
            *compiler,
            "-std=c++17",
            "-DUSE_CPU",
            "-DUSE_MPI",
            '-DSPONGE_VERSION_STR="test"',
            "-ffunction-sections",
            "-fdata-sections",
            "-w",
            *_thread_flags(),
            *fake_mpi_flags,
            f"-I{FAKE_MPI_INCLUDE}",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(MPI_INIT_PROBE_SOURCE),
            str(REPOSITORY_ROOT / "SPONGE" / "control.cpp"),
            *_dead_strip_flags(),
            "-o",
            str(executable),
        ]
    compile_result = subprocess.run(
        command,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile MPI-init probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )
    return executable


@pytest.fixture(scope="module")
def fatal_cleanup_probe(tmp_path_factory):
    return _build_probe(tmp_path_factory, "fatal_cleanup_probe")


@pytest.fixture(scope="module")
def fatal_cleanup_mpi_probe(tmp_path_factory):
    return _build_probe(
        tmp_path_factory,
        "fatal_cleanup_mpi_probe",
        ("-DUSE_MPI", f"-I{FAKE_MPI_INCLUDE}"),
    )


@pytest.fixture(scope="module")
def mpi_init_success_probe(tmp_path_factory):
    return _build_mpi_init_probe(tmp_path_factory, "mpi_init_success_probe")


@pytest.fixture(scope="module")
def mpi_init_failure_probe(tmp_path_factory):
    return _build_mpi_init_probe(
        tmp_path_factory,
        "mpi_init_failure_probe",
        ("-DFAKE_MPI_INIT_STATUS=17", "-DFAKE_MPI_INITIALIZED=0"),
    )


@pytest.fixture(scope="module")
def mpi_init_insufficient_probe(tmp_path_factory):
    return _build_mpi_init_probe(
        tmp_path_factory,
        "mpi_init_insufficient_probe",
        ("-DFAKE_MPI_PROVIDED=2", "-DFAKE_MPI_INITIALIZED=1"),
    )


def _run_fatal_mode(
    executable,
    mode,
    marker_path,
    expected_error="spongeErrorValueErrorCommand",
):
    result = subprocess.run(
        [str(executable), mode, str(marker_path)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode != 0
    assert expected_error in result.stdout + result.stderr
    return result


def test_registered_callback_runs_on_fatal_from_another_controller(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "fatal.marker"
    _run_fatal_mode(fatal_cleanup_probe, "fatal", marker)
    assert marker.read_bytes() == b"F"


def test_registered_callback_runs_before_mpi_abort(
    fatal_cleanup_mpi_probe, tmp_path
):
    marker = tmp_path / "mpi-fatal.marker"
    result = _run_fatal_mode(fatal_cleanup_mpi_probe, "fatal", marker)
    expected_code = (
        SPONGE_VALUE_ERROR
        if sys.platform == "win32"
        else SPONGE_VALUE_ERROR & 0xFF
    )
    assert result.returncode == expected_code
    assert "FAKE_MPI_ABORT_RETURNED" in result.stdout + result.stderr
    assert marker.read_bytes() == b"F"


def test_mpi_thread_initialization_accepts_required_support(
    mpi_init_success_probe,
):
    result = subprocess.run(
        [str(mpi_init_success_probe)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode == 0


def test_failed_mpi_init_never_calls_mpi_abort(mpi_init_failure_probe):
    result = subprocess.run(
        [str(mpi_init_failure_probe)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    expected_code = (
        SPONGE_SIMULATION_BREAKDOWN
        if sys.platform == "win32"
        else SPONGE_SIMULATION_BREAKDOWN & 0xFF
    )
    diagnostic = result.stdout + result.stderr
    assert result.returncode == expected_code
    assert "MPI_Init_thread failed with status 17" in diagnostic
    assert "FAKE_MPI_ABORT_RETURNED" not in diagnostic


def test_insufficient_mpi_thread_support_aborts_then_hard_exits(
    mpi_init_insufficient_probe,
):
    result = subprocess.run(
        [str(mpi_init_insufficient_probe)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    expected_code = (
        SPONGE_SIMULATION_BREAKDOWN
        if sys.platform == "win32"
        else SPONGE_SIMULATION_BREAKDOWN & 0xFF
    )
    diagnostic = result.stdout + result.stderr
    assert result.returncode == expected_code
    assert "requires MPI_THREAD_MULTIPLE" in diagnostic
    assert "FAKE_MPI_ABORT_RETURNED" in diagnostic


def test_unregistered_callback_does_not_run_on_fatal(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "unregistered.marker"
    _run_fatal_mode(fatal_cleanup_probe, "unregister", marker)
    assert not marker.exists()


def test_fatal_callbacks_are_lifo_exactly_once_under_reentry(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "reentrant.marker"
    _run_fatal_mode(fatal_cleanup_probe, "reentrant_lifo", marker)
    assert marker.read_bytes() == b"BA"


def test_unregistering_a_middle_callback_preserves_lifo_order(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "ordered.marker"
    _run_fatal_mode(fatal_cleanup_probe, "preserve_order", marker)
    assert marker.read_bytes() == b"DCA"


def test_full_registry_is_reported_without_overwrite(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "capacity.marker"
    result = subprocess.run(
        [str(fatal_cleanup_probe), "capacity", str(marker)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert result.returncode != 0
    assert "spongeErrorValueErrorCommand" in result.stdout + result.stderr
    assert marker.read_bytes() == b"Y" + b"X" * 31


def test_concurrent_duplicate_registration_is_atomic(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "concurrent-registry.marker"
    _run_fatal_mode(fatal_cleanup_probe, "concurrent_registry", marker)
    assert marker.read_bytes() == b"C"


def test_only_one_concurrent_fatal_caller_owns_the_entire_sequence(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "concurrent-fatal.marker"
    result = _run_fatal_mode(fatal_cleanup_probe, "concurrent_fatal", marker)
    diagnostic = result.stdout + result.stderr
    assert diagnostic.count("spongeErrorValueErrorCommand") == 1
    assert "fatal sequence owner" in diagnostic
    assert "spongeErrorMissingCommand" not in diagnostic
    assert "secondary fatal caller" not in diagnostic
    assert marker.read_bytes() == b"T"


def test_callback_runs_outside_registry_lock_and_self_unregister_is_safe(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "callback-operations.marker"
    _run_fatal_mode(fatal_cleanup_probe, "callback_operations", marker)
    assert marker.read_bytes() == b"O"


def test_throwing_callback_does_not_skip_older_cleanup(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "callback-throw.marker"
    _run_fatal_mode(fatal_cleanup_probe, "callback_throw", marker)
    assert marker.read_bytes() == b"BA"


def test_register_unregister_and_fatal_race_has_lifetime_guarantee(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "fatal-registry-race.marker"
    _run_fatal_mode(fatal_cleanup_probe, "fatal_registry_race", marker)
    assert marker.read_bytes() == b"TG"


def test_fatal_termination_bypasses_atexit_handlers(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "hard-exit.marker"
    _run_fatal_mode(fatal_cleanup_probe, "hard_exit", marker)
    assert marker.read_bytes() == b"A"


@pytest.mark.parametrize(
    ("mode", "expected"),
    [
        ("direct_allocation_failure", b"D"),
        ("formatted_allocation_failure", b"M"),
    ],
)
def test_cleanup_and_base_diagnostic_survive_rejected_cpp_new(
    fatal_cleanup_probe, tmp_path, mode, expected
):
    marker = tmp_path / f"{mode}.marker"
    result = subprocess.run(
        [str(fatal_cleanup_probe), mode, str(marker)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    expected_code = (
        SPONGE_MALLOC_ERROR
        if sys.platform == "win32"
        else SPONGE_MALLOC_ERROR & 0xFF
    )
    assert result.returncode == expected_code
    diagnostic = result.stdout + result.stderr
    assert "spongeErrorMallocFailed" in diagnostic
    assert "fatal_cleanup_probe" in diagnostic
    assert marker.read_bytes() == expected


@pytest.mark.parametrize(
    ("mode", "expected_marker", "expected_detail", "expected_error"),
    [
        (
            "direct_snapshot_mutation",
            b"S",
            "ORIGINAL_EXTRA",
            "spongeErrorValueErrorCommand",
        ),
        (
            "formatted_snapshot_mutation",
            b"S",
            "Formatted source: ORIGINAL_FORMATTED",
            "spongeErrorValueErrorCommand",
        ),
        (
            "direct_snapshot_free",
            b"F",
            "ORIGINAL_EXTRA",
            "spongeErrorValueErrorCommand",
        ),
        (
            "formatted_snapshot_free",
            b"F",
            "Formatted source: ORIGINAL_FORMATTED",
            "spongeErrorValueErrorCommand",
        ),
        (
            "direct_malloc_snapshot_free",
            b"F",
            "ORIGINAL_EXTRA",
            "spongeErrorMallocFailed",
        ),
        (
            "formatted_malloc_snapshot_free",
            b"F",
            "Formatted source: ORIGINAL_FORMATTED",
            "spongeErrorMallocFailed",
        ),
    ],
)
def test_diagnostic_sources_are_snapshotted_before_cleanup(
    fatal_cleanup_probe,
    tmp_path,
    mode,
    expected_marker,
    expected_detail,
    expected_error,
):
    marker = tmp_path / f"{mode}.marker"
    result = _run_fatal_mode(fatal_cleanup_probe, mode, marker, expected_error)
    diagnostic = result.stdout + result.stderr
    assert "ORIGINAL_ERROR_BY" in diagnostic
    assert expected_detail in diagnostic
    assert "DESTROYED" not in diagnostic
    assert marker.read_bytes() == expected_marker


def test_nonmalloc_capture_failure_cleans_up_and_uses_fixed_snapshot(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "capture-allocation-failure.marker"
    result = _run_fatal_mode(
        fatal_cleanup_probe, "capture_allocation_failure", marker
    )
    diagnostic = result.stdout + result.stderr
    assert "spongeErrorValueErrorCommand" in diagnostic
    assert "UNAVAILABLE_ERROR_BY" in diagnostic
    assert "UNAVAILABLE_EXTRA" in diagnostic
    assert "DESTROYED" not in diagnostic
    assert marker.read_bytes() == b"C"


def test_fixed_snapshot_marks_truncation_when_heap_is_unavailable(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "capture-allocation-truncation.marker"
    result = _run_fatal_mode(
        fatal_cleanup_probe, "capture_allocation_truncation", marker
    )
    diagnostic = result.stdout + result.stderr
    assert "TRUNCATED_ERROR_BY" in diagnostic
    assert (
        "[diagnostic truncated: dynamic allocation unavailable]" in diagnostic
    )
    assert "FALLBACK_TAIL_SHOULD_BE_TRUNCATED" not in diagnostic
    assert marker.read_bytes() == b"T"


@pytest.mark.skipif(sys.platform == "win32", reason="requires POSIX pipes")
def test_fixed_snapshot_cleanup_precedes_blocking_stderr_write(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "blocked-stderr.marker"
    read_fd, write_fd = os.pipe()
    process = None
    try:
        os.set_blocking(write_fd, False)
        filler = b"P" * 65536
        filled = 0
        while True:
            try:
                filled += os.write(write_fd, filler)
            except BlockingIOError:
                break
        assert filled != 0
        os.set_blocking(write_fd, True)

        process = subprocess.Popen(
            [
                str(fatal_cleanup_probe),
                "capture_allocation_failure",
                str(marker),
            ],
            stdout=subprocess.DEVNULL,
            stderr=write_fd,
        )
        os.close(write_fd)
        write_fd = None

        deadline = time.monotonic() + 5
        while (
            not marker.exists()
            and process.poll() is None
            and time.monotonic() < deadline
        ):
            time.sleep(0.01)
        assert marker.read_bytes() == b"C"
        assert process.poll() is None

        reader = os.fdopen(read_fd, "rb")
        read_fd = None
        with reader:
            captured_stderr = reader.read()
        assert process.wait(timeout=30) != 0
        assert b"spongeErrorValueErrorCommand" in captured_stderr
    finally:
        if write_fd is not None:
            os.close(write_fd)
        if process is not None and process.poll() is None:
            process.kill()
            process.wait(timeout=30)
        if read_fd is not None:
            os.close(read_fd)


def test_formatted_fatal_diagnostic_is_not_fixed_buffer_truncated(
    fatal_cleanup_probe, tmp_path
):
    marker = tmp_path / "long-formatted.marker"
    result = _run_fatal_mode(fatal_cleanup_probe, "long_formatted", marker)
    assert "LONG_DIAGNOSTIC_END" in result.stdout + result.stderr
    assert marker.read_bytes() == b"L"
