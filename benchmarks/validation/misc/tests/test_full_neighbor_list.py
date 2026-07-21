import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("full_neighbor_list_probe.cpp")


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = (
        shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the full-list probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.append(
        REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include"
    )
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (
            candidate / "fftw3.h"
        ).is_file():
            return candidate
    pytest.skip("OpenMP and FFTW headers are required for the full-list probe")


def _openmp_flags(openmp_library_dir):
    if sys.platform == "darwin":
        return [
            "-Xpreprocessor",
            "-fopenmp",
            f"-L{openmp_library_dir}",
            f"-Wl,-rpath,{openmp_library_dir}",
            "-lomp",
        ]
    return ["-fopenmp"]


def test_full_neighbor_list_dd_semantics_under_sanitizers(tmp_path):
    executable = tmp_path / "full_neighbor_list_probe"
    dependency_include = _dependency_include()
    openmp_library_dir = dependency_include.parent / "lib"
    if sys.platform == "darwin":
        # Keep the pixi directory itself out of the executable rpath: it also
        # contains a compiler-version-specific ASan runtime.  A private libomp
        # copy lets Apple clang load its matching sanitizer runtime.
        shutil.copy2(openmp_library_dir / "libomp.dylib", tmp_path)
        openmp_library_dir = tmp_path
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-O1",
            "-g",
            "-DUSE_CPU",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
            *_openmp_flags(openmp_library_dir),
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{dependency_include}",
            str(PROBE_SOURCE),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile full-neighbor-list probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )

    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        env={
            **os.environ,
            "ASAN_OPTIONS": "halt_on_error=1",
            "UBSAN_OPTIONS": "halt_on_error=1:print_stacktrace=1",
        },
        timeout=120,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr
    assert "full-neighbor-list probe passed" in run_result.stdout
