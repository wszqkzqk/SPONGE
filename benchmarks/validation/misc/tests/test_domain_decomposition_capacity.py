import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name(
    "domain_decomposition_capacity_probe.cpp"
)
DD_SOURCE = (
    REPOSITORY_ROOT
    / "SPONGE"
    / "Domain_decomposition"
    / "Domain_decomposition.cpp"
)
FAKE_MPI_INCLUDE = Path(__file__).with_name("fake_mpi")


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
        pytest.skip(
            "a C++17 compiler is required for the domain-capacity probe"
        )
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
    pytest.skip(
        "OpenMP and FFTW headers are required for the domain-capacity probe"
    )


def test_domain_corner_storage_has_no_compile_time_rank_limit(tmp_path):
    executable = tmp_path / "domain_decomposition_capacity_probe"
    dead_code_flags = (
        ["-Wl,-dead_strip"]
        if sys.platform == "darwin"
        else ["-Wl,--gc-sections"]
    )
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
            "-DUSE_MPI",
            "-DNO_GLOBAL_CONTROLLER",
            "-ffunction-sections",
            "-fdata-sections",
            "-w",
            f"-I{FAKE_MPI_INCLUDE}",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(PROBE_SOURCE),
            str(DD_SOURCE),
            *dead_code_flags,
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile domain-capacity probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )

    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr
