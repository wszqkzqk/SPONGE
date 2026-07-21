import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE = Path(__file__).with_name("reaxff_input_probe.cpp")


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = shutil.which("c++") or shutil.which("clang++")
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the ReaxFF parser probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.extend(
        (
            REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include",
            REPOSITORY_ROOT
            / ".pixi"
            / "envs"
            / "dev-cpu-mpi"
            / "include",
        )
    )
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (
            candidate / "fftw3.h"
        ).is_file():
            return candidate
    pytest.skip("OpenMP and FFTW headers are required for the parser probe")


def test_reaxff_shared_input_parser_contract(tmp_path):
    executable = tmp_path / "reaxff_input_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
            "-ffast-math",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(PROBE),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile ReaxFF parser probe\n"
        f"stdout:\n{compile_result.stdout}\n"
        f"stderr:\n{compile_result.stderr}"
    )
    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, (
        "ReaxFF parser probe failed\n"
        f"stdout:\n{run_result.stdout}\n"
        f"stderr:\n{run_result.stderr}"
    )
