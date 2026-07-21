import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("ri_workspace_probe.cpp")


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = shutil.which("c++") or shutil.which("clang++") or shutil.which(
        "g++"
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the RI workspace probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.append(REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include")
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (candidate / "fftw3.h").is_file():
            return candidate
    pytest.skip("OpenMP and FFTW headers are required for the RI workspace probe")


def test_ri_dynamic_workspace_covers_high_angular_momentum_under_sanitizers(
    tmp_path,
):
    executable = tmp_path / "ri_workspace_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-O1",
            "-g",
            "-DUSE_CPU",
            "-fsanitize=address,undefined",
            "-fno-omit-frame-pointer",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
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
        f"failed to compile RI workspace probe\nstdout:\n"
        f"{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
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
    assert "energy_workers=" in run_result.stdout
