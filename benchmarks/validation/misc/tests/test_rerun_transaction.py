import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("rerun_transaction_probe.cpp")


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
            "a C++17 compiler is required for the rerun transaction probe"
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
    pytest.skip("OpenMP and FFTW headers are required for the rerun probe")


def test_semantically_invalid_rerun_box_does_not_publish_staged_frame(tmp_path):
    executable = tmp_path / "rerun_transaction_probe"
    transaction_marker = tmp_path / "transaction.marker"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
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
        "failed to compile rerun transaction probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )

    run_result = subprocess.run(
        [str(executable), str(transaction_marker)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    diagnostic = run_result.stdout + run_result.stderr
    assert run_result.returncode != 0
    assert "rerun box deformation matrix is not representable" in diagnostic
    assert transaction_marker.read_bytes() == b"P"
