import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("ecp_error_policy_probe.cpp")


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
        pytest.skip("a C++17 compiler is required for the ECP policy probe")
    return [compiler]


def test_ecp_error_policy_cast_cancellation_and_transaction(tmp_path):
    executable = tmp_path / "ecp_error_policy_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-O3",
            "-ffast-math",
            "-Wall",
            "-Wextra",
            "-Werror",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
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
        "failed to compile ECP error-policy probe\n"
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
    assert "storage_error=" in run_result.stdout
    assert "cancelled_observable=0" in run_result.stdout
    assert "amplified_error=" in run_result.stdout
    assert "fallback_threshold_bracket=[1e6,1.01e6]" in run_result.stdout
