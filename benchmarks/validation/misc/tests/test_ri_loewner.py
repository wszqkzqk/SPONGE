import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("ri_loewner_probe.cpp")


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
        pytest.skip("a C++17 compiler is required for the RI Loewner probe")
    return [compiler]


def test_ri_inverse_sqrt_response_handles_degenerate_and_truncated_spectra(
    tmp_path,
):
    executable = tmp_path / "ri_loewner_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-O2",
            "-ffast-math",
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
        f"failed to compile RI Loewner probe\nstdout:\n{compile_result.stdout}"
        f"\nstderr:\n{compile_result.stderr}"
    )
    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr
