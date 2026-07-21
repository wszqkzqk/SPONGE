import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("xponge_load_transaction_probe.cpp")


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
        pytest.skip("a C++17 compiler is required for the load transaction probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.append(REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include")
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (candidate / "fftw3.h").is_file():
            return candidate
    pytest.skip("OpenMP and FFTW headers are required for the load transaction probe")


def test_system_input_transaction_publishes_only_complete_state(tmp_path):
    executable = tmp_path / "xponge_load_transaction_probe"
    input_directory = tmp_path / "native_inputs"
    input_directory.mkdir()
    (input_directory / "mass_first.txt").write_text("2\n12\n13\n")
    (input_directory / "charge_first.txt").write_text("2\n0.1\n-0.1\n")
    (input_directory / "coordinate_first.txt").write_text(
        "2 1.5 3\n0 0 0\n1 0 0\n10 11 12\n90 90 90\n"
    )
    (input_directory / "mass_second.txt").write_text("3\n2\n3\n4\n")
    (input_directory / "charge_second.txt").write_text("3\n-1\n0\n1\n")
    (input_directory / "coordinate_second.txt").write_text(
        "3 8.5 9\n0 0 0\n1 0 0\n2 0 0\n"
        "21 22 23\n80 90 100\n"
    )
    (input_directory / "gromacs_first.top").write_text(
        """
[ defaults ]
1 2 yes 0.5 0.833333
[ atomtypes ]
A A 12.0 0.0 A 0.3 0.4184
B B 14.0 0.0 A 0.4 0.8368
[ moleculetype ]
MOL 1
[ atoms ]
1 A 1 MOL A1 1 -0.1 12.0
2 B 1 MOL B1 1 0.1 14.0
[ bonds ]
1 2 1 0.1 100.0
[ pairs ]
1 2 1
[ system ]
first transaction probe
[ molecules ]
MOL 1
""".strip()
        + "\n"
    )
    (input_directory / "gromacs_first.gro").write_text(
        "first t=1.25\n2\n"
        "    1MOL     A1    1   0.000   0.000   0.000\n"
        "    1MOL     B1    2   0.100   0.000   0.000\n"
        "   5.00000   5.00000   5.00000\n"
    )
    (input_directory / "gromacs_second.top").write_text(
        """
[ defaults ]
1 2 yes 1.0 1.0
[ atomtypes ]
C C 20.0 0.0 A 0.2 0.2
[ moleculetype ]
ONE 0
[ atoms ]
1 C 1 ONE C1 1 0.5 20.0
[ system ]
second transaction probe
[ molecules ]
ONE 1
""".strip()
        + "\n"
    )
    (input_directory / "gromacs_second.gro").write_text(
        "second t=4.25\n1\n"
        "    1ONE     C1    1   0.200   0.300   0.400\n"
        "   6.00000   7.00000   8.00000\n"
    )
    amber_reference = (
        REPOSITORY_ROOT
        / "benchmarks"
        / "comparison"
        / "reference"
        / "amber"
        / "statics"
        / "alanine_dipeptide_tip4pew"
    )
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(PROBE_SOURCE),
            str(REPOSITORY_ROOT / "SPONGE" / "xponge" / "xponge.cpp"),
            str(REPOSITORY_ROOT / "SPONGE" / "common.cpp"),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile load transaction probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )
    run_result = subprocess.run(
        [
            str(executable),
            str(input_directory),
            str(amber_reference / "system.parm7"),
            str(amber_reference / "system.rst7"),
            str(amber_reference / "system_minimized.rst7"),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr
