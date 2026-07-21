import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("reaxff_atom_identity_probe.cpp")
REAXFF_STATIC_ROOT = (
    REPOSITORY_ROOT
    / "benchmarks"
    / "comparison"
    / "tests"
    / "lammps"
    / "statics"
    / "reaxff"
    / "sponge"
)


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
            "a C++17 compiler is required for the ReaxFF identity probe"
        )
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.append(
        REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include"
    )
    candidates.append(
        REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu-mpi" / "include"
    )
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (
            candidate / "fftw3.h"
        ).is_file():
            return candidate
    pytest.skip("OpenMP and FFTW headers are required for the identity probe")


def _write_counted_lines(path, values):
    path.write_text(
        str(len(values)) + "\n" + "\n".join(values) + "\n",
        encoding="utf-8",
    )


def _write_two_atom_reaxff_case(case_dir):
    case_dir.mkdir()
    shutil.copy2(REAXFF_STATIC_ROOT / "ffield.reax.cho", case_dir)
    _write_counted_lines(case_dir / "mass.txt", ["1.008", "1.008"])
    _write_counted_lines(case_dir / "type.txt", ["H", "H"])
    (case_dir / "coordinate.txt").write_text(
        "2\n5 12 12\n17 12 12\n24 24 24\n90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            (
                "mode = 'NVE'",
                "dt = 0.0",
                "step_limit = 0",
                "cutoff = 4.0",
                "skin = 1.0",
                "coordinate_in_file = 'coordinate.txt'",
                "mass_in_file = 'mass.txt'",
                "write_restart_file_interval = 0",
                "PM.MPI_size = 0",
                "DOM_DEC.split_nx = 2",
                "DOM_DEC.split_ny = 1",
                "DOM_DEC.split_nz = 1",
                "",
                "[REAXFF]",
                "in_file = 'ffield.reax.cho'",
                "type_in_file = 'type.txt'",
            )
        )
        + "\n",
        encoding="utf-8",
    )


def test_reaxff_global_id_gather_scatter_and_history_probe(tmp_path):
    executable = tmp_path / "reaxff_atom_identity_probe"
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
        "failed to compile ReaxFF identity probe\n"
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

    multi_pp_result = subprocess.run(
        [str(executable), "--multi-pp"],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    multi_pp_output = multi_pp_result.stdout + multi_pp_result.stderr
    assert multi_pp_result.returncode != 0, multi_pp_output
    assert "REAXFF requires one PP rank" in multi_pp_output
    assert "EEQ is a globally coupled linear system" in multi_pp_output
    assert "Refusing this configuration" in multi_pp_output


def test_reaxff_dd_multi_pp_fails_before_rank_local_eeq(tmp_path, mpi_np):
    if mpi_np != 2:
        pytest.skip("the ReaxFF DD fail-fast contract requires --mpi 2")

    case_dir = tmp_path / "reaxff_multi_pp"
    _write_two_atom_reaxff_case(case_dir)
    command = [
        "mpirun",
        "--oversubscribe",
        "-np",
        str(mpi_np),
        os.environ.get("SPONGE_BIN", "SPONGE"),
        "-mdin",
        "mdin.spg.toml",
    ]
    result = subprocess.run(
        command,
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "REAXFF requires one PP rank" in output
    assert "EEQ is a globally coupled linear system" in output
    assert "Refusing this configuration" in output
