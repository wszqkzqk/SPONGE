import json
import os
import shlex
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("voronoi_detector_probe.cpp")


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
            "a C++17 compiler is required for the Voronoi detector probe"
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
    pytest.skip("OpenMP and FFTW headers are required for the detector probe")


def test_source_interface_recrossing_and_strict_detector_contract(tmp_path):
    executable = tmp_path / "voronoi_detector_probe"
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
            "-O3",
            "-march=native",
            "-ffast-math",
            "-ffunction-sections",
            "-fdata-sections",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(PROBE_SOURCE),
            str(REPOSITORY_ROOT / "SPONGE" / "common.cpp"),
            *dead_code_flags,
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=180,
    )
    assert compile_result.returncode == 0, (
        "failed to compile Voronoi detector probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )

    run_result = subprocess.run(
        [str(executable), str(tmp_path)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr


def test_full_engine_terminal_hit_exports_committed_state_once(tmp_path):
    case_dir = tmp_path / "full_engine_terminal_hit"
    case_dir.mkdir()

    (case_dir / "mass.txt").write_text("2\n12\n12\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("2\n0\n0\n", encoding="utf-8")
    (case_dir / "coordinate.txt").write_text(
        "2 11.25\n5 5 5\n6.4 5 5\n1000 1000 1000\n90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "velocity.txt").write_text(
        "2\n0 0 0\n1.2 0 0\n", encoding="utf-8"
    )
    (case_dir / "cv.toml").write_text(
        '[distance]\nCV_type = "distance"\natom = [0, 1]\n\n'
        '[voronoi_detector]\nCV = ["distance"]\n'
        'milestone_file = "milestones.txt"\n'
        'source_interface = "S_0_1"\n',
        encoding="utf-8",
    )
    (case_dir / "milestones.txt").write_text(
        "3\nM_0 1\nM_1 2\nM_2 3\n2\nS_0_1 0 1\nS_1_2 1 2\n",
        encoding="utf-8",
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 9,
        "dt": 0.01,
        "pbc": False,
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "velocity_in_file": "velocity.txt",
        "cv_in_file": "cv.toml",
        "mdout": "mdout.txt",
        "crd": "trajectory.dat",
        "rst": "periodic",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 2,
        "max_restart_export_count": 4,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )

    run_result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    output = run_result.stdout + "\n" + run_result.stderr
    assert run_result.returncode == 0, output
    event_log = (case_dir / "mdinfo.txt").read_text(encoding="utf-8")
    assert (
        "VORONOI_HIT source=S_0_1 from=1 destination=S_1_2 to=2 "
        "completed_steps=5 hit_time_ps=0.05 source_recrossings=1 "
        "artifact=voronoi_hit_S_1_2"
    ) in event_log
    assert event_log.count("VORONOI_HIT ") == 1

    hit_coordinate = case_dir / "voronoi_hit_S_1_2_coordinate.txt"
    hit_velocity = case_dir / "voronoi_hit_S_1_2_velocity.txt"
    assert hit_coordinate.is_file()
    assert hit_velocity.is_file()
    assert sorted(path.name for path in case_dir.glob("*voronoi_hit*")) == [
        hit_coordinate.name,
        hit_velocity.name,
    ]

    coordinate_lines = hit_coordinate.read_text(encoding="utf-8").splitlines()
    atom_count, hit_time, hit_step = coordinate_lines[0].split()
    assert int(atom_count) == 2
    assert int(hit_step) == 5
    assert float(hit_time) == pytest.approx(11.30, abs=1.0e-12)

    velocity_header = hit_velocity.read_text(encoding="utf-8").splitlines()[0]
    velocity_atom_count, velocity_time, velocity_step = velocity_header.split()
    assert int(velocity_atom_count) == 2
    assert int(velocity_step) == int(hit_step)
    assert float(velocity_time) == pytest.approx(float(hit_time), abs=1.0e-12)

    hit_coordinates = tuple(
        float(value) for line in coordinate_lines[1:3] for value in line.split()
    )
    expected_x = 6.4 + 5 * (0.01 * 20.455) * 1.2
    assert hit_coordinates == pytest.approx(
        (5.0, 5.0, 5.0, expected_x, 5.0, 5.0), abs=2.0e-6
    )

    # Steps 1 and 3 have already emitted rotating periodic restarts.  Step 5
    # is itself the next ordinary restart boundary, so neither a third restart
    # nor a sixth trajectory frame may appear after the committed x_5 hit.
    assert (case_dir / "periodic_coordinate.txt").is_file()
    assert (case_dir / "1_periodic_coordinate.txt").is_file()
    assert not (case_dir / "2_periodic_coordinate.txt").exists()
    assert not (case_dir / "2_periodic_velocity.txt").exists()

    trajectory = (case_dir / "trajectory.dat").read_bytes()
    assert len(trajectory) == 5 * 2 * 3 * 4
    final_frame = struct.unpack("=6f", trajectory[-6 * 4 :])
    assert final_frame == pytest.approx(hit_coordinates, abs=2.0e-6)
