import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from Xponge.analysis import MdoutReader


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
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
EEQ_PAIR_PROBE = Path(__file__).with_name("reaxff_eeq_pair_probe.cpp")
TRANSACTION_PROBE = Path(__file__).with_name("reaxff_transaction_probe.cpp")


def _write_counted_lines(path, values):
    path.write_text(
        str(len(values)) + "\n" + "\n".join(map(str, values)) + "\n",
        encoding="utf-8",
    )


def _write_case(case_dir, coordinates, atom_types):
    case_dir.mkdir()
    shutil.copy2(
        REAXFF_STATIC_ROOT / "ffield.reax.cho",
        case_dir / "ffield.reax.cho",
    )
    mass_by_type = {"H": "1.008", "C": "12.011", "O": "15.999"}
    _write_counted_lines(
        case_dir / "mass.txt", [mass_by_type[t] for t in atom_types]
    )
    _write_counted_lines(case_dir / "type.txt", atom_types)

    coordinate_lines = [str(len(coordinates))]
    coordinate_lines.extend(
        " ".join(f"{component:.9g}" for component in coordinate)
        for coordinate in coordinates
    )
    coordinate_lines.extend(("25 25 25", "90 90 90"))
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            (
                "mode = 'NVE'",
                "dt = 0.0",
                "step_limit = 0",
                "cutoff = 4.0",
                "skin = 2.0",
                "coordinate_in_file = 'coordinate.txt'",
                "mass_in_file = 'mass.txt'",
                "frc = 'frc.dat'",
                "mdout = 'mdout.txt'",
                "print_zeroth_frame = true",
                "write_information_interval = 1",
                "write_mdout_interval = 1",
                "write_restart_file_interval = 0",
                "",
                "[REAXFF]",
                "in_file = 'ffield.reax.cho'",
                "type_in_file = 'type.txt'",
            )
        )
        + "\n",
        encoding="utf-8",
    )


def _run_case(case_dir):
    return subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


def _combined_output(result):
    return result.stdout + "\n" + result.stderr


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = (
        shutil.which("c++")
        or shutil.which("clang++")
        or shutil.which("g++")
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the ReaxFF EEQ probe")
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
    pytest.skip("OpenMP and FFTW headers are required for the ReaxFF EEQ probe")


def test_reaxff_eeq_tiny_nonzero_pair_force_probe(tmp_path):
    executable = tmp_path / "reaxff_eeq_pair_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
            "-ffast-math",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(EEQ_PAIR_PROBE),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile ReaxFF EEQ pair probe\n"
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
        "ReaxFF EEQ pair probe failed\n"
        f"stdout:\n{run_result.stdout}\n"
        f"stderr:\n{run_result.stderr}"
    )


def test_reaxff_evaluation_transaction_probe(tmp_path):
    executable = tmp_path / "reaxff_transaction_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
            "-ffast-math",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(TRANSACTION_PROBE),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile ReaxFF transaction probe\n"
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
        "ReaxFF transaction probe failed\n"
        f"stdout:\n{run_result.stdout}\n"
        f"stderr:\n{run_result.stderr}"
    )


def _replace_first_general_parameter(force_field, token):
    lines = force_field.read_text(encoding="utf-8").splitlines()
    assert int(lines[1].split()[0]) > 0
    comment = ""
    if "!" in lines[2]:
        comment = " !" + lines[2].split("!", 1)[1]
    lines[2] = token + comment
    force_field.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _disable_hydrogen_sigma_bond_order(force_field):
    lines = force_field.read_text(encoding="utf-8").splitlines()
    for line_index, line in enumerate(lines):
        fields = line.split()
        if len(fields) >= 9 and fields[0] == "H":
            fields[1] = "0.0000"
            lines[line_index] = "  " + "  ".join(fields)
            force_field.write_text(
                "\n".join(lines) + "\n", encoding="utf-8"
            )
            return
    raise AssertionError("failed to find the hydrogen atom parameter entry")


@pytest.mark.parametrize("bad_token", ["nan", "inf", "1e999", "1junk"])
def test_reaxff_rejects_nonfinite_or_partial_parameter_tokens(
    tmp_path, bad_token
):
    case_dir = tmp_path / f"bad_parameter_{re.sub('[^a-z0-9]', '_', bad_token)}"
    _write_case(
        case_dir,
        np.asarray(((12.0, 12.5, 12.5), (13.0, 12.5, 12.5))),
        ("H", "H"),
    )
    _replace_first_general_parameter(
        case_dir / "ffield.reax.cho", bad_token
    )
    result = _run_case(case_dir)
    output = _combined_output(result)
    assert result.returncode != 0, output
    assert "failed to parse general parameter at index 1" in output


def test_reaxff_active_exact_overlap_is_fatal(tmp_path):
    case_dir = tmp_path / "active_overlap"
    _write_case(
        case_dir,
        np.asarray(((12.5, 12.5, 12.5), (12.5, 12.5, 12.5))),
        ("H", "H"),
    )
    result = _run_case(case_dir)
    output = _combined_output(result)
    assert result.returncode != 0, output
    assert "overlap exactly" in output
    assert "radial force direction is undefined" in output


def test_reaxff_shielded_vdw_exact_overlap_has_zero_force_limit(tmp_path):
    case_dir = tmp_path / "shielded_vdw_overlap"
    _write_case(
        case_dir,
        np.asarray(((12.5, 12.5, 12.5), (12.5, 12.5, 12.5))),
        ("H", "H"),
    )
    # Disable only the H-H bond-order radius so bond order is exactly
    # inactive.  The active shielded vdW interaction then owns the overlap.
    _disable_hydrogen_sigma_bond_order(case_dir / "ffield.reax.cho")
    result = _run_case(case_dir)
    output = _combined_output(result)
    assert result.returncode == 0, output
    forces = np.fromfile(case_dir / "frc.dat", dtype=np.float32)
    assert forces.size >= 6
    assert np.all(np.isfinite(forces))
    assert np.max(np.abs(forces[-6:])) == 0.0
    mdout = MdoutReader(str(case_dir / "mdout.txt"))
    assert np.isfinite(mdout.REAXFF_VDW[0])
    assert abs(mdout.REAXFF_VDW[0]) > 1.0e-6


def test_reaxff_nonzero_sub_threshold_pair_is_evaluated(tmp_path):
    case_dir = tmp_path / "short_nonzero_pair"
    _write_case(
        case_dir,
        np.asarray(((12.5, 12.5, 12.5), (12.505, 12.5, 12.5))),
        ("H", "H"),
    )
    result = _run_case(case_dir)
    output = _combined_output(result)
    assert result.returncode == 0, output
    forces = np.fromfile(case_dir / "frc.dat", dtype=np.float32)
    assert forces.size >= 6
    assert np.all(np.isfinite(forces))
    mdout = MdoutReader(str(case_dir / "mdout.txt"))
    assert np.isfinite(mdout.REAXFF_BOND[0])
    assert abs(mdout.REAXFF_BOND[0]) > 1.0e-3


def test_reaxff_input_paths_preserve_spaces(tmp_path):
    case_dir = tmp_path / "case with spaces"
    _write_case(
        case_dir,
        np.asarray(((12.0, 12.5, 12.5), (13.0, 12.5, 12.5))),
        ("H", "H"),
    )
    input_dir = case_dir / "input files"
    input_dir.mkdir()
    (case_dir / "ffield.reax.cho").rename(input_dir / "ffield reax.cho")
    (case_dir / "type.txt").rename(input_dir / "atom types.txt")
    mdin = (case_dir / "mdin.spg.toml").read_text(encoding="utf-8")
    mdin = mdin.replace(
        "in_file = 'ffield.reax.cho'",
        "in_file = 'input files/ffield reax.cho'",
    ).replace(
        "type_in_file = 'type.txt'",
        "type_in_file = 'input files/atom types.txt'",
    )
    (case_dir / "mdin.spg.toml").write_text(mdin, encoding="utf-8")

    result = _run_case(case_dir)
    assert result.returncode == 0, _combined_output(result)


def test_reaxff_undercoordinated_bond_order_correction_is_stable(tmp_path):
    case_dir = tmp_path / "undercoordinated_carbon_pair"
    _write_case(
        case_dir,
        np.asarray(((11.0, 12.5, 12.5), (14.0, 12.5, 12.5))),
        ("C", "C"),
    )
    result = _run_case(case_dir)
    output = _combined_output(result)
    assert result.returncode == 0, output
    forces = np.fromfile(case_dir / "frc.dat", dtype=np.float32)
    assert forces.size >= 6
    assert np.all(np.isfinite(forces))
    mdout = MdoutReader(str(case_dir / "mdout.txt"))
    assert np.isfinite(mdout.REAXFF_BOND[0])


def test_reaxff_collinear_active_valence_angle_is_fatal(tmp_path):
    case_dir = tmp_path / "collinear_angle"
    _write_case(
        case_dir,
        np.asarray(
            (
                (12.5, 12.5, 12.5),
                (13.5, 12.5, 12.5),
                (11.5, 12.5, 12.5),
            )
        ),
        ("O", "H", "H"),
    )
    result = _run_case(case_dir)
    output = _combined_output(result)
    assert result.returncode != 0, output
    assert "undefined zero-arm/collinear ReaxFF valence-angle" in output


def _disable_all_valence_angle_entries(force_field):
    lines = force_field.read_text(encoding="utf-8").splitlines()
    count_line = next(
        i for i, line in enumerate(lines) if "Nr of angles" in line
    )
    count = int(lines[count_line].split()[0])
    for line_index in range(count_line + 1, count_line + count + 1):
        fields = lines[line_index].split()
        assert len(fields) >= 10
        fields[4] = "0.0000"
        lines[line_index] = "  " + "  ".join(fields)
    force_field.write_text("\n".join(lines) + "\n", encoding="utf-8")


def test_reaxff_collinear_active_torsion_is_fatal(tmp_path):
    case_dir = tmp_path / "collinear_torsion"
    _write_case(
        case_dir,
        np.asarray(
            (
                (10.0, 12.5, 12.5),
                (11.0, 12.5, 12.5),
                (12.4, 12.5, 12.5),
                (13.4, 12.5, 12.5),
            )
        ),
        ("H", "C", "C", "H"),
    )
    _disable_all_valence_angle_entries(case_dir / "ffield.reax.cho")
    result = _run_case(case_dir)
    output = _combined_output(result)
    assert result.returncode != 0, output
    assert "undefined zero-bond/collinear ReaxFF torsion" in output
