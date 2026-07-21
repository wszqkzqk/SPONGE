import json
import math
import os
import struct
import subprocess
import textwrap
from pathlib import Path

import pytest


def _write_counted_values(path, values):
    Path(path).write_text(
        str(len(values))
        + "\n"
        + "\n".join(f"{value:.12g}" for value in values)
        + "\n",
        encoding="utf-8",
    )


def _write_native_case(
    case_dir,
    bonds,
    coordinates,
    *,
    bond_text=None,
    box_length=40.0,
    pbc=True,
):
    case_dir.mkdir()
    atom_count = len(coordinates)
    _write_counted_values(case_dir / "mass.txt", [12.0] * atom_count)
    _write_counted_values(case_dir / "charge.txt", [0.0] * atom_count)
    coordinate_lines = [str(atom_count)]
    coordinate_lines.extend(
        " ".join(f"{value:.12g}" for value in coordinate)
        for coordinate in coordinates
    )
    if isinstance(box_length, (tuple, list)):
        box_lengths = box_length
    else:
        box_lengths = (box_length, box_length, box_length)
    coordinate_lines.extend(
        (" ".join(f"{value:.12g}" for value in box_lengths), "90 90 90")
    )
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )
    bond_lines = [str(len(bonds))]
    bond_lines.extend(
        f"{atom_i} {atom_j} {force_constant:.12g} {distance:.12g}"
        for atom_i, atom_j, force_constant, distance in bonds
    )
    (case_dir / "bond.txt").write_text(
        "\n".join(bond_lines) + "\n" if bond_text is None else bond_text,
        encoding="utf-8",
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 0,
        "dt": 0,
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "bond_in_file": "bond.txt",
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_pressure": True,
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
        "pbc": pbc,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )


def _run_case(case_dir, *, check=True):
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if check and result.returncode != 0:
        raise AssertionError(
            f"SPONGE failed with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    return result


def _read_mdout(case_dir):
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    return {
        name: float(value)
        for name, value in zip(lines[0].split(), lines[-1].split())
    }


def _read_forces(case_dir, atom_count):
    raw = (case_dir / "frc.dat").read_bytes()
    values = struct.unpack(f"={len(raw) // 4}f", raw)
    return values[-3 * atom_count :]


@pytest.mark.parametrize(
    ("force_constant", "equilibrium_distance"),
    [(0.0, 1.25), (2.0, 0.0)],
    ids=["zero_force_constant", "zero_equilibrium_distance"],
)
def test_exact_zero_harmonic_bond_has_zero_energy_force_and_virial(
    tmp_path, force_constant, equilibrium_distance
):
    case_dir = tmp_path / f"bond_zero_{force_constant}_{equilibrium_distance}"
    _write_native_case(
        case_dir,
        [(0, 1, force_constant, equilibrium_distance)],
        [(10.0, 10.0, 10.0), (10.0, 10.0, 10.0)],
    )

    _run_case(case_dir)
    mdout = _read_mdout(case_dir)
    forces = _read_forces(case_dir, 2)
    assert all(math.isfinite(value) for value in mdout.values())
    assert all(math.isfinite(value) for value in forces)
    assert mdout["bond"] == pytest.approx(0.0, abs=1.0e-7)
    assert max(abs(value) for value in forces) == pytest.approx(0.0, abs=1.0e-7)


def test_regular_harmonic_bond_energy_and_force_are_unchanged(tmp_path):
    case_dir = tmp_path / "bond_regular"
    _write_native_case(
        case_dir,
        [(0, 1, 2.0, 1.0)],
        [(10.0, 10.0, 10.0), (10.5, 10.0, 10.0)],
    )

    _run_case(case_dir)
    assert _read_mdout(case_dir)["bond"] == pytest.approx(0.5, abs=0.01)
    assert _read_forces(case_dir, 2) == pytest.approx(
        (-2.0, 0.0, 0.0, 2.0, 0.0, 0.0), abs=1.0e-6
    )


def test_nopbc_bond_uses_direct_displacement_with_an_ordinary_box(tmp_path):
    case_dir = tmp_path / "bond_nopbc_direct_displacement"
    _write_native_case(
        case_dir,
        [(0, 1, 2.0, 1.0)],
        [(5.0, 10.0, 10.0), (35.0, 10.0, 10.0)],
        box_length=40.0,
        pbc=False,
    )

    result = _run_case(case_dir)
    output = result.stdout + result.stderr
    # Direct distance is 30 A: E = k * (r-r0)^2 = 2 * 29^2.  A periodic
    # minimum image would incorrectly use 10 A and produce 162 instead.
    assert _read_mdout(case_dir)["bond"] == pytest.approx(1682.0, abs=0.01)
    assert "may be inaccurate" not in output


def test_nopbc_direct_displacement_does_not_form_fractional_images(tmp_path):
    case_dir = tmp_path / "bond_nopbc_fractional_overflow"
    _write_native_case(
        case_dir,
        [(0, 1, 1.0e-20, 0.0)],
        [(0.0, 0.0, 0.0), (1.0e10, 0.0, 0.0)],
        # Every cell component, reciprocal, pair product, and the volume are
        # representable normal floats.  Only the irrelevant NOPBC fractional
        # image calculation (dr.x * rcell.a11 = 1e40) would overflow.
        box_length=(1.0e-30, 1.0e15, 1.0e15),
        pbc=False,
    )

    result = _run_case(case_dir)
    output = result.stdout + result.stderr
    assert all(math.isfinite(value) for value in _read_mdout(case_dir).values())
    assert _read_mdout(case_dir)["bond"] == pytest.approx(1.0, abs=0.01)
    assert all(math.isfinite(value) for value in _read_forces(case_dir, 2))
    assert "non-finite geometry/energy/force/virial" not in output


def test_undefined_zero_distance_bond_fails_with_global_context(tmp_path):
    case_dir = tmp_path / "bond_undefined_zero"
    _write_native_case(
        case_dir,
        [(0, 1, 2.0, 1.0)],
        [(10.0, 10.0, 10.0), (10.0, 10.0, 10.0)],
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "bond term 0 (global atoms 0 1)" in output
    assert "undefined zero-distance geometry" in output


def test_bond_rejects_nonfinite_single_term_energy_with_global_context(
    tmp_path,
):
    case_dir = tmp_path / "bond_single_term_overflow"
    _write_native_case(
        case_dir,
        [(0, 1, 3.0e38, 0.0)],
        [(10.0, 10.0, 10.0), (12.0, 10.0, 10.0)],
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "bond term 0 (global atoms 0 1)" in output
    assert "non-finite geometry/energy/force/virial" in output


def test_bond_rejects_finite_force_accumulation_overflow(tmp_path):
    case_dir = tmp_path / "bond_accumulator_overflow"
    _write_native_case(
        case_dir,
        [(0, 1, 1.0e38, 0.0), (0, 1, 1.0e38, 0.0)],
        [(10.0, 10.0, 10.0), (11.0, 10.0, 10.0)],
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "bond term " in output
    assert "global atoms 0 1" in output
    assert "non-finite force/energy/virial accumulator" in output


@pytest.mark.parametrize(
    ("case_name", "bond_text", "message"),
    [
        ("negative_count", "-1\n", "negative interaction count"),
        (
            "count_over_int",
            "2147483648\n",
            "not a strict signed integer in range",
        ),
        (
            "max_count_truncated",
            "2147483647\n",
            "interaction 0 is missing atom A",
        ),
        ("count_suffix", "1junk\n", "not a strict signed integer in range"),
        (
            "atom_over_int",
            "1\n2147483648 1 2.0 1.0\n",
            "not a strict signed integer in range",
        ),
        (
            "atom_suffix",
            "1\n0junk 1 2.0 1.0\n",
            "not a strict signed integer in range",
        ),
        (
            "float_suffix",
            "1\n0 1 2.0junk 1.0\n",
            "invalid or outside the finite double range",
        ),
        (
            "nonfinite",
            "1\n0 1 nan 1.0\n",
            "non-finite parameter",
        ),
        (
            "double_overflow",
            "1\n0 1 1e309 1.0\n",
            "invalid or outside the finite double range",
        ),
        (
            "float_overflow",
            "1\n0 1 3.5e38 1.0\n",
            "outside the finite float range",
        ),
        (
            "float_underflow",
            "1\n0 1 1e-50 1.0\n",
            "underflows the finite float range",
        ),
        (
            "truncated",
            "1\n0 1 2.0\n",
            "interaction 0 is missing equilibrium distance",
        ),
        (
            "trailing",
            "1\n0 1 2.0 1.0\nunexpected\n",
            "trailing data beginning with 'unexpected'",
        ),
    ],
)
def test_native_bond_parser_is_strict_and_contextual(
    tmp_path, case_name, bond_text, message
):
    case_dir = tmp_path / case_name
    _write_native_case(
        case_dir,
        [],
        [(10.0, 10.0, 10.0), (11.0, 10.0, 10.0)],
        bond_text=bond_text,
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert message in output
    assert "Input file: bond.txt" in output


@pytest.mark.parametrize(
    ("case_name", "bond_text", "message"),
    [
        ("atom_outside_system", "1\n0 2 2.0 1.0\n", "outside [0, 2)"),
        ("self_bond", "1\n1 1 0.0 0.0\n", "repeats global atom 1"),
        (
            "negative_equilibrium_distance",
            "1\n0 1 2.0 -1.0\n",
            "negative equilibrium distance",
        ),
    ],
)
def test_bond_ir_validation_precedes_allocation_and_connectivity(
    tmp_path, case_name, bond_text, message
):
    case_dir = tmp_path / case_name
    _write_native_case(
        case_dir,
        [],
        [(10.0, 10.0, 10.0), (11.0, 10.0, 10.0)],
        bond_text=bond_text,
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "invalid bond data" in output
    assert message in output


def test_native_bond_rejects_smallest_positive_float_subnormal(tmp_path):
    case_dir = tmp_path / "smallest_float_subnormal"
    _write_native_case(
        case_dir,
        [],
        [(10.0, 10.0, 10.0), (10.0, 10.0, 10.0)],
        bond_text="1\n0 1 1.401298464324817e-45 0.0\n",
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "force constant parameter is a subnormal float" in output
    assert "consistent FTZ behavior" in output
    assert "Input file: bond.txt" in output


def _gro_atom(atom_name, atom_number, coordinate):
    x, y, z = coordinate
    return (
        f"{1:5d}{'MOL':<5}{atom_name:>5}{atom_number:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )


def test_gromacs_vsite1_retained_bond_fails_at_undefined_zero_distance(
    tmp_path,
):
    case_dir = tmp_path / "gromacs_vsite1_retained_bond"
    case_dir.mkdir()
    topology = """
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.0 0.0
        X X 0.0 0.0 V 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A 1 0.0 12.0
        2 X 1 MOL VS 1 0.0 0.0

        [ virtual_sites1 ]
        2 1

        [ bonds ]
        2 1 1 0.20 836.8

        [ exclusions ]
        1 2

        [ system ]
        retained vsite1 bond

        [ molecules ]
        MOL 1
    """
    (case_dir / "topol.top").write_text(
        textwrap.dedent(topology).strip() + "\n", encoding="utf-8"
    )
    gro_lines = [
        "retained vsite1 bond",
        "2",
        _gro_atom("A", 1, (1.0, 1.0, 1.0)),
        _gro_atom("VS", 2, (2.0, 1.0, 1.0)),
        "   4.00000   4.00000   4.00000",
    ]
    (case_dir / "conf.gro").write_text(
        "\n".join(gro_lines) + "\n", encoding="utf-8"
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 0,
        "dt": 0,
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "gromacs_top": "topol.top",
        "gromacs_gro": "conf.gro",
        "mdout": "mdout.txt",
        "print_pressure": True,
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_restart_file_interval": 0,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "bond term 0 (global atoms 1 0)" in output
    assert "undefined zero-distance geometry" in output
