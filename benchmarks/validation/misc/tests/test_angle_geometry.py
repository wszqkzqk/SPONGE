import json
import math
import os
import struct
import subprocess
import textwrap
from pathlib import Path

import pytest


def _float32(value):
    return struct.unpack("=f", struct.pack("=f", value))[0]


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
    interactions,
    coordinates,
    *,
    interaction_text=None,
    kind="angle",
    print_pressure=False,
    box_length=40.0,
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
    coordinate_lines.extend(
        (f"{box_length:.12g} {box_length:.12g} {box_length:.12g}", "90 90 90")
    )
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )

    if kind == "angle":
        filename = "angle.txt"
        command = "angle_in_file"
        lines = [str(len(interactions))]
        lines.extend(
            f"{a} {b} {c} {force_constant:.12g} {theta0:.12g}"
            for a, b, c, force_constant, theta0 in interactions
        )
    else:
        filename = "urey_bradley.txt"
        command = "urey_bradley_in_file"
        lines = [str(len(interactions))]
        lines.extend(
            f"{a} {b} {c} {angle_k:.12g} {theta0:.12g} {bond_k:.12g} {r0:.12g}"
            for a, b, c, angle_k, theta0, bond_k, r0 in interactions
        )
    (case_dir / filename).write_text(
        "\n".join(lines) + "\n"
        if interaction_text is None
        else interaction_text,
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
        command: filename,
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_pressure": print_pressure,
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
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


def _read_forces(case_dir, atom_count=3):
    raw = (case_dir / "frc.dat").read_bytes()
    values = struct.unpack(f"={len(raw) // 4}f", raw)
    return values[-3 * atom_count :]


def _subtract(first, second):
    return tuple(a - b for a, b in zip(first, second))


def _dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def _cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def _angle_energy(coordinates, force_constant, theta0):
    u = _subtract(coordinates[0], coordinates[1])
    v = _subtract(coordinates[2], coordinates[1])
    normal = _cross(u, v)
    theta = math.atan2(math.sqrt(_dot(normal, normal)), _dot(u, v))
    return force_constant * (theta - theta0) ** 2


def _finite_difference_forces(coordinates, force_constant, theta0, epsilon):
    forces = []
    for atom in range(3):
        for axis in range(3):
            plus = [list(coordinate) for coordinate in coordinates]
            minus = [list(coordinate) for coordinate in coordinates]
            plus[atom][axis] += epsilon
            minus[atom][axis] -= epsilon
            forces.append(
                -(
                    _angle_energy(plus, force_constant, theta0)
                    - _angle_energy(minus, force_constant, theta0)
                )
                / (2.0 * epsilon)
            )
    return forces


def test_angle_energy_force_and_virial_are_finite_and_match_finite_difference(
    tmp_path,
):
    coordinates = [(1.1, 0.2, 0.3), (0.0, 0.0, 0.0), (-0.4, 1.3, 0.2)]
    force_constant = 2.3
    theta0 = 0.8
    case_dir = tmp_path / "ordinary_angle"
    _write_native_case(
        case_dir,
        [(0, 1, 2, force_constant, theta0)],
        coordinates,
        print_pressure=True,
    )

    _run_case(case_dir)
    mdout = _read_mdout(case_dir)
    forces = _read_forces(case_dir)
    expected_forces = _finite_difference_forces(
        coordinates, force_constant, theta0, 1.0e-5
    )
    assert all(math.isfinite(value) for value in mdout.values())
    assert all(math.isfinite(value) for value in forces)
    assert mdout["angle"] == pytest.approx(
        _angle_energy(coordinates, force_constant, theta0), abs=0.011
    )
    assert forces == pytest.approx(expected_forces, rel=3.0e-3, abs=3.0e-3)


@pytest.mark.parametrize(
    "coordinates",
    [
        [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 1.0e-5, 0.0)],
        [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (-1.0, 1.0e-5, 0.0)],
    ],
    ids=["near_zero", "near_pi"],
)
def test_near_collinear_angle_is_not_empirically_clamped(tmp_path, coordinates):
    theta0 = 0.3 if coordinates[2][0] > 0.0 else 2.8
    case_dir = tmp_path / (
        "near_zero" if coordinates[2][0] > 0.0 else "near_pi"
    )
    _write_native_case(case_dir, [(0, 1, 2, 1.0, theta0)], coordinates)

    _run_case(case_dir)
    forces = _read_forces(case_dir)
    expected = _finite_difference_forces(coordinates, 1.0, theta0, 1.0e-7)
    assert all(math.isfinite(value) for value in forces)
    assert forces == pytest.approx(expected, rel=5.0e-3, abs=5.0e-3)


def test_inactive_angle_skips_exact_zero_arm(tmp_path):
    case_dir = tmp_path / "inactive_zero_arm"
    _write_native_case(
        case_dir,
        [(0, 1, 2, 0.0, 1.0)],
        [(1.0, 2.0, 3.0), (1.0, 2.0, 3.0), (2.0, 2.0, 3.0)],
        print_pressure=True,
    )

    _run_case(case_dir)
    assert _read_mdout(case_dir)["angle"] == pytest.approx(0.0, abs=1.0e-7)
    assert max(abs(value) for value in _read_forces(case_dir)) == pytest.approx(
        0.0, abs=1.0e-7
    )


def test_active_zero_arm_fails_with_global_term_and_atoms(tmp_path):
    case_dir = tmp_path / "active_zero_arm"
    _write_native_case(
        case_dir,
        [(0, 1, 2, 1.0, 1.0)],
        [(1.0, 2.0, 3.0), (1.0, 2.0, 3.0), (2.0, 2.0, 3.0)],
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "angle term 0 (global atoms 0 1 2)" in output
    assert "undefined zero-arm/collinear geometry" in output


@pytest.mark.parametrize(
    ("coordinates", "theta0"),
    [
        ([(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (2.0, 0.0, 0.0)], 0.0),
        (
            [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (-2.0, 0.0, 0.0)],
            _float32(math.pi),
        ),
    ],
    ids=["theta_zero", "theta_pi"],
)
def test_exact_collinear_equilibrium_has_zero_energy_force_and_virial(
    tmp_path, coordinates, theta0
):
    case_dir = tmp_path / f"collinear_equilibrium_{theta0}"
    _write_native_case(
        case_dir,
        [(0, 1, 2, 2.0, theta0)],
        coordinates,
        print_pressure=True,
    )

    _run_case(case_dir)
    assert _read_mdout(case_dir)["angle"] == pytest.approx(0.0, abs=1.0e-7)
    assert max(abs(value) for value in _read_forces(case_dir)) == pytest.approx(
        0.0, abs=1.0e-7
    )


@pytest.mark.parametrize(
    ("coordinates", "theta0"),
    [
        ([(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (2.0, 0.0, 0.0)], 0.1),
        ([(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (-2.0, 0.0, 0.0)], 3.0),
    ],
    ids=["theta_zero", "theta_pi"],
)
def test_exact_collinear_off_equilibrium_is_rejected(
    tmp_path, coordinates, theta0
):
    case_dir = tmp_path / f"collinear_off_equilibrium_{theta0}"
    _write_native_case(case_dir, [(0, 1, 2, 2.0, theta0)], coordinates)

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "angle term 0 (global atoms 0 1 2)" in output
    assert "undefined zero-arm/collinear geometry" in output


@pytest.mark.parametrize(
    ("case_name", "angle_text", "message"),
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
            "1\n2147483648 1 2 1.0 1.0\n",
            "not a strict signed integer in range",
        ),
        (
            "float_suffix",
            "1\n0 1 2 1.0junk 1.0\n",
            "invalid or outside the finite double range",
        ),
        ("nonfinite", "1\n0 1 2 nan 1.0\n", "non-finite parameter"),
        (
            "double_overflow",
            "1\n0 1 2 1e309 1.0\n",
            "invalid or outside the finite double range",
        ),
        (
            "float_overflow",
            "1\n0 1 2 3.5e38 1.0\n",
            "outside the finite float range",
        ),
        (
            "float_underflow",
            "1\n0 1 2 1e-50 1.0\n",
            "underflows the finite float range",
        ),
        (
            "truncated",
            "1\n0 1 2 1.0\n",
            "interaction 0 is missing equilibrium angle",
        ),
        (
            "trailing",
            "1\n0 1 2 1.0 1.0\nunexpected\n",
            "trailing data beginning with 'unexpected'",
        ),
    ],
)
def test_native_angle_parser_is_strict_and_contextual(
    tmp_path, case_name, angle_text, message
):
    case_dir = tmp_path / case_name
    _write_native_case(
        case_dir,
        [],
        [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        interaction_text=angle_text,
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert message in output
    assert "Input file: angle.txt" in output


@pytest.mark.parametrize(
    ("case_name", "angle_text", "message"),
    [
        ("atom_outside", "1\n0 1 3 1.0 1.0\n", "outside [0, 3)"),
        ("theta_negative", "1\n0 1 2 1.0 -0.1\n", "outside [0, pi]"),
        ("theta_above_pi", "1\n0 1 2 1.0 3.2\n", "outside [0, pi]"),
        (
            "repeat_a_b",
            "1\n0 0 2 0.0 1.0\n",
            "requires three distinct atoms",
        ),
        (
            "repeat_a_c",
            "1\n0 1 0 0.0 1.0\n",
            "requires three distinct atoms",
        ),
        (
            "repeat_b_c",
            "1\n0 1 1 0.0 1.0\n",
            "requires three distinct atoms",
        ),
    ],
)
def test_angle_ir_validation_precedes_allocation(
    tmp_path, case_name, angle_text, message
):
    case_dir = tmp_path / case_name
    _write_native_case(
        case_dir,
        [],
        [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        interaction_text=angle_text,
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "invalid angle data" in output
    assert message in output


def test_native_angle_rejects_smallest_positive_float_subnormal(tmp_path):
    case_dir = tmp_path / "smallest_angle_float_subnormal"
    _write_native_case(
        case_dir,
        [],
        [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (2.0, 0.0, 0.0)],
        interaction_text="1\n0 1 2 1.401298464324817e-45 0.0\n",
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "force constant parameter is a subnormal float" in output
    assert "consistent FTZ behavior" in output
    assert "Input file: angle.txt" in output


def test_angle_rejects_finite_energy_accumulation_overflow(tmp_path):
    case_dir = tmp_path / "angle_accumulator_overflow"
    coordinates = [
        (10.0, 0.0, 0.0),
        (0.0, 0.0, 0.0),
        (0.0, 10.0, 0.0),
    ]
    _write_native_case(
        case_dir,
        [(0, 1, 2, 1.0e38, 0.0), (0, 1, 2, 1.0e38, 0.0)],
        coordinates,
        box_length=40.0,
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "angle term " in output
    assert "global atoms 0 1 2" in output
    assert "non-finite force/energy/virial accumulator" in output


def test_urey_bradley_inner_angle_uses_contextual_runtime_failure(tmp_path):
    case_dir = tmp_path / "urey_zero_arm"
    _write_native_case(
        case_dir,
        [(0, 1, 2, 1.0, 1.0, 0.0, 0.0)],
        [(1.0, 2.0, 3.0), (1.0, 2.0, 3.0), (2.0, 2.0, 3.0)],
        kind="urey_bradley",
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "urey_bradley_angle angle term 0 (global atoms 0 1 2)" in output
    assert "undefined zero-arm/collinear geometry" in output


def test_urey_bradley_reused_angle_path_matches_finite_difference(tmp_path):
    coordinates = [(1.1, 0.2, 0.3), (0.0, 0.0, 0.0), (-0.4, 1.3, 0.2)]
    force_constant = 2.3
    theta0 = 0.8
    case_dir = tmp_path / "urey_ordinary_angle"
    _write_native_case(
        case_dir,
        [(0, 1, 2, force_constant, theta0, 0.0, 0.0)],
        coordinates,
        kind="urey_bradley",
        print_pressure=True,
    )

    _run_case(case_dir)
    mdout = _read_mdout(case_dir)
    forces = _read_forces(case_dir)
    expected_forces = _finite_difference_forces(
        coordinates, force_constant, theta0, 1.0e-5
    )
    assert all(math.isfinite(value) for value in mdout.values())
    assert mdout["urey_bradley"] == pytest.approx(
        _angle_energy(coordinates, force_constant, theta0), abs=0.011
    )
    assert forces == pytest.approx(expected_forces, rel=3.0e-3, abs=3.0e-3)


@pytest.mark.parametrize(
    ("case_name", "urey_text", "message", "validation"),
    [
        ("negative_count", "-1\n", "negative interaction count", False),
        (
            "max_count_truncated",
            "2147483647\n",
            "interaction 0 is missing atom A",
            False,
        ),
        (
            "truncated",
            "1\n0 1 2 1.0 1.0 1.0\n",
            "interaction 0 is missing equilibrium distance",
            False,
        ),
        (
            "trailing",
            "1\n0 1 2 1.0 1.0 1.0 1.0\nextra\n",
            "trailing data beginning with 'extra'",
            False,
        ),
        (
            "underflow",
            "1\n0 1 2 1e-50 1.0 1.0 1.0\n",
            "underflows the finite float range",
            False,
        ),
        (
            "subnormal",
            "1\n0 1 2 1.40129846e-45 1.0 1.0 1.0\n",
            "angle force constant parameter is a subnormal float",
            False,
        ),
        (
            "atom_outside",
            "1\n0 1 3 1.0 1.0 1.0 1.0\n",
            "outside [0, 3)",
            True,
        ),
        (
            "repeat_atom",
            "1\n0 1 0 0.0 1.0 0.0 0.0\n",
            "requires three distinct atoms",
            True,
        ),
        (
            "theta_outside",
            "1\n0 1 2 1.0 3.2 1.0 1.0\n",
            "outside [0, pi]",
            True,
        ),
        (
            "negative_r0",
            "1\n0 1 2 1.0 1.0 1.0 -1.0\n",
            "negative equilibrium distance",
            True,
        ),
    ],
)
def test_native_urey_bradley_parser_and_ir_validation_are_strict(
    tmp_path, case_name, urey_text, message, validation
):
    case_dir = tmp_path / f"urey_{case_name}"
    _write_native_case(
        case_dir,
        [],
        [(1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 1.0, 0.0)],
        interaction_text=urey_text,
        kind="urey_bradley",
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert message in output
    if validation:
        assert "invalid Urey-Bradley data" in output
    else:
        assert "Input file: urey_bradley.txt" in output


def _gro_atom(atom_name, atom_number, coordinate):
    x, y, z = coordinate
    return (
        f"{1:5d}{'MOL':<5}{atom_name:>5}{atom_number:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )


def test_gromacs_vsite1_retained_angle_fails_for_runtime_zero_arm(tmp_path):
    case_dir = tmp_path / "gromacs_vsite1_retained_angle"
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
        1 A 1 MOL SRC 1 0.0 12.0
        2 X 1 MOL VS 1 0.0 0.0
        3 A 1 MOL END 1 0.0 12.0

        [ virtual_sites1 ]
        2 1

        [ angles ]
        2 1 3 1 60.0 8.368

        [ exclusions ]
        1 2 3
        2 3

        [ system ]
        retained vsite1 angle

        [ molecules ]
        MOL 1
    """
    (case_dir / "topol.top").write_text(
        textwrap.dedent(topology).strip() + "\n", encoding="utf-8"
    )
    (case_dir / "conf.gro").write_text(
        "\n".join(
            [
                "retained vsite1 angle",
                "3",
                _gro_atom("SRC", 1, (1.0, 1.0, 1.0)),
                _gro_atom("VS", 2, (2.0, 1.0, 1.0)),
                _gro_atom("END", 3, (1.0, 2.0, 1.0)),
                "   4.00000   4.00000   4.00000",
            ]
        )
        + "\n",
        encoding="utf-8",
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
    assert "angle term 0 (global atoms 1 0 2)" in output
    assert "undefined zero-arm/collinear geometry" in output
