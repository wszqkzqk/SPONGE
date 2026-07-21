import json
import math
import os
import struct
import subprocess
from pathlib import Path

import pytest

BOX = (40.0, 40.0, 40.0)
PROPER_K = 1.7
PROPER_PHASE = 0.61
PROPER_MULTIPLICITY = -3
IMPROPER_K = 1.3
IMPROPER_PHASE = 0.41


def _write_counted_values(path, values):
    Path(path).write_text(
        str(len(values))
        + "\n"
        + "\n".join(str(value) for value in values)
        + "\n",
        encoding="utf-8",
    )


def _write_case(
    case_dir,
    coordinates,
    kind,
    *,
    phase=None,
    torsion_text=None,
    print_pressure=False,
):
    case_dir.mkdir(parents=True, exist_ok=True)
    _write_counted_values(case_dir / "mass.txt", [12.0] * len(coordinates))
    _write_counted_values(case_dir / "charge.txt", [0.0] * len(coordinates))
    coordinate_lines = [str(len(coordinates))]
    coordinate_lines.extend(
        f"{x:.12g} {y:.12g} {z:.12g}" for x, y, z in coordinates
    )
    coordinate_lines.extend(
        (" ".join(str(value) for value in BOX), "90.0 90.0 90.0")
    )
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )

    if kind == "proper":
        filename = "dihedral.txt"
        command = "dihedral_in_file"
        if torsion_text is None:
            phase = PROPER_PHASE if phase is None else phase
            torsion_text = (
                f"1\n0 1 2 3 {PROPER_MULTIPLICITY} {PROPER_K} {phase:.12g}\n"
            )
    else:
        filename = "improper_dihedral.txt"
        command = "improper_dihedral_in_file"
        if torsion_text is None:
            phase = IMPROPER_PHASE if phase is None else phase
            torsion_text = f"1\n0 1 2 3 {IMPROPER_K} {phase:.12g}\n"
    (case_dir / filename).write_text(torsion_text, encoding="utf-8")

    mdin = {
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
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
    }
    if print_pressure:
        mdin["print_pressure"] = True
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(f"{key} = {json.dumps(value)}" for key, value in mdin.items())
        + "\n",
        encoding="utf-8",
    )


def _run_case(case_dir, *, check=True, mpi_np=None):
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]
    if mpi_np is not None:
        command = [
            "mpirun",
            "--oversubscribe",
            "-np",
            str(mpi_np),
            *command,
        ]
    result = subprocess.run(
        command,
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
        for name, value in zip(lines[0].split(), lines[1].split())
    }


def _read_forces(case_dir, atom_count=4):
    raw = (case_dir / "frc.dat").read_bytes()
    values = struct.unpack(f"={len(raw) // 4}f", raw)
    return values[-3 * atom_count :]


def _subtract(first, second):
    return tuple(a - b for a, b in zip(first, second))


def _cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def _dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def _minimum_image(displacement, box=BOX):
    return tuple(
        value - math.floor(value / length + 0.5) * length
        for value, length in zip(displacement, box)
    )


def _torsion_phi(coordinates):
    drij = _minimum_image(_subtract(coordinates[0], coordinates[1]))
    drkj = _minimum_image(_subtract(coordinates[2], coordinates[1]))
    drkl = _minimum_image(_subtract(coordinates[2], coordinates[3]))
    normal_1 = _cross(drij, drkj)
    normal_2 = _cross(drkl, drkj)
    normal_product = math.sqrt(
        _dot(normal_1, normal_1) * _dot(normal_2, normal_2)
    )
    central_length = math.sqrt(_dot(drkj, drkj))
    cosine = -_dot(normal_1, normal_2) / normal_product
    sine = _dot(_cross(normal_2, normal_1), drkj) / (
        normal_product * central_length
    )
    return math.atan2(sine, cosine)


def _energy(kind, coordinates, *, phase=None):
    phi = _torsion_phi(coordinates)
    if kind == "proper":
        phase = PROPER_PHASE if phase is None else phase
        return PROPER_K * (1.0 + math.cos(PROPER_MULTIPLICITY * phi - phase))
    phase = IMPROPER_PHASE if phase is None else phase
    delta = math.remainder(phi - phase, 2.0 * math.pi)
    return IMPROPER_K * delta * delta


GENERAL_COORDINATES = [
    (1.0, 2.2, 0.3),
    (2.1, 0.8, 1.9),
    (4.3, 3.1, 1.1),
    (5.8, 1.7, 4.2),
]
NEAR_PLANAR_COORDINATES = [
    (1.0, 2.0, 0.0),
    (1.0, 1.0, 0.0),
    (3.0, 1.0, 0.0),
    (3.0, 0.0, 0.00001),
]
# Both plane normals are nonzero, but their squared lengths are about 1e-44.
# They therefore cannot be classified or inverted in float even though the
# final roughly 1e12 Cartesian force components are finite floats.
SMALL_NONZERO_NORMAL_COORDINATES = [
    (0.0, 1.0e-12, 0.0),
    (0.0, 0.0, 0.0),
    (1.0e-10, 0.0, 0.0),
    (1.0e-10, 1.0e-12, 1.0e-12),
]
VIRIAL_OVERFLOW_COORDINATES = [
    (3.3858608141483106, 6.64555505564911, -1.3542055366957726),
    (3.572704134714767, -4.932186749113466, 7.857499102449379),
    (3.563509518537904, -5.426970300787671, -0.858875941343161),
    (7.096057488161737, -5.691484491058793, -4.901344798020586),
]


@pytest.mark.parametrize("kind", ["proper", "improper"])
@pytest.mark.parametrize(
    ("geometry", "coordinates", "epsilon"),
    [
        ("general", GENERAL_COORDINATES, 1.0e-5),
        ("near_planar", NEAR_PLANAR_COORDINATES, 1.0e-6),
        (
            "small_nonzero_normals",
            SMALL_NONZERO_NORMAL_COORDINATES,
            1.0e-14,
        ),
    ],
)
def test_torsion_energy_and_force_match_finite_difference(
    tmp_path, kind, geometry, coordinates, epsilon
):
    case_dir = tmp_path / f"{kind}_{geometry}"
    _write_case(case_dir, coordinates, kind)
    _run_case(case_dir)

    term_name = "dihedral" if kind == "proper" else "improper_dihedral"
    mdout = _read_mdout(case_dir)
    assert all(math.isfinite(value) for value in mdout.values())
    assert mdout[term_name] == pytest.approx(
        _energy(kind, coordinates), abs=0.02
    )
    forces = _read_forces(case_dir)
    assert all(math.isfinite(value) for value in forces)
    for atom in range(4):
        for axis in range(3):
            plus = [list(coordinate) for coordinate in coordinates]
            minus = [list(coordinate) for coordinate in coordinates]
            plus[atom][axis] += epsilon
            minus[atom][axis] -= epsilon
            numerical_force = -(_energy(kind, plus) - _energy(kind, minus)) / (
                2.0 * epsilon
            )
            assert forces[3 * atom + axis] == pytest.approx(
                numerical_force, rel=5.0e-3, abs=5.0e-3
            )


@pytest.mark.parametrize("kind", ["proper", "improper"])
@pytest.mark.parametrize(
    ("degeneracy", "coordinates"),
    [
        (
            "collinear",
            [
                (0.0, 0.0, 0.0),
                (1.0, 0.0, 0.0),
                (2.0, 0.0, 0.0),
                (3.0, 1.0, 0.0),
            ],
        ),
    ],
)
def test_degenerate_torsion_fails_with_term_and_atom_context(
    tmp_path, kind, degeneracy, coordinates
):
    case_dir = tmp_path / f"{kind}_{degeneracy}"
    _write_case(case_dir, coordinates, kind)
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"{kind} dihedral term 0" in output
    assert "global atoms 0 1 2 3" in output
    assert (
        "undefined, non-finite, or unrepresentable torsion "
        "geometry/energy/force/virial" in output
    )


@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_torsion_rejects_final_force_outside_float_range(tmp_path, kind):
    coordinates = [
        tuple(0.1 * value for value in atom)
        for atom in VIRIAL_OVERFLOW_COORDINATES
    ]
    phi = _torsion_phi(coordinates)
    if kind == "proper":
        row = f"0 1 2 3 1 1e38 {phi + math.pi / 2.0:.12g}"
    else:
        row = f"0 1 2 3 1e38 {phi - 1.0:.12g}"

    case_dir = tmp_path / f"{kind}_final_force_outside_float_range"
    _write_case(case_dir, coordinates, kind, torsion_text=f"1\n{row}\n")
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"{kind} dihedral term 0" in output
    assert "global atoms 0 1 2 3" in output
    assert "unrepresentable torsion geometry/energy/force/virial" in output
    force_path = case_dir / "frc.dat"
    if force_path.exists():
        assert all(math.isfinite(value) for value in _read_forces(case_dir))


@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_torsion_rejects_nonfinite_virial_before_accumulation(tmp_path, kind):
    phi = _torsion_phi(VIRIAL_OVERFLOW_COORDINATES)
    if kind == "proper":
        torsion_text = f"1\n0 1 2 3 1 2e38 {phi + math.pi / 2.0:.12g}\n"
    else:
        torsion_text = f"1\n0 1 2 3 1e38 {phi - 1.0:.12g}\n"

    force_only_dir = tmp_path / f"{kind}_large_finite_force"
    _write_case(
        force_only_dir,
        VIRIAL_OVERFLOW_COORDINATES,
        kind,
        torsion_text=torsion_text,
    )
    _run_case(force_only_dir)
    assert all(math.isfinite(value) for value in _read_forces(force_only_dir))

    pressure_dir = tmp_path / f"{kind}_virial_overflow"
    _write_case(
        pressure_dir,
        VIRIAL_OVERFLOW_COORDINATES,
        kind,
        torsion_text=torsion_text,
        print_pressure=True,
    )
    result = _run_case(pressure_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"{kind} dihedral term 0" in output
    assert (
        "undefined, non-finite, or unrepresentable torsion "
        "geometry/energy/force/virial" in output
    )


@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_torsion_atom_accumulator_overflow_fails_without_storing_inf(
    tmp_path, kind
):
    coordinates = [
        tuple(0.1 * value for value in atom)
        for atom in VIRIAL_OVERFLOW_COORDINATES
    ]
    phi = _torsion_phi(coordinates)
    if kind == "proper":
        row = f"0 1 2 3 1 5e37 {phi + math.pi / 2.0:.12g}"
    else:
        row = f"0 1 2 3 2.5e37 {phi - 1.0:.12g}"

    single_dir = tmp_path / f"{kind}_single_large_finite_term"
    _write_case(
        single_dir,
        coordinates,
        kind,
        torsion_text=f"1\n{row}\n",
    )
    _run_case(single_dir)
    single_forces = _read_forces(single_dir)
    assert all(math.isfinite(value) for value in single_forces)
    assert max(abs(value) for value in single_forces) > 2.0e38
    assert all(
        math.isfinite(value) for value in _read_mdout(single_dir).values()
    )

    duplicate_dir = tmp_path / f"{kind}_atom_accumulator_overflow"
    _write_case(
        duplicate_dir,
        coordinates,
        kind,
        torsion_text=f"2\n{row}\n{row}\n",
    )
    result = _run_case(duplicate_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"{kind} dihedral term" in output
    assert "global atoms 0 1 2 3" in output
    reason = output[output.index("Reason:") :]
    assert "accumulator" in reason
    force_path = duplicate_dir / "frc.dat"
    if force_path.exists():
        assert all(
            math.isfinite(value) for value in _read_forces(duplicate_dir)
        )


def test_torsion_global_potential_reduction_overflow_fails(tmp_path):
    coordinates = GENERAL_COORDINATES + [
        tuple(
            value + (10.0 if axis == 0 else 0.0)
            for axis, value in enumerate(atom)
        )
        for atom in GENERAL_COORDINATES
    ]
    torsion_text = "2\n0 1 2 3 0 1e38 0\n4 5 6 7 0 1e38 0\n"
    case_dir = tmp_path / "global_potential_reduction_overflow"
    _write_case(case_dir, coordinates, "proper", torsion_text=torsion_text)
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "global potential-energy reduction is not finite" in output
    assert "input contribution is non-finite" in output
    assert "finite reduction overflowed" in output
    assert all(math.isfinite(value) for value in _read_forces(case_dir, 8))


@pytest.mark.parametrize(
    ("kind", "force_constant"),
    [("proper", 9.0e37), ("improper", 4.5e37)],
)
def test_torsion_global_stress_reduction_overflow_fails(
    tmp_path, kind, force_constant
):
    first_coordinates = VIRIAL_OVERFLOW_COORDINATES
    coordinates = first_coordinates + [
        tuple(
            value + (12.0 if axis == 0 else 0.0)
            for axis, value in enumerate(atom)
        )
        for atom in first_coordinates
    ]
    phi = _torsion_phi(first_coordinates)
    phase = phi + math.pi / 2.0 if kind == "proper" else phi - 1.0

    def row(offset):
        atoms = [offset + atom for atom in range(4)]
        if kind == "proper":
            return (
                f"{atoms[0]} {atoms[1]} {atoms[2]} {atoms[3]} 1 "
                f"{force_constant:.12g} {phase:.12g}"
            )
        return (
            f"{atoms[0]} {atoms[1]} {atoms[2]} {atoms[3]} "
            f"{force_constant:.12g} {phase:.12g}"
        )

    single_dir = tmp_path / f"{kind}_single_large_finite_virial"
    _write_case(
        single_dir,
        first_coordinates,
        kind,
        torsion_text=f"1\n{row(0)}\n",
        print_pressure=True,
    )
    _run_case(single_dir)
    assert all(
        math.isfinite(value) for value in _read_mdout(single_dir).values()
    )

    aggregate_dir = tmp_path / f"{kind}_global_stress_reduction_overflow"
    _write_case(
        aggregate_dir,
        coordinates,
        kind,
        torsion_text=f"2\n{row(0)}\n{row(4)}\n",
        print_pressure=True,
    )
    result = _run_case(aggregate_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "global stress reduction is not finite" in output
    assert "input virial or kinetic contribution is non-finite" in output
    assert "finite reduction overflowed" in output


@pytest.mark.parametrize(
    ("speed", "message"),
    [
        (1.0e20, "global kinetic-energy reduction"),
        (1.0e18, "temperature cannot be represented"),
    ],
)
def test_kinetic_energy_and_temperature_overflow_fail_fast(
    tmp_path, speed, message
):
    case_dir = tmp_path / message.replace(" ", "_")
    _write_case(case_dir, GENERAL_COORDINATES, "proper")
    velocity_lines = [
        str(len(GENERAL_COORDINATES)),
        f"{speed:.12g} 0 0",
        f"{-speed:.12g} 0 0",
        "0 0 0",
        "0 0 0",
    ]
    (case_dir / "velocity.txt").write_text(
        "\n".join(velocity_lines) + "\n", encoding="utf-8"
    )
    with (case_dir / "mdin.spg.toml").open("a", encoding="utf-8") as handle:
        handle.write('velocity_in_file = "velocity.txt"\n')
    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert message in output


@pytest.mark.parametrize(
    ("masses", "source", "message"),
    [
        (
            [-1.0, 12.0, 12.0, 12.0],
            "MD_INFORMATION::Read_Mass",
            "non-finite or negative mass",
        ),
        (
            ["nan", 12.0, 12.0, 12.0],
            "mass_in_file mass entry 0",
            "strict finite decimal",
        ),
            (
                [1.0e-39, 12.0, 12.0, 12.0],
                "mass_in_file mass entry 0",
                "finite zero or normal float",
            ),
        (
            [1.0e38, 12.0, 12.0, 12.0],
            "MD_INFORMATION::Read_Mass",
            "reciprocal mass of atom 0",
        ),
        (
            [8.507059173023462e37] * 4,
            "MD_INFORMATION::Read_Mass",
            "system total mass",
        ),
    ],
)
def test_invalid_mass_is_rejected_at_the_common_runtime_boundary(
    tmp_path, masses, source, message
):
    case_dir = tmp_path / message.replace(" ", "_")
    _write_case(case_dir, GENERAL_COORDINATES, "proper")
    _write_counted_values(case_dir / "mass.txt", masses)
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert source in output
    assert message in output


def _continuous_coordinates(coordinates):
    drij = _minimum_image(_subtract(coordinates[0], coordinates[1]))
    drkj = _minimum_image(_subtract(coordinates[2], coordinates[1]))
    drkl = _minimum_image(_subtract(coordinates[2], coordinates[3]))
    atom_l = _subtract(drkj, drkl)
    return (drij, (0.0, 0.0, 0.0), drkj, atom_l)


def _reference_stress(coordinates, forces):
    tensor = [[0.0] * 3 for _ in range(3)]
    for atom, coordinate in enumerate(coordinates):
        for row in range(3):
            for column in range(3):
                tensor[row][column] += (
                    coordinate[row] * forces[3 * atom + column]
                )
    scale = 6.946827162543585e4 / math.prod(BOX)
    return {
        "Pxx": tensor[0][0] * scale,
        "Pyy": tensor[1][1] * scale,
        "Pzz": tensor[2][2] * scale,
        "Pxy": 0.5 * (tensor[0][1] + tensor[1][0]) * scale,
        "Pxz": 0.5 * (tensor[0][2] + tensor[2][0]) * scale,
        "Pyz": 0.5 * (tensor[1][2] + tensor[2][1]) * scale,
    }


@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_periodic_images_preserve_torsion_energy_force_and_stress(
    tmp_path, kind, mpi_np
):
    coordinates = [
        (18.8, 19.9, 18.7),
        (19.5, 20.4, 19.0),
        (20.3, 20.0, 19.6),
        (20.9, 20.7, 20.1),
    ]
    shifts = [(1, 0, -1), (0, -1, 0), (-1, 1, 0), (0, 0, 1)]
    imaged = [
        tuple(
            value + shift * length
            for value, shift, length in zip(atom, image, BOX)
        )
        for atom, image in zip(coordinates, shifts)
    ]

    launch_modes = [("direct", None)]
    if mpi_np is not None:
        launch_modes.append((f"mpi_{mpi_np}", mpi_np))
    all_results = {}
    for launch_name, launch_np in launch_modes:
        results = []
        for name, case_coordinates in (
            ("base", coordinates),
            ("imaged", imaged),
        ):
            case_dir = tmp_path / f"{kind}_{launch_name}_{name}"
            _write_case(case_dir, case_coordinates, kind, print_pressure=True)
            _run_case(case_dir, mpi_np=launch_np)
            mdout = _read_mdout(case_dir)
            forces = _read_forces(case_dir)
            expected_stress = _reference_stress(
                _continuous_coordinates(case_coordinates), forces
            )
            for component, expected in expected_stress.items():
                assert mdout[component] == pytest.approx(expected, abs=0.02)
            expected_pressure = (
                sum(expected_stress[key] for key in ("Pxx", "Pyy", "Pzz")) / 3.0
            )
            assert mdout["pressure"] == pytest.approx(
                expected_pressure, abs=0.02
            )
            term_name = "dihedral" if kind == "proper" else "improper_dihedral"
            assert mdout[term_name] == pytest.approx(
                _energy(kind, case_coordinates), abs=0.02
            )
            results.append((mdout, forces))

        term_name = "dihedral" if kind == "proper" else "improper_dihedral"
        assert results[1][0][term_name] == pytest.approx(
            results[0][0][term_name], abs=3.0e-5
        )
        assert results[1][1] == pytest.approx(
            results[0][1], rel=5.0e-5, abs=5.0e-5
        )
        for component in (
            "pressure",
            "Pxx",
            "Pyy",
            "Pzz",
            "Pxy",
            "Pxz",
            "Pyz",
        ):
            assert results[1][0][component] == pytest.approx(
                results[0][0][component], abs=0.02
            )
        all_results[launch_name] = results

    if mpi_np is not None:
        direct_results = all_results["direct"]
        mpi_results = all_results[f"mpi_{mpi_np}"]
        for direct, distributed in zip(direct_results, mpi_results):
            assert distributed[1] == pytest.approx(
                direct[1], rel=5.0e-5, abs=5.0e-5
            )
            for component in (
                term_name,
                "pressure",
                "Pxx",
                "Pyy",
                "Pzz",
                "Pxy",
                "Pxz",
                "Pyz",
            ):
                assert distributed[0][component] == pytest.approx(
                    direct[0][component], abs=0.02
                )


def test_improper_phase_is_periodic_across_multiple_turns(tmp_path):
    phases = (IMPROPER_PHASE, IMPROPER_PHASE + 12.0 * math.pi)
    results = []
    for index, phase in enumerate(phases):
        case_dir = tmp_path / f"phase_{index}"
        _write_case(case_dir, GENERAL_COORDINATES, "improper", phase=phase)
        _run_case(case_dir)
        results.append((_read_mdout(case_dir), _read_forces(case_dir)))
    assert results[1][0]["improper_dihedral"] == pytest.approx(
        results[0][0]["improper_dihedral"], abs=3.0e-5
    )
    assert results[1][1] == pytest.approx(results[0][1], rel=2.0e-5, abs=2.0e-5)


@pytest.mark.parametrize("kind", ["proper", "improper"])
@pytest.mark.parametrize(
    ("bad_data", "message"),
    [
        ("-1\n", "negative interaction count"),
        (None, "outside the system"),
        ("repeat", "repeats atom"),
        ("nonfinite", "non-finite parameter"),
        ("trailing", "trailing data"),
    ],
)
def test_native_torsion_input_validation(tmp_path, kind, bad_data, message):
    assert_source_context = bad_data in ("-1\n", "nonfinite", "trailing")
    if bad_data is None:
        row = "0 1 2 4 1 1.0 0.0" if kind == "proper" else "0 1 2 4 1.0 0.0"
        bad_data = f"1\n{row}\n"
    elif bad_data == "repeat":
        row = "0 1 1 3 1 1.0 0.0" if kind == "proper" else "0 1 1 3 1.0 0.0"
        bad_data = f"1\n{row}\n"
    elif bad_data == "nonfinite":
        row = "0 1 2 3 1 nan 0.0" if kind == "proper" else "0 1 2 3 nan 0.0"
        bad_data = f"1\n{row}\n"
    elif bad_data == "trailing":
        row = "0 1 2 3 1 1.0 0.0" if kind == "proper" else "0 1 2 3 1.0 0.0"
        bad_data = f"1\n{row}\nunexpected\n"

    case_dir = tmp_path / f"{kind}_{message.replace(' ', '_')}"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text=bad_data,
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert message in output
    if assert_source_context:
        filename = (
            "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
        )
        assert f"Input file: {filename}" in output


@pytest.mark.parametrize("kind", ["proper", "improper"])
@pytest.mark.parametrize(
    ("field", "token"),
    [
        ("interaction count", "1junk"),
        ("atom A", "0junk"),
        ("force constant", "1.0junk"),
        ("phase", "0.0junk"),
    ],
)
def test_native_torsion_requires_complete_numeric_tokens(
    tmp_path, kind, field, token
):
    if field == "interaction count":
        torsion_text = f"{token}\n"
    elif kind == "proper":
        values = ["0", "1", "2", "3", "1", "1.0", "0.0"]
        field_index = {"atom A": 0, "force constant": 5, "phase": 6}[field]
        values[field_index] = token
        torsion_text = f"1\n{' '.join(values)}\n"
    else:
        values = ["0", "1", "2", "3", "1.0", "0.0"]
        field_index = {"atom A": 0, "force constant": 4, "phase": 5}[field]
        values[field_index] = token
        torsion_text = f"1\n{' '.join(values)}\n"

    case_dir = tmp_path / f"{kind}_{field.replace(' ', '_')}_suffix"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text=torsion_text,
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"{field} token '{token}'" in output
    if field in ("interaction count", "atom A"):
        assert "strict signed integer in range" in output
    else:
        assert "invalid or outside the finite double range" in output
    filename = "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
    assert f"Input file: {filename}" in output


@pytest.mark.parametrize("kind", ["proper", "improper"])
@pytest.mark.parametrize("empty", [True, False])
def test_native_torsion_reports_truncated_input_with_field_context(
    tmp_path, kind, empty
):
    if empty:
        torsion_text = ""
        missing_field = "interaction count"
    elif kind == "proper":
        torsion_text = "1\n0 1 2 3 1 1.0\n"
        missing_field = "phase"
    else:
        torsion_text = "1\n0 1 2 3 1.0\n"
        missing_field = "phase"

    case_dir = tmp_path / f"{kind}_truncated_{missing_field.replace(' ', '_')}"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text=torsion_text,
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"missing {missing_field}" in output
    if not empty:
        assert "interaction 0" in output
    filename = "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
    assert f"Input file: {filename}" in output


def test_native_proper_rejects_inexact_float_multiplicity(tmp_path):
    case_dir = tmp_path / "inexact_multiplicity"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        "proper",
        torsion_text="1\n0 1 2 3 16777217 1.0 0.0\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert (
        "16777217 cannot be represented exactly by the force kernel" in output
    )


def test_native_proper_rejects_nonzero_token_below_double_range(tmp_path):
    case_dir = tmp_path / "nonzero_token_below_double_range"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        "proper",
        torsion_text="1\n0 1 2 3 1 1e-999 0.37\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "force constant token '1e-999'" in output
    assert "outside the finite double range" in output
    assert "Input file: dihedral.txt" in output


def test_native_proper_rejects_multiplicity_outside_signed_int(tmp_path):
    case_dir = tmp_path / "multiplicity_outside_signed_int"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        "proper",
        torsion_text="1\n0 1 2 3 2147483648 1.0 0.37\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "multiplicity token '2147483648'" in output
    assert "strict signed integer in range" in output
    assert "Input file: dihedral.txt" in output


@pytest.mark.parametrize(
    ("phase", "coefficient"),
    [("1e-8", "sine"), ("1.5707963267948966", "cosine")],
)
def test_native_proper_rejects_nonzero_derived_coefficient_underflow(
    tmp_path, phase, coefficient
):
    case_dir = tmp_path / f"derived_{coefficient}_coefficient_underflow"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        "proper",
        torsion_text=f"1\n0 1 2 3 1 1.17549435e-38 {phase}\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"nonzero {coefficient} coefficient underflows" in output
    assert "Input file: dihedral.txt" in output


@pytest.mark.parametrize(
    ("phase", "coefficient"), [("0.37", "cosine"), ("1e-4", "sine")]
)
def test_native_proper_rejects_subnormal_derived_coefficient(
    tmp_path, phase, coefficient
):
    case_dir = tmp_path / f"derived_{coefficient}_coefficient_subnormal"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        "proper",
        torsion_text=f"1\n0 1 2 3 1 1.17549435e-38 {phase}\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"{coefficient} coefficient is a subnormal float" in output
    assert "consistent FTZ behavior" in output
    assert "Input file: dihedral.txt" in output


@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_native_torsion_rejects_interaction_count_outside_signed_int(
    tmp_path, kind
):
    case_dir = tmp_path / f"{kind}_count_outside_signed_int"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text="2147483648\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "interaction count token '2147483648'" in output
    assert "strict signed integer in range" in output
    filename = "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
    assert f"Input file: {filename}" in output


@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_native_torsion_max_count_is_parsed_without_eager_allocation(
    tmp_path, kind
):
    case_dir = tmp_path / f"{kind}_max_count_truncated"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text="2147483647\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "interaction 0 is missing atom A" in output
    filename = "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
    assert f"Input file: {filename}" in output


@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_native_torsion_rejects_atom_index_outside_signed_int(tmp_path, kind):
    case_dir = tmp_path / f"{kind}_atom_outside_signed_int"
    if kind == "proper":
        row = "2147483648 1 2 3 1 1.0 0.0"
    else:
        row = "2147483648 1 2 3 1.0 0.0"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text=f"1\n{row}\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "atom A token '2147483648'" in output
    assert "strict signed integer in range" in output
    filename = "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
    assert f"Input file: {filename}" in output


@pytest.mark.parametrize(
    ("force_constant", "expected_error"),
    [
        ("1e-999", "outside the finite double range"),
        ("1e-46", "nonzero force constant parameter underflows"),
        ("1.40129846e-45", "force constant parameter is a subnormal float"),
        (
            "3.5e38",
            "force constant parameter is outside the finite float range",
        ),
    ],
)
@pytest.mark.parametrize("kind", ["proper", "improper"])
def test_native_torsion_rejects_unrepresentable_force_constant(
    tmp_path, force_constant, expected_error, kind
):
    case_dir = tmp_path / f"{kind}_unrepresentable_{force_constant}"
    if kind == "proper":
        torsion_text = f"1\n0 1 2 3 1 {force_constant} 0.0\n"
    else:
        torsion_text = f"1\n0 1 2 3 {force_constant} 0.0\n"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text=torsion_text,
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert expected_error in output
    filename = "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
    assert f"Input file: {filename}" in output


@pytest.mark.parametrize("kind", ["proper", "improper"])
@pytest.mark.parametrize(
    ("phase", "expected_error"),
    [
        ("1e-999", "outside the finite double range"),
        ("1e-46", "nonzero phase parameter underflows"),
        ("1.40129846e-45", "phase parameter is a subnormal float"),
        ("3.5e38", "phase parameter is outside the finite float range"),
        ("inf", "non-finite parameter"),
    ],
)
def test_native_torsion_rejects_unrepresentable_phase(
    tmp_path, kind, phase, expected_error
):
    case_dir = tmp_path / f"{kind}_unrepresentable_phase_{phase}"
    if kind == "proper":
        torsion_text = f"1\n0 1 2 3 1 1.0 {phase}\n"
    else:
        torsion_text = f"1\n0 1 2 3 1.0 {phase}\n"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        kind,
        torsion_text=torsion_text,
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert expected_error in output
    filename = "dihedral.txt" if kind == "proper" else "improper_dihedral.txt"
    assert f"Input file: {filename}" in output


def test_native_improper_rejects_negative_force_constant(tmp_path):
    case_dir = tmp_path / "negative_force_constant"
    _write_case(
        case_dir,
        GENERAL_COORDINATES,
        "improper",
        torsion_text="1\n0 1 2 3 -1.0 0.0\n",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "negative harmonic force constant" in output
