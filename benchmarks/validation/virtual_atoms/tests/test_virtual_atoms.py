import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest


def _write_values(path, values):
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
    virtual_records,
    bonds,
    *,
    step_limit=1,
    dom_update_interval=None,
    cv_text=None,
):
    case_dir.mkdir(parents=True, exist_ok=True)
    atom_numbers = len(coordinates)
    _write_values(
        case_dir / "mass.txt",
        [
            0.0 if index in {record[1] for record in virtual_records} else 12.0
            for index in range(atom_numbers)
        ],
    )
    _write_values(case_dir / "charge.txt", [0.0] * atom_numbers)

    coordinate_lines = [str(atom_numbers)]
    coordinate_lines.extend(
        f"{x:.9f} {y:.9f} {z:.9f}" for x, y, z in coordinates
    )
    coordinate_lines.extend(("40.0 40.0 40.0", "90.0 90.0 90.0"))
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )

    (case_dir / "virtual.txt").write_text(
        "\n".join(
            " ".join(str(value) for value in record)
            for record in virtual_records
        )
        + "\n",
        encoding="utf-8",
    )
    bond_lines = [str(len(bonds))]
    bond_lines.extend(" ".join(str(value) for value in bond) for bond in bonds)
    (case_dir / "bond.txt").write_text(
        "\n".join(bond_lines) + "\n", encoding="utf-8"
    )

    mdin = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": step_limit,
        "dt": 0,
        "cutoff": 8.0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "bond_in_file": "bond.txt",
        "virtual_atom_in_file": "virtual.txt",
        "mdout": "mdout.txt",
        "crd": "crd.dat",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 1,
    }
    if cv_text is not None:
        (case_dir / "cv.txt").write_text(cv_text, encoding="utf-8")
        mdin["cv_in_file"] = "cv.txt"
    mdin_lines = [f"{key} = {json.dumps(value)}" for key, value in mdin.items()]
    if dom_update_interval is not None:
        mdin_lines.extend(
            ("", "[DOM_DEC]", f"update_interval = {dom_update_interval}")
        )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(mdin_lines) + "\n", encoding="utf-8"
    )


def _run_case(case_dir, *, omp_threads=4, check=True, mpi_np=None):
    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(omp_threads)
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]
    if mpi_np is not None:
        command = ["mpirun", "--oversubscribe", "-np", str(mpi_np), *command]
    result = subprocess.run(
        command,
        cwd=case_dir,
        env=env,
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


def test_type3_nonzero_distance_rejects_zero_direction(tmp_path):
    case_dir = tmp_path / "type3_zero_direction"
    _write_case(
        case_dir,
        [(0.1, 0.2, 0.3), (0.1, 0.2, 0.3), (0.4, 0.5, 0.6)],
        [(3, 2, 0, 1, 1, 0.2, 0.0)],
        [],
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "zero construction direction" in output


def test_type3_repeated_second_and_third_source_with_zero_k_is_valid(tmp_path):
    coordinates = np.asarray(
        [[0.3, 0.4, 0.5], [1.3, 0.4, 0.5], [8.0, 8.0, 8.0]],
        dtype=np.float64,
    )
    case_dir = tmp_path / "type3_repeated_source_zero_k"
    _write_case(case_dir, coordinates, [(3, 2, 0, 1, 1, 0.2, 0.0)], [])
    _run_case(case_dir)
    observed = _read_coordinates(case_dir, len(coordinates))[2]
    _assert_periodic_coordinate(observed, [0.5, 0.4, 0.5])


def test_type3_rejects_overflowing_construction_geometry(tmp_path):
    case_dir = tmp_path / "type3_overflowing_direction"
    _write_case(
        case_dir,
        [(0.1, 0.2, 0.3), (1.1, 0.2, 0.3), (1.1, 1.2, 0.3), (8, 8, 8)],
        [(3, 3, 0, 1, 2, 0.2, 3.0e38)],
        [],
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "type-3 virtual atom 3" in output
    assert "non-finite or unrepresentable data" in output


@pytest.mark.parametrize(
    "record,error_fragment",
    [
        ("1 2 0 1 0.5 trailing", "unexpected trailing data"),
        ("1 2 0 1 nan", "non-finite virtual-atom parameter"),
        (
            "1 2 0 1 1e9999",
            "virtual-atom parameter is outside the supported finite range",
        ),
        (
            "1 2 0 1 1e-9999",
            "virtual-atom parameter is outside the supported finite range",
        ),
        (
            "1 2147483648 0 1 0.5",
            "virtual-atom integer is outside the supported int range",
        ),
        (
            "1 2 0 2147483648 0.5",
            "virtual-atom integer is outside the supported int range",
        ),
        (
            "1 2 0 1 1e-46",
            "virtual-atom parameter is outside the supported finite float range",
        ),
        (
            "1 2 0 1 1.40129846e-45",
            "virtual-atom parameter is a subnormal float",
        ),
        ("5 3 0 1 2 nan", "non-finite virtual-atom parameter"),
        ("5 3 0 1 2 -0.1", "negative type-5 virtual-atom distance"),
        ("5 3 0 1 2", "invalid type-5 virtual-atom record"),
        ("not a record", "invalid virtual-atom record"),
    ],
)
def test_native_virtual_atom_parser_is_strict_and_source_aware(
    tmp_path, record, error_fragment
):
    case_dir = tmp_path / error_fragment.replace(" ", "_")
    _write_case(
        case_dir,
        [(0.0, 0.0, 0.0)] * 4,
        [(1, 2, 0, 1, 0.5)],
        [],
    )
    (case_dir / "virtual.txt").write_text(record + "\n", encoding="utf-8")
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output
    assert "virtual.txt:1" in output


def test_native_virtual_atom_parser_accepts_long_valid_records(tmp_path):
    case_dir = tmp_path / "long_virtual_atom_record"
    _write_case(
        case_dir,
        [(0.0, 0.0, 0.0)] * 3,
        [(1, 2, 0, 1, 0.5)],
        [],
    )
    separator = " " * 2048
    (case_dir / "virtual.txt").write_text(
        separator.join(("1", "2", "0", "1", "0.5")) + "\n",
        encoding="utf-8",
    )

    result = _run_case(case_dir, check=False)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def _read_potential(case_dir):
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    headers = lines[0].split()
    values = lines[1].split()
    return float(dict(zip(headers, values))["eff_pot"])


def _read_forces(case_dir, atom_numbers):
    raw = np.fromfile(case_dir / "frc.dat", dtype=np.float32)
    assert raw.size >= 3 * atom_numbers
    return raw[-3 * atom_numbers :].reshape(atom_numbers, 3).astype(np.float64)


def _read_coordinates(case_dir, atom_numbers):
    raw = np.fromfile(case_dir / "crd.dat", dtype=np.float32)
    assert raw.size == 3 * atom_numbers
    return raw.reshape(atom_numbers, 3).astype(np.float64)


def _type5_position(oxygen, hydrogen_1, hydrogen_2, distance):
    oxygen = np.asarray(oxygen, dtype=np.float64)
    if distance == 0.0:
        return oxygen.copy()
    oh1 = np.asarray(hydrogen_1, dtype=np.float64) - oxygen
    oh2 = np.asarray(hydrogen_2, dtype=np.float64) - oxygen
    bisector = oh1 / np.linalg.norm(oh1) + oh2 / np.linalg.norm(oh2)
    return oxygen + distance * bisector / np.linalg.norm(bisector)


def _type5_parent_forces(oxygen, hydrogen_1, hydrogen_2, distance, site_force):
    site_force = np.asarray(site_force, dtype=np.float64)
    if distance == 0.0:
        return np.asarray([site_force, np.zeros(3), np.zeros(3)])
    oh1 = np.asarray(hydrogen_1, dtype=np.float64) - np.asarray(
        oxygen, dtype=np.float64
    )
    oh2 = np.asarray(hydrogen_2, dtype=np.float64) - np.asarray(
        oxygen, dtype=np.float64
    )
    unit_oh1 = oh1 / np.linalg.norm(oh1)
    unit_oh2 = oh2 / np.linalg.norm(oh2)
    bisector = unit_oh1 + unit_oh2
    unit_bisector = bisector / np.linalg.norm(bisector)
    q = (
        distance
        / np.linalg.norm(bisector)
        * (site_force - unit_bisector * np.dot(unit_bisector, site_force))
    )
    force_h1 = (q - unit_oh1 * np.dot(unit_oh1, q)) / np.linalg.norm(oh1)
    force_h2 = (q - unit_oh2 * np.dot(unit_oh2, q)) / np.linalg.norm(oh2)
    return np.asarray([site_force - force_h1 - force_h2, force_h1, force_h2])


def _assert_periodic_coordinate(actual, expected, box_length=40.0, atol=5.0e-6):
    difference = np.asarray(actual) - np.asarray(expected)
    difference -= box_length * np.rint(difference / box_length)
    np.testing.assert_allclose(difference, 0.0, rtol=0.0, atol=atol)


def _run_energy_force(tmp_path, name, coordinates, virtual_records, bonds):
    case_dir = tmp_path / name
    _write_case(case_dir, coordinates, virtual_records, bonds)
    _run_case(case_dir)
    return _read_potential(case_dir), _read_forces(case_dir, len(coordinates))


def test_type0_reflects_once_about_the_requested_plane(tmp_path):
    coordinates = np.asarray(
        [[0.7, -0.2, 1.1], [8.0, 8.0, 8.0]], dtype=np.float64
    )
    case_dir = tmp_path / "type0_reflection"
    _write_case(case_dir, coordinates, [(0, 1, 0, 2.25)], [])
    _run_case(case_dir)
    observed = _read_coordinates(case_dir, len(coordinates))[1]
    _assert_periodic_coordinate(observed, [0.7, -0.2, 2.0 * 2.25 - 1.1])


def test_type3_force_accumulator_overflow_fails_before_storing_infinity(
    tmp_path,
):
    parent_coordinates = [[1.0, 1.0, 1.0], [2.0, 1.0, 1.0], [2.0, 2.0, 1.0]]
    virtual_count = 4
    coordinates = [*parent_coordinates]
    coordinates.extend([[8.0, 8.0, 8.0]] * virtual_count)
    coordinates.extend([[2.0, 1.0, 1.0]] * virtual_count)
    virtual_records = [
        (3, 3 + index, 0, 1, 2, 0.0, 0.0) for index in range(virtual_count)
    ]
    bonds = [
        (3 + index, 3 + virtual_count + index, 5.0e37, 0.0)
        for index in range(virtual_count)
    ]
    case_dir = tmp_path / "type3_force_accumulator_overflow"
    _write_case(case_dir, coordinates, virtual_records, bonds)
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "type-3 virtual atom" in output
    assert "overflowing force accumulator" in output


def test_type1_force_redistribution_matches_finite_difference(tmp_path):
    coordinates = np.asarray(
        [
            [0.1, -0.2, 0.3],
            [2.2, 0.7, -0.4],
            [9.0, 9.0, 9.0],
            [3.1, -0.8, 1.4],
        ],
        dtype=np.float64,
    )
    virtual_records = [(1, 2, 0, 1, 0.23)]
    bonds = [(2, 3, 1.7, 0.4)]
    _, forces = _run_energy_force(
        tmp_path, "type1_base", coordinates, virtual_records, bonds
    )

    epsilon = 2.0e-3
    for atom, axis in ((0, 0), (0, 2), (1, 0), (1, 2)):
        displaced_plus = coordinates.copy()
        displaced_minus = coordinates.copy()
        displaced_plus[atom, axis] += epsilon
        displaced_minus[atom, axis] -= epsilon
        energy_plus, _ = _run_energy_force(
            tmp_path,
            f"type1_plus_{atom}_{axis}",
            displaced_plus,
            virtual_records,
            bonds,
        )
        energy_minus, _ = _run_energy_force(
            tmp_path,
            f"type1_minus_{atom}_{axis}",
            displaced_minus,
            virtual_records,
            bonds,
        )
        numerical_force = -(energy_plus - energy_minus) / (2.0 * epsilon)
        assert forces[atom, axis] == pytest.approx(
            numerical_force, rel=3.0e-3, abs=3.0e-3
        )


def _bond_force(position_i, position_j, k, r0):
    displacement = np.asarray(position_i) - np.asarray(position_j)
    distance = np.linalg.norm(displacement)
    return -2.0 * k * (distance - r0) * displacement / distance


def test_type2_shared_parents_are_race_free(tmp_path):
    parent_coordinates = np.asarray(
        [[0.2, -0.1, 0.3], [1.8, 0.4, -0.2], [-0.5, 1.6, 0.7]],
        dtype=np.float64,
    )
    a = 0.21
    b = 0.34
    weights = np.asarray([1.0 - a - b, a, b])
    virtual_position = weights @ parent_coordinates
    site_count = 64
    virtual_start = 3
    probe_start = virtual_start + site_count
    coordinates = [*parent_coordinates.tolist()]
    coordinates.extend([[8.0, 8.0, 8.0]] * site_count)
    probes = []
    for index in range(site_count):
        probes.append(
            [
                2.7 + 0.003 * index,
                -0.9 + 0.002 * (index % 7),
                1.2 - 0.001 * (index % 11),
            ]
        )
    coordinates.extend(probes)
    virtual_records = [
        (2, virtual_start + index, 0, 1, 2, a, b) for index in range(site_count)
    ]
    bonds = [
        (virtual_start + index, probe_start + index, 0.8 + 0.01 * index, 0.5)
        for index in range(site_count)
    ]

    case_dir = tmp_path / "type2_shared"
    _write_case(case_dir, coordinates, virtual_records, bonds)
    expected_parent_force = np.zeros((3, 3), dtype=np.float64)
    expected_probe_forces = []
    for probe, (_, _, k, r0) in zip(probes, bonds):
        force_on_site = _bond_force(virtual_position, probe, k, r0)
        expected_parent_force += weights[:, None] * force_on_site
        expected_probe_forces.append(-force_on_site)

    observed_runs = []
    for _ in range(4):
        _run_case(case_dir, omp_threads=8)
        observed_runs.append(_read_forces(case_dir, len(coordinates)))

    for forces in observed_runs:
        np.testing.assert_allclose(
            forces[:3], expected_parent_force, rtol=4.0e-5, atol=4.0e-5
        )
        np.testing.assert_allclose(
            forces[probe_start:],
            expected_probe_forces,
            rtol=4.0e-5,
            atol=4.0e-5,
        )
    for forces in observed_runs[1:]:
        # Atomic accumulation order is not bitwise deterministic, but the
        # round-off spread must stay at the scale of single-precision sums.
        np.testing.assert_allclose(
            forces, observed_runs[0], rtol=0.0, atol=1.0e-4
        )


def test_type5_coordinate_matches_amber_flexible_bisector(tmp_path):
    coordinates = np.asarray(
        [
            [1.2, -0.7, 2.1],
            [2.31, -0.19, 2.44],
            [0.83, 0.51, 1.72],
            [8.0, 8.0, 8.0],
        ],
        dtype=np.float64,
    )
    distance = 0.1546
    case_dir = tmp_path / "type5_flexible_coordinate"
    _write_case(case_dir, coordinates, [(5, 3, 0, 1, 2, distance)], [])
    _run_case(case_dir)
    observed = _read_coordinates(case_dir, len(coordinates))[3]
    expected = np.asarray([1.26363197, -0.55911780, 2.09791679])
    _assert_periodic_coordinate(observed, expected)


def test_type5_force_redistribution_matches_oracle_and_all_finite_differences(
    tmp_path,
):
    coordinates = np.asarray(
        [
            [1.2, -0.7, 2.1],
            [2.31, -0.19, 2.44],
            [0.83, 0.51, 1.72],
            [8.0, 8.0, 8.0],
            [3.4, -1.2, 0.6],
        ],
        dtype=np.float64,
    )
    distance = 0.1546
    virtual_records = [(5, 3, 0, 1, 2, distance)]
    bonds = [(3, 4, 1.37, 0.42)]
    _, forces = _run_energy_force(
        tmp_path, "type5_flexible_base", coordinates, virtual_records, bonds
    )

    site_position = _type5_position(*coordinates[:3], distance)
    site_force = _bond_force(site_position, coordinates[4], 1.37, 0.42)
    expected_parent_forces = _type5_parent_forces(
        *coordinates[:3], distance, site_force
    )
    np.testing.assert_allclose(
        forces[:3], expected_parent_forces, rtol=5.0e-5, atol=5.0e-5
    )
    np.testing.assert_allclose(forces[3], 0.0, rtol=0.0, atol=1.0e-7)
    np.testing.assert_allclose(forces[4], -site_force, rtol=5.0e-5, atol=5.0e-5)

    epsilon = 2.0e-3
    for atom in range(3):
        for axis in range(3):
            displaced_plus = coordinates.copy()
            displaced_minus = coordinates.copy()
            displaced_plus[atom, axis] += epsilon
            displaced_minus[atom, axis] -= epsilon
            energy_plus, _ = _run_energy_force(
                tmp_path,
                f"type5_plus_{atom}_{axis}",
                displaced_plus,
                virtual_records,
                bonds,
            )
            energy_minus, _ = _run_energy_force(
                tmp_path,
                f"type5_minus_{atom}_{axis}",
                displaced_minus,
                virtual_records,
                bonds,
            )
            numerical_force = -(energy_plus - energy_minus) / (2.0 * epsilon)
            assert forces[atom, axis] == pytest.approx(
                numerical_force, rel=4.0e-3, abs=4.0e-3
            )


def test_type5_coordinate_and_force_use_minimum_image_across_pbc(tmp_path):
    oxygen = np.asarray([39.6, 0.3, 39.5], dtype=np.float64)
    oh1 = np.asarray([1.0, 0.5, 0.3], dtype=np.float64)
    oh2 = np.asarray([-0.4, -1.1, 0.6], dtype=np.float64)
    parents_unwrapped = np.asarray([oxygen, oxygen + oh1, oxygen + oh2])
    parents_wrapped = np.mod(parents_unwrapped, 40.0)
    distance = 0.1546
    site_unwrapped = _type5_position(*parents_unwrapped, distance)
    probe_unwrapped = site_unwrapped + np.asarray([0.7, -0.4, 0.6])
    coordinates = np.vstack(
        [
            parents_wrapped,
            np.asarray([8.0, 8.0, 8.0]),
            np.mod(probe_unwrapped, 40.0),
        ]
    )
    virtual_records = [(5, 3, 0, 1, 2, distance)]
    bonds = [(3, 4, 1.37, 0.42)]
    _, forces = _run_energy_force(
        tmp_path, "type5_cross_pbc", coordinates, virtual_records, bonds
    )

    observed_coordinates = _read_coordinates(
        tmp_path / "type5_cross_pbc", len(coordinates)
    )
    _assert_periodic_coordinate(observed_coordinates[3], site_unwrapped)
    site_force = _bond_force(site_unwrapped, probe_unwrapped, 1.37, 0.42)
    expected_parent_forces = _type5_parent_forces(
        *parents_unwrapped, distance, site_force
    )
    np.testing.assert_allclose(
        forces[:3], expected_parent_forces, rtol=7.0e-5, atol=7.0e-5
    )
    np.testing.assert_allclose(forces[3], 0.0, rtol=0.0, atol=1.0e-7)
    np.testing.assert_allclose(forces[4], -site_force, rtol=7.0e-5, atol=7.0e-5)


def test_type5_near_singular_but_finite_geometry_is_not_threshold_rejected(
    tmp_path,
):
    coordinates = np.asarray(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [-1.0, 5.0e-7, 0.0],
            [8.0, 8.0, 8.0],
            [0.6, 0.4, 0.2],
        ],
        dtype=np.float64,
    )
    distance = 0.125
    case_dir = tmp_path / "type5_near_singular_finite"
    _write_case(
        case_dir,
        coordinates,
        [(5, 3, 0, 1, 2, distance)],
        [(3, 4, 1.0e-6, 0.3)],
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    observed = _read_coordinates(case_dir, len(coordinates))[3]
    expected = _type5_position(*coordinates[:3], distance)
    _assert_periodic_coordinate(observed, expected, atol=2.0e-5)
    assert np.isfinite(_read_forces(case_dir, len(coordinates))).all()


def test_type5_zero_distance_is_exact_without_normalizing_sources(tmp_path):
    coordinates = np.asarray(
        [
            [0.3, -0.2, 0.7],
            [0.3, -0.2, 0.7],
            [0.3, -0.2, 0.7],
            [8.0, 8.0, 8.0],
            [2.1, 0.9, -0.4],
        ],
        dtype=np.float64,
    )
    virtual_records = [(5, 3, 0, 1, 2, 0.0)]
    bonds = [(3, 4, 1.23, 0.51)]
    _, forces = _run_energy_force(
        tmp_path, "type5_zero_distance", coordinates, virtual_records, bonds
    )
    observed_coordinates = _read_coordinates(
        tmp_path / "type5_zero_distance", len(coordinates)
    )
    _assert_periodic_coordinate(observed_coordinates[3], coordinates[0])
    site_force = _bond_force(coordinates[0], coordinates[4], 1.23, 0.51)
    np.testing.assert_allclose(forces[0], site_force, rtol=4.0e-5, atol=4.0e-5)
    np.testing.assert_allclose(forces[1:4], 0.0, rtol=0.0, atol=1.0e-7)
    np.testing.assert_allclose(forces[4], -site_force, rtol=4.0e-5, atol=4.0e-5)


def test_type5_shared_parents_are_race_free(tmp_path):
    parents = np.asarray(
        [[0.2, -0.1, 0.3], [1.7, 0.45, -0.15], [-0.4, 1.8, 0.82]],
        dtype=np.float64,
    )
    distance = 0.125
    site_position = _type5_position(*parents, distance)
    site_count = 64
    virtual_start = 3
    probe_start = virtual_start + site_count
    probes = np.asarray(
        [
            [
                2.9 + 0.004 * index,
                -0.8 + 0.003 * (index % 7),
                1.3 - 0.002 * (index % 11),
            ]
            for index in range(site_count)
        ],
        dtype=np.float64,
    )
    coordinates = np.vstack([parents, np.full((site_count, 3), 8.0), probes])
    virtual_records = [
        (5, virtual_start + index, 0, 1, 2, distance)
        for index in range(site_count)
    ]
    bonds = [
        (virtual_start + index, probe_start + index, 0.9 + 0.01 * index, 0.47)
        for index in range(site_count)
    ]
    expected_parent_force = np.zeros((3, 3), dtype=np.float64)
    expected_probe_forces = []
    for probe, (_, _, k, r0) in zip(probes, bonds):
        site_force = _bond_force(site_position, probe, k, r0)
        expected_parent_force += _type5_parent_forces(
            *parents, distance, site_force
        )
        expected_probe_forces.append(-site_force)

    case_dir = tmp_path / "type5_shared"
    _write_case(case_dir, coordinates, virtual_records, bonds)
    observed_runs = []
    for _ in range(5):
        _run_case(case_dir, omp_threads=8)
        observed_runs.append(_read_forces(case_dir, len(coordinates)))

    for forces in observed_runs:
        np.testing.assert_allclose(
            forces[:3], expected_parent_force, rtol=6.0e-5, atol=6.0e-5
        )
        np.testing.assert_allclose(
            forces[virtual_start:probe_start], 0.0, rtol=0.0, atol=1.0e-7
        )
        np.testing.assert_allclose(
            forces[probe_start:],
            expected_probe_forces,
            rtol=6.0e-5,
            atol=6.0e-5,
        )
    for forces in observed_runs[1:]:
        np.testing.assert_allclose(
            forces, observed_runs[0], rtol=0.0, atol=3.0e-4
        )


@pytest.mark.parametrize(
    "parents",
    [
        [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]],
        [[0.0, 0.0, 0.0], [1.0e-30, 0.0, 0.0], [0.0, 1.0, 0.0]],
    ],
    ids=["zero_oh", "canceling_bisector", "unrepresentable_norm"],
)
def test_type5_singular_geometry_fails_fast(tmp_path, parents):
    case_dir = tmp_path / "type5_singular"
    _write_case(
        case_dir,
        [*parents, [8.0, 8.0, 8.0]],
        [(5, 3, 0, 1, 2, 0.125)],
        [],
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "type-5 virtual atom" in output
    assert "singular O-H or bisector geometry" in output


def test_type5_nested_dependency_uses_reverse_force_order_and_local_copy(
    tmp_path,
):
    coordinates = np.asarray(
        [
            [0.1, 0.2, -0.3],
            [1.4, -0.2, 0.5],
            [-0.4, 1.3, 0.8],
            [8.0, 8.0, 8.0],
            [9.0, 9.0, 9.0],
            [2.8, -0.7, 1.3],
        ],
        dtype=np.float64,
    )
    distance = 0.17
    # The child record deliberately precedes its type-5 parent.
    virtual_records = [
        (1, 3, 4, 2, 0.37),
        (5, 4, 0, 1, 2, distance),
    ]
    bonds = [(3, 5, 1.3, 0.6)]
    energy, forces = _run_energy_force(
        tmp_path, "type5_nested_reverse", coordinates, virtual_records, bonds
    )
    parent_site = _type5_position(*coordinates[:3], distance)
    child_site = parent_site + 0.37 * (coordinates[2] - parent_site)
    expected_energy = (
        1.3 * (np.linalg.norm(child_site - coordinates[5]) - 0.6) ** 2
    )
    assert energy == pytest.approx(expected_energy, rel=3.0e-5, abs=3.0e-5)
    assert np.isfinite(forces).all()
    np.testing.assert_allclose(forces[3:5], 0.0, rtol=0.0, atol=1.0e-7)


def test_type3_initialization_and_reverse_order_dependency(tmp_path):
    coordinates = np.asarray(
        [
            [0.1, 0.2, -0.3],
            [1.4, -0.2, 0.5],
            [2.0, 1.1, -0.4],
            [8.0, 8.0, 8.0],
            [9.0, 9.0, 9.0],
            [2.8, -0.7, 1.3],
        ],
        dtype=np.float64,
    )
    # The child record deliberately precedes the parent record.
    virtual_records = [
        (1, 3, 4, 2, 0.37),
        (3, 4, 0, 1, 2, 1.2, 0.31),
    ]
    bonds = [(3, 5, 1.3, 0.6)]
    energy, forces = _run_energy_force(
        tmp_path,
        "type3_reverse_dependency",
        coordinates,
        virtual_records,
        bonds,
    )

    direction = (coordinates[1] - coordinates[0]) + 0.31 * (
        coordinates[2] - coordinates[1]
    )
    parent_site = coordinates[0] + 1.2 * direction / np.linalg.norm(direction)
    child_site = parent_site + 0.37 * (coordinates[2] - parent_site)
    expected_energy = (
        1.3 * (np.linalg.norm(child_site - coordinates[5]) - 0.6) ** 2
    )
    assert math.isfinite(energy)
    assert energy == pytest.approx(expected_energy, rel=2.0e-5, abs=2.0e-5)
    assert np.isfinite(forces).all()


@pytest.mark.parametrize("distance", [1.2, 0.0])
def test_type3_force_redistribution_matches_finite_difference(
    tmp_path, distance
):
    coordinates = np.asarray(
        [
            [0.1, 0.2, -0.3],
            [1.4, -0.2, 0.5],
            [2.0, 1.1, -0.4],
            [8.0, 8.0, 8.0],
            [2.8, -0.7, 1.3],
        ],
        dtype=np.float64,
    )
    virtual_records = [(3, 3, 0, 1, 2, distance, 0.31)]
    bonds = [(3, 4, 1.3, 0.6)]
    _, forces = _run_energy_force(
        tmp_path,
        f"type3_d_{distance}_base",
        coordinates,
        virtual_records,
        bonds,
    )

    epsilon = 2.0e-3
    for atom in range(3):
        for axis in range(3):
            displaced_plus = coordinates.copy()
            displaced_minus = coordinates.copy()
            displaced_plus[atom, axis] += epsilon
            displaced_minus[atom, axis] -= epsilon
            energy_plus, _ = _run_energy_force(
                tmp_path,
                f"type3_d_{distance}_plus_{atom}_{axis}",
                displaced_plus,
                virtual_records,
                bonds,
            )
            energy_minus, _ = _run_energy_force(
                tmp_path,
                f"type3_d_{distance}_minus_{atom}_{axis}",
                displaced_minus,
                virtual_records,
                bonds,
            )
            numerical_force = -(energy_plus - energy_minus) / (2.0 * epsilon)
            assert forces[atom, axis] == pytest.approx(
                numerical_force, rel=3.0e-3, abs=3.0e-3
            )


@pytest.mark.parametrize(
    ("records", "message"),
    [
        (
            [(1, 2, 3, 0, 0.5), (1, 3, 2, 1, 0.5)],
            "dependency cycle",
        ),
        (
            [(1, 2, 0, 1, 0.5), (1, 2, 0, 1, 0.25)],
            "target of more than one",
        ),
        (
            [(5, 2, 3, 0, 1, 0.125), (1, 3, 2, 1, 0.5)],
            "dependency cycle",
        ),
        ([(1, 2, 0, 4, 0.5)], "source atom index 4 outside"),
        ([(1, 2, 2, 1, 0.5)], "uses its target atom as a source"),
        ([(5, 3, 0, 1, 1, 0.125)], "repeated type-5 source atoms"),
        ([(0, 2, 0, 2.0e38)], "doubled value"),
    ],
)
def test_invalid_virtual_atom_graph_fails_before_initialization(
    tmp_path, records, message
):
    case_dir = tmp_path / message.replace(" ", "_")
    coordinates = [[0.0, 0.0, 0.0]] * 4
    _write_case(case_dir, coordinates, records, [])
    result = _run_case(case_dir, check=False)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr


def test_nested_cv_centers_support_more_than_one_device_warp(tmp_path):
    source_count = 65
    coordinates = np.asarray(
        [[0.2 + 0.05 * index, 0.1, 0.2] for index in range(source_count)]
        + [[7.3, 0.1, 0.2]],
        dtype=np.float64,
    )
    source_indices = "\n".join(str(index) for index in range(source_count))
    equal_weight = 1.0 / source_count
    weights = "\n".join(f"{equal_weight:.17g}" for _ in range(source_count))
    cv_text = "\n".join(
        (
            "many",
            "{",
            "    vatom_type = center",
            "    atom_in_file = sources.list",
            "    weight_in_file = weights.list",
            "}",
            "nested",
            "{",
            "    vatom_type = center",
            f"    atom = many {source_count}",
            "    weight = 0.25 0.75",
            "}",
            "cx",
            "{",
            "    CV_type = position_x",
            "    atom = nested",
            "}",
            "print",
            "{",
            "    CV = cx",
            "}",
            "",
        )
    )
    case_dir = tmp_path / "nested_cv_more_than_warp"
    _write_case(case_dir, coordinates, [], [], cv_text=cv_text)
    (case_dir / "sources.list").write_text(
        source_indices + "\n", encoding="utf-8"
    )
    (case_dir / "weights.list").write_text(weights + "\n", encoding="utf-8")
    _run_case(case_dir)
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    row = dict(zip(lines[0].split(), lines[-1].split()))
    many_x = np.mean(coordinates[:source_count, 0])
    expected_x = 0.25 * many_x + 0.75 * coordinates[source_count, 0]
    assert float(row["cx"]) == pytest.approx(expected_x, rel=2.0e-6, abs=2.0e-6)


@pytest.mark.parametrize(
    "virtual_type",
    [0, 1, 2, 3, 5],
    ids=["type0", "type1", "type2", "type3", "type5"],
)
def test_virtual_atom_local_layout_survives_every_step_dd_refresh(
    tmp_path, mpi_np, virtual_type
):
    source_count = {0: 1, 1: 2, 2: 3, 3: 3, 5: 3}[virtual_type]
    parents = [
        [19.7, 1.0, 1.0],
        [19.9, 1.2, 1.0],
        [19.6, 0.8, 1.2],
    ][:source_count]
    target = source_count
    probe = target + 1
    coordinates = np.asarray(
        [*parents, [8.0, 8.0, 8.0], [20.4, 1.0, 1.0]], dtype=np.float64
    )
    records = {
        0: (0, target, 0, 0.6),
        1: (1, target, 0, 1, 0.37),
        2: (2, target, 0, 1, 2, 0.21, 0.34),
        3: (3, target, 0, 1, 2, 0.2, 0.31),
        5: (5, target, 0, 1, 2, 0.125),
    }
    case_dir = tmp_path / f"dd_refresh_{virtual_type}"
    _write_case(
        case_dir,
        coordinates,
        [records[virtual_type]],
        [(target, probe, 1.3, 0.4)],
        step_limit=3,
        dom_update_interval=1,
    )
    result = _run_case(case_dir, mpi_np=mpi_np, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    forces = _read_forces(case_dir, len(coordinates))
    assert np.isfinite(forces).all()
    np.testing.assert_allclose(forces[target], 0.0, rtol=0.0, atol=1.0e-7)
