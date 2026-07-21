import json
import math
import os
import subprocess
from pathlib import Path

import numpy as np
import pytest


def _write_counted_values(path, values):
    Path(path).write_text(
        str(len(values))
        + "\n"
        + "\n".join(str(value) for value in values)
        + "\n",
        encoding="utf-8",
    )


def _write_case(
    case_dir, coordinates, grid, *, print_pressure=False, cmap_terms=None
):
    case_dir.mkdir(parents=True, exist_ok=True)
    atom_count = len(coordinates)
    _write_counted_values(case_dir / "mass.txt", [12.0] * atom_count)
    _write_counted_values(case_dir / "charge.txt", [0.0] * atom_count)

    coordinate_lines = [str(atom_count)]
    coordinate_lines.extend(
        f"{x:.10f} {y:.10f} {z:.10f}" for x, y, z in coordinates
    )
    coordinate_lines.extend(("40.0 40.0 40.0", "90.0 90.0 90.0"))
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )

    resolution = grid.shape[0]
    if cmap_terms is None:
        cmap_terms = [(0, 1, 2, 3, 4, 0)]
    cmap_lines = [str(len(cmap_terms)), "1", str(resolution)]
    cmap_lines.extend(f"{value:.12g}" for value in grid.reshape(-1))
    cmap_lines.extend(" ".join(str(value) for value in term) for term in cmap_terms)
    (case_dir / "cmap.txt").write_text(
        "\n".join(cmap_lines) + "\n", encoding="utf-8"
    )

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
        "cmap_in_file": "cmap.txt",
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


def _read_cmap_energy(case_dir):
    return _read_mdout_row(case_dir)["cmap"]


def _read_mdout_row(case_dir):
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    return {
        key: float(value)
        for key, value in zip(lines[0].split(), lines[1].split())
    }


def _read_forces(case_dir, atom_count):
    raw = np.fromfile(case_dir / "frc.dat", dtype=np.float32)
    assert raw.size >= 3 * atom_count
    return raw[-3 * atom_count :].reshape(atom_count, 3).astype(np.float64)


def _minimum_image(displacement, box):
    return displacement - np.floor(displacement / box + 0.5) * box


def _continuous_cmap_coordinates(coordinates, box):
    coordinates = np.asarray(coordinates, dtype=np.float64)
    box = np.asarray(box, dtype=np.float64)
    drij = _minimum_image(coordinates[0] - coordinates[1], box)
    drkj = _minimum_image(coordinates[2] - coordinates[1], box)
    drkl = _minimum_image(coordinates[2] - coordinates[3], box)
    drml = _minimum_image(coordinates[4] - coordinates[3], box)
    atom_l = drkj - drkl
    return np.asarray(
        [drij, np.zeros(3), drkj, atom_l, atom_l + drml],
        dtype=np.float64,
    )


def _reference_stress(coordinates, forces, box):
    coordinates = np.asarray(coordinates, dtype=np.float64)
    forces = np.asarray(forces, dtype=np.float64)
    box = np.asarray(box, dtype=np.float64)
    scale = 6.946827162543585e4 / float(np.prod(box))
    tensor = coordinates.T @ forces
    return {
        "Pxx": tensor[0, 0] * scale,
        "Pyy": tensor[1, 1] * scale,
        "Pzz": tensor[2, 2] * scale,
        "Pxy": 0.5 * (tensor[0, 1] + tensor[1, 0]) * scale,
        "Pxz": 0.5 * (tensor[0, 2] + tensor[2, 0]) * scale,
        "Pyz": 0.5 * (tensor[1, 2] + tensor[2, 1]) * scale,
    }


def _natural_spline_second_derivatives(values):
    values = np.asarray(values, dtype=np.float64)
    count = len(values)
    second = np.zeros(count, dtype=np.float64)
    if count <= 2:
        return second
    work = np.zeros(count, dtype=np.float64)
    for index in range(1, count - 1):
        pivot = 0.5 * second[index - 1] + 2.0
        second[index] = -0.5 / pivot
        second_difference = (
            values[index + 1] - 2.0 * values[index] + values[index - 1]
        )
        work[index] = (3.0 * second_difference - 0.5 * work[index - 1]) / pivot
    for index in range(count - 2, -1, -1):
        second[index] = second[index] * second[index + 1] + work[index]
    return second


def _natural_spline_interpolate(values, second, coordinate):
    values = np.asarray(values, dtype=np.float64)
    second = np.asarray(second, dtype=np.float64)
    left = int(coordinate)
    lower_weight = left + 1 - coordinate
    upper_weight = coordinate - left
    value = (
        lower_weight * values[left]
        + upper_weight * values[left + 1]
        + (
            (lower_weight**3 - lower_weight) * second[left]
            + (upper_weight**3 - upper_weight) * second[left + 1]
        )
        / 6.0
    )
    derivative = (
        values[left + 1]
        - values[left]
        - (3.0 * lower_weight**2 - 1.0) * second[left] / 6.0
        + (3.0 * upper_weight**2 - 1.0) * second[left + 1] / 6.0
    )
    return value, derivative


def _gromacs_tiled_spline_samples(values):
    values = np.asarray(values, dtype=np.float64)
    count = len(values)
    indices = (np.arange(2 * count) + count - count // 2) % count
    tiled = values[indices]
    second = _natural_spline_second_derivatives(tiled)
    coordinates = np.arange(count, dtype=np.float64) + 0.5 * count
    samples = [
        _natural_spline_interpolate(tiled, second, coordinate)
        for coordinate in coordinates
    ]
    return np.asarray(samples, dtype=np.float64)


def _gromacs_cmap_knot_data(grid):
    grid = np.asarray(grid, dtype=np.float64)
    resolution = grid.shape[0]
    extended_resolution = 2 * resolution
    indices = (
        np.arange(extended_resolution) + resolution - resolution // 2
    ) % resolution
    extended = grid[np.ix_(indices, indices)]
    coordinates = np.arange(resolution, dtype=np.float64) + 0.5 * resolution

    values_at_psi = np.empty(
        (extended_resolution, resolution), dtype=np.float64
    )
    dpsi_at_psi = np.empty_like(values_at_psi)
    for phi in range(extended_resolution):
        second = _natural_spline_second_derivatives(extended[phi])
        for psi, coordinate in enumerate(coordinates):
            values_at_psi[phi, psi], dpsi_at_psi[phi, psi] = (
                _natural_spline_interpolate(extended[phi], second, coordinate)
            )

    dphi = np.empty_like(grid)
    dpsi = np.empty_like(grid)
    mixed = np.empty_like(grid)
    for psi in range(resolution):
        line = values_at_psi[:, psi]
        second = _natural_spline_second_derivatives(line)
        for phi, coordinate in enumerate(coordinates):
            _, dphi[phi, psi] = _natural_spline_interpolate(
                line, second, coordinate
            )

        line = dpsi_at_psi[:, psi]
        second = _natural_spline_second_derivatives(line)
        for phi, coordinate in enumerate(coordinates):
            dpsi[phi, psi], mixed[phi, psi] = _natural_spline_interpolate(
                line, second, coordinate
            )
    return grid, dphi, dpsi, mixed


def _hermite(value_0, value_1, derivative_0, derivative_1, fraction):
    fraction_2 = fraction * fraction
    fraction_3 = fraction_2 * fraction
    return (
        (2.0 * fraction_3 - 3.0 * fraction_2 + 1.0) * value_0
        + (fraction_3 - 2.0 * fraction_2 + fraction) * derivative_0
        + (-2.0 * fraction_3 + 3.0 * fraction_2) * value_1
        + (fraction_3 - fraction_2) * derivative_1
    )


def _signed_dihedral(atom_i, atom_j, atom_k, atom_l):
    drij = atom_i - atom_j
    drkj = atom_k - atom_j
    drkl = atom_k - atom_l
    normal_1 = np.cross(drij, drkj)
    normal_2 = np.cross(drkl, drkj)
    inverse_norms = 1.0 / (np.linalg.norm(normal_1) * np.linalg.norm(normal_2))
    cosine = -np.dot(normal_1, normal_2) * inverse_norms
    sine = (
        np.dot(np.cross(normal_2, normal_1), drkj)
        * inverse_norms
        / np.linalg.norm(drkj)
    )
    return math.atan2(sine, cosine)


def _reference_energy(coordinates, grid):
    phi = _signed_dihedral(*coordinates[:4])
    psi = _signed_dihedral(*coordinates[1:])
    resolution = grid.shape[0]
    phi_grid = (phi + math.pi) * resolution / (2.0 * math.pi)
    psi_grid = (psi + math.pi) * resolution / (2.0 * math.pi)
    values, dphi, dpsi, mixed = _gromacs_cmap_knot_data(grid)
    phi_floor = math.floor(phi_grid)
    psi_floor = math.floor(psi_grid)
    phi_0 = phi_floor % resolution
    phi_1 = (phi_0 + 1) % resolution
    psi_0 = psi_floor % resolution
    psi_1 = (psi_0 + 1) % resolution
    phi_fraction = phi_grid - phi_floor
    psi_fraction = psi_grid - psi_floor

    value_at_phi_0 = _hermite(
        values[phi_0, psi_0],
        values[phi_0, psi_1],
        dpsi[phi_0, psi_0],
        dpsi[phi_0, psi_1],
        psi_fraction,
    )
    value_at_phi_1 = _hermite(
        values[phi_1, psi_0],
        values[phi_1, psi_1],
        dpsi[phi_1, psi_0],
        dpsi[phi_1, psi_1],
        psi_fraction,
    )
    dphi_at_phi_0 = _hermite(
        dphi[phi_0, psi_0],
        dphi[phi_0, psi_1],
        mixed[phi_0, psi_0],
        mixed[phi_0, psi_1],
        psi_fraction,
    )
    dphi_at_phi_1 = _hermite(
        dphi[phi_1, psi_0],
        dphi[phi_1, psi_1],
        mixed[phi_1, psi_0],
        mixed[phi_1, psi_1],
        psi_fraction,
    )
    return _hermite(
        value_at_phi_0,
        value_at_phi_1,
        dphi_at_phi_0,
        dphi_at_phi_1,
        phi_fraction,
    )


def _asymmetric_grid(resolution):
    angles = -math.pi + np.arange(resolution) * 2.0 * math.pi / resolution
    phi, psi = np.meshgrid(angles, angles, indexing="ij")
    return (
        0.37
        + 0.83 * np.sin(phi)
        - 0.29 * np.cos(psi)
        + 0.41 * np.sin(phi - 2.0 * psi)
        + 0.17 * np.cos(2.0 * phi + psi)
    )


def test_gromacs_tiled_natural_spline_reference():
    values = np.asarray([0.0, 1.0, 4.0, -2.0, 3.0])
    samples = _gromacs_tiled_spline_samples(values)
    # setup_cmap uses xmin=-360 degrees.  With an odd-sized map the central
    # output coordinates therefore fall halfway between extended samples;
    # this regression intentionally covers that non-obvious behavior.
    np.testing.assert_allclose(
        samples,
        [
            [-0.399417314095449, 1.055216426193119],
            [3.466009988901221, 4.188401775804662],
            [0.660377358490566, -8.058823529411764],
            [0.142480577136515, 7.046892341842398],
            [2.144700332963374, -4.378745837957824],
        ],
        rtol=0.0,
        atol=1.0e-12,
    )

    indices = (np.arange(10) + 5 - 5 // 2) % 5
    tiled = values[indices]
    second = _natural_spline_second_derivatives(tiled)
    quarter_value, _ = _natural_spline_interpolate(tiled, second, 2.75)
    assert quarter_value == pytest.approx(0.085790094339623, abs=1.0e-12)


@pytest.mark.parametrize(
    "coordinates",
    [
        np.asarray(
            [
                [0.3, 0.4, 0.2],
                [1.0, 0.9, 0.5],
                [1.8, 0.5, 1.1],
                [2.4, 1.2, 1.6],
                [3.0, 0.8, 2.5],
            ],
            dtype=np.float64,
        ),
        # The first torsion is very close to zero.  The old acos clamp moved
        # it by about 1.4e-3 rad and made its reported energy and force
        # inconsistent.
        np.asarray(
            [
                [0.0, 1.0, 0.0],
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [1.0, 1.0, 1.0e-5],
                [1.8, 1.4, 0.9],
            ],
            dtype=np.float64,
        ),
    ],
    ids=["generic", "near-planar"],
)
def test_periodic_cmap_energy_and_force_match_reference(tmp_path, coordinates):
    grid = _asymmetric_grid(5)
    case_dir = tmp_path / "base"
    _write_case(case_dir, coordinates, grid)
    _run_case(case_dir)

    assert _read_cmap_energy(case_dir) == pytest.approx(
        _reference_energy(coordinates, grid), rel=2.0e-5, abs=2.0e-5
    )
    forces = _read_forces(case_dir, len(coordinates))
    assert np.isfinite(forces).all()

    epsilon = 2.0e-3
    for atom, axis in ((0, 2), (1, 0), (2, 1), (3, 2), (4, 0)):
        plus = coordinates.copy()
        minus = coordinates.copy()
        plus[atom, axis] += epsilon
        minus[atom, axis] -= epsilon
        plus_dir = tmp_path / f"plus_{atom}_{axis}"
        minus_dir = tmp_path / f"minus_{atom}_{axis}"
        _write_case(plus_dir, plus, grid)
        _write_case(minus_dir, minus, grid)
        _run_case(plus_dir)
        _run_case(minus_dir)
        numerical_force = -(
            _read_cmap_energy(plus_dir) - _read_cmap_energy(minus_dir)
        ) / (2.0 * epsilon)
        assert forces[atom, axis] == pytest.approx(
            numerical_force, rel=4.0e-3, abs=4.0e-3
        )


def test_periodic_images_preserve_cmap_energy_force_and_stress(
    tmp_path, mpi_np
):
    coordinates = np.asarray(
        [
            [0.3, 0.4, 0.2],
            [1.0, 0.9, 0.5],
            [1.8, 0.5, 1.1],
            [2.4, 1.2, 1.6],
            [3.0, 0.8, 2.5],
        ],
        dtype=np.float64,
    )
    # Keep the five atoms close in minimum-image space while straddling the
    # midpoint in every Cartesian direction.  A two-rank decomposition can
    # therefore split the CMAP term regardless of which box axis it chooses.
    coordinates += np.asarray([18.5, 19.5, 18.5], dtype=np.float64)
    box = np.asarray([40.0, 40.0, 40.0], dtype=np.float64)
    image_shifts = np.asarray(
        [
            [1, 0, -1],
            [0, -1, 0],
            [-1, 1, 0],
            [0, 0, 1],
            [1, -1, 1],
        ],
        dtype=np.float64,
    )
    imaged_coordinates = coordinates + image_shifts * box
    grid = _asymmetric_grid(5)

    results = []
    for name, case_coordinates in (
        ("base", coordinates),
        ("periodic_images", imaged_coordinates),
    ):
        case_dir = tmp_path / name
        _write_case(case_dir, case_coordinates, grid, print_pressure=True)
        _run_case(case_dir, mpi_np=mpi_np)
        mdout = _read_mdout_row(case_dir)
        forces = _read_forces(case_dir, len(case_coordinates))
        continuous_coordinates = _continuous_cmap_coordinates(
            case_coordinates, box
        )
        assert mdout["cmap"] == pytest.approx(
            _reference_energy(continuous_coordinates, grid),
            rel=3.0e-5,
            abs=3.0e-5,
        )
        expected_stress = _reference_stress(continuous_coordinates, forces, box)
        for component, expected in expected_stress.items():
            assert mdout[component] == pytest.approx(expected, abs=0.02)
        expected_pressure = (
            sum(
                expected_stress[component]
                for component in ("Pxx", "Pyy", "Pzz")
            )
            / 3.0
        )
        assert mdout["pressure"] == pytest.approx(expected_pressure, abs=0.02)
        results.append((mdout, forces))

    base_mdout, base_forces = results[0]
    imaged_mdout, imaged_forces = results[1]
    assert imaged_mdout["cmap"] == pytest.approx(
        base_mdout["cmap"], rel=3.0e-5, abs=3.0e-5
    )
    np.testing.assert_allclose(
        imaged_forces, base_forces, rtol=5.0e-5, atol=5.0e-5
    )
    for component in ("pressure", "Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz"):
        assert imaged_mdout[component] == pytest.approx(
            base_mdout[component], abs=0.02
        )


@pytest.mark.parametrize(
    ("bad_cmap", "message"),
    [
        ("-1\n0\n", "negative CMAP interaction count"),
        ("0\n-1\n", "negative CMAP type count"),
        ("2147483648\n", "strict signed integer"),
        (
            "2147483647\n0\n",
            "truncated while reading interaction 0 atom A",
        ),
        (
            "0\n2147483647\n",
            "truncated while reading resolution for CMAP type 0",
        ),
        ("0\n1\n1suffix\n", "strict signed integer"),
        ("0\n1\n", "truncated while reading resolution"),
        ("0\n1\n0\n", "non-positive resolution"),
        ("0\n1\n10000\n", "truncated while reading grid value 0"),
        ("0\n1\n1\n1.0suffix\n", "invalid or outside"),
        ("0\n1\n1\n1e39\n", "outside the finite float range"),
        ("0\n1\n1\n1e-46\n", "cannot be represented as a finite float"),
        ("0\n1\n1\nnan\n", "outside the finite float range"),
        ("0\n1\n1\n0.0\ntrailing\n", "trailing data"),
        (
            "1\n1\n1\n0.0\n0 1 0 3 4 0\n",
            "repeats atom 0 within its first torsion",
        ),
        (
            "1\n1\n1\n0.0\n0 1 2 3 2 0\n",
            "repeats atom 2 within its second torsion",
        ),
    ],
)
def test_invalid_native_cmap_fails_before_initialization(
    tmp_path, bad_cmap, message
):
    coordinates = np.asarray(
        [[0.0, 0.0, 0.0]] * 5,
        dtype=np.float64,
    )
    case_dir = tmp_path / message.replace(" ", "_")
    _write_case(case_dir, coordinates, np.zeros((1, 1)))
    (case_dir / "cmap.txt").write_text(bad_cmap, encoding="utf-8")
    result = _run_case(case_dir, check=False)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr


def test_native_cmap_rejects_representable_float_subnormal(tmp_path):
    coordinates = np.zeros((5, 3), dtype=np.float64)
    case_dir = tmp_path / "representable_subnormal"
    _write_case(case_dir, coordinates, np.zeros((1, 1)))
    (case_dir / "cmap.txt").write_text(
        "0\n1\n1\n1.401298464324817e-45\n", encoding="utf-8"
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "cmap_in_file grid value 0 token" in output
    assert "is a subnormal float" in output
    assert "consistent FTZ behavior" in output
    assert "Input file: cmap.txt" in output


def test_cmap_rejects_subnormal_derived_interpolation_coefficient(tmp_path):
    coordinates = np.zeros((5, 3), dtype=np.float64)
    case_dir = tmp_path / "derived_interpolation_subnormal"
    _write_case(case_dir, coordinates, np.zeros((2, 2)))
    # Both inputs are normal float32 values (FLT_MIN and its next larger
    # neighbor), but their finite difference is the smallest subnormal.  A
    # CPU build that preserves it and a GPU FTZ build would otherwise create
    # different interpolation maps from the same input.
    (case_dir / "cmap.txt").write_text(
        "1\n1\n2\n"
        "1.1754943508222875e-38\n"
        "1.175494490952134e-38\n"
        "1.1754943508222875e-38\n"
        "1.175494490952134e-38\n"
        "0 1 2 3 4 0\n",
        encoding="utf-8",
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    preserved_subnormal = "subnormal interpolation coefficient" in output
    flushed_underflow = (
        "interpolation coefficient outside the finite float range" in output
    )
    assert preserved_subnormal or flushed_underflow, output
    if preserved_subnormal:
        assert "consistent FTZ behavior" in output


def test_collinear_cmap_geometry_fails_instead_of_producing_nan(tmp_path):
    coordinates = np.asarray(
        [[float(atom), 0.0, 0.0] for atom in range(5)], dtype=np.float64
    )
    case_dir = tmp_path / "collinear"
    _write_case(case_dir, coordinates, np.zeros((1, 1), dtype=np.float64))
    result = _run_case(case_dir, check=False)
    assert result.returncode != 0
    assert "undefined or non-finite torsion geometry" in (
        result.stdout + result.stderr
    )


def test_nonfinite_cmap_energy_fails_before_accumulation(tmp_path):
    # Every interpolated coefficient remains representable as float, but this
    # cell of the bicubic polynomial overshoots FLT_MAX.  Checking only input
    # knots, coefficients, or Cartesian forces would let an infinite energy
    # reach the atom-energy accumulator.
    grid = 2.0e38 * np.asarray(
        [
            [1.0, 1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0, 1.0],
            [1.0, 1.0, 1.0, 1.0],
        ],
        dtype=np.float64,
    )
    coordinates = np.asarray(
        [
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 0.7071067812, -0.7071067812],
            [1.71153568, 0.21025822, -1.20395534],
        ],
        dtype=np.float64,
    )
    case_dir = tmp_path / "nonfinite_energy"
    _write_case(case_dir, coordinates, grid)
    result = _run_case(case_dir, check=False)
    assert result.returncode != 0
    assert "undefined or non-finite torsion geometry/energy" in (
        result.stdout + result.stderr
    )


def test_cmap_atom_energy_accumulation_overflow_fails(tmp_path):
    coordinates = np.asarray(
        [
            [0.3, 0.4, 0.2],
            [1.0, 0.9, 0.5],
            [1.8, 0.5, 1.1],
            [2.4, 1.2, 1.6],
            [3.0, 0.8, 2.5],
        ],
        dtype=np.float64,
    )
    grid = np.full((1, 1), 2.0e38, dtype=np.float64)
    case_dir = tmp_path / "atom_energy_accumulation_overflow"
    _write_case(
        case_dir,
        coordinates,
        grid,
        cmap_terms=[(0, 1, 2, 3, 4, 0), (0, 1, 2, 3, 4, 0)],
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "cmap CMAP term" in output
    assert "accumulator" in output


def test_cmap_total_energy_reduction_overflow_fails(tmp_path):
    first_coordinates = np.asarray(
        [
            [0.3, 0.4, 0.2],
            [1.0, 0.9, 0.5],
            [1.8, 0.5, 1.1],
            [2.4, 1.2, 1.6],
            [3.0, 0.8, 2.5],
        ],
        dtype=np.float64,
    )
    coordinates = np.concatenate(
        [first_coordinates, first_coordinates + np.asarray([5.0, 0.0, 0.0])]
    )
    grid = np.full((1, 1), 2.0e38, dtype=np.float64)
    case_dir = tmp_path / "total_energy_reduction_overflow"
    _write_case(
        case_dir,
        coordinates,
        grid,
        cmap_terms=[(0, 1, 2, 3, 4, 0), (5, 6, 7, 8, 9, 0)],
    )
    result = _run_case(case_dir, check=False)
    assert result.returncode != 0
    output = result.stdout + result.stderr
    assert "global potential-energy reduction is not finite" in output
    assert "system accumulator" in output
