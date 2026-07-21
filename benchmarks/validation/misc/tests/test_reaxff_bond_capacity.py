import os
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

from benchmarks.utils import Extractor


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
ATOM_COUNT = 3
EXPECTED_PAIR_COUNT = 2
INITIAL_CAPACITY = 1
GROWN_CAPACITY = 2
DEFAULT_OVERFLOW_ATOM_COUNT = 66
DEFAULT_OVERFLOW_PAIR_COUNT = (
    DEFAULT_OVERFLOW_ATOM_COUNT * (DEFAULT_OVERFLOW_ATOM_COUNT - 1) // 2
)
DEFAULT_INITIAL_CAPACITY = DEFAULT_OVERFLOW_ATOM_COUNT * 32
DEFAULT_GROWN_CAPACITY = DEFAULT_INITIAL_CAPACITY * 3 // 2


def _write_counted_lines(path, values):
    path.write_text(
        str(len(values)) + "\n" + "\n".join(map(str, values)) + "\n",
        encoding="utf-8",
    )


def _water_coordinates():
    angle = np.deg2rad(160.0)
    center = np.array((12.5, 12.5, 12.5), dtype=np.float64)
    bond_length = 1.1
    return np.asarray(
        (
            center,
            center + np.array((bond_length, 0.0, 0.0)),
            center
            + bond_length * np.array((np.cos(angle), np.sin(angle), 0.0)),
        )
    )


def _set_hydrogen_eeq_hardness(force_field, hardness):
    """Keep the synthetic dense case's QEq Hessian positive definite."""
    lines = force_field.read_text(encoding="utf-8").splitlines()
    for line_index, line in enumerate(lines):
        fields = line.split()
        if len(fields) >= 9 and fields[0] == "H":
            eeq_fields = lines[line_index + 1].split()
            assert len(eeq_fields) >= 7
            eeq_fields[6] = f"{hardness:.8g}"
            lines[line_index + 1] = "  " + "  ".join(eeq_fields)
            force_field.write_text(
                "\n".join(lines) + "\n", encoding="utf-8"
            )
            return
    raise AssertionError("failed to find the hydrogen EEQ parameter entry")


def _write_case(
    case_dir,
    initial_capacity=None,
    *,
    coordinates=None,
    masses=("15.999", "1.008", "1.008"),
    types=("O", "H", "H"),
):
    case_dir.mkdir()
    force_field_path = case_dir / "ffield.reax.cho"
    shutil.copy2(REAXFF_STATIC_ROOT / "ffield.reax.cho", force_field_path)
    _write_counted_lines(case_dir / "mass.txt", masses)
    _write_counted_lines(case_dir / "type.txt", types)

    if coordinates is None:
        coordinates = _water_coordinates()
    coordinate_lines = [str(len(coordinates))]
    coordinate_lines.extend(
        " ".join(f"{component:.8f}" for component in coordinate)
        for coordinate in coordinates
    )
    coordinate_lines.extend(("25 25 25", "90 90 90"))
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )

    settings = [
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
    ]
    if initial_capacity is not None:
        settings.append(f"initial_bond_capacity = {initial_capacity}")
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(settings) + "\n",
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


def _last_potential(case_dir):
    mdout = (case_dir / "mdout.txt").read_text(encoding="utf-8")
    return Extractor.extract_sponge_potential(case_dir), mdout


def test_reaxff_bond_storage_grows_without_pair_truncation(tmp_path):
    reference_dir = tmp_path / "reference"
    expanded_dir = tmp_path / "expanded"
    _write_case(reference_dir)
    _write_case(expanded_dir, INITIAL_CAPACITY)

    reference_result = _run_case(reference_dir)
    result = _run_case(expanded_dir)
    reference_output = reference_result.stdout + "\n" + reference_result.stderr
    output = result.stdout + "\n" + result.stderr
    assert reference_result.returncode == 0, reference_output
    assert result.returncode == 0, output

    expansion = (
        f"Expanding sparse bond storage from {INITIAL_CAPACITY} to "
        f"{GROWN_CAPACITY} for {EXPECTED_PAIR_COUNT} bonds"
    )
    assert output.count(expansion) == 1, output
    assert "results may be incorrect" not in output
    assert "num_pairs" not in output

    reference_forces = np.fromfile(reference_dir / "frc.dat", dtype=np.float32)
    expanded_forces = np.fromfile(expanded_dir / "frc.dat", dtype=np.float32)
    assert reference_forces.shape == expanded_forces.shape
    assert reference_forces.size >= 3 * ATOM_COUNT
    assert reference_forces.size % (3 * ATOM_COUNT) == 0
    assert np.all(np.isfinite(reference_forces))
    assert np.all(np.isfinite(expanded_forces))
    np.testing.assert_allclose(
        expanded_forces, reference_forces, rtol=1e-6, atol=1e-5
    )

    reference_potential, reference_mdout = _last_potential(reference_dir)
    expanded_potential, expanded_mdout = _last_potential(expanded_dir)
    assert expanded_potential == reference_potential
    for mdout in (reference_mdout, expanded_mdout):
        assert (
            re.search(r"(^|[^a-z])(?:nan|[+-]?inf)([^a-z]|$)", mdout.lower())
            is None
        )


def test_reaxff_default_capacity_grows_for_a_truly_dense_pair_set(tmp_path):
    # 66 fully connected atoms have 2145 unique pairs, which is greater than
    # the historical fixed atom_count * 32 capacity (2112).  Keeping every
    # separation below one angstrom makes all H-H pairs exceed the ReaxFF
    # bond-order cutoff, so this exercises the old production overflow path
    # without relying on the explicit small-capacity test hook above.
    coordinates = []
    spacing = 0.15
    origin = np.array((12.2, 12.2, 12.3), dtype=np.float64)
    for x in range(5):
        for y in range(5):
            for z in range(3):
                coordinates.append(origin + spacing * np.array((x, y, z)))
    coordinates = np.asarray(coordinates[:DEFAULT_OVERFLOW_ATOM_COUNT])

    case_dir = tmp_path / "default_capacity_dense"
    _write_case(
        case_dir,
        coordinates=coordinates,
        masses=("1.008",) * DEFAULT_OVERFLOW_ATOM_COUNT,
        types=("H",) * DEFAULT_OVERFLOW_ATOM_COUNT,
    )
    # With 65 almost-overlapping neighbors per atom, the stock H QEq
    # off-diagonal couplings overwhelm its diagonal hardness and the charge
    # Hessian is correctly rejected as non-positive-definite/non-convergent.
    # Increase only eta in this synthetic capacity case so the charge solve is
    # well posed while the H-H bond order (and all 2145 stored pairs) is
    # unchanged.
    _set_hydrogen_eeq_hardness(case_dir / "ffield.reax.cho", 2000.0)
    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    expansion = (
        f"Expanding sparse bond storage from {DEFAULT_INITIAL_CAPACITY} to "
        f"{DEFAULT_GROWN_CAPACITY} for {DEFAULT_OVERFLOW_PAIR_COUNT} bonds"
    )
    assert output.count(expansion) == 1, output
    assert "results may be incorrect" not in output

    forces = np.fromfile(case_dir / "frc.dat", dtype=np.float32)
    assert forces.size >= 3 * DEFAULT_OVERFLOW_ATOM_COUNT
    assert forces.size % (3 * DEFAULT_OVERFLOW_ATOM_COUNT) == 0
    assert np.all(np.isfinite(forces))


def test_reaxff_initial_bond_capacity_must_be_positive(tmp_path):
    case_dir = tmp_path / "invalid_capacity"
    _write_case(case_dir, 0)
    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "REAXFF.initial_bond_capacity must be positive" in output
