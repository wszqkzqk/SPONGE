import json
import math
import os
import re
import struct
import subprocess
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


def _write_coordinates(path, coordinates, box_length=80.0):
    lines = [str(len(coordinates))]
    lines.extend(
        " ".join(f"{value:.12g}" for value in xyz) for xyz in coordinates
    )
    lines.extend((f"{box_length} {box_length} {box_length}", "90 90 90"))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_mdin(case_dir, **overrides):
    mdin = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 1,
        "dt": 0,
        "cutoff": 8.0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "LJ_in_file": "lj.txt",
        "mdout": "mdout.txt",
        "crd": "crd.dat",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
    }
    mdin.update(overrides)
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}"
            for key, value in mdin.items()
            if value is not None
        )
        + "\n",
        encoding="utf-8",
    )


def _run_case(case_dir, *, check=True, mpi_np=None):
    command = [
        os.environ.get("SPONGE_BIN", "SPONGE"),
        "-mdin",
        "mdin.spg.toml",
    ]
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


def _read_last_mdout_row(case_dir):
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    return {
        name: float(value)
        for name, value in zip(lines[0].split(), lines[-1].split())
    }


def _read_last_forces(case_dir, atom_count):
    raw = (case_dir / "frc.dat").read_bytes()
    values = struct.unpack(f"={len(raw) // 4}f", raw)
    return values[-3 * atom_count :]


def _float32(value):
    return struct.unpack("=f", struct.pack("=f", value))[0]


def test_zero_parameter_probe_overlap_is_finite_with_solvent_dispatch(tmp_path):
    case_dir = tmp_path / "tip4p_d0_solvent"
    case_dir.mkdir()

    # A zero-parameter solute probe overlaps the first O/M pair. On GPU, ten
    # trailing four-point waters select solvent dispatch, while the triangular
    # probe-water pair is owned by the general kernel. On CPU, the general
    # kernel owns the complete system.
    coordinates = [(10.0, 10.0, 10.0)]
    masses = [12.0]
    charges = [0.0]
    lj_types = [0]
    virtual_records = []
    exclusion_rows = ["0"]
    for water_index in range(10):
        oxygen = (10.0 + 4.0 * water_index, 10.0, 10.0)
        hydrogen_1 = (oxygen[0] + 0.9572, oxygen[1], oxygen[2])
        hydrogen_2 = (oxygen[0] - 0.239987, oxygen[1] + 0.926627, oxygen[2])
        oxygen_index = len(coordinates)
        coordinates.extend((oxygen, hydrogen_1, hydrogen_2, (0.0, 0.0, 0.0)))
        masses.extend((15.999, 1.008, 1.008, 0.0))
        charges.extend((0.0, 0.52422, 0.52422, -1.04844))
        lj_types.extend((1, 0, 0, 0))
        virtual_records.append(
            f"5 {oxygen_index + 3} {oxygen_index} {oxygen_index + 1} "
            f"{oxygen_index + 2} 0"
        )
        exclusion_rows.extend(
            (
                f"3 {oxygen_index + 1} {oxygen_index + 2} {oxygen_index + 3}",
                f"2 {oxygen_index + 2} {oxygen_index + 3}",
                f"1 {oxygen_index + 3}",
                "0",
            )
        )

    _write_counted_values(case_dir / "mass.txt", masses)
    _write_counted_values(case_dir / "charge.txt", charges)
    _write_coordinates(case_dir / "coordinate.txt", coordinates)
    (case_dir / "residue.txt").write_text(
        "41 11\n1\n" + "4\n" * 10, encoding="utf-8"
    )
    (case_dir / "virtual.txt").write_text(
        "\n".join(virtual_records) + "\n", encoding="utf-8"
    )
    (case_dir / "exclude.txt").write_text(
        "41 60\n" + "\n".join(exclusion_rows) + "\n", encoding="utf-8"
    )
    # Pair order is 0-0, 0-1, 1-1. Only O-O has a nonzero LJ interaction;
    # every coincident probe/O/M and O/M pair is exactly inactive.
    (case_dir / "lj.txt").write_text(
        "41 2\n0 0 582000\n0 0 595\n"
        + "\n".join(str(value) for value in lj_types)
        + "\n",
        encoding="utf-8",
    )
    _write_mdin(
        case_dir,
        residue_in_file="residue.txt",
        virtual_atom_in_file="virtual.txt",
        exclude_in_file="exclude.txt",
        solvent_LJ=True,
    )

    result = _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, len(coordinates))
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    # CPU intentionally disables this GPU-only optimization. If it is enabled,
    # prove that this case really selected all ten four-point waters.
    output = result.stdout + result.stderr
    if "optimized hard-LJ solvent dispatch is inactive" not in output:
        assert "the number of solvent atoms is 40" in output
        assert "the solvent is 4-point" in output


def test_ten_coincident_zero_parameter_waters_are_finite_in_solvent_kernel(
    tmp_path,
):
    case_dir = tmp_path / "zero_d0_solvent_kernel"
    case_dir.mkdir()

    # With no leading solute atoms, GPU solvent dispatch owns every pair in
    # this system. All ten four-point residues use the same coordinates, so
    # inter-residue pairs include exact overlaps, but every LJ and Coulomb
    # coefficient is exactly zero. CPU builds exercise the equivalent general
    # kernel and the output check below proves GPU builds selected the intended
    # optimized path.
    water_coordinates = (
        (10.0, 10.0, 10.0),
        (10.9572, 10.0, 10.0),
        (9.760013, 10.926627, 10.0),
        (10.0, 10.0, 10.0),
    )
    coordinates = list(water_coordinates) * 10
    _write_counted_values(
        case_dir / "mass.txt", [15.999, 1.008, 1.008, 0.0] * 10
    )
    _write_counted_values(case_dir / "charge.txt", [0.0] * 40)
    _write_coordinates(case_dir / "coordinate.txt", coordinates)
    (case_dir / "residue.txt").write_text(
        "40 10\n" + "4\n" * 10, encoding="utf-8"
    )
    (case_dir / "lj.txt").write_text(
        "40 1\n0\n0\n" + "0\n" * 40, encoding="utf-8"
    )
    _write_mdin(case_dir, residue_in_file="residue.txt", solvent_LJ=True)

    result = _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 40)
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    assert abs(values["eff_pot"]) < 1.0e-7
    assert max(abs(value) for value in forces) < 1.0e-7
    output = result.stdout + result.stderr
    if "optimized hard-LJ solvent dispatch is inactive" not in output:
        assert "the number of solvent atoms is 40" in output
        assert "the solvent is 4-point" in output


def test_active_coincident_waters_fail_in_owning_hard_nonbond_kernel(tmp_path):
    case_dir = tmp_path / "active_d0_solvent_kernel"
    case_dir.mkdir()
    water_coordinates = (
        (10.0, 10.0, 10.0),
        (10.9572, 10.0, 10.0),
        (9.760013, 10.926627, 10.0),
        (10.0, 10.0, 10.0),
    )
    coordinates = list(water_coordinates) * 10
    _write_counted_values(
        case_dir / "mass.txt", [15.999, 1.008, 1.008, 0.0] * 10
    )
    _write_counted_values(case_dir / "charge.txt", [0.0] * 40)
    _write_coordinates(case_dir / "coordinate.txt", coordinates)
    (case_dir / "residue.txt").write_text(
        "40 10\n" + "4\n" * 10, encoding="utf-8"
    )
    # Type 1 is oxygen and has an active repulsive pair. All other pairs are
    # inactive, so the exact O/O overlap across residues is the sole failure.
    atom_types = [1, 0, 0, 0] * 10
    (case_dir / "lj.txt").write_text(
        "40 2\n0 0 582000\n0 0 0\n"
        + "\n".join(str(value) for value in atom_types)
        + "\n",
        encoding="utf-8",
    )
    _write_mdin(case_dir, residue_in_file="residue.txt", solvent_LJ=True)

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert re.search(r"global atoms \d+ \d+ overlap exactly", output)
    assert "LJ component" in output
    # Dispatch equivalence is covered by test_solvent_dispatch. A fatal owner
    # uses _Exit (CPU) or a device trap (GPU), so buffered initialization logs
    # are deliberately not an oracle at this failure boundary.


@pytest.mark.parametrize("use_sits", [False, True], ids=["standard", "sits"])
def test_zero_parameter_overlap_is_finite_in_softcore_and_selective_kernels(
    tmp_path, use_sits
):
    case_dir = tmp_path / f"zero_softcore_{'sits' if use_sits else 'standard'}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    # A and B each have one atom type. All four coefficient tables are exact
    # zero, followed by the (A-type, B-type) pair for each atom.
    (case_dir / "lj_soft.txt").write_text(
        "2 1 1\n0\n0\n0\n0\n0 0\n0 0\n", encoding="utf-8"
    )
    overrides = {
        "LJ_soft_core_in_file": "lj_soft.txt",
        "lambda_lj": 0.5,
    }
    if use_sits:
        # Observation mode activates selective nonbond kernels without an
        # adaptive bias update. Marking one of two atoms as the selected system
        # keeps the selective path enabled.
        overrides.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    _write_mdin(case_dir, **overrides)

    _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 2)
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    assert abs(values["eff_pot"]) < 1.0e-7
    assert max(abs(value) for value in forces) < 1.0e-7


@pytest.mark.parametrize("use_sits", [False, True], ids=["standard", "sits"])
@pytest.mark.parametrize("active_state", ["A", "B"], ids=["state-a", "state-b"])
@pytest.mark.parametrize(
    "active_component", ["repulsive", "attractive"], ids=["r12", "r6"]
)
def test_active_softcore_endpoint_is_not_skipped_at_overlap(
    tmp_path, use_sits, active_state, active_component
):
    case_dir = (
        tmp_path
        / f"active_{active_state}_{active_component}_{'sits' if use_sits else 'standard'}"
    )
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    coefficient_index = 0 if active_component == "repulsive" else 1
    if active_state == "B":
        coefficient_index += 2
    coefficients = [0.0, 0.0, 0.0, 0.0]
    coefficients[coefficient_index] = 10000.0
    (case_dir / "lj_soft.txt").write_text(
        "2 1 1\n"
        + "\n".join(str(value) for value in coefficients)
        + "\n0 0\n0 0\n",
        encoding="utf-8",
    )
    overrides = {
        "LJ_soft_core_in_file": "lj_soft.txt",
        "lambda_lj": 0.5,
    }
    if use_sits:
        overrides.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    _write_mdin(case_dir, **overrides)

    _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 2)
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    # At lambda=1/2 the active endpoint is softened to a finite distance. Its
    # energy remains nonzero whether only r^-12 or only r^-6 is present; an
    # activity test that only inspects the repulsive coefficient misses the
    # attractive-only transition and evaluates it at the singular r=0.
    if active_component == "repulsive":
        assert 0.1 < values["eff_pot"] < 0.3
    else:
        assert -100.0 < values["eff_pot"] < -1.0
    assert max(abs(value) for value in forces) < 1.0e-7


@pytest.mark.parametrize("use_sits", [False, True], ids=["standard", "sits"])
def test_two_active_attractive_softcore_states_remain_hard(tmp_path, use_sits):
    case_dir = tmp_path / f"two_active_r6_{'sits' if use_sits else 'standard'}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (22.0, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    # Both states contain the same attractive-only -64/r^6 potential. At r=2
    # the unsoftened short-range energy is exactly -1 for every lambda.
    (case_dir / "lj_soft.txt").write_text(
        "2 1 1\n0\n64\n0\n64\n0 0\n0 0\n", encoding="utf-8"
    )
    overrides = {
        "LJ_soft_core_in_file": "lj_soft.txt",
        "lambda_lj": 0.5,
    }
    if use_sits:
        overrides.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    _write_mdin(case_dir, **overrides)

    _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 2)
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    assert values["LJ_soft_short"] == pytest.approx(-1.0, abs=1.0e-6)
    assert max(abs(value) for value in forces) > 1.0e-3


@pytest.mark.parametrize("use_sits", [False, True], ids=["standard", "sits"])
@pytest.mark.parametrize("interaction", ["lj-b-only", "coulomb-only"])
def test_nonzero_hard_pair_components_are_not_skipped(
    tmp_path, use_sits, interaction
):
    case_dir = (
        tmp_path / f"active_{interaction}_{'sits' if use_sits else 'standard'}"
    )
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    charges = [0.0, 0.0]
    pair_b = 0.0
    if interaction == "lj-b-only":
        # At r=2 this attractive-only term has short-range energy -1 exactly.
        pair_b = 64.0
    else:
        charges = [1.0, -1.0]
    _write_counted_values(case_dir / "charge.txt", charges)
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (22.0, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text(
        f"2 1\n0\n{pair_b}\n0\n0\n", encoding="utf-8"
    )
    overrides = {"PM.print_detail": True}
    if use_sits:
        overrides.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    _write_mdin(case_dir, **overrides)

    _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 2)
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    assert max(abs(value) for value in forces) > 1.0e-3
    if interaction == "lj-b-only":
        assert values["LJ_short"] == pytest.approx(-1.0, abs=1.0e-6)
    else:
        assert abs(values["PM_direct"]) > 0.1


@pytest.mark.parametrize("use_sits", [False, True], ids=["standard", "sits"])
@pytest.mark.parametrize(
    ("active_components", "expected_component"),
    [
        ("lj", "LJ component"),
        ("coulomb", "Coulomb component"),
        ("both", "LJ and Coulomb components"),
    ],
)
def test_active_hard_overlap_fails_with_global_pair_and_component(
    tmp_path, use_sits, active_components, expected_component
):
    case_dir = tmp_path / f"hard_overlap_{active_components}_{use_sits}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    charges = (
        [1.0, -1.0] if active_components in ("coulomb", "both") else [0.0, 0.0]
    )
    _write_counted_values(case_dir / "charge.txt", charges)
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0, 20.0, 20.0)],
        box_length=40.0,
    )
    pair_a = 12000.0 if active_components in ("lj", "both") else 0.0
    (case_dir / "lj.txt").write_text(
        f"2 1\n{pair_a}\n0\n0\n0\n", encoding="utf-8"
    )
    overrides = {}
    if use_sits:
        overrides.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    _write_mdin(case_dir, **overrides)

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "global atoms 0 1 overlap exactly" in output
    assert expected_component in output


def test_nopbc_zero_parameter_overlap_is_exactly_inactive(tmp_path):
    case_dir = tmp_path / "nopbc_zero_overlap"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0, 20.0, 20.0)],
        box_length=1000.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    _write_mdin(case_dir, pbc=False, cutoff=100.0)

    _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 2)
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    assert abs(values["eff_pot"]) < 1.0e-7
    assert max(abs(value) for value in forces) < 1.0e-7


@pytest.mark.parametrize(
    ("active_component", "charges", "pair_a", "expected_component"),
    [
        ("lj", [0.0, 0.0], 12000.0, "LJ component"),
        ("coulomb", [1.0, -1.0], 0.0, "Coulomb component"),
    ],
)
def test_nopbc_active_overlap_fails_at_owning_component(
    tmp_path, active_component, charges, pair_a, expected_component
):
    case_dir = tmp_path / f"nopbc_active_{active_component}_overlap"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", charges)
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0, 20.0, 20.0)],
        box_length=1000.0,
    )
    (case_dir / "lj.txt").write_text(
        f"2 1\n{pair_a}\n0\n0\n0\n", encoding="utf-8"
    )
    _write_mdin(case_dir, pbc=False, cutoff=100.0)

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "global atoms 0 1 overlap exactly" in output
    assert expected_component in output


def test_nopbc_gb_pair_energy_is_not_truncated_by_neighbor_cutoff(tmp_path):
    case_dir = tmp_path / "nopbc_gb_full_pair"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [10.0, -10.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(5.0, 10.0, 10.0), (35.0, 10.0, 10.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    # Zero scale factors remove descreening geometry from this analytic case,
    # so the effective radius is exactly radius-offset. The atoms remain 30 A
    # apart, beyond the unrelated 8 A periodic neighbor cutoff.
    (case_dir / "gb.txt").write_text("2\n1.5 0\n1.5 0\n", encoding="utf-8")
    _write_mdin(case_dir, pbc=False, cutoff=8.0, gb_in_file="gb.txt")

    _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    effective_radius = 1.5 - 0.09
    epsilon = 78.5
    charge_square = 100.0
    self_energy = 2.0 * charge_square * (0.5 / epsilon - 0.5) / effective_radius
    pair_energy = (
        -charge_square
        * (1.0 / epsilon - 1.0)
        / math.sqrt(
            30.0**2
            + effective_radius**2
            * math.exp(-0.25 * 30.0**2 / effective_radius**2)
        )
    )
    # GB registers two-decimal output. Omitting the beyond-cutoff pair would
    # differ by about 3.29 kcal/mol and cannot pass this comparison.
    assert values["gb"] == pytest.approx(self_energy + pair_energy, abs=1.1e-2)


def test_nopbc_gb_overlap_fails_before_force_or_energy_publication(tmp_path):
    case_dir = tmp_path / "nopbc_gb_overlap"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    # Keep ordinary Coulomb exactly inactive so GB owns this failure.
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(10.0, 10.0, 10.0), (10.0, 10.0, 10.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    (case_dir / "gb.txt").write_text("2\n1.5 0.8\n1.5 0.8\n", encoding="utf-8")
    _write_mdin(case_dir, pbc=False, cutoff=8.0, gb_in_file="gb.txt")

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "global atoms 0 1 overlap exactly" in output
    assert "generalized-Born pairwise descreening is undefined" in output


def test_nopbc_gb_rejects_nonpositive_effective_radius_before_force(tmp_path):
    case_dir = tmp_path / "nopbc_gb_invalid_effective_radius"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(10.0, 10.0, 10.0), (10.001, 10.0, 10.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    # These individually finite parameters make the pairwise-descreening sum
    # exceed 1/self_radius at this geometry. Consuming the resulting negative
    # effective radius would otherwise produce NaN pair energies.
    (case_dir / "gb.txt").write_text("2\n1.5 2.0\n1.5 2.0\n", encoding="utf-8")
    _write_mdin(case_dir, pbc=False, cutoff=8.0, gb_in_file="gb.txt")

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "atom 0 has invalid effective Born radius" in output
    assert "pairwise-descreening denominator" in output


@pytest.mark.parametrize("use_sits", [False, True], ids=["standard", "sits"])
def test_two_active_softcore_states_fail_as_hard_pair_at_overlap(
    tmp_path, use_sits
):
    case_dir = tmp_path / f"soft_both_active_overlap_{use_sits}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    # Both endpoints contain an r^-12 coefficient, so endpoint activity does
    # not change and this is intentionally a hard pair.
    (case_dir / "lj_soft.txt").write_text(
        "2 1 1\n12000\n0\n12000\n0\n0 0\n0 0\n", encoding="utf-8"
    )
    overrides = {
        "LJ_soft_core_in_file": "lj_soft.txt",
        "lambda_lj": 0.5,
    }
    if use_sits:
        overrides.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    _write_mdin(case_dir, **overrides)

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "global atoms 0 1 overlap exactly" in output
    assert "LJ component" in output


@pytest.mark.parametrize("use_sits", [False, True], ids=["standard", "sits"])
def test_softcore_force_api_does_not_guess_charge_endpoints_at_overlap(
    tmp_path, use_sits
):
    case_dir = tmp_path / f"soft_charge_without_endpoints_{use_sits}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [1.0, -1.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj_soft.txt").write_text(
        # The LJ activity transition legitimately softens LJ. It must not be
        # borrowed as evidence that the charged component changes endpoints.
        "2 1 1\n12000\n0\n0\n0\n0 0\n0 0\n",
        encoding="utf-8",
    )
    overrides = {
        "LJ_in_file": None,
        "LJ_soft_core_in_file": "lj_soft.txt",
        "lambda_lj": 0.5,
    }
    if use_sits:
        overrides.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    _write_mdin(case_dir, **overrides)

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "global atoms 0 1 overlap exactly" in output
    assert "Coulomb component" in output


@pytest.mark.parametrize("separation", [0.0, 1.0e-7, 1.0e-4])
def test_pme_excluded_opposite_charges_are_finite_at_and_near_overlap(
    tmp_path, separation
):
    case_dir = tmp_path / f"pme_excluded_{separation:.0e}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [1.0, -1.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0 + separation, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    # Store the pair once. A single PP rank consumes this triangular exclusion
    # directly, so correction energy and force are not double counted.
    (case_dir / "exclude.txt").write_text("2 1\n1 1\n0\n", encoding="utf-8")
    _write_mdin(
        case_dir,
        exclude_in_file="exclude.txt",
        **{"PM.print_detail": True},
    )

    _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 2)
    assert all(math.isfinite(value) for value in values.values())
    assert all(math.isfinite(value) for value in forces)
    # With q1 + q2 = 0 and the only pair excluded, PME self/reciprocal and the
    # exclusion correction represent a zero physical electrostatic energy.
    assert abs(values["eff_pot"]) < 2.0e-5
    assert max(abs(value) for value in forces) < 2.0e-4


@pytest.mark.parametrize("separation", [0.0, 1.0e-4, 0.7, 0.75])
def test_pme_excluded_radial_kernel_matches_analytic_result(
    tmp_path, separation
):
    case_dir = tmp_path / f"pme_radial_{separation:.0e}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [1.0, -1.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(20.0, 20.0, 20.0), (20.0 + separation, 20.0, 20.0)],
        box_length=40.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    (case_dir / "exclude.txt").write_text("2 1\n1 1\n0\n", encoding="utf-8")
    _write_mdin(
        case_dir,
        exclude_in_file="exclude.txt",
        **{
            "PM.print_detail": True,
            "PME.calculate_reciprocal_part": False,
        },
    )

    result = _run_case(case_dir)
    values = _read_last_mdout_row(case_dir)
    forces = _read_last_forces(case_dir, 2)
    beta_match = re.search(r"\bbeta:\s*([0-9.eE+-]+)", result.stdout)
    assert beta_match is not None
    beta = float(beta_match.group(1))
    distance = _float32(20.0 + separation) - _float32(20.0)
    two_over_sqrt_pi = 2.0 / math.sqrt(math.pi)
    x = beta * distance
    if distance == 0.0:
        expected_energy = two_over_sqrt_pi * beta
        expected_force_x = 0.0
    elif abs(x) <= 0.25:
        x2 = x * x
        energy_series = 1.0 + x2 * (
            -1.0 / 3.0
            + x2
            * (
                1.0 / 10.0
                + x2 * (-1.0 / 42.0 + x2 * (1.0 / 216.0 - x2 / 1320.0))
            )
        )
        force_series = 2.0 / 3.0 + x2 * (
            -2.0 / 5.0 + x2 * (1.0 / 7.0 + x2 * (-1.0 / 27.0 + x2 / 132.0))
        )
        expected_energy = two_over_sqrt_pi * beta * energy_series
        # q_i q_j = -1 and the displacement points in +x.
        expected_force_x = -two_over_sqrt_pi * beta**3 * force_series * distance
    else:
        erf_x = math.erf(x)
        expected_energy = erf_x / distance
        expected_force_x = -(
            erf_x - two_over_sqrt_pi * x * math.exp(-(x * x))
        ) / (distance * distance)

    assert values["eff_pot"] == pytest.approx(
        expected_energy, rel=5.0e-6, abs=2.0e-7
    )
    assert forces[0] == pytest.approx(expected_force_x, rel=1.0e-5, abs=2.0e-7)
    assert forces[3] == pytest.approx(-expected_force_x, rel=1.0e-5, abs=2.0e-7)
    assert max(abs(value) for value in forces[1:3] + forces[4:6]) < 1.0e-7


def test_pme_excluded_pair_outside_spatial_halo_matches_pp_mpi(
    tmp_path, mpi_np
):
    def write_case(case_dir, *, distributed, rank_count=None):
        case_dir.mkdir()
        _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
        _write_counted_values(case_dir / "charge.txt", [1.0, -1.0])
        # The minimum-image distance is 30 A in an 80 A box, well beyond the
        # 8 A cutoff plus its 2 A skin. With an x-only split, both atoms remain
        # at least 15 A from a domain boundary for np=2 and 11.67 A for np=3,
        # so neither can reach its exclusion partner through the spatial halo.
        _write_coordinates(
            case_dir / "coordinate.txt",
            [(15.0, 40.0, 40.0), (65.0, 40.0, 40.0)],
            box_length=80.0,
        )
        (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
        (case_dir / "exclude.txt").write_text("2 1\n1 1\n0\n", encoding="utf-8")
        overrides = {
            "exclude_in_file": "exclude.txt",
            "PM.print_detail": True,
            # Compare only the exclusion correction on both sides. A separate
            # PM rank would add reciprocal/self terms to the direct baseline.
            "PM.MPI_size": 0,
        }
        if distributed:
            # Keep every MPI rank in the PP communicator so the comparison
            # specifically exercises distributed exclusion ownership.
            overrides["DOM_DEC.split_nx"] = rank_count
            overrides["DOM_DEC.split_ny"] = 1
            overrides["DOM_DEC.split_nz"] = 1
        _write_mdin(case_dir, **overrides)

    direct_dir = tmp_path / "pme_excluded_halo_direct"
    write_case(direct_dir, distributed=False)
    _run_case(direct_dir)
    direct_energy = _read_last_mdout_row(direct_dir)["eff_pot"]
    direct_forces = _read_last_forces(direct_dir, 2)
    assert math.isfinite(direct_energy)
    assert all(math.isfinite(value) for value in direct_forces)

    if mpi_np is None:
        return

    mpi_dir = tmp_path / f"pme_excluded_halo_mpi_{mpi_np}"
    write_case(mpi_dir, distributed=True, rank_count=mpi_np)
    _run_case(mpi_dir, mpi_np=mpi_np)
    mpi_energy = _read_last_mdout_row(mpi_dir)["eff_pot"]
    mpi_forces = _read_last_forces(mpi_dir, 2)
    assert math.isfinite(mpi_energy)
    assert all(math.isfinite(value) for value in mpi_forces)
    assert mpi_energy == pytest.approx(direct_energy, rel=3.0e-5, abs=3.0e-6)
    assert mpi_forces == pytest.approx(direct_forces, rel=5.0e-5, abs=5.0e-6)


def test_pme_excluded_mpi_accepts_empty_csr_and_empty_pp_rank(tmp_path, mpi_np):
    if mpi_np is None:
        return

    case_dir = tmp_path / f"pme_excluded_empty_csr_mpi_{mpi_np}"
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [12.0, 12.0])
    _write_counted_values(case_dir / "charge.txt", [0.0, 0.0])
    _write_coordinates(
        case_dir / "coordinate.txt",
        [(15.0, 40.0, 40.0), (65.0, 40.0, 40.0)],
        box_length=80.0,
    )
    (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
    _write_mdin(
        case_dir,
        **{
            "PM.MPI_size": 0,
            "DOM_DEC.split_nx": mpi_np,
            "DOM_DEC.split_ny": 1,
            "DOM_DEC.split_nz": 1,
        },
    )

    result = _run_case(case_dir, check=False, mpi_np=mpi_np)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
