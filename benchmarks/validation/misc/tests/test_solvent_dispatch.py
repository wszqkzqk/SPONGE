import json
import math
import os
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


def _write_coordinates(path, coordinates, box_length=40.0):
    lines = [str(len(coordinates))]
    lines.extend(
        " ".join(f"{value:.12g}" for value in coordinate)
        for coordinate in coordinates
    )
    lines.extend((f"{box_length} {box_length} {box_length}", "90 90 90"))
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _write_water_case(case_dir, *, interaction, optimized, rank_count=None):
    case_dir.mkdir()
    coordinates = []
    for water in range(10):
        oxygen_x = 4.5 + 4.0 * (water % 5)
        oxygen_y = 8.0 + 4.0 * (water // 5)
        coordinates.extend(
            (
                (oxygen_x, oxygen_y, 10.0),
                (oxygen_x + 0.9572, oxygen_y, 10.0),
                (oxygen_x - 0.239987, oxygen_y + 0.926627, 10.0),
            )
        )

    atom_count = len(coordinates)
    _write_counted_values(case_dir / "mass.txt", [15.999, 1.008, 1.008] * 10)
    _write_counted_values(case_dir / "charge.txt", [0.0] * atom_count)
    _write_coordinates(case_dir / "coordinate.txt", coordinates)
    (case_dir / "residue.txt").write_text(
        f"{atom_count} 10\n" + "3\n" * 10, encoding="utf-8"
    )

    # Pair order for two types is 0-0, 0-1, 1-1. Hydrogens use inactive type
    # 0 and oxygens use type 1, leaving a clean O-O dispatch comparison.
    atom_types = ["1", "0", "0"] * 10
    settings = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 1,
        "dt": 0,
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "residue_in_file": "residue.txt",
        "solvent_LJ": optimized,
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_pressure": True,
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
    }
    if interaction == "hard":
        (case_dir / "lj.txt").write_text(
            "30 2\n0 0 582000\n0 0 595\n" + "\n".join(atom_types) + "\n",
            encoding="utf-8",
        )
        settings["LJ_in_file"] = "lj.txt"
    else:
        # State A is exactly inactive and state B has an O-O potential. The
        # old optimized path discarded state B and lambda by converting this
        # record to hard coordinates and always selecting the A table.
        (case_dir / "lj_soft.txt").write_text(
            "30 2 2\n"
            "0 0 0\n"
            "0 0 0\n"
            "0 0 582000\n"
            "0 0 595\n"
            + "\n".join(f"{atom_type} {atom_type}" for atom_type in atom_types)
            + "\n",
            encoding="utf-8",
        )
        settings["LJ_soft_core_in_file"] = "lj_soft.txt"
        settings["lambda_lj"] = 0.4

    if rank_count is not None:
        settings.update(
            {
                "DOM_DEC.split_nx": rank_count,
                "DOM_DEC.split_ny": 1,
                "DOM_DEC.split_nz": 1,
            }
        )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )


def _run_case(case_dir, mpi_np=None):
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
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    return result


def _read_mdout(case_dir):
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    return {
        key: float(value)
        for key, value in zip(lines[0].split(), lines[-1].split())
    }


def _read_forces(case_dir, atom_count=30):
    raw = (case_dir / "frc.dat").read_bytes()
    values = struct.unpack(f"={len(raw) // 4}f", raw)
    return values[-3 * atom_count :]


def _assert_same_observables(reference_dir, candidate_dir):
    reference_mdout = _read_mdout(reference_dir)
    candidate_mdout = _read_mdout(candidate_dir)
    assert candidate_mdout.keys() == reference_mdout.keys()
    for key, reference in reference_mdout.items():
        if key in {"step", "time"}:
            continue
        assert math.isfinite(reference)
        assert math.isfinite(candidate_mdout[key])
        assert candidate_mdout[key] == pytest.approx(
            reference, rel=2.0e-5, abs=2.0e-5
        )

    reference_forces = _read_forces(reference_dir)
    candidate_forces = _read_forces(candidate_dir)
    assert all(math.isfinite(value) for value in reference_forces)
    assert all(math.isfinite(value) for value in candidate_forces)
    assert candidate_forces == pytest.approx(
        reference_forces, rel=3.0e-5, abs=3.0e-6
    )


@pytest.mark.parametrize("interaction", ["hard", "soft"])
def test_solvent_dispatch_matches_general_across_soft_state_and_dd(
    tmp_path, mpi_np, interaction
):
    direct_general = tmp_path / f"{interaction}_direct_general"
    direct_requested = tmp_path / f"{interaction}_direct_requested"
    _write_water_case(direct_general, interaction=interaction, optimized=False)
    _write_water_case(direct_requested, interaction=interaction, optimized=True)
    _run_case(direct_general)
    requested_result = _run_case(direct_requested)
    _assert_same_observables(direct_general, direct_requested)

    direct_forces = _read_forces(direct_general)
    assert max(abs(value) for value in direct_forces) > 1.0e-4
    if interaction == "soft":
        # Soft-core solvent remains fully supported, but the lossy optimized
        # hard-LJ conversion is forbidden; every atom stays in general
        # dispatch so A/B tables and lambda are preserved.
        output = requested_result.stdout + requested_result.stderr
        assert "interactions remain in the general nonbond dispatch" in output

    if mpi_np is None:
        return

    mpi_general = tmp_path / f"{interaction}_mpi_general"
    mpi_requested = tmp_path / f"{interaction}_mpi_requested"
    _write_water_case(
        mpi_general,
        interaction=interaction,
        optimized=False,
        rank_count=mpi_np,
    )
    _write_water_case(
        mpi_requested,
        interaction=interaction,
        optimized=True,
        rank_count=mpi_np,
    )
    _run_case(mpi_general, mpi_np=mpi_np)
    _run_case(mpi_requested, mpi_np=mpi_np)
    _assert_same_observables(mpi_general, mpi_requested)
    _assert_same_observables(direct_general, mpi_requested)
