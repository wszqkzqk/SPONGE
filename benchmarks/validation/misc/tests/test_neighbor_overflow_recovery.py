import json
import os
import subprocess

import numpy as np
import pytest


def _write_case(case_dir, *, constrained):
    case_dir.mkdir()
    (case_dir / "mass.txt").write_text("3\n12\n12\n12\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("3\n0\n0\n0\n", encoding="utf-8")
    (case_dir / "coordinate.txt").write_text(
        "3\n9.0 9.0 9.0\n9.7 9.0 9.0\n9.2 9.6 9.0\n"
        "24 24 24\n90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "lj.txt").write_text(
        "3 1\n1\n0\n0\n0\n0\n", encoding="utf-8"
    )
    settings = {
        "mode": "nve",
        "step_limit": 0,
        "dt": 0.0,
        "cutoff": 4.0,
        "skin": 1.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "LJ_in_file": "lj.txt",
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_restart_file_interval": 0,
    }
    if constrained:
        settings.update(
            {
                "neighbor_list.max_neighbor_numbers": 1,
                "neighbor_list.max_atom_in_grid_numbers": 1,
            }
        )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(f"{key} = {json.dumps(value)}" for key, value in settings.items())
        + "\n",
        encoding="utf-8",
    )


def _run(case_dir):
    return subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


def _last_mdout(case_dir):
    rows = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    return {name: float(value) for name, value in zip(rows[0].split(), rows[-1].split())}


def test_neighbor_storage_overflow_grows_and_retries_same_state(tmp_path):
    reference_dir = tmp_path / "reference"
    constrained_dir = tmp_path / "constrained"
    _write_case(reference_dir, constrained=False)
    _write_case(constrained_dir, constrained=True)

    reference = _run(reference_dir)
    constrained = _run(constrained_dir)
    assert reference.returncode == 0, reference.stdout + reference.stderr
    output = constrained.stdout + constrained.stderr
    assert constrained.returncode == 0, output
    assert "Neighbor-list capacity was insufficient for the current state" in output

    reference_force = np.fromfile(reference_dir / "frc.dat", dtype=np.float32)
    constrained_force = np.fromfile(constrained_dir / "frc.dat", dtype=np.float32)
    np.testing.assert_allclose(
        constrained_force, reference_force, rtol=1.0e-6, atol=1.0e-6
    )
    assert _last_mdout(constrained_dir)["potential"] == pytest.approx(
        _last_mdout(reference_dir)["potential"], abs=1.0e-6
    )
