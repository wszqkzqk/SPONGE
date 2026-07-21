import ctypes
import json
import math
import os
import struct
import subprocess
import threading

import pytest

BOX_LENGTH = 10.0
PAIR_DISTANCE = 2.0
CROSS_FACTOR = 0.4
AMD_ALPHA = 2.0
AMD_THRESHOLD = 1.0
GAMD_K = 0.4
GAMD_THRESHOLD = 1.2
EMPIRICAL_SMOOTHING_ENERGY = 0.5
EMPIRICAL_THRESHOLD = 1.0
EMPIRICAL_T_LOW = 300.0
EMPIRICAL_T_HIGH = 600.0
REFERENCE_TEMPERATURE = 300.0
PRESSURE_CONVERSION = 6.946827162543585e4
REPULSIVE_COEFFICIENT = 4096.0
SOFT_CORE_LAMBDA = 0.25
SOFT_CORE_ALPHA = 0.5
SOFT_CORE_SIGMA = 1.0


def _write_counted_values(path, values):
    path.write_text(
        str(len(values))
        + "\n"
        + "\n".join(f"{value:.12g}" for value in values)
        + "\n",
        encoding="utf-8",
    )


def _write_coordinates(path, distance):
    center = 0.5 * BOX_LENGTH
    half_distance = 0.5 * distance
    path.write_text(
        "\n".join(
            (
                "2",
                f"{center - half_distance:.12g} {center:.12g} {center:.12g}",
                f"{center + half_distance:.12g} {center:.12g} {center:.12g}",
                f"{BOX_LENGTH:.12g} {BOX_LENGTH:.12g} {BOX_LENGTH:.12g}",
                "90 90 90",
            )
        )
        + "\n",
        encoding="utf-8",
    )


def _write_case(
    case_dir,
    distance,
    potential_kind,
    sits_atom_numbers=1,
    sits_overrides=None,
    extra_files=None,
):
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", (12.0, 12.0))
    _write_counted_values(case_dir / "charge.txt", (0.0, 0.0))
    _write_coordinates(case_dir / "coordinate.txt", distance)

    settings = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 0,
        "dt": 0,
        "cutoff": 4.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_pressure": True,
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
        "SITS.mode": "amd",
        "SITS.atom_numbers": sits_atom_numbers,
        "SITS.cross_enhance_factor": CROSS_FACTOR,
        "SITS.pe_a": AMD_ALPHA,
        "SITS.pe_b": AMD_THRESHOLD,
    }
    if potential_kind == "hard":
        # Native input stores U(r) = A / r^12; B=0 also makes the LJ long-range
        # correction exactly zero, so the pair is the whole Hamiltonian.
        (case_dir / "lj.txt").write_text(
            f"2 1\n{REPULSIVE_COEFFICIENT:.12g}\n0\n0\n0\n",
            encoding="utf-8",
        )
        settings["LJ_in_file"] = "lj.txt"
    else:
        # State A contains the same repulsive pair and state B is inactive.
        # This forces the selective soft-core branch rather than merely using
        # the soft-core module with two hard active endpoints.
        (case_dir / "lj_soft.txt").write_text(
            f"2 1 1\n{REPULSIVE_COEFFICIENT:.12g}\n0\n0\n0\n0 0\n0 0\n",
            encoding="utf-8",
        )
        settings.update(
            {
                "LJ_soft_core_in_file": "lj_soft.txt",
                "lambda_lj": SOFT_CORE_LAMBDA,
                "soft_core_alpha": SOFT_CORE_ALPHA,
                "soft_core_sigma": SOFT_CORE_SIGMA,
            }
        )

    for key, value in (sits_overrides or {}).items():
        if value is None:
            settings.pop(key, None)
        else:
            settings[key] = value
    for name, content in (extra_files or {}).items():
        (case_dir / name).write_text(content, encoding="utf-8")

    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )


def _run_state(
    case_dir,
    distance,
    potential_kind,
    sits_atom_numbers=1,
    sits_overrides=None,
):
    _write_case(
        case_dir,
        distance,
        potential_kind,
        sits_atom_numbers,
        sits_overrides=sits_overrides,
    )
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"SPONGE failed with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    mdout_lines = (
        (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    )
    mdout = {
        name: float(value)
        for name, value in zip(mdout_lines[0].split(), mdout_lines[-1].split())
    }
    raw_forces = (case_dir / "frc.dat").read_bytes()
    force_values = struct.unpack(f"={len(raw_forces) // 4}f", raw_forces)[-6:]
    forces = (force_values[:3], force_values[3:])
    return mdout, forces


def _run_sponge_raw(case_dir, **kwargs):
    return subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
        **kwargs,
    )


def _sits_iteration_persistence_overrides(
    *, trajectory_name="nk_traj.dat", restart_name="nk_rest.txt"
):
    return {
        "SITS.mode": "iteration",
        "SITS.k_numbers": 2,
        "SITS.T": "300/600",
        "SITS.record_interval": 1,
        "SITS.update_interval": 1,
        "SITS.pe_a": 1.0,
        "SITS.pe_b": 0.0,
        "SITS.nk_traj_file": trajectory_name,
        "SITS.nk_rest_file": restart_name,
    }


def _run_rerun_states(case_dir, distances, sits_overrides=None):
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", (12.0, 12.0))
    _write_counted_values(case_dir / "charge.txt", (0.0, 0.0))
    _write_coordinates(case_dir / "coordinate.txt", distances[0])
    (case_dir / "lj.txt").write_text(
        f"2 1\n{REPULSIVE_COEFFICIENT:.12g}\n0\n0\n0\n",
        encoding="utf-8",
    )

    center = 0.5 * BOX_LENGTH
    trajectory_values = []
    for distance in distances:
        half_distance = 0.5 * distance
        trajectory_values.extend(
            (
                center - half_distance,
                center,
                center,
                center + half_distance,
                center,
                center,
            )
        )
    (case_dir / "trajectory.dat").write_bytes(
        struct.pack(f"={len(trajectory_values)}f", *trajectory_values)
    )
    (case_dir / "box.txt").write_text(
        "".join(
            f"{BOX_LENGTH} {BOX_LENGTH} {BOX_LENGTH} 90 90 90\n"
            for _ in distances
        ),
        encoding="utf-8",
    )

    settings = {
        "md_name": case_dir.name,
        "mode": "rerun",
        "cutoff": 4.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "LJ_in_file": "lj.txt",
        "crd": "trajectory.dat",
        "box": "box.txt",
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
        "SITS.mode": "amd",
        "SITS.atom_numbers": "ALL",
        "SITS.pe_a": AMD_ALPHA,
        "SITS.pe_b": AMD_THRESHOLD,
        "SITS.fb_interval": 2,
    }
    settings.update(sits_overrides or {})
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"SPONGE failed with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )

    mdout_lines = (
        (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    )
    names = mdout_lines[0].split()
    mdout = [
        {name: float(value) for name, value in zip(names, line.split())}
        for line in mdout_lines[1:]
    ]
    raw_forces = (case_dir / "frc.dat").read_bytes()
    force_values = struct.unpack(f"={len(raw_forces) // 4}f", raw_forces)
    forces = [
        (
            force_values[offset : offset + 3],
            force_values[offset + 3 : offset + 6],
        )
        for offset in range(0, len(force_values), 6)
    ]
    return mdout, forces


def _c_rand_pair():
    libc = ctypes.CDLL(None)
    libc.srand(0)
    rand_max = 2_147_483_647
    proposal = 2.0 * libc.rand() / rand_max - 1.0
    acceptance_draw = libc.rand() / rand_max
    return proposal, acceptance_draw


def _run_mc_state(
    case_dir,
    target_pressure,
    initial_ratio,
    *,
    amd_alpha=0.01,
    amd_threshold=5.0,
    settings_overrides=None,
    extra_files=None,
):
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", (12.0, 12.0))
    _write_counted_values(case_dir / "charge.txt", (0.0, 0.0))
    _write_coordinates(case_dir / "coordinate.txt", PAIR_DISTANCE)
    (case_dir / "velocity.txt").write_text(
        "2\n0 0 0\n0 0 0\n", encoding="utf-8"
    )
    (case_dir / "lj.txt").write_text(
        f"2 1\n{REPULSIVE_COEFFICIENT:.12g}\n0\n0\n0\n",
        encoding="utf-8",
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "npt",
        "step_limit": 0,
        "dt": 0.001,
        "cutoff": 4.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "velocity_in_file": "velocity.txt",
        "LJ_in_file": "lj.txt",
        "thermostat": "berendsen_thermostat",
        "thermostat_tau": 1.0,
        "barostat": "monte_carlo_barostat",
        "target_temperature": 300.0,
        "target_pressure": target_pressure,
        "monte_carlo_barostat_initial_ratio": initial_ratio,
        "monte_carlo_barostat_update_interval": 1,
        "monte_carlo_barostat_check_interval": 100,
        "monte_carlo_barostat_couple_dimension": "XYZ",
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
        "SITS.mode": "amd",
        "SITS.atom_numbers": "ALL",
        "SITS.pe_a": amd_alpha,
        "SITS.pe_b": amd_threshold,
    }
    settings.update(settings_overrides or {})
    for name, contents in (extra_files or {}).items():
        (case_dir / name).write_text(contents, encoding="utf-8")
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"SPONGE failed with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    mdout_lines = (
        (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    )
    mdout = {
        name: float(value)
        for name, value in zip(mdout_lines[0].split(), mdout_lines[-1].split())
    }
    raw_forces = (case_dir / "frc.dat").read_bytes()
    force_values = struct.unpack(f"={len(raw_forces) // 4}f", raw_forces)[-6:]
    return mdout, (force_values[:3], force_values[3:])


def _run_coulomb_state(case_dir, dispatch):
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", (12.0, 12.0))
    _write_counted_values(case_dir / "charge.txt", (0.2, -0.3))
    _write_coordinates(case_dir / "coordinate.txt", PAIR_DISTANCE)
    settings = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 0,
        "dt": 0,
        "cutoff": 4.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
    }
    if dispatch in ("standard", "sits-hard"):
        (case_dir / "lj.txt").write_text("2 1\n0\n0\n0\n0\n", encoding="utf-8")
        settings["LJ_in_file"] = "lj.txt"
    else:
        (case_dir / "lj_soft.txt").write_text(
            "2 1 1\n0\n0\n0\n0\n0 0\n0 0\n", encoding="utf-8"
        )
        settings.update(
            {
                "LJ_soft_core_in_file": "lj_soft.txt",
                "lambda_lj": SOFT_CORE_LAMBDA,
            }
        )
    if dispatch != "standard":
        settings.update({"SITS.mode": "observation", "SITS.atom_numbers": 1})
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    if result.returncode != 0:
        raise AssertionError(
            f"SPONGE failed with code {result.returncode}\n"
            f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
        )
    mdout_lines = (
        (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    )
    mdout = {
        name: float(value)
        for name, value in zip(mdout_lines[0].split(), mdout_lines[-1].split())
    }
    raw_forces = (case_dir / "frc.dat").read_bytes()
    forces = struct.unpack(f"={len(raw_forces) // 4}f", raw_forces)[-6:]
    return mdout, forces


def _pair_energy(distance, potential_kind):
    if potential_kind == "hard":
        return REPULSIVE_COEFFICIENT / distance**12
    softened_r6 = (
        distance**6 + SOFT_CORE_ALPHA * SOFT_CORE_LAMBDA * SOFT_CORE_SIGMA**6
    )
    return (1.0 - SOFT_CORE_LAMBDA) * REPULSIVE_COEFFICIENT / softened_r6**2


def _amd_bias(enhancing_energy):
    if enhancing_energy >= AMD_THRESHOLD:
        return 0.0
    delta = AMD_THRESHOLD - enhancing_energy
    return delta * delta / (AMD_ALPHA + delta)


def _amd_factor(enhancing_energy):
    if enhancing_energy >= AMD_THRESHOLD:
        return 1.0
    delta = AMD_THRESHOLD - enhancing_energy
    return 1.0 - delta * (2.0 * AMD_ALPHA + delta) / ((AMD_ALPHA + delta) ** 2)


def _gamd_bias(enhancing_energy):
    if enhancing_energy >= GAMD_THRESHOLD:
        return 0.0
    delta = GAMD_THRESHOLD - enhancing_energy
    return 0.5 * GAMD_K * delta**2


def _gamd_factor(enhancing_energy):
    if enhancing_energy >= GAMD_THRESHOLD:
        return 1.0
    return 1.0 - GAMD_K * (GAMD_THRESHOLD - enhancing_energy)


def _empirical_bias(enhancing_energy):
    displacement = enhancing_energy - EMPIRICAL_THRESHOLD
    low_temperature_factor = REFERENCE_TEMPERATURE / EMPIRICAL_T_LOW
    high_temperature_factor = REFERENCE_TEMPERATURE / EMPIRICAL_T_HIGH
    if displacement <= 0.0:
        return (high_temperature_factor - 1.0) * displacement
    return (low_temperature_factor - 1.0) * displacement - (
        low_temperature_factor - high_temperature_factor
    ) * EMPIRICAL_SMOOTHING_ENERGY * math.log1p(
        displacement / EMPIRICAL_SMOOTHING_ENERGY
    )


def _empirical_factor(enhancing_energy):
    displacement = enhancing_energy - EMPIRICAL_THRESHOLD
    low_temperature_factor = REFERENCE_TEMPERATURE / EMPIRICAL_T_LOW
    high_temperature_factor = REFERENCE_TEMPERATURE / EMPIRICAL_T_HIGH
    if displacement <= 0.0:
        return high_temperature_factor
    return low_temperature_factor - (
        (low_temperature_factor - high_temperature_factor)
        * EMPIRICAL_SMOOTHING_ENERGY
        / (displacement + EMPIRICAL_SMOOTHING_ENERGY)
    )


def _mode_overrides(mode, fb_interval=1):
    if mode == "gamd":
        return {
            "SITS.mode": mode,
            "SITS.pe_a": GAMD_K,
            "SITS.pe_b": GAMD_THRESHOLD,
            "SITS.fb_interval": fb_interval,
        }
    if mode == "empirical":
        return {
            "SITS.mode": mode,
            "SITS.pe_a": EMPIRICAL_SMOOTHING_ENERGY,
            "SITS.pe_b": EMPIRICAL_THRESHOLD,
            "SITS.T_low": EMPIRICAL_T_LOW,
            "SITS.T_high": EMPIRICAL_T_HIGH,
            "SITS.fb_interval": fb_interval,
        }
    raise AssertionError(f"unsupported test mode: {mode}")


def _hamiltonian(distance, potential_kind):
    pair_energy = _pair_energy(distance, potential_kind)
    return pair_energy + _amd_bias(CROSS_FACTOR * pair_energy)


def _reported_hamiltonian(mdout):
    return mdout["eff_pot"]


@pytest.mark.parametrize("potential_kind", ("hard", "soft"))
def test_sits_pair_force_and_virial_are_hamiltonian_derivatives(
    tmp_path, potential_kind
):
    coordinate_delta = 1.0e-3
    log_scale_delta = 1.0e-3
    base_mdout, base_forces = _run_state(
        tmp_path / f"{potential_kind}_base", PAIR_DISTANCE, potential_kind
    )
    coordinate_plus_mdout, _ = _run_state(
        tmp_path / f"{potential_kind}_coordinate_plus",
        PAIR_DISTANCE + coordinate_delta,
        potential_kind,
    )
    coordinate_minus_mdout, _ = _run_state(
        tmp_path / f"{potential_kind}_coordinate_minus",
        PAIR_DISTANCE - coordinate_delta,
        potential_kind,
    )
    scale_plus_mdout, _ = _run_state(
        tmp_path / f"{potential_kind}_scale_plus",
        PAIR_DISTANCE * math.exp(log_scale_delta),
        potential_kind,
    )
    scale_minus_mdout, _ = _run_state(
        tmp_path / f"{potential_kind}_scale_minus",
        PAIR_DISTANCE * math.exp(-log_scale_delta),
        potential_kind,
    )

    # The two owned endpoints of a local-local pair must receive exact
    # opposite enhanced forces.  The old hard kernel only wrote endpoint i.
    for component in range(3):
        assert base_forces[0][component] + base_forces[1][
            component
        ] == pytest.approx(0.0, abs=2.0e-5)
    assert base_forces[0][1:] == pytest.approx((0.0, 0.0), abs=1.0e-7)
    assert base_forces[1][1:] == pytest.approx((0.0, 0.0), abs=1.0e-7)

    force_from_energy = (
        _reported_hamiltonian(coordinate_plus_mdout)
        - _reported_hamiltonian(coordinate_minus_mdout)
    ) / (2.0 * coordinate_delta)
    assert base_forces[0][0] == pytest.approx(
        force_from_energy, rel=3.0e-4, abs=2.0e-4
    )
    assert base_forces[1][0] == pytest.approx(
        -force_from_energy, rel=3.0e-4, abs=2.0e-4
    )

    d_hamiltonian_d_log_scale = (
        _reported_hamiltonian(scale_plus_mdout)
        - _reported_hamiltonian(scale_minus_mdout)
    ) / (2.0 * log_scale_delta)
    expected_pxx = (
        -d_hamiltonian_d_log_scale / BOX_LENGTH**3 * PRESSURE_CONVERSION
    )
    assert base_mdout["Pxx"] == pytest.approx(expected_pxx, abs=0.06)
    assert base_mdout["Pyy"] == pytest.approx(0.0, abs=0.01)
    assert base_mdout["Pzz"] == pytest.approx(0.0, abs=0.01)

    pair_energy = _pair_energy(PAIR_DISTANCE, potential_kind)
    enhancing_energy = CROSS_FACTOR * pair_energy
    # The physical-potential diagnostic is intentionally printed with two
    # decimals; eff_pot below carries enough precision for derivative checks.
    assert base_mdout["potential"] == pytest.approx(pair_energy, abs=5.1e-3)
    assert base_mdout["eff_pot"] == pytest.approx(
        _hamiltonian(PAIR_DISTANCE, potential_kind), abs=3.0e-6
    )
    assert _reported_hamiltonian(base_mdout) == pytest.approx(
        _hamiltonian(PAIR_DISTANCE, potential_kind), abs=3.0e-6
    )
    assert base_mdout["SITS_AA_kAB"] == pytest.approx(
        enhancing_energy, abs=0.011
    )
    assert base_mdout["SITS_bias"] == pytest.approx(
        _amd_bias(enhancing_energy), abs=1.1e-4
    )
    assert base_mdout["SITS_fb"] == pytest.approx(
        _amd_factor(enhancing_energy), abs=1.1e-4
    )


@pytest.mark.parametrize("potential_kind", ("hard", "soft"))
def test_sits_all_mode_uses_current_reduced_potential(tmp_path, potential_kind):
    coordinate_delta = 1.0e-3
    distance = 2.1
    base_mdout, base_forces = _run_state(
        tmp_path / f"{potential_kind}_all_base",
        distance,
        potential_kind,
        "ALL",
    )
    plus_mdout, _ = _run_state(
        tmp_path / f"{potential_kind}_all_plus",
        distance + coordinate_delta,
        potential_kind,
        "ALL",
    )
    minus_mdout, _ = _run_state(
        tmp_path / f"{potential_kind}_all_minus",
        distance - coordinate_delta,
        potential_kind,
        "ALL",
    )

    pair_energy = _pair_energy(distance, potential_kind)
    expected_effective = pair_energy + _amd_bias(pair_energy)
    assert base_mdout["potential"] == pytest.approx(pair_energy, abs=5.1e-3)
    assert base_mdout["eff_pot"] == pytest.approx(
        expected_effective, abs=3.0e-6
    )
    assert base_mdout["SITS_AA_kAB"] == pytest.approx(pair_energy, abs=0.011)
    assert base_mdout["SITS_bias"] == pytest.approx(
        _amd_bias(pair_energy), abs=1.1e-4
    )
    assert base_mdout["SITS_fb"] == pytest.approx(
        _amd_factor(pair_energy), abs=1.1e-4
    )

    force_from_energy = (plus_mdout["eff_pot"] - minus_mdout["eff_pot"]) / (
        2.0 * coordinate_delta
    )
    assert base_forces[0][0] == pytest.approx(
        force_from_energy, rel=5.0e-4, abs=2.0e-4
    )
    assert base_forces[1][0] == pytest.approx(
        -force_from_energy, rel=5.0e-4, abs=2.0e-4
    )


@pytest.mark.parametrize(
    "mode,distance,bias_function,factor_function",
    (
        ("gamd", 2.05, _gamd_bias, _gamd_factor),
        ("empirical", 1.95, _empirical_bias, _empirical_factor),
        ("empirical", 2.05, _empirical_bias, _empirical_factor),
    ),
    ids=(
        "gamd-below-threshold",
        "empirical-above-threshold",
        "empirical-below-threshold",
    ),
)
def test_sits_exact_gamd_and_empirical_forces_are_bias_derivatives(
    tmp_path, mode, distance, bias_function, factor_function
):
    coordinate_delta = 1.0e-3
    overrides = _mode_overrides(mode)
    base_mdout, base_forces = _run_state(
        tmp_path / f"{mode}_exact_base",
        distance,
        "hard",
        "ALL",
        overrides,
    )
    plus_mdout, _ = _run_state(
        tmp_path / f"{mode}_exact_plus",
        distance + coordinate_delta,
        "hard",
        "ALL",
        overrides,
    )
    minus_mdout, _ = _run_state(
        tmp_path / f"{mode}_exact_minus",
        distance - coordinate_delta,
        "hard",
        "ALL",
        overrides,
    )

    energy = _pair_energy(distance, "hard")
    expected_bias = bias_function(energy)
    expected_factor = factor_function(energy)
    assert base_mdout["SITS_bias"] == pytest.approx(expected_bias, abs=1.1e-4)
    assert base_mdout["SITS_fb"] == pytest.approx(expected_factor, abs=1.1e-4)
    assert base_mdout["eff_pot"] == pytest.approx(
        energy + expected_bias, abs=4.0e-6
    )

    force_from_energy = (plus_mdout["eff_pot"] - minus_mdout["eff_pot"]) / (
        2.0 * coordinate_delta
    )
    physical_force = -12.0 * REPULSIVE_COEFFICIENT / distance**13
    assert base_forces[0][0] == pytest.approx(
        force_from_energy, rel=6.0e-4, abs=4.0e-4
    )
    assert base_forces[1][0] == pytest.approx(
        -force_from_energy, rel=6.0e-4, abs=4.0e-4
    )
    assert base_forces[0][0] == pytest.approx(
        expected_factor * physical_force, rel=5.0e-4, abs=3.0e-4
    )


@pytest.mark.parametrize(
    "mode,threshold,bias_function,factor_function",
    (
        ("gamd", GAMD_THRESHOLD, _gamd_bias, _gamd_factor),
        ("empirical", EMPIRICAL_THRESHOLD, _empirical_bias, _empirical_factor),
    ),
)
def test_sits_gamd_and_empirical_are_continuous_across_threshold(
    tmp_path, mode, threshold, bias_function, factor_function
):
    energy_offset = 0.02
    energies = (threshold - energy_offset, threshold, threshold + energy_offset)
    distances = tuple(
        (REPULSIVE_COEFFICIENT / energy) ** (1.0 / 12.0) for energy in energies
    )
    mdout, _ = _run_rerun_states(
        tmp_path / f"{mode}_threshold",
        distances,
        _mode_overrides(mode),
    )

    for frame, energy in zip(mdout, energies):
        assert frame["SITS_bias"] == pytest.approx(
            bias_function(energy), abs=1.1e-4
        )
        assert frame["SITS_fb"] == pytest.approx(
            factor_function(energy), abs=1.1e-4
        )
        assert frame["eff_pot"] == pytest.approx(
            energy + bias_function(energy), abs=5.0e-6
        )

    assert mdout[1]["SITS_bias"] == pytest.approx(0.0, abs=1.0e-7)
    assert abs(mdout[0]["SITS_bias"] - mdout[1]["SITS_bias"]) < 0.011
    assert abs(mdout[2]["SITS_bias"] - mdout[1]["SITS_bias"]) < 0.011


@pytest.mark.parametrize(
    "mode,distances,bias_function,factor_function",
    (
        ("gamd", (2.05, 1.95, 2.0), _gamd_bias, _gamd_factor),
        ("empirical", (1.95, 2.05, 2.0), _empirical_bias, _empirical_factor),
    ),
)
def test_sits_gamd_and_empirical_fb_interval_uses_tangent_hamiltonian(
    tmp_path, mode, distances, bias_function, factor_function
):
    mdout, forces = _run_rerun_states(
        tmp_path / f"{mode}_fb_interval",
        distances,
        _mode_overrides(mode, fb_interval=2),
    )
    reference_energy = _pair_energy(distances[0], "hard")
    reference_bias = bias_function(reference_energy)
    reference_factor = factor_function(reference_energy)
    middle_energy = _pair_energy(distances[1], "hard")
    middle_bias = reference_bias + (reference_factor - 1.0) * (
        middle_energy - reference_energy
    )

    assert mdout[1]["SITS_bias"] == pytest.approx(middle_bias, abs=1.1e-4)
    assert mdout[1]["SITS_fb"] == pytest.approx(reference_factor, abs=1.1e-4)
    assert mdout[1]["eff_pot"] == pytest.approx(
        middle_energy + middle_bias, abs=5.0e-6
    )
    physical_force = -12.0 * REPULSIVE_COEFFICIENT / distances[1] ** 13
    assert forces[1][0][0] == pytest.approx(
        reference_factor * physical_force, rel=5.0e-4, abs=4.0e-4
    )
    assert forces[1][1][0] == pytest.approx(
        -reference_factor * physical_force, rel=5.0e-4, abs=4.0e-4
    )

    final_energy = _pair_energy(distances[2], "hard")
    assert mdout[2]["SITS_bias"] == pytest.approx(
        bias_function(final_energy), abs=1.1e-4
    )
    assert mdout[2]["SITS_fb"] == pytest.approx(
        factor_function(final_energy), abs=1.1e-4
    )
    assert mdout[2]["eff_pot"] == pytest.approx(
        final_energy + bias_function(final_energy), abs=5.0e-6
    )


def test_sits_fb_interval_uses_one_consistent_tangent_hamiltonian(tmp_path):
    distances = (2.1, 1.95, 2.2)
    mdout, forces = _run_rerun_states(tmp_path / "fb_interval", distances)
    assert len(mdout) == len(distances)
    assert len(forces) == len(distances)

    reference_energy = _pair_energy(distances[0], "hard")
    reference_bias = _amd_bias(reference_energy)
    reference_factor = _amd_factor(reference_energy)
    middle_energy = _pair_energy(distances[1], "hard")
    middle_bias = reference_bias + (reference_factor - 1.0) * (
        middle_energy - reference_energy
    )

    # Step 1 is deliberately between feedback updates.  Its cached force
    # factor and reported bias must describe the same tangent Hamiltonian.
    assert mdout[1]["potential"] == pytest.approx(middle_energy, abs=5.1e-3)
    assert mdout[1]["SITS_bias"] == pytest.approx(middle_bias, abs=1.1e-4)
    assert mdout[1]["SITS_fb"] == pytest.approx(reference_factor, abs=1.1e-4)
    assert mdout[1]["eff_pot"] == pytest.approx(
        middle_energy + middle_bias, abs=4.0e-6
    )
    physical_force = -12.0 * REPULSIVE_COEFFICIENT / distances[1] ** 13
    assert forces[1][0][0] == pytest.approx(
        reference_factor * physical_force, rel=4.0e-4, abs=3.0e-4
    )
    assert forces[1][1][0] == pytest.approx(
        -reference_factor * physical_force, rel=4.0e-4, abs=3.0e-4
    )

    # Step 2 is an exact feedback-update step again, not another tangent
    # extrapolation from the stale step-0 state.
    final_energy = _pair_energy(distances[2], "hard")
    assert mdout[2]["SITS_bias"] == pytest.approx(
        _amd_bias(final_energy), abs=1.1e-4
    )
    assert mdout[2]["SITS_fb"] == pytest.approx(
        _amd_factor(final_energy), abs=1.1e-4
    )
    assert mdout[2]["eff_pot"] == pytest.approx(
        final_energy + _amd_bias(final_energy), abs=4.0e-6
    )


def test_mc_barostat_acceptance_uses_the_sits_effective_hamiltonian(tmp_path):
    initial_ratio = 0.05
    amd_alpha = 0.01
    amd_threshold = 5.0
    proposal, acceptance_draw = _c_rand_pair()
    scale = 1.0 + proposal * initial_ratio
    assert scale > 0.0
    assert 0.0 < acceptance_draw < 1.0

    old_energy = _pair_energy(PAIR_DISTANCE, "hard")
    new_energy = _pair_energy(PAIR_DISTANCE * scale, "hard")

    def bias(energy):
        delta = amd_threshold - energy
        return 0.0 if delta <= 0.0 else delta * delta / (amd_alpha + delta)

    physical_delta = new_energy - old_energy
    effective_delta = (
        new_energy + bias(new_energy) - old_energy - bias(old_energy)
    )
    # Tune P*DeltaV - NkT*log(V'/V) so the one deterministic MC draw lies
    # strictly between the physical- and effective-Hamiltonian decisions.
    temperature = 300.0
    thermal_energy = 0.00198716 * temperature
    volume_ratio = scale**3
    delta_volume = BOX_LENGTH**3 * (volume_ratio - 1.0)
    decision_threshold = -thermal_energy * math.log(acceptance_draw)
    desired_extra = decision_threshold - 0.5 * (
        physical_delta + effective_delta
    )
    molecule_count = 2
    target_pressure = (
        (
            desired_extra
            + molecule_count * thermal_energy * math.log(volume_ratio)
        )
        / delta_volume
        * PRESSURE_CONVERSION
    )

    physical_accept = physical_delta + desired_extra < decision_threshold
    effective_accept = effective_delta + desired_extra < decision_threshold
    assert physical_accept != effective_accept

    mdout, _ = _run_mc_state(
        tmp_path / "mc_effective_acceptance",
        target_pressure,
        initial_ratio,
        amd_alpha=amd_alpha,
        amd_threshold=amd_threshold,
    )
    old_density = 24.0 * 1.0e24 / 6.023e23 / BOX_LENGTH**3
    expected_density = (
        old_density / volume_ratio if effective_accept else old_density
    )
    assert mdout["density"] == pytest.approx(expected_density, abs=5.1e-5)
    expected_energy = new_energy if effective_accept else old_energy
    assert mdout["potential"] == pytest.approx(expected_energy, abs=5.1e-3)
    assert mdout["eff_pot"] == pytest.approx(
        expected_energy + bias(expected_energy), abs=5.0e-6
    )


def test_rejected_mc_trial_restores_all_sits_and_energy_diagnostics(tmp_path):
    initial_ratio = 0.05
    amd_alpha = 0.01
    amd_threshold = 5.0
    proposal, _ = _c_rand_pair()
    scale = 1.0 + proposal * initial_ratio
    delta_volume = BOX_LENGTH**3 * (scale**3 - 1.0)
    # Make P*DeltaV overwhelmingly positive on either libc rand sequence, so
    # this proposal is rejected independently of the pair-energy change.
    rejecting_pressure = math.copysign(1.0e8, delta_volume)
    mdout, forces = _run_mc_state(
        tmp_path / "mc_reject_restore",
        rejecting_pressure,
        initial_ratio,
        amd_alpha=amd_alpha,
        amd_threshold=amd_threshold,
    )

    old_energy = _pair_energy(PAIR_DISTANCE, "hard")
    old_delta = amd_threshold - old_energy
    old_bias = old_delta * old_delta / (amd_alpha + old_delta)
    old_factor = (
        1.0
        - old_delta
        * (2.0 * amd_alpha + old_delta)
        / (amd_alpha + old_delta) ** 2
    )
    old_density = 24.0 * 1.0e24 / 6.023e23 / BOX_LENGTH**3
    assert mdout["density"] == pytest.approx(old_density, abs=1.0e-4)
    assert mdout["potential"] == pytest.approx(old_energy, abs=5.1e-3)
    assert mdout["LJ"] == pytest.approx(old_energy, abs=5.1e-3)
    assert mdout["SITS_AA_kAB"] == pytest.approx(old_energy, abs=0.011)
    assert mdout["SITS_bias"] == pytest.approx(old_bias, abs=1.1e-4)
    assert mdout["SITS_fb"] == pytest.approx(old_factor, abs=1.1e-4)
    assert mdout["eff_pot"] == pytest.approx(old_energy + old_bias, abs=5.0e-6)

    physical_force = -12.0 * REPULSIVE_COEFFICIENT / PAIR_DISTANCE**13
    assert forces[0][0] == pytest.approx(
        old_factor * physical_force, rel=3.0e-3, abs=3.0e-5
    )
    assert forces[1][0] == pytest.approx(
        -old_factor * physical_force, rel=3.0e-3, abs=3.0e-5
    )


def test_sits_nk_persistence_replaces_restart_after_complete_record(tmp_path):
    case_dir = tmp_path / "nk_persistence_success"
    _write_case(
        case_dir,
        PAIR_DISTANCE,
        "hard",
        sits_overrides=_sits_iteration_persistence_overrides(),
    )
    restart_path = case_dir / "nk_rest.txt"
    restart_path.write_text("previous complete restart\n", encoding="utf-8")

    result = _run_sponge_raw(case_dir)
    output = result.stdout + result.stderr
    assert result.returncode == 0, output

    trajectory_values = struct.unpack(
        "=2f", (case_dir / "nk_traj.dat").read_bytes()
    )
    restart_values = [
        float(value)
        for value in restart_path.read_text(encoding="utf-8").split()
    ]
    assert all(
        math.isfinite(value) and value > 0.0 for value in trajectory_values
    )
    assert restart_values == pytest.approx(trajectory_values, rel=2.0e-6)
    assert not list(case_dir.glob("nk_rest.txt.sponge-tmp.*"))


def test_sits_nk_restart_replace_failure_preserves_existing_target(tmp_path):
    case_dir = tmp_path / "nk_restart_replace_failure"
    _write_case(
        case_dir,
        PAIR_DISTANCE,
        "hard",
        sits_overrides=_sits_iteration_persistence_overrides(
            restart_name="nk_restart_target"
        ),
    )
    restart_target = case_dir / "nk_restart_target"
    restart_target.mkdir()
    sentinel = restart_target / "previous-restart-marker"
    sentinel.write_text("preserve me\n", encoding="utf-8")

    result = _run_sponge_raw(case_dir)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "failed to atomically replace SITS restart" in output
    assert "the previous restart was preserved" in output
    assert sentinel.read_text(encoding="utf-8") == "preserve me\n"
    assert not list(case_dir.glob("nk_restart_target.sponge-tmp.*"))
    # The trajectory record precedes the restart transaction and must already
    # have been written and flushed completely.
    assert len((case_dir / "nk_traj.dat").read_bytes()) == 2 * 4


@pytest.mark.skipif(not hasattr(os, "mkfifo"), reason="requires a POSIX FIFO")
def test_sits_nk_trajectory_flush_failure_is_fatal(tmp_path):
    case_dir = tmp_path / "nk_trajectory_flush_failure"
    _write_case(
        case_dir,
        PAIR_DISTANCE,
        "hard",
        sits_overrides={
            **_sits_iteration_persistence_overrides(
                trajectory_name="nk_traj.pipe"
            ),
            # Delay persistence until the second sample so the FIFO reader is
            # guaranteed to have closed after unblocking fopen during setup.
            "step_limit": 1,
            "SITS.update_interval": 2,
        },
    )
    trajectory_fifo = case_dir / "nk_traj.pipe"
    os.mkfifo(trajectory_fifo)
    restart_path = case_dir / "nk_rest.txt"
    restart_path.write_text("previous complete restart\n", encoding="utf-8")

    reader_errors = []

    def open_and_close_fifo_reader():
        try:
            descriptor = os.open(trajectory_fifo, os.O_RDONLY)
            os.close(descriptor)
        except OSError as error:
            reader_errors.append(error)

    reader = threading.Thread(target=open_and_close_fifo_reader, daemon=True)
    reader.start()
    # Python ignores SIGPIPE.  Preserve that disposition in SPONGE so the
    # stdio call returns EPIPE and exercises the checked fatal path.
    result = _run_sponge_raw(case_dir, restore_signals=False)
    reader.join(timeout=5)
    output = result.stdout + result.stderr

    assert not reader.is_alive()
    assert not reader_errors
    assert result.returncode != 0, output
    assert "failed to flush the complete SITS Nk trajectory record" in output
    assert "nk_traj.pipe" in output
    assert restart_path.read_text(encoding="utf-8") == (
        "previous complete restart\n"
    )
    assert not list(case_dir.glob("nk_rest.txt.sponge-tmp.*"))


def test_mc_trial_does_not_record_sits_statistics_twice(tmp_path):
    case_dir = tmp_path / "mc_single_sits_sample"
    _run_mc_state(
        case_dir,
        target_pressure=1.0,
        initial_ratio=0.0,
        settings_overrides=_sits_iteration_persistence_overrides(),
    )
    # One dynamics force evaluation records exactly one two-float Nk frame.
    # The accepted zero-volume MC trial is another force evaluation at the
    # same step, but must not be counted as a second statistical sample.
    assert (case_dir / "nk_traj.dat").stat().st_size == 2 * 4
    nk_values = [
        float(value)
        for value in (case_dir / "nk_rest.txt")
        .read_text(encoding="utf-8")
        .split()
    ]
    assert len(nk_values) == 2
    assert all(math.isfinite(value) and value > 0.0 for value in nk_values)


@pytest.mark.parametrize(
    "name,initial_ratio,overrides,expected_error",
    (
        (
            "invalid_ratio",
            1.0,
            {},
            "initial_ratio must be finite",
        ),
        (
            "zero_update_interval",
            0.01,
            {"monte_carlo_barostat_update_interval": 0},
            "update_interval must be positive",
        ),
        (
            "zero_check_interval",
            0.01,
            {"monte_carlo_barostat_check_interval": 0},
            "check_interval must be positive",
        ),
        (
            "unknown_coupling",
            0.01,
            {"monte_carlo_barostat_couple_dimension": "invalid"},
            "couple_dimension must be one of",
        ),
    ),
)
def test_mc_barostat_rejects_invalid_transaction_inputs(
    tmp_path, name, initial_ratio, overrides, expected_error
):
    with pytest.raises(AssertionError, match=expected_error):
        _run_mc_state(
            tmp_path / name,
            target_pressure=1.0,
            initial_ratio=initial_ratio,
            settings_overrides=overrides,
        )


def test_sits_selective_dispatch_preserves_direct_coulomb_energy(tmp_path):
    standard_mdout, standard_forces = _run_coulomb_state(
        tmp_path / "standard", "standard"
    )
    assert abs(standard_mdout["eff_pot"]) > 1.0e-3
    assert max(abs(value) for value in standard_forces) > 1.0e-3

    for dispatch in ("sits-hard", "sits-soft"):
        mdout, forces = _run_coulomb_state(tmp_path / dispatch, dispatch)
        assert mdout["eff_pot"] == pytest.approx(
            standard_mdout["eff_pot"], abs=2.0e-6
        )
        assert forces == pytest.approx(standard_forces, abs=2.0e-5)


@pytest.mark.parametrize(
    "name,overrides,extra_files,expected_error",
    (
        (
            "zero_fb_interval",
            {"SITS.fb_interval": 0},
            {},
            "SITS_fb_interval must be positive",
        ),
        (
            "zero_reference_temperature",
            {"target_temperature": 0.0},
            {},
            "target_temperature must be finite, positive, and normal",
        ),
        (
            "subnormal_reference_temperature",
            {"target_temperature": 1.0e-45},
            {},
            "target_temperature must be finite, positive, and normal",
        ),
        (
            "zero_amd_alpha",
            {"SITS.pe_a": 0.0},
            {},
            "SITS_pe_a must be positive and finite",
        ),
        (
            "negative_gamd_force_factor",
            {
                "SITS.mode": "gamd",
                "SITS.atom_numbers": "ALL",
                "SITS.pe_a": 2.0,
                "SITS.pe_b": 2.0,
            },
            {},
            "GaMD force factor is negative",
        ),
        (
            "conflicting_selection",
            {"SITS.atom_in_file": "selected.txt"},
            {"selected.txt": "0\n"},
            "exactly one of SITS_atom_in_file and SITS_atom_numbers",
        ),
        (
            "selection_out_of_range",
            {"SITS.atom_numbers": None, "SITS.atom_in_file": "selected.txt"},
            {"selected.txt": "2\n"},
            "contains atom index 2 outside [0, 2)",
        ),
        (
            "selection_bad_token",
            {"SITS.atom_numbers": None, "SITS.atom_in_file": "selected.txt"},
            {"selected.txt": "not-an-index\n"},
            "contains a non-integer token",
        ),
        (
            "one_temperature",
            {
                "SITS.mode": "iteration",
                "SITS.k_numbers": 1,
                "SITS.T": "300",
            },
            {},
            "SITS_k_numbers must be at least 2",
        ),
        (
            "bad_temperature_list",
            {
                "SITS.mode": "iteration",
                "SITS.k_numbers": 2,
                "SITS.T": "300/not-a-temperature",
            },
            {},
            "SITS_T entry 1 must be a positive finite temperature",
        ),
        (
            "zero_record_interval",
            {
                "SITS.mode": "iteration",
                "SITS.k_numbers": 2,
                "SITS.T": "300/400",
                "SITS.record_interval": 0,
            },
            {},
            "SITS_record_interval must be positive",
        ),
        (
            "zero_update_interval",
            {
                "SITS.mode": "iteration",
                "SITS.k_numbers": 2,
                "SITS.T": "300/400",
                "SITS.update_interval": 0,
            },
            {},
            "SITS_update_interval must be positive",
        ),
        (
            "truncated_nk",
            {
                "SITS.mode": "production",
                "SITS.k_numbers": 2,
                "SITS.T": "300/400",
                "SITS.nk_in_file": "nk.txt",
            },
            {"nk.txt": "1\n"},
            "SITS_nk_in_file entry 1 must be a positive finite number",
        ),
    ),
)
def test_sits_rejects_invalid_input_contracts(
    tmp_path, name, overrides, extra_files, expected_error
):
    case_dir = tmp_path / name
    _write_case(
        case_dir,
        PAIR_DISTANCE,
        "hard",
        sits_overrides=overrides,
        extra_files=extra_files,
    )
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode != 0
    assert expected_error in result.stdout + result.stderr
