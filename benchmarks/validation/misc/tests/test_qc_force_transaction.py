import json
import math
import os
import re
import struct
import subprocess

import pytest

SCF_ITERATION = re.compile(
    r"Step\s+(-?\d+)\s+\|\s+SCF Iter\s+(\d+)\s+\|\s+"
    r"E\(Ha\)=([-+0-9.eE]+)"
)
SCF_CONVERGENCE = re.compile(
    r"SCF Iter\s+(\d+).*?dE\(Ha\)=([-+0-9.eE]+).*?"
    r"dP\(rms\)=([-+0-9.eE]+).*?stable=(\d+)/(\d+)"
)
SCF_FIXED_POINT = re.compile(
    r"SCF Iter\s+(\d+).*?dE\(Ha\)=([-+0-9.eE]+).*?"
    r"dP\(rms\)=([-+0-9.eE]+).*?map=(\w+).*?"
    r"shift=([-+0-9.eE]+).*?stable=(\d+)/(\d+)"
)
SCF_FOCK_BUILD = re.compile(r"SCF Iter\s+(\d+).*?map=(\w+).*?fock=(\w+)")
SCF_ENSEMBLE = re.compile(
    r"Step\s+(-?\d+)\s+\|\s+SCF Ensemble\s+(\d+)\s+\|\s+"
    r"E\(Ha\)=([-+0-9.eE]+)\s+\|\s+phase=([a-z-]+)\s+\|\s+"
    r"x=([-+0-9.eE]+)\s+\|\s+(?:g|global-gap)\(Ha\)="
    r"([-+0-9.eE]+)\s+\|\s+"
    r"comm=([-+0-9.eE]+)\s+\|\s+active-gap=([-+0-9.eE]+)\s+\|\s+"
    r"active=(\d+)"
)


def _write_counted_values(path, values):
    path.write_text(
        str(len(values))
        + "\n"
        + "\n".join(f"{value:.12g}" for value in values)
        + "\n",
        encoding="utf-8",
    )


def _write_qc_case(
    case_dir,
    *,
    atoms,
    coordinates,
    charge=0,
    multiplicity=1,
    mc=False,
    scf_print=False,
    max_scf_iter=None,
    step_limit=0,
    dt=None,
    velocities=None,
    print_pressure=False,
    box=(20.0, 20.0, 20.0),
    extra_settings=None,
    raw_settings=None,
):
    case_dir.mkdir()
    _write_counted_values(case_dir / "mass.txt", [mass for _, mass in atoms])
    _write_counted_values(case_dir / "charge.txt", [0.0] * len(atoms))
    coordinate_lines = [str(len(atoms))]
    coordinate_lines.extend(
        " ".join(f"{value:.12g}" for value in coordinate)
        for coordinate in coordinates
    )
    coordinate_lines.extend(
        (" ".join(f"{value:.12g}" for value in box), "90 90 90")
    )
    (case_dir / "coordinate.txt").write_text(
        "\n".join(coordinate_lines) + "\n", encoding="utf-8"
    )
    qc_lines = [f"{len(atoms)} {charge} {multiplicity}"]
    qc_lines.extend(
        f"{index} {symbol}" for index, (symbol, _) in enumerate(atoms)
    )
    (case_dir / "qc_type.txt").write_text(
        "\n".join(qc_lines) + "\n", encoding="utf-8"
    )
    if velocities is not None:
        velocity_lines = [str(len(atoms))]
        velocity_lines.extend(
            " ".join(f"{value:.12g}" for value in velocity)
            for velocity in velocities
        )
        (case_dir / "velocity.txt").write_text(
            "\n".join(velocity_lines) + "\n", encoding="utf-8"
        )

    settings = {
        "md_name": case_dir.name,
        "mode": "npt" if mc else "nve",
        "step_limit": step_limit,
        "dt": (1.0e-8 if mc else 0.0) if dt is None else dt,
        "cutoff": 6.0,
        "skin": 0.5,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "qc_type_in_file": "qc_type.txt",
        "qc_model_chemistry": "HF/sto-3g",
        "qc_scf_print_iter": 1 if scf_print else 0,
        "mdout": "mdout.txt",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
        "dont_check_input": 1,
        "print_pressure": print_pressure,
    }
    if velocities is not None:
        settings["velocity_in_file"] = "velocity.txt"
    if max_scf_iter is not None:
        settings["qc_scf_max_iter"] = max_scf_iter
    if mc:
        settings.update(
            {
                "thermostat": "berendsen_thermostat",
                "thermostat_tau": 1.0,
                "barostat": "monte_carlo_barostat",
                "target_temperature": 300.0,
                "target_pressure": 1.0,
                "monte_carlo_barostat_initial_ratio": 0.0,
                "monte_carlo_barostat_update_interval": 1,
                "monte_carlo_barostat_check_interval": 100,
                "monte_carlo_barostat_couple_dimension": "XYZ",
            }
        )
    if extra_settings is not None:
        settings.update(extra_settings)
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + (
            "\n"
            + "\n".join(
                f"{key} = {value}" for key, value in raw_settings.items()
            )
            if raw_settings
            else ""
        )
        + "\n",
        encoding="utf-8",
    )


def _run_case(case_dir, *, mpi_np=None, check=True):
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
        timeout=180,
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


def _scf_evaluations(output):
    evaluations = []
    for step, iteration, energy in SCF_ITERATION.findall(output):
        assert int(step) == 0
        iteration = int(iteration)
        if iteration == 1:
            evaluations.append([])
        assert evaluations, output
        assert iteration == len(evaluations[-1]) + 1
        evaluations[-1].append(float(energy))
    return evaluations


def _final_scf_energy_by_step(output):
    final = {}
    for step, _iteration, energy in SCF_ITERATION.findall(output):
        final[int(step)] = float(energy)
    return final


@pytest.mark.parametrize(
    "command,value,expected_message",
    [
        ("qc_scf_energy_tol", 0.0, "must be > 0"),
        ("qc_scf_energy_tol", -1.0e-6, "must be > 0"),
        ("qc_scf_energy_tol", "nan", "is not a float"),
        ("qc_scf_energy_tol", "+inf", "is not a float"),
        ("qc_scf_energy_tol", "-inf", "is not a float"),
        ("qc_scf_density_tol", 0.0, "must be > 0"),
        ("qc_scf_density_tol", -1.0e-6, "must be > 0"),
        ("qc_scf_density_tol", "nan", "is not a float"),
        ("qc_scf_density_tol", "+inf", "is not a float"),
        ("qc_scf_density_tol", "-inf", "is not a float"),
        ("qc_diis_reg", -1.0e-10, "must be finite and >= 0"),
        ("qc_diis_reg", "nan", "is not a float"),
        ("qc_diis_reg", "+inf", "is not a float"),
        ("qc_diis_reg", "-inf", "is not a float"),
        ("qc_level_shift", -1.0e-6, "must be finite and >= 0"),
        ("qc_level_shift", "nan", "is not a float"),
        ("qc_level_shift", "+inf", "is not a float"),
        ("qc_level_shift", "-inf", "is not a float"),
    ],
)
def test_scf_floating_controls_reject_invalid_values(
    tmp_path, command, value, expected_message
):
    case_dir = tmp_path / f"invalid_{command}_{str(value).replace('-', 'neg')}"
    is_nonfinite_token = isinstance(value, str)
    _write_qc_case(
        case_dir,
        atoms=[("He", 4.0)],
        coordinates=[(10.0, 10.0, 10.0)],
        extra_settings=None if is_nonfinite_token else {command: value},
        raw_settings={command: value} if is_nonfinite_token else None,
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert command in output
    assert expected_message in output


@pytest.mark.parametrize(
    "command,value",
    [
        ("qc_scf_energy_tol", 1.0e-300),
        ("qc_scf_energy_tol", 1.7976931348623157e308),
        ("qc_scf_density_tol", 1.0e-300),
        ("qc_scf_density_tol", 1.7976931348623157e308),
        ("qc_diis_reg", 0.0),
        ("qc_diis_reg", 1.0e-300),
        ("qc_diis_reg", 1.7976931348623157e308),
        ("qc_level_shift", 0.0),
    ],
)
def test_scf_floating_controls_accept_finite_boundary_values(
    tmp_path, command, value
):
    case_dir = tmp_path / f"finite_{command}_{value}"
    settings = {command: value}
    if command == "qc_diis_reg":
        settings["qc_diis"] = 0
    _write_qc_case(
        case_dir,
        atoms=[("He", 4.0)],
        coordinates=[(10.0, 10.0, 10.0)],
        scf_print=command == "qc_level_shift",
        extra_settings=settings,
    )

    result = _run_case(case_dir)
    if command == "qc_level_shift":
        shifts = [
            float(shift)
            for *_prefix, shift, _stable, _need in SCF_FIXED_POINT.findall(
                result.stdout + result.stderr
            )
        ]
        assert shifts
        assert all(shift == 0.0 for shift in shifts)


@pytest.mark.parametrize(
    "command,value",
    [
        ("qc_dft_radial_points", 0),
        ("qc_dft_radial_points", -1),
        ("qc_dft_angular_points", 0),
        ("qc_dft_angular_points", -1),
    ],
)
def test_dft_grid_counts_must_be_positive(tmp_path, command, value):
    case_dir = tmp_path / f"dft_grid_{command}_{value}"
    _write_qc_case(
        case_dir,
        atoms=[("He", 4.0)],
        coordinates=[(10.0, 10.0, 10.0)],
        extra_settings={
            "qc_model_chemistry": "LDA/sto-3g",
            command: value,
        },
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert f"{command} must be positive" in output


@pytest.mark.parametrize(
    "command",
    ["qc_dft_radial_points", "qc_dft_angular_points"],
)
def test_dft_grid_counts_reject_values_outside_int_range(tmp_path, command):
    case_dir = tmp_path / f"dft_grid_unrepresentable_{command}"
    _write_qc_case(
        case_dir,
        atoms=[("He", 4.0)],
        coordinates=[(10.0, 10.0, 10.0)],
        extra_settings={
            "qc_model_chemistry": "LDA/sto-3g",
            command: 2**31,
        },
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert f"{command} must be an exactly representable integer" in output


def test_dft_grid_product_respects_kernel_index_range(tmp_path):
    case_dir = tmp_path / "dft_grid_kernel_index_overflow"
    _write_qc_case(
        case_dir,
        atoms=[("He", 4.0)],
        coordinates=[(10.0, 10.0, 10.0)],
        extra_settings={
            "qc_model_chemistry": "LDA/sto-3g",
            "qc_dft_radial_points": 2**31 // 3 + 1,
            "qc_dft_angular_points": 1,
        },
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert (
        "DFT grid dimensions exceed the supported kernel index range" in output
    )


def test_dft_grid_accepts_positive_counts_below_old_clamp(tmp_path):
    case_dir = tmp_path / "dft_grid_small_positive"
    _write_qc_case(
        case_dir,
        atoms=[("He", 4.0)],
        coordinates=[(10.0, 10.0, 10.0)],
        extra_settings={
            "qc_model_chemistry": "LDA/sto-3g",
            "qc_dft_radial_points": 1,
            "qc_dft_angular_points": 1,
            "qc_need_gradient": 0,
        },
    )

    result = _run_case(case_dir)
    output = result.stdout + result.stderr
    assert "DFT grid: radial=1 angular=1" in output


def test_dft_grid_update_preserves_live_storage_for_scf(tmp_path):
    case_dir = tmp_path / "dft_grid_storage_lifetime"
    _write_qc_case(
        case_dir,
        atoms=[("He", 4.0)],
        coordinates=[(10.0, 10.0, 10.0)],
        extra_settings={"qc_model_chemistry": "PBE/6-31g"},
    )

    _run_case(case_dir)
    # The CPU backend aliases device pointers to host vectors.  Updating the
    # grid must republish those aliases after committing staged vector storage.
    # A dangling grid pointer drives this one-electron-density calculation to
    # energies hundreds of Hartree away from the reference value.
    assert _read_mdout(case_dir)["QC"] == pytest.approx(-1810.02, abs=3.2)


def test_mc_old_trial_and_commit_restore_one_accepted_scf_snapshot(tmp_path):
    case_dir = tmp_path / "qc_mc_transaction"
    _write_qc_case(
        case_dir,
        atoms=[("H", 1.0), ("H", 1.0)],
        coordinates=[(9.63, 10.0, 10.0), (10.37, 10.0, 10.0)],
        mc=True,
        scf_print=True,
    )
    result = _run_case(case_dir)

    evaluations = _scf_evaluations(result.stdout + result.stderr)
    assert len(evaluations) == 3
    assert len(evaluations[0]) >= 2
    assert evaluations[1] == pytest.approx(evaluations[0], abs=1.0e-12)
    assert evaluations[2] == pytest.approx(evaluations[0], abs=1.0e-12)

    mdout = _read_mdout(case_dir)
    forces = _read_forces(case_dir, 2)
    # QC reports kcal/mol. The total Hamiltonian must contain that converted
    # energy exactly once, rather than zero times, once per SCF call, or in Ha.
    assert abs(mdout["QC"]) > 100.0
    assert mdout["QC"] == pytest.approx(-700.777, abs=0.05)
    assert mdout["potential"] == pytest.approx(mdout["QC"], abs=0.02)
    assert mdout["eff_pot"] == pytest.approx(mdout["QC"], abs=5.1e-3)
    assert all(math.isfinite(value) for value in forces)


def test_qc_scf_nonconvergence_is_fatal(tmp_path):
    case_dir = tmp_path / "qc_nonconverged"
    _write_qc_case(
        case_dir,
        atoms=[("H", 1.0), ("H", 1.0)],
        coordinates=[(9.63, 10.0, 10.0), (10.37, 10.0, 10.0)],
        max_scf_iter=1,
        scf_print=True,
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    evaluations = _scf_evaluations(output)
    assert len(evaluations) == 1
    assert len(evaluations[0]) == 1
    assert "SCF failed to converge at MD step 0 within 1 iteration" in output
    assert "last finite energy" in output


@pytest.mark.parametrize(
    "fock_backend,backend_settings",
    [
        pytest.param(
            "ri-stored",
            {"qc_density_fit": 1, "qc_density_fitting_mode": "stored"},
            id="ri-stored",
        ),
        pytest.param(
            "ri-direct",
            {"qc_density_fit": 1, "qc_density_fitting_mode": "direct"},
            id="ri-direct",
        ),
        pytest.param(
            "direct-eri",
            {"qc_density_fit": 0},
            id="direct-eri",
        ),
    ],
)
def test_qc_uks_scf_confirms_unshifted_physical_fixed_point(
    tmp_path, fock_backend, backend_settings
):
    case_dir = tmp_path / f"qc_uks_physical_scf_{fock_backend}"
    _write_qc_case(
        case_dir,
        atoms=[("N", 14.0), ("O", 16.0)],
        coordinates=[
            (-0.2, 0.17, -0.094),
            (0.1714285714, 0.7271428571, 1.0042857143),
        ],
        charge=0,
        multiplicity=2,
        scf_print=True,
        extra_settings={
            "qc_model_chemistry": "LDA/def2-svp",
            "qc_restricted": 0,
            "qc_scf_energy_tol": 1.0e-6,
            "qc_scf_density_tol": 1.0e-6,
            "qc_scf_max_iter": 300,
            **backend_settings,
        },
    )
    result = _run_case(case_dir)
    output = result.stdout + result.stderr
    assert "SCF failed to converge" not in output

    records = [
        (int(iteration), float(delta_e), float(delta_p), int(stable), int(need))
        for iteration, delta_e, delta_p, stable, need in SCF_CONVERGENCE.findall(
            output
        )
    ]
    assert records

    fixed_point_records = [
        (
            int(iteration),
            float(delta_e),
            float(delta_p),
            mapping,
            float(shift),
            int(stable),
            int(need),
        )
        for iteration, delta_e, delta_p, mapping, shift, stable, need in (
            SCF_FIXED_POINT.findall(output)
        )
    ]
    assert len(fixed_point_records) == len(records)

    fock_build_records = [
        (int(iteration), mapping, mode)
        for iteration, mapping, mode in SCF_FOCK_BUILD.findall(output)
    ]
    assert len(fock_build_records) == len(records)

    ensemble_records = [
        (
            int(step),
            int(iteration),
            float(energy),
            phase,
            float(fraction),
            float(global_derivative),
            float(commutator),
            float(active_gap),
            int(active_count),
        )
        for (
            step,
            iteration,
            energy,
            phase,
            fraction,
            global_derivative,
            commutator,
            active_gap,
            active_count,
        ) in SCF_ENSEMBLE.findall(output)
    ]

    # A failed trial is never a legal handoff to ordinary SCF.  The explicit
    # recovery iteration rebuilds raw E/F from the last committed density and
    # deliberately publishes no diagonalized P_new; the very next observable
    # SCF map must therefore be unshifted physical, not accelerated.
    for recovery in (
        record
        for record in ensemble_records
        if record[3] == "recover-committed"
    ):
        subsequent_maps = [
            record for record in fixed_point_records if record[0] > recovery[1]
        ]
        assert subsequent_maps
        assert subsequent_maps[0][0] == recovery[1] + 1
        assert subsequent_maps[0][3] == "physical"
        assert subsequent_maps[0][4] == 0.0

    ensemble_certified = bool(ensemble_records) and (
        ensemble_records[-1][1] > fixed_point_records[-1][0]
    )
    if ensemble_certified:
        # A fractional solution is certified by two consecutive raw builds at
        # exactly the same committed density.  Both builds must independently
        # satisfy the global oracle, fixed-hull oracle, and commutator tests;
        # the repeat additionally proves that all four observables reproduce.
        assert len(ensemble_records) >= 2
        first, repeat = ensemble_records[-2:]
        assert first[0] == repeat[0] == 0
        assert first[3] == "verify-kkt"
        assert repeat[3] == "verify-kkt-repeat"
        for record in (first, repeat):
            assert abs(record[5]) <= 1.0e-6
            assert record[6] <= 1.0e-6
            assert record[7] <= 1.0e-6
            assert record[8] > 0
        assert abs(repeat[2] - first[2]) <= 1.0e-6
        assert abs(repeat[5] - first[5]) <= 1.0e-6
        assert abs(repeat[6] - first[6]) <= 1.0e-6
        assert abs(repeat[7] - first[7]) <= 1.0e-6

        # The raw physical map that exposed the non-idempotent solution must
        # itself be unshifted and unextrapolated.  Subsequent ensemble records
        # carry the stronger two-build certificate above.
        assert fixed_point_records[-1][3] == "physical"
        assert fixed_point_records[-1][4] == 0.0
    else:
        # For an ordinary idempotent solution, one unshifted, non-extrapolated
        # F[P] -> P' map is the physical fixed-point oracle.  Requiring a
        # second physical-to-physical energy comparison is neither stronger
        # nor part of the stopping criterion: the final density residual is
        # the certificate and P is retained with its own raw E and F.  This
        # branch also covers a failed ensemble search which transactionally
        # restored its last committed density before returning to an ordinary
        # physical map.
        final = fixed_point_records[-1]
        assert final[2] < 1.0e-6
        assert final[3] == "physical"
        assert final[4] == 0.0
        assert final[5:] == (1, 1)
        assert any(
            record[3] == "accelerated" and record[5:] == (1, 1)
            for record in fixed_point_records[:-1]
        )

    if fock_backend == "direct-eri":
        # qc_density_fit=0 reaches the four-center direct-ERI implementation.
        # It may use incremental accumulation during candidate generation, but
        # the physical map which either certifies an ordinary solution or
        # launches the ensemble solver must be a fresh full F[P] build.
        assert any(
            mode == "incremental" for _, _, mode in fock_build_records[:-1]
        )
        assert fock_build_records[-1][1:] == ("physical", "full")
    else:
        assert fock_build_records[-1][1:] == ("physical", "ri")


@pytest.mark.parametrize(
    "case_name,atoms,coordinates,expected_stage,expected_location",
    [
        pytest.param(
            "conversion-overflow",
            [("He", 4.0), ("He", 4.0)],
            [(0.0, 0.0, 0.0), (1.9e38, 0.0, 0.0)],
            "Angstrom-to-Bohr conversion",
            "local atom 1 (global atom 1), axis x",
            id="intermediate-conversion",
        ),
        pytest.param(
            "pair-overflow",
            [("He", 4.0), ("He", 4.0), ("He", 4.0)],
            [
                (0.0, 0.0, 0.0),
                (1.7e38, 0.0, 0.0),
                (-1.7e38, 0.0, 0.0),
            ],
            "pair displacement",
            "local atoms 1/2 (global atoms 1/2), axis x",
            id="pair-displacement",
        ),
    ],
)
def test_qc_nonfinite_geometry_is_fatal_at_the_producing_stage(
    tmp_path,
    case_name,
    atoms,
    coordinates,
    expected_stage,
    expected_location,
):
    case_dir = tmp_path / f"qc_nonfinite_geometry_{case_name}"
    _write_qc_case(
        case_dir,
        atoms=atoms,
        coordinates=coordinates,
        extra_settings={"pbc": False},
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert f"non-finite QC geometry during {expected_stage}" in output
    assert expected_location in output
    # Geometry must fail where it is produced.  Letting it flow into integral
    # kernels and later reporting a non-finite SCF energy loses the responsible
    # atom, axis, and transformation stage.
    assert "non-finite SCF energy" not in output
    assert "SCF failed to converge" not in output


def test_qc_nonfinite_raw_geometry_from_md_motion_is_fatal(tmp_path):
    case_dir = tmp_path / "qc_nonfinite_raw_geometry"
    _write_qc_case(
        case_dir,
        atoms=[("H", 1.0), ("H", 1.0)],
        coordinates=[(9.63, 10.0, 10.0), (10.37, 10.0, 10.0)],
        velocities=[(1.0e18, 0.0, 0.0), (0.0, 0.0, 0.0)],
        step_limit=1,
        dt=1.0e21,
    )
    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "QUANTUM_CHEMISTRY::Update_Coordinates_From_MD" in output
    assert "non-finite QC geometry during raw-coordinate input" in output
    assert "local atom 0 (global atom 0), axis x" in output
    assert "non-finite SCF energy" not in output
    assert "SCF failed to converge" not in output


def test_qc_coordinate_generation_applies_a_real_second_step_motion(tmp_path):
    case_dir = tmp_path / "qc_coordinate_generation_motion"
    _write_qc_case(
        case_dir,
        atoms=[("H", 1.0), ("H", 1.0)],
        coordinates=[(9.63, 10.0, 10.0), (10.37, 10.0, 10.0)],
        velocities=[(-1.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
        step_limit=1,
        dt=1.0e-4,
        scf_print=True,
    )
    result = _run_case(case_dir)
    final_energy = _final_scf_energy_by_step(result.stdout + result.stderr)
    assert set(final_energy) == {0, 1}
    assert all(math.isfinite(value) for value in final_energy.values())
    assert abs(final_energy[1] - final_energy[0]) > 1.0e-7


def test_qc_dd_global_coordinates_gradient_mapping_and_single_energy(
    tmp_path, mpi_np
):
    if mpi_np is None or mpi_np < 2:
        pytest.skip("run with --mpi 2 to exercise a reordered DD layout")

    atoms = [("F", 19.0), ("H", 1.0), ("He", 4.0)]
    coordinates = [
        (9.45, 9.85, 10.10),
        (10.35, 10.20, 9.90),
        (8.70, 10.75, 10.35),
    ]
    serial_dir = tmp_path / "qc_serial"
    mpi_dir = tmp_path / "qc_mpi"
    _write_qc_case(serial_dir, atoms=atoms, coordinates=coordinates)
    _write_qc_case(mpi_dir, atoms=atoms, coordinates=coordinates)
    _run_case(serial_dir)
    _run_case(mpi_dir, mpi_np=mpi_np)

    serial_mdout = _read_mdout(serial_dir)
    mpi_mdout = _read_mdout(mpi_dir)
    assert serial_mdout["potential"] == pytest.approx(
        serial_mdout["QC"], abs=0.02
    )
    assert mpi_mdout["potential"] == pytest.approx(mpi_mdout["QC"], abs=0.02)
    assert mpi_mdout["QC"] == pytest.approx(serial_mdout["QC"], abs=5.1e-3)
    assert _read_forces(mpi_dir, len(atoms)) == pytest.approx(
        _read_forces(serial_dir, len(atoms)), abs=2.0e-3
    )


def test_qc_virial_uses_the_same_unwrapped_geometry_as_scf(tmp_path):
    atoms = [("H", 1.0), ("H", 1.0)]
    zero_velocities = [(0.0, 0.0, 0.0)] * 2
    wrapped_dir = tmp_path / "qc_wrapped_virial"
    centered_dir = tmp_path / "qc_centered_virial"
    _write_qc_case(
        wrapped_dir,
        atoms=atoms,
        coordinates=[(0.1, 10.0, 10.0), (19.9, 10.0, 10.0)],
        velocities=zero_velocities,
        print_pressure=True,
    )
    _write_qc_case(
        centered_dir,
        atoms=atoms,
        coordinates=[(10.1, 10.0, 10.0), (9.9, 10.0, 10.0)],
        velocities=zero_velocities,
        print_pressure=True,
    )
    _run_case(wrapped_dir)
    _run_case(centered_dir)

    wrapped_mdout = _read_mdout(wrapped_dir)
    centered_mdout = _read_mdout(centered_dir)
    assert wrapped_mdout["QC"] == pytest.approx(
        centered_mdout["QC"], abs=5.1e-3
    )
    assert wrapped_mdout["pressure"] == pytest.approx(
        centered_mdout["pressure"], abs=0.02
    )
    assert _read_forces(wrapped_dir, 2) == pytest.approx(
        _read_forces(centered_dir, 2), abs=2.0e-3
    )


def test_qc_nopbc_uses_direct_displacement_for_energy_and_gradient(tmp_path):
    atoms = [("H", 1.0), ("H", 1.0)]
    coordinates = [(0.1, 10.0, 10.0), (19.9, 10.0, 10.0)]
    nopbc_dir = tmp_path / "qc_nopbc_direct"
    direct_pbc_dir = tmp_path / "qc_pbc_large_box_direct"
    wrapped_pbc_dir = tmp_path / "qc_pbc_small_box_wrapped"
    _write_qc_case(
        nopbc_dir,
        atoms=atoms,
        coordinates=coordinates,
        extra_settings={"pbc": False},
    )
    _write_qc_case(
        direct_pbc_dir,
        atoms=atoms,
        coordinates=coordinates,
        box=(40.0, 40.0, 40.0),
        extra_settings={"pbc": True},
    )
    _write_qc_case(
        wrapped_pbc_dir,
        atoms=atoms,
        coordinates=coordinates,
        extra_settings={"pbc": True},
    )

    _run_case(nopbc_dir)
    _run_case(direct_pbc_dir)
    _run_case(wrapped_pbc_dir)

    nopbc_mdout = _read_mdout(nopbc_dir)
    direct_pbc_mdout = _read_mdout(direct_pbc_dir)
    wrapped_pbc_mdout = _read_mdout(wrapped_pbc_dir)
    assert nopbc_mdout["QC"] == pytest.approx(
        direct_pbc_mdout["QC"], abs=5.1e-3
    )
    assert _read_forces(nopbc_dir, 2) == pytest.approx(
        _read_forces(direct_pbc_dir, 2), abs=2.0e-3
    )
    assert abs(nopbc_mdout["QC"] - wrapped_pbc_mdout["QC"]) > 100.0


@pytest.mark.parametrize("periodic", [False, True], ids=["nopbc", "pbc"])
def test_qc_exact_nuclear_overlap_is_rejected_without_distance_floor(
    tmp_path, periodic
):
    case_dir = tmp_path / f"qc_nuclear_overlap_{periodic}"
    _write_qc_case(
        case_dir,
        atoms=[("H", 1.0), ("H", 1.0)],
        coordinates=[(10.0, 10.0, 10.0), (10.0, 10.0, 10.0)],
        extra_settings={"pbc": periodic},
    )

    result = _run_case(case_dir, check=False)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "QC nuclei 0 (global atom 0) and 1 (global atom 1)" in output
    boundary = "periodic" if periodic else "non-periodic"
    assert f"overlap exactly under {boundary} boundary geometry" in output
