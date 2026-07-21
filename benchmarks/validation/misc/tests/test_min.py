import math
import os
import struct
import subprocess
from pathlib import Path

import pytest

from benchmarks.utils import Outputer, Runner


def write_minimization_mdin(case_dir, *, step_limit=100000, overrides=None):
    settings = {
        "md_name": '"validation tip3p minimization bad coordinate"',
        "mode": '"minimization"',
        "step_limit": str(step_limit),
        "cutoff": "8.0",
        "default_in_file_prefix": '"tip3p"',
        "coordinate_in_file": '"bad_coordinate.txt"',
        "minimization_dynamic_dt": "1",
        "minimization_max_move": "0.05",
        "print_zeroth_frame": "1",
        "write_mdout_interval": "1000",
        "write_information_interval": "1000",
    }
    settings.update(overrides or {})
    mdin = "".join(f"{key} = {value}\n" for key, value in settings.items())
    Path(case_dir, "mdin.spg.toml").write_text(mdin, encoding="utf-8")


def run_sponge_raw(case_dir, *, mpi_np=None):
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]
    if mpi_np is not None:
        command = ["mpirun", "--oversubscribe", "-np", str(mpi_np)] + command
    return subprocess.run(
        command,
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


def parse_potential_by_step(mdout_path):
    lines = Path(mdout_path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError(f"Invalid mdout file: {mdout_path}")

    headers = lines[0].split()
    if "step" not in headers or "potential" not in headers:
        raise ValueError(f"Missing step/potential columns in {mdout_path}")
    step_idx = headers.index("step")
    pot_idx = headers.index("potential")

    values = {}
    for line in lines[1:]:
        fields = line.split()
        if len(fields) <= max(step_idx, pot_idx):
            continue
        try:
            step = int(fields[step_idx])
            potential = float(fields[pot_idx])
        except ValueError:
            continue
        values[step] = potential
    return values


def test_tip3p_bad_coordinate_minimization_runs(
    statics_path, outputs_path, mpi_np
):
    step_limit = 100000
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p",
        mpi_np=mpi_np,
        run_name="tip3p_min_bad_coordinate",
    )
    write_minimization_mdin(case_dir, step_limit=step_limit)
    Runner.run_sponge(case_dir, timeout=1200, mpi_np=mpi_np)

    potential_by_step = parse_potential_by_step(case_dir / "mdout.txt")
    assert all(math.isfinite(value) for value in potential_by_step.values())
    assert potential_by_step[1000] < potential_by_step[0]
    milestones = [
        ("Initial", 0),
        ("1/4", step_limit // 4),
        ("1/2", step_limit // 2),
        ("3/4", (step_limit * 3) // 4),
        ("Final", step_limit),
    ]
    milestone_rows = []
    for stage, step in milestones:
        if step not in potential_by_step:
            raise AssertionError(
                f"Missing step {step} in mdout.txt; "
                "check write_mdout_interval setting."
            )
        milestone_rows.append(
            [stage, str(step), f"{potential_by_step[step]:.6e}"]
        )
    Outputer.print_table(
        ["Stage", "Step", "Potential"],
        milestone_rows,
        title="Misc Validation: TIP3P Minimization Energy Milestones",
    )

    Outputer.print_table(
        ["Metric", "Value"],
        [
            ["Case", "tip3p"],
            ["InputCoordinate", "bad_coordinate.txt"],
            ["Mode", "minimization"],
            ["StepLimit", str(step_limit)],
            ["Status", "PASS"],
        ],
        title="Misc Validation: TIP3P Minimization",
    )

    final_potential = potential_by_step[step_limit]
    assert final_potential < -4000.0


def test_tip3p_extreme_finite_force_has_finite_first_adam_step(
    statics_path, outputs_path, mpi_np
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p",
        mpi_np=mpi_np,
        run_name="tip3p_min_extreme_force_first_step",
    )
    write_minimization_mdin(
        case_dir,
        step_limit=1,
        overrides={
            "write_mdout_interval": "1",
            "write_information_interval": "1",
        },
    )
    Runner.run_sponge(case_dir, timeout=120, mpi_np=mpi_np)

    potential_by_step = parse_potential_by_step(case_dir / "mdout.txt")
    assert set(potential_by_step) == {0, 1}
    assert math.isfinite(potential_by_step[0])
    assert math.isfinite(potential_by_step[1])
    assert potential_by_step[1] < potential_by_step[0]


def test_adam_keeps_physical_force_trajectory_unmodified(
    statics_path, outputs_path, mpi_np
):
    force_payloads = {}
    for dynamic_dt in (0, 1):
        case_dir = Outputer.prepare_output_case(
            statics_path=statics_path,
            outputs_path=outputs_path,
            case_name="tip3p",
            mpi_np=mpi_np,
            run_name=f"tip3p_min_force_contract_{dynamic_dt}",
        )
        write_minimization_mdin(
            case_dir,
            step_limit=0,
            overrides={
                "minimization_dynamic_dt": str(dynamic_dt),
                "frc": '"force.bin"',
                "write_trajectory_interval": "1",
                "write_mdout_interval": "1",
                "write_information_interval": "1",
            },
        )
        Runner.run_sponge(
            case_dir,
            timeout=120,
            mpi_np=mpi_np,
            env={**os.environ, "OMP_NUM_THREADS": "1"},
        )
        payload = (case_dir / "force.bin").read_bytes()
        assert payload
        assert len(payload) % (3 * struct.calcsize("=f")) == 0
        values = struct.unpack(f"={len(payload) // 4}f", payload)
        assert all(math.isfinite(value) for value in values)
        force_payloads[dynamic_dt] = values

    # Optimizer selection happens after force assembly.  It may change the
    # coordinate move, but never the physical force exposed to trajectories.
    assert force_payloads[1] == pytest.approx(
        force_payloads[0], rel=1.0e-6, abs=1.0e-6
    )


def test_adam_coordinate_update_is_not_implicitly_mass_preconditioned(
    tmp_path, mpi_np
):
    case_dir = tmp_path / "adam_mass_contract"
    case_dir.mkdir()
    (case_dir / "mass.txt").write_text("2\n1\n12\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("2\n0\n0\n", encoding="utf-8")
    (case_dir / "coordinate.txt").write_text(
        "2\n19 20 20\n21 20 20\n40 40 40\n90 90 90\n",
        encoding="utf-8",
    )
    # U(r) = 4096/r^12 gives equal-and-opposite, nonzero forces at r=2.
    (case_dir / "lj.txt").write_text("2 1\n4096\n0\n0\n0\n", encoding="utf-8")
    (case_dir / "mdin.spg.toml").write_text(
        """md_name = "adam mass contract"
mode = "minimization"
step_limit = 1
cutoff = 4.0
mass_in_file = "mass.txt"
charge_in_file = "charge.txt"
coordinate_in_file = "coordinate.txt"
LJ_in_file = "lj.txt"
PM.MPI_size = 0
minimization_dynamic_dt = 1
minimization_learning_rate = 0.01
minimization_max_move = 0.1
crd = "coordinate.bin"
print_zeroth_frame = 1
write_mdout_interval = 1
write_information_interval = 1
write_trajectory_interval = 1
write_restart_file_interval = 0
""",
        encoding="utf-8",
    )
    Runner.run_sponge(case_dir, timeout=120, mpi_np=mpi_np)

    payload = (case_dir / "coordinate.bin").read_bytes()
    coordinates = struct.unpack("=6f", payload[-6 * 4 :])
    first_displacement = coordinates[0] - 19.0
    second_displacement = coordinates[3] - 21.0
    assert first_displacement * second_displacement < 0.0
    # The theoretical moves are equal; their reconstructed displacements can
    # differ by one float ULP because they are added to different coordinates.
    assert abs(first_displacement) == pytest.approx(
        abs(second_displacement), rel=1.0e-5, abs=2.5e-6
    )
    assert abs(first_displacement) > 0.009


def write_minimization_numeric_update_case(
    case_dir,
    *,
    coordinate="10 20 30",
    velocity="0 0 0",
    dynamic_dt,
    max_move,
    dt="0.001",
    learning_rate=None,
    constant_force=False,
):
    case_dir.mkdir()
    (case_dir / "mass.txt").write_text("1\n1\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("1\n0\n", encoding="utf-8")
    (case_dir / "coordinate.txt").write_text(
        f"1\n{coordinate}\n1000 1000 1000\n90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "velocity.txt").write_text(f"1\n{velocity}\n", encoding="utf-8")
    settings = [
        'md_name = "minimization numeric update contract"',
        'mode = "minimization"',
        "step_limit = 0",
        f"dt = {dt}",
        "cutoff = 101.0",
        "pbc = false",
        "PM.MPI_size = 0",
        'mass_in_file = "mass.txt"',
        'charge_in_file = "charge.txt"',
        'coordinate_in_file = "coordinate.txt"',
        'velocity_in_file = "velocity.txt"',
        f"minimization_dynamic_dt = {dynamic_dt}",
        f"minimization_max_move = {max_move}",
        "minimization_momentum_keep = 1.0",
        "write_mdout_interval = 1",
        "write_information_interval = 1",
        "write_trajectory_interval = 0",
        "write_restart_file_interval = 0",
    ]
    if learning_rate is not None:
        settings.append(f"minimization_learning_rate = {learning_rate}")
    if constant_force:
        (case_dir / "listed_force.txt").write_text(
            """[[[ constant_push ]]]
[[ parameters ]]
int atom_i, float magnitude_i
[[ potential ]]
E = -magnitude_i * r_i.x;
[[ end ]]
""",
            encoding="utf-8",
        )
        (case_dir / "constant_push.txt").write_text(
            "1\n0 1\n", encoding="utf-8"
        )
        settings.extend(
            (
                'listed_forces_in_file = "listed_force.txt"',
                'constant_push_in_file = "constant_push.txt"',
            )
        )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(settings) + "\n", encoding="utf-8"
    )


def write_minimization_mpi_contract_case(
    case_dir,
    *,
    coordinates,
    velocities,
    dynamic_dt,
    step_limit,
    dt=1.0,
    force_atom=None,
    learning_rate=0.001,
    explicit_two_x_domains=True,
):
    atom_numbers = len(coordinates)
    assert len(velocities) == atom_numbers
    case_dir.mkdir()
    (case_dir / "mass.txt").write_text(
        f"{atom_numbers}\n" + "1\n" * atom_numbers, encoding="utf-8"
    )
    (case_dir / "charge.txt").write_text(
        f"{atom_numbers}\n" + "0\n" * atom_numbers, encoding="utf-8"
    )
    (case_dir / "residue.txt").write_text(
        f"{atom_numbers} {atom_numbers}\n" + "1\n" * atom_numbers,
        encoding="utf-8",
    )
    (case_dir / "coordinate.txt").write_text(
        f"{atom_numbers}\n"
        + "\n".join(
            " ".join(str(value) for value in crd) for crd in coordinates
        )
        + "\n40 40 40\n90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "velocity.txt").write_text(
        f"{atom_numbers}\n"
        + "\n".join(" ".join(str(value) for value in vel) for vel in velocities)
        + "\n",
        encoding="utf-8",
    )

    settings = [
        'md_name = "minimization MPI state contract"',
        'mode = "minimization"',
        f"step_limit = {step_limit}",
        f"dt = {dt}",
        "cutoff = 4.0",
        "PM.MPI_size = 0",
        'mass_in_file = "mass.txt"',
        'charge_in_file = "charge.txt"',
        'residue_in_file = "residue.txt"',
        'coordinate_in_file = "coordinate.txt"',
        'velocity_in_file = "velocity.txt"',
        f"minimization_dynamic_dt = {dynamic_dt}",
        "minimization_max_move = 0.0",
        "minimization_momentum_keep = 1.0",
        f"minimization_learning_rate = {learning_rate}",
        'crd = "coordinate.bin"',
        "print_zeroth_frame = 1",
        "write_mdout_interval = 1",
        "write_information_interval = 1",
        "write_trajectory_interval = 1",
        "write_restart_file_interval = 0",
        "",
        "[DOM_DEC]",
        "update_interval = 1",
    ]
    if explicit_two_x_domains:
        settings.extend(("split_nx = 2", "split_ny = 1", "split_nz = 1"))
    if force_atom is not None:
        (case_dir / "listed_force.txt").write_text(
            """[[[ constant_push ]]]
[[ parameters ]]
int atom_i, float magnitude_i
[[ potential ]]
E = -magnitude_i * r_i.x;
[[ end ]]
""",
            encoding="utf-8",
        )
        (case_dir / "constant_push.txt").write_text(
            f"1\n{force_atom} 1\n", encoding="utf-8"
        )
        settings[settings.index("")] = (
            'listed_forces_in_file = "listed_force.txt"'
        )
        settings.insert(
            settings.index("[DOM_DEC]"),
            'constant_push_in_file = "constant_push.txt"',
        )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(settings) + "\n", encoding="utf-8"
    )


def read_coordinate_frames(path, atom_numbers):
    payload = path.read_bytes()
    values = struct.unpack(f"={len(payload) // 4}f", payload)
    frame_width = 3 * atom_numbers
    assert len(values) % frame_width == 0
    return [
        values[offset : offset + frame_width]
        for offset in range(0, len(values), frame_width)
    ]


@pytest.mark.parametrize("max_move", ["0.0", "0.1"], ids=["uncapped", "capped"])
def test_classic_minimization_rejects_overflow_at_coordinate_update(
    tmp_path, max_move
):
    case_dir = tmp_path / f"classic_update_overflow_{max_move}"
    write_minimization_numeric_update_case(
        case_dir,
        velocity="3e38 0 0",
        dynamic_dt=0,
        max_move=max_move,
        dt="1.0",
    )

    result = run_sponge_raw(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "minimization coordinate update" in output
    assert "finite zero or normal float" in output
    assert "global atom 0 at step 0" in output


def test_mpi_minimization_reports_invalid_nonroot_global_atom(tmp_path, mpi_np):
    if mpi_np != 2:
        pytest.skip("requires the explicit two-PP x-domain layout")
    case_dir = tmp_path / "mpi_nonroot_invalid_atom"
    write_minimization_mpi_contract_case(
        case_dir,
        coordinates=((5, 5, 5), (25, 5, 5)),
        velocities=((0, 0, 0), ("3e38", 0, 0)),
        dynamic_dt=0,
        step_limit=0,
        dt=2.0,
    )

    result = run_sponge_raw(case_dir, mpi_np=mpi_np)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "minimization coordinate update" in output
    assert "global atom 1 at step 0" in output


def test_mpi_adam_state_follows_atom_across_domain(tmp_path, mpi_np):
    if mpi_np != 2:
        pytest.skip("requires the explicit two-PP x-domain layout")
    coordinates = ((19, 5, 5), (5, 10, 5), (30, 15, 5))
    velocities = ((0, 0, 0),) * 3
    serial_dir = tmp_path / "adam_migration_serial"
    mpi_dir = tmp_path / "adam_migration_mpi"
    write_minimization_mpi_contract_case(
        serial_dir,
        coordinates=coordinates,
        velocities=velocities,
        dynamic_dt=1,
        step_limit=3,
        force_atom=0,
        learning_rate=2.0,
        explicit_two_x_domains=False,
    )
    write_minimization_mpi_contract_case(
        mpi_dir,
        coordinates=coordinates,
        velocities=velocities,
        dynamic_dt=1,
        step_limit=3,
        force_atom=0,
        learning_rate=2.0,
    )

    serial_result = run_sponge_raw(serial_dir)
    mpi_result = run_sponge_raw(mpi_dir, mpi_np=mpi_np)
    assert serial_result.returncode == 0, (
        serial_result.stdout + serial_result.stderr
    )
    assert mpi_result.returncode == 0, mpi_result.stdout + mpi_result.stderr
    serial_frames = read_coordinate_frames(serial_dir / "coordinate.bin", 3)
    mpi_frames = read_coordinate_frames(mpi_dir / "coordinate.bin", 3)
    assert len(serial_frames) == len(mpi_frames) == 4
    # The first minimization move crosses the explicitly fixed x=20 domain
    # boundary, and update_interval=1 forces immediate ownership migration.
    assert serial_frames[0][0] > 20.0
    assert mpi_frames[0][0] > 20.0
    for serial_frame, mpi_frame in zip(serial_frames, mpi_frames):
        assert mpi_frame == pytest.approx(serial_frame, rel=2.0e-6, abs=2.0e-6)


def test_adam_rejects_coordinate_addition_overflow_at_its_source(tmp_path):
    case_dir = tmp_path / "adam_coordinate_addition_overflow"
    write_minimization_numeric_update_case(
        case_dir,
        coordinate="1e37 20 30",
        dynamic_dt=1,
        max_move="0.0",
        learning_rate="3.4e38",
        constant_force=True,
    )

    result = run_sponge_raw(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "minimization coordinate update" in output
    assert "finite zero or normal float" in output
    assert "global atom 0 at step 0" in output


@pytest.mark.parametrize(
    ("overrides", "message"),
    [
        (
            {"minimization_dynamic_dt": "2"},
            "minimization_dynamic_dt must be 0 or 1",
        ),
        (
            {"minimization_max_move": "-0.01"},
            "minimization_max_move must be finite and non-negative",
        ),
        (
            {"minimization_max_move": "3.5e38"},
            "minimization_max_move must be finite and non-negative",
        ),
        (
            {"minimization_max_move": "1e-45"},
            "minimization_max_move must be zero or a normal float",
        ),
        (
            {"minimization_beta1": "1.0"},
            "Adam minimization requires 0 <= beta1,beta2 < 1",
        ),
        (
            {"minimization_beta2": "-0.01"},
            "Adam minimization requires 0 <= beta1,beta2 < 1",
        ),
        (
            {"minimization_epsilon": "0.0"},
            "epsilon > 0",
        ),
        (
            {"minimization_epsilon": "1e-45"},
            "positive values must be normal floats",
        ),
        (
            {"minimization_learning_rate": "0.0"},
            "learning_rate > 0",
        ),
        (
            {"minimization_learning_rate": "3.5e38"},
            "learning_rate > 0",
        ),
        (
            {"minimization_learning_rate": "1e-45"},
            "positive values must be normal floats",
        ),
        (
            {"dt": "0.0"},
            "minimization dt must be a finite positive normal float",
        ),
        (
            {"dt": "-0.001"},
            "minimization dt must be a finite positive normal float",
        ),
        (
            {"dt": "3.5e38"},
            "minimization dt must be a finite positive normal float",
        ),
        (
            {"dt": "1e-45"},
            "minimization dt must be a finite positive normal float",
        ),
        (
            {
                "minimization_dynamic_dt": "0",
                "minimization_momentum_keep": "-0.01",
            },
            "minimization_momentum_keep must be finite and within [0, 1]",
        ),
        (
            {
                "minimization_dynamic_dt": "0",
                "minimization_momentum_keep": "1e-45",
            },
            "any nonzero value must be a normal float",
        ),
    ],
    ids=[
        "dynamic_dt_not_boolean",
        "negative_max_move",
        "overflowed_max_move",
        "subnormal_max_move",
        "beta1_endpoint",
        "negative_beta2",
        "zero_epsilon",
        "subnormal_epsilon",
        "zero_learning_rate",
        "overflowed_learning_rate",
        "subnormal_learning_rate",
        "zero_dt",
        "negative_dt",
        "overflowed_dt",
        "subnormal_dt",
        "negative_momentum_keep",
        "subnormal_momentum_keep",
    ],
)
def test_minimization_rejects_invalid_optimizer_parameters(
    statics_path, outputs_path, mpi_np, overrides, message
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p",
        mpi_np=mpi_np,
        run_name="tip3p_min_invalid_" + next(iter(overrides)),
    )
    write_minimization_mdin(case_dir, step_limit=0, overrides=overrides)
    result = run_sponge_raw(case_dir, mpi_np=mpi_np)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert message in output
