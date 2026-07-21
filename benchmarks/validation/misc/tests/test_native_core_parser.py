import json
import math
import os
import struct
import subprocess

import pytest


def _valid_files(*, atom_count=2, include_optional=False):
    coordinate_lines = [f"{atom_count} 1.25 7"]
    coordinate_lines.extend(
        f"{atom + 1}.0 2.0 3.0" for atom in range(atom_count)
    )
    coordinate_lines.extend(("40.0 40.0 40.0", "90.0 90.0 90.0"))
    files = {
        "mass.txt": f"{atom_count}\n" + "12.0\n" * atom_count,
        "charge.txt": f"{atom_count}\n" + "0.0\n" * atom_count,
        "coordinate.txt": "\n".join(coordinate_lines) + "\n",
    }
    if include_optional:
        velocity_lines = [f"{atom_count} 1.25 7"]
        velocity_lines.extend("0.1 0.2 0.3" for _ in range(atom_count))
        exclusion_lines = [
            f"{atom_count} {1 if atom_count >= 2 else 0}",
            "1 1" if atom_count >= 2 else "0",
        ]
        exclusion_lines.extend("0" for _ in range(1, atom_count))
        files.update(
            {
                "velocity.txt": "\n".join(velocity_lines) + "\n",
                "residue.txt": f"{atom_count} 1\n{atom_count}\n",
                "exclude.txt": "\n".join(exclusion_lines) + "\n",
            }
        )
    return files


def _write_case(
    case_dir,
    *,
    atom_count=2,
    overrides=None,
    include_optional=False,
    settings_overrides=None,
):
    case_dir.mkdir()
    files = _valid_files(
        atom_count=atom_count, include_optional=include_optional
    )
    if overrides:
        files.update(overrides)
    for name, contents in files.items():
        (case_dir / name).write_text(contents, encoding="utf-8")

    settings = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 0,
        "dt": 0,
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    if include_optional or "velocity.txt" in files:
        settings["velocity_in_file"] = "velocity.txt"
    if include_optional or "residue.txt" in files:
        settings["residue_in_file"] = "residue.txt"
    if include_optional or "exclude.txt" in files:
        settings["exclude_in_file"] = "exclude.txt"
    if settings_overrides:
        settings.update(settings_overrides)
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
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


def _assert_rejected(case_dir, filename, *messages):
    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert filename in output
    for message in messages:
        assert message in output


def _coordinate_text(box_lengths, box_angles, coordinates=None):
    if coordinates is None:
        coordinates = ((1.0, 2.0, 3.0), (2.0, 2.0, 3.0))
    lines = [str(len(coordinates))]
    lines.extend(" ".join(str(value) for value in atom) for atom in coordinates)
    lines.append(" ".join(str(value) for value in box_lengths))
    lines.append(" ".join(str(value) for value in box_angles))
    return "\n".join(lines) + "\n"


def _write_rerun_case(
    case_dir,
    *,
    trajectory_values,
    box_text,
    velocity_values=None,
    extra_files=None,
    settings_overrides=None,
):
    settings = {
        "mode": "rerun",
        "pbc": False,
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "crd": "trajectory.dat",
        "box": "box.txt",
        "rerun_need_box_update": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    if velocity_values is not None:
        settings["vel"] = "velocity.dat"
    if settings_overrides:
        settings.update(settings_overrides)
    _write_case(
        case_dir,
        overrides=extra_files,
        settings_overrides=settings,
    )
    (case_dir / "trajectory.dat").write_bytes(
        struct.pack(f"={len(trajectory_values)}f", *trajectory_values)
    )
    (case_dir / "box.txt").write_text(box_text, encoding="utf-8")
    if velocity_values is not None:
        (case_dir / "velocity.dat").write_bytes(
            struct.pack(f"={len(velocity_values)}f", *velocity_values)
        )


def test_valid_native_core_files_and_restart_headers_run(tmp_path):
    case_dir = tmp_path / "valid_native_core"
    _write_case(case_dir, include_optional=True)

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_missing_optional_velocity_residue_and_exclusion_keep_defaults(
    tmp_path,
):
    case_dir = tmp_path / "native_core_defaults"
    _write_case(case_dir)

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_mixed_gromacs_and_amber_source_selection_is_rejected(tmp_path):
    case_dir = tmp_path / "mixed_external_sources"
    _write_case(
        case_dir,
        settings_overrides={
            "gromacs_top": "unused.top",
            "amber_rst7": "unused.rst7",
        },
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "GROMACS and AMBER input sources cannot be selected together" in output


def test_toml_scalar_paths_workspace_and_task_name_preserve_full_value(
    tmp_path,
):
    launch_dir = tmp_path / "launcher"
    launch_dir.mkdir()
    workspace_name = "workspace # with spaces"
    workspace = launch_dir / workspace_name
    _write_case(workspace)

    renamed_inputs = {
        "mass.txt": "mass values # full.txt",
        "charge.txt": "charge values # full.txt",
        "coordinate.txt": "coordinate values # full.txt",
    }
    for old_name, new_name in renamed_inputs.items():
        (workspace / old_name).rename(workspace / new_name)

    mdin_source = workspace / "mdin.spg.toml"
    contents = mdin_source.read_text(encoding="utf-8")
    contents = contents.replace(
        f"md_name = {json.dumps(workspace.name)}",
        'md_name = "task # name with spaces"',
    )
    for command, old_name in (
        ("mass_in_file", "mass.txt"),
        ("charge_in_file", "charge.txt"),
        ("coordinate_in_file", "coordinate.txt"),
    ):
        contents = contents.replace(
            f"{command} = {json.dumps(old_name)}",
            f"{command} = {json.dumps(renamed_inputs[old_name])}",
        )
    output_name = "md output # full.txt"
    contents += (
        f"workspace = {json.dumps(workspace_name)}\n"
        f"mdout = {json.dumps(output_name)}\n"
    )
    (launch_dir / "mdin.spg.toml").write_text(contents, encoding="utf-8")
    mdin_source.unlink()

    result = _run_case(launch_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "task # name with spaces" in output
    assert (workspace / output_name).is_file()


@pytest.mark.parametrize(
    "settings",
    [
        {"write_mdout_interval": 0},
        {"write_information_interval": 0},
    ],
)
def test_zero_output_interval_cleanly_disables_output(tmp_path, settings):
    case_dir = tmp_path / next(iter(settings))
    # Run beyond step zero so the historical short-circuit cannot hide a
    # modulo-by-zero in the force-evaluation output predicate.
    _write_case(
        case_dir,
        settings_overrides={"step_limit": 1, **settings},
    )

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


@pytest.mark.parametrize(
    ("settings", "message"),
    [
        ({"step_limit": -1}, "step/frame limit must be nonnegative"),
        (
            {"write_mdout_interval": -1},
            "output intervals must be nonnegative",
        ),
        (
            {"write_trajectory_interval": -1},
            "output intervals must be nonnegative",
        ),
        (
            {"write_restart_file_interval": -1},
            "output intervals must be nonnegative",
        ),
    ],
)
def test_negative_step_or_output_limits_are_rejected(
    tmp_path, settings, message
):
    case_dir = tmp_path / next(iter(settings))
    _write_case(case_dir, settings_overrides=settings)

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert message in output


@pytest.mark.parametrize(
    ("name", "settings"),
    [
        ("zero_cutoff", {"cutoff": 0.0}),
        ("negative_cutoff", {"cutoff": -1.0}),
        ("overflowed_cutoff", {"cutoff": 1.0e39}),
        ("subnormal_cutoff", {"cutoff": 1.0e-45}),
        ("negative_skin", {"skin": -0.1}),
        ("overflowed_skin", {"skin": 1.0e39}),
        ("subnormal_skin", {"skin": 1.0e-45}),
    ],
)
def test_nonbond_geometry_rejects_invalid_cutoff_or_skin(
    tmp_path, name, settings
):
    case_dir = tmp_path / name
    _write_case(case_dir, settings_overrides=settings)

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "cutoff must be finite and positive" in output
    assert "skin must be finite and nonnegative" in output


@pytest.mark.parametrize(
    ("coordinate_header", "velocity_header"),
    [("2", None), ("2 .5", "2 .5"), ("2 1. 7", "2 1. 7")],
)
def test_valid_native_header_forms_run(
    tmp_path, coordinate_header, velocity_header
):
    case_dir = tmp_path / f"valid_header_{len(coordinate_header)}"
    overrides = {
        "mass.txt": "2\n.5\n1.\n",
        "charge.txt": "2\n+1e-1\n-2.E-1\n",
        "coordinate.txt": (
            f"{coordinate_header}\n"
            ".5 2. 3e0\n"
            "+2.0 2.0 3.0\n"
            "4e1 40. 40\n"
            "90 9e1 +90.0\n"
        ),
    }
    if velocity_header is not None:
        overrides["velocity.txt"] = (
            f"{velocity_header}\n.1 .2 .3\n-.1 -.2 -.3\n"
        )
    _write_case(case_dir, overrides=overrides)

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_missing_native_core_file_reports_input_path(tmp_path):
    case_dir = tmp_path / "missing_mass_file"
    _write_case(case_dir)
    (case_dir / "mass.txt").unlink()

    _assert_rejected(case_dir, "mass.txt", "failed to open mass_in_file")


@pytest.mark.parametrize(
    ("name", "box_lengths", "box_angles", "messages"),
    [
        (
            "zero_length",
            (0, 40, 40),
            (90, 90, 90),
            ("box length x", "strictly positive"),
        ),
        (
            "negative_length",
            (40, -1, 40),
            (90, 90, 90),
            ("box length y", "strictly positive"),
        ),
        (
            "zero_angle",
            (40, 40, 40),
            (90, 90, 0),
            ("box angle gamma", "strictly between 0 and 180"),
        ),
        (
            "straight_angle",
            (40, 40, 40),
            (180, 90, 90),
            ("box angle alpha", "strictly between 0 and 180"),
        ),
        (
            "degenerate_angles",
            (40, 40, 40),
            (10, 10, 30),
            ("triclinic Gram determinant", "non-positive"),
        ),
        (
            "sin_gamma_ftz",
            (40, 40, 40),
            (90, 90, "1e-37"),
            ("sin(gamma)", "flushed to zero"),
        ),
        (
            "volume_overflow",
            ("1e13", "1e13", "1e13"),
            (90, 90, 90),
            ("cell volume", "finite float range"),
        ),
        (
            "volume_underflow",
            ("1e-20", "1e-20", "1e-20"),
            (90, 90, 90),
            ("cell volume", "normal float"),
        ),
        (
            "reciprocal_underflow",
            ("1e38", 1, 1),
            (90, 90, 90),
            ("reciprocal cell a11", "normal float"),
        ),
    ],
)
def test_invalid_periodic_cell_is_rejected_at_shared_builder(
    tmp_path, name, box_lengths, box_angles, messages
):
    case_dir = tmp_path / f"invalid_cell_{name}"
    coordinate_text = _coordinate_text(box_lengths, box_angles)
    _write_case(case_dir, overrides={"coordinate.txt": coordinate_text})

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "spongeErrorBadFileFormat" in output
    assert "periodic_box_condition_information::Initial" in output
    assert "initial periodic-cell construction" in output
    for message in messages:
        assert message in output


def test_subnormal_box_length_is_rejected_by_native_float_contract(tmp_path):
    case_dir = tmp_path / "subnormal_box_length"
    coordinate_text = _coordinate_text(("1e-40", 1, 1), (90, 90, 90))
    _write_case(case_dir, overrides={"coordinate.txt": coordinate_text})

    _assert_rejected(
        case_dir,
        "coordinate.txt",
        "box length x",
        "subnormal float",
        "consistent FTZ behavior",
    )


def _write_triclinic_bond_case(case_dir, image_count):
    cell_a = (20.0, 0.0, 0.0)
    cell_b = (3.0, 25.0, 0.0)
    cell_c = (5.0, 6.0, 30.0)

    def norm(vector):
        return math.sqrt(sum(value * value for value in vector))

    def angle(first, second):
        dot = sum(a * b for a, b in zip(first, second))
        return math.degrees(math.acos(dot / norm(first) / norm(second)))

    box_lengths = (norm(cell_a), norm(cell_b), norm(cell_c))
    box_angles = (
        angle(cell_b, cell_c),
        angle(cell_a, cell_c),
        angle(cell_a, cell_b),
    )
    origin = (1.0, 2.0, 3.0)
    displacement = (1.2, -0.7, 0.9)
    second = tuple(
        origin[axis] + displacement[axis] + image_count * cell_c[axis]
        for axis in range(3)
    )
    coordinate_text = _coordinate_text(
        box_lengths, box_angles, (origin, second)
    )
    _write_case(
        case_dir,
        overrides={
            "coordinate.txt": coordinate_text,
            "bond.txt": "1\n0 1 1.3 0.5\n",
        },
    )
    with (case_dir / "mdin.spg.toml").open("a", encoding="utf-8") as mdin:
        mdin.write(
            'bond_in_file = "bond.txt"\n'
            'mdout = "mdout.txt"\n'
            'frc = "frc.dat"\n'
            "print_zeroth_frame = true\n"
            "write_mdout_interval = 1\n"
        )


def _read_mdout_terms(case_dir):
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    return dict(zip(lines[0].split(), map(float, lines[-1].split())))


def _read_last_forces(case_dir, atom_count):
    raw = (case_dir / "frc.dat").read_bytes()
    values = struct.unpack(f"={len(raw) // 4}f", raw)
    return values[-3 * atom_count :]


def test_triclinic_cell_and_shared_inverse_are_periodic_image_invariant(
    tmp_path,
):
    base_dir = tmp_path / "triclinic_base"
    imaged_dir = tmp_path / "triclinic_many_c_images"
    _write_triclinic_bond_case(base_dir, 0)
    # Seventeen c-vector images make the historical extra-a33 division in
    # LTMatrix3::inv(a31) cross an x-image boundary, so this is also a direct
    # regression for the shared triclinic inverse.
    _write_triclinic_bond_case(imaged_dir, 17)

    base_result = _run_case(base_dir)
    imaged_result = _run_case(imaged_dir)
    assert base_result.returncode == 0, base_result.stdout + base_result.stderr
    assert imaged_result.returncode == 0, (
        imaged_result.stdout + imaged_result.stderr
    )
    assert _read_mdout_terms(imaged_dir)["bond"] == pytest.approx(
        _read_mdout_terms(base_dir)["bond"], abs=2.0e-4
    )
    assert _read_last_forces(imaged_dir, 2) == pytest.approx(
        _read_last_forces(base_dir, 2), abs=2.0e-4
    )


def test_nopbc_box_length_cv_reads_the_validated_reference_cell(tmp_path):
    case_dir = tmp_path / "nopbc_reference_box_cv"
    cv_text = (
        "box_y\n"
        "{\n"
        "    CV_type = box_length_y\n"
        "}\n"
        "print\n"
        "{\n"
        "    CV = box_y\n"
        "}\n"
    )
    _write_case(
        case_dir,
        overrides={"cv.txt": cv_text},
        settings_overrides={
            "pbc": False,
            "cv_in_file": "cv.txt",
            "mdout": "mdout.txt",
            "print_zeroth_frame": True,
            "write_mdout_interval": 1,
        },
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    values = _read_mdout_terms(case_dir)
    assert math.isfinite(values["box_y"])
    assert values["box_y"] == pytest.approx(40.0, abs=1.0e-4)


def test_rerun_identical_triclinic_record_is_not_replayed_from_roundoff(
    tmp_path,
):
    case_dir = tmp_path / "rerun_identical_triclinic_record"
    frame = (1000.0, 1000.0, 1000.0, 0.0, 0.0, 0.0)
    triclinic_box = (
        "0.8950731158256531 15.028543472290039 3.018704891204834 "
        "104.5488052368164 107.60084533691406 29.174039840698242"
    )
    cv_text = (
        "scaled\n"
        "{\n"
        "    CV_type = scaled_position_x\n"
        "    atom = 0\n"
        "}\n"
        "print\n"
        "{\n"
        "    CV = scaled\n"
        "}\n"
    )
    _write_rerun_case(
        case_dir,
        trajectory_values=frame * 3,
        box_text=(
            "40 40 40 90 90 90\n"
            f"{triclinic_box}\n"
            f"{triclinic_box}\n"
        ),
        extra_files={"cv.txt": cv_text},
        settings_overrides={
            "cv_in_file": "cv.txt",
            "mdout": "mdout.txt",
        },
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    scaled_index = lines[0].split().index("scaled")
    scaled_values = [float(line.split()[scaled_index]) for line in lines[1:]]
    assert len(scaled_values) == 3
    assert scaled_values[1] == pytest.approx(-474.1510, abs=2.0e-4)
    # Reconstructing public box angles from the float cell changes gamma by
    # one ulp for this record.  The input record must remain authoritative, so
    # the identical third frame cannot trigger a second deformation.
    assert scaled_values[2] == scaled_values[1]


@pytest.mark.parametrize(
    ("name", "trajectory_values", "box_text", "velocity_values", "messages"),
    [
        (
            "empty_first_frame",
            (),
            "",
            None,
            ("no complete selected frame",),
        ),
        (
            "partial_second_coordinate",
            (
                1.0,
                2.0,
                3.0,
                2.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
            ),
            "40 40 40 90 90 90\n40 40 40 90 90 90\n",
            None,
            ("truncated binary rerun frame 1", "coordinate 1/2 atoms"),
        ),
        (
            "box_ends_before_coordinate",
            (
                1.0,
                2.0,
                3.0,
                2.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
            ),
            "40 40 40 90 90 90\n",
            None,
            ("out of sync at frame 1", "coordinate=complete", "box=clean EOF"),
        ),
        (
            "coordinate_ends_before_box",
            (1.0, 2.0, 3.0, 2.0, 2.0, 3.0),
            "40 40 40 90 90 90\n40 40 40 90 90 90\n",
            None,
            ("out of sync at frame 1", "coordinate=clean EOF", "box=complete"),
        ),
        (
            "partial_second_box",
            (
                1.0,
                2.0,
                3.0,
                2.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
            ),
            "40 40 40 90 90 90\n40 40 40\n",
            None,
            ("rerun box line 2", "exactly six"),
        ),
        (
            "velocity_ends_before_coordinate",
            (
                1.0,
                2.0,
                3.0,
                2.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
            ),
            "40 40 40 90 90 90\n40 40 40 90 90 90\n",
            (0.1, 0.2, 0.3, 0.4, 0.5, 0.6),
            ("out of sync at frame 1", "velocity=clean EOF"),
        ),
        (
            "partial_second_velocity",
            (
                1.0,
                2.0,
                3.0,
                2.0,
                2.0,
                3.0,
                4.0,
                5.0,
                6.0,
                7.0,
                8.0,
                9.0,
            ),
            "40 40 40 90 90 90\n40 40 40 90 90 90\n",
            (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9),
            ("truncated binary rerun frame 1", "velocity 1/2 atoms"),
        ),
    ],
)
def test_rerun_frame_streams_are_transactional(
    tmp_path,
    name,
    trajectory_values,
    box_text,
    velocity_values,
    messages,
):
    case_dir = tmp_path / name
    _write_rerun_case(
        case_dir,
        trajectory_values=trajectory_values,
        box_text=box_text,
        velocity_values=velocity_values,
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "spongeErrorBadFileFormat" in output
    for message in messages:
        assert message in output


def test_invalid_runtime_rerun_cell_is_a_simulation_breakdown(tmp_path):
    case_dir = tmp_path / "invalid_runtime_rerun_cell"
    _write_case(case_dir)
    trajectory_values = (
        1.0,
        2.0,
        3.0,
        2.0,
        2.0,
        3.0,
        1.0,
        2.0,
        3.0,
        2.0,
        2.0,
        3.0,
    )
    (case_dir / "trajectory.dat").write_bytes(
        struct.pack(f"={len(trajectory_values)}f", *trajectory_values)
    )
    (case_dir / "box.txt").write_text(
        "40 40 40 90 90 90\n0 40 40 90 90 90\n", encoding="utf-8"
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "rerun",
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "crd": "trajectory.dat",
        "box": "box.txt",
        "rerun_need_box_update": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "spongeErrorSimulationBrokenDown" in output
    assert "periodic_box_condition_information::Get_Cell" in output
    assert "box length x" in output
    assert "runtime periodic-cell rebuild" in output


def test_rerun_paths_preserve_spaces_and_hash_characters(tmp_path):
    case_dir = tmp_path / "rerun_paths_with_spaces"
    _write_case(case_dir)
    trajectory_name = "trajectory # one.dat"
    box_name = "box frames # one.txt"
    trajectory_values = (1.0, 2.0, 3.0, 2.0, 2.0, 3.0)
    (case_dir / trajectory_name).write_bytes(
        struct.pack(f"={len(trajectory_values)}f", *trajectory_values)
    )
    (case_dir / box_name).write_text("40 40 40 90 90 90\n", encoding="utf-8")
    settings = {
        "md_name": case_dir.name,
        "mode": "rerun",
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "crd": trajectory_name,
        "box": box_name,
        "rerun_need_box_update": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_nopbc_rerun_publishes_every_representable_reference_box_change(
    tmp_path,
):
    case_dir = tmp_path / "nopbc_rerun_small_reference_box_change"
    _write_case(case_dir)
    # Keep the printed CV within its fixed-width field while making the small
    # reciprocal-cell change visible at four decimal places.
    coordinate_x = 3.6e3
    trajectory_values = (
        coordinate_x,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        coordinate_x,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        coordinate_x,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    )
    (case_dir / "trajectory.dat").write_bytes(
        struct.pack(f"={len(trajectory_values)}f", *trajectory_values)
    )
    length_bits = struct.unpack("=I", struct.pack("=f", 40.0))[0]
    # Twenty-six ulps at this magnitude are still below the legacy 1e-4
    # threshold, but make the stale reciprocal scale externally observable.
    next_length = struct.unpack("=f", struct.pack("=I", length_bits + 26))[0]
    assert 0.0 < next_length - 40.0 < 1.0e-4
    (case_dir / "box.txt").write_text(
        "40 40 40 90 90 90\n"
        f"{next_length:.9g} 40 40 90 90 90\n"
        f"{next_length:.9g} 40 40 90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "cv.txt").write_text(
        "scaled\n"
        "{\n"
        "    CV_type = scaled_position_x\n"
        "    atom = 0\n"
        "}\n"
        "print\n"
        "{\n"
        "    CV = scaled\n"
        "}\n",
        encoding="utf-8",
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "rerun",
        "pbc": False,
        "cutoff": 8.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "crd": "trajectory.dat",
        "box": "box.txt",
        "cv_in_file": "cv.txt",
        "rerun_need_box_update": 1,
        "mdout": "mdout.txt",
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    lines = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    names = lines[0].split()
    scaled_index = names.index("scaled")
    scaled_values = [float(line.split()[scaled_index]) for line in lines[1:]]
    assert len(scaled_values) == 3
    assert scaled_values[0] == pytest.approx(coordinate_x / 40.0, abs=1.1e-4)
    assert scaled_values[1] == pytest.approx(
        coordinate_x / next_length, abs=1.1e-4
    )
    assert scaled_values[1] < scaled_values[0]
    # An unchanged frame must not replay the previous deformation matrix.
    assert scaled_values[2] == scaled_values[1]


def test_unrepresentable_runtime_box_update_is_a_simulation_breakdown(
    tmp_path,
):
    case_dir = tmp_path / "unrepresentable_runtime_box_update"
    _write_case(
        case_dir,
        overrides={"coordinate.txt": _coordinate_text((1, 1, 1), (90, 90, 90))},
    )
    trajectory_values = (
        0.1,
        0.2,
        0.3,
        0.2,
        0.2,
        0.3,
        0.1,
        0.2,
        0.3,
        0.2,
        0.2,
        0.3,
    )
    (case_dir / "trajectory.dat").write_bytes(
        struct.pack(f"={len(trajectory_values)}f", *trajectory_values)
    )
    # Both cells are individually valid. Their x-length ratio divided by the
    # rerun timestep cannot be represented by LTMatrix3's float-valued g, so
    # rerun must reject the deformation before mutating the active cell.
    (case_dir / "box.txt").write_text(
        "1 1 1 90 90 90\n1e37 1 1 90 90 90\n", encoding="utf-8"
    )
    settings = {
        "md_name": case_dir.name,
        "mode": "rerun",
        "cutoff": 0.1,
        "skin": 0.1,
        "rerun_need_box_update": 1,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "crd": "trajectory.dat",
        "box": "box.txt",
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "spongeErrorSimulationBrokenDown" in output
    assert "RERUN_information::Iteration" in output
    assert "deformation matrix" in output
    assert "finite zero-or-normal" in output


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("negative_count", "-1\n", ("negative atom count",)),
        (
            "count_over_int",
            "2147483648\n",
            ("atom count", "not a strict signed integer in range"),
        ),
        (
            "triple_count_overflow",
            "715827883\n",
            ("cannot safely represent all 3 * atom count",),
        ),
        (
            "large_truncated_count_is_streamed",
            "700000000\n",
            ("truncated", "mass entry 0"),
        ),
        (
            "bad_count_suffix",
            "2junk\n12\n1\n",
            ("atom count", "not a strict signed integer in range"),
        ),
        ("nan", "2\nnan\n1\n", ("mass entry 0", "strict finite decimal")),
        ("inf", "2\ninf\n1\n", ("mass entry 0", "strict finite decimal")),
        (
            "double_overflow",
            "2\n1e309\n1\n",
            ("mass entry 0", "finite double range"),
        ),
        (
            "float_overflow",
            "2\n3.5e38\n1\n",
            ("mass entry 0", "finite float range"),
        ),
        (
            "float_underflow",
            "2\n1e-50\n1\n",
            ("nonzero mass entry 0", "underflows the finite float range"),
        ),
        (
            "float_subnormal",
            "2\n1e-40\n1\n",
            ("mass entry 0", "subnormal float", "consistent FTZ behavior"),
        ),
        (
            "hex_real",
            "2\n0x1p0\n1\n",
            ("mass entry 0", "strict finite decimal"),
        ),
        ("truncated", "2\n12\n", ("truncated", "mass entry 1")),
        ("trailing", "2\n12\n1\nextra\n", ("trailing data", "extra")),
    ],
)
def test_mass_parser_rejects_invalid_input(tmp_path, name, contents, messages):
    case_dir = tmp_path / f"mass_{name}"
    _write_case(case_dir, overrides={"mass.txt": contents})

    _assert_rejected(case_dir, "mass.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("nan", "2\n0\nnan\n", ("charge entry 1", "strict finite decimal")),
        (
            "float_overflow",
            "2\n0\n-3.5e38\n",
            ("charge entry 1", "finite float range"),
        ),
        (
            "float_underflow",
            "2\n0\n-1e-50\n",
            ("nonzero charge entry 1", "underflows the finite float range"),
        ),
        (
            "float_subnormal",
            "2\n0\n-1e-40\n",
            ("charge entry 1", "subnormal float", "consistent FTZ behavior"),
        ),
        ("truncated", "2\n0\n", ("truncated", "charge entry 1")),
        ("trailing", "2\n0\n0\n3\n", ("trailing data",)),
    ],
)
def test_charge_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"charge_{name}"
    _write_case(case_dir, overrides={"charge.txt": contents})

    _assert_rejected(case_dir, "charge.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        (
            "bad_count_suffix",
            "2junk\n1 2 3\n2 2 3\n40 40 40\n90 90 90\n",
            ("header atom count", "not a strict signed integer in range"),
        ),
        (
            "long_header_with_payload",
            "2 0 " + " " * 700 + "1 2 3\n2 2 3\n40 40 40\n90 90 90\n",
            ("header has unexpected field 4",),
        ),
        (
            "header_trailing_field",
            "2 0 7 extra\n1 2 3\n2 2 3\n40 40 40\n90 90 90\n",
            ("header has unexpected field 4", "extra"),
        ),
        (
            "nonfinite_start_time",
            "2 inf\n1 2 3\n2 2 3\n40 40 40\n90 90 90\n",
            ("header start time", "strict finite decimal"),
        ),
        (
            "double_overflow_start_time",
            "2 1e309\n1 2 3\n2 2 3\n40 40 40\n90 90 90\n",
            ("header start time", "finite double range"),
        ),
        (
            "bad_step",
            "2 0 7junk\n1 2 3\n2 2 3\n40 40 40\n90 90 90\n",
            ("header step", "not a strict signed integer in range"),
        ),
        (
            "nonfinite_coordinate",
            "2\n1 nan 3\n2 2 3\n40 40 40\n90 90 90\n",
            ("coordinate y entry 0", "strict finite decimal"),
        ),
        (
            "coordinate_float_overflow",
            "2\n1 2 3\n2 3.5e38 3\n40 40 40\n90 90 90\n",
            ("coordinate y entry 1", "finite float range"),
        ),
        (
            "coordinate_float_underflow",
            "2\n1 2 3\n2 1e-50 3\n40 40 40\n90 90 90\n",
            ("nonzero coordinate y entry 1", "underflows"),
        ),
        (
            "coordinate_float_subnormal",
            "2\n1 2 3\n2 1e-40 3\n40 40 40\n90 90 90\n",
            (
                "coordinate y entry 1",
                "subnormal float",
                "consistent FTZ behavior",
            ),
        ),
        (
            "nonfinite_box",
            "2\n1 2 3\n2 2 3\n40 inf 40\n90 90 90\n",
            ("box length y", "strict finite decimal"),
        ),
        (
            "truncated",
            "2\n1 2 3\n2 2 3\n40 40 40\n90 90\n",
            ("truncated", "box angle gamma"),
        ),
        (
            "trailing",
            "2\n1 2 3\n2 2 3\n40 40 40\n90 90 90\nextra\n",
            ("trailing data", "extra"),
        ),
    ],
)
def test_coordinate_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"coordinate_{name}"
    _write_case(case_dir, overrides={"coordinate.txt": contents})

    _assert_rejected(case_dir, "coordinate.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        (
            "bad_header",
            "2junk\n0 0 0\n0 0 0\n",
            ("header atom count", "not a strict signed integer in range"),
        ),
        (
            "nonfinite_header_time",
            "2 nan\n0 0 0\n0 0 0\n",
            ("header start time", "strict finite decimal"),
        ),
        (
            "header_trailing_field",
            "2 0 7 extra\n0 0 0\n0 0 0\n",
            ("header has unexpected field 4",),
        ),
        (
            "nonfinite_velocity",
            "2\n0 0 0\n0 0 inf\n",
            ("velocity z entry 1", "strict finite decimal"),
        ),
        (
            "velocity_underflow",
            "2\n0 0 0\n0 0 1e-50\n",
            ("nonzero velocity z entry 1", "underflows"),
        ),
        (
            "velocity_subnormal",
            "2\n0 0 0\n0 0 1e-40\n",
            (
                "velocity z entry 1",
                "subnormal float",
                "consistent FTZ behavior",
            ),
        ),
        ("truncated", "2\n0 0 0\n0 0\n", ("truncated", "velocity z entry 1")),
        ("trailing", "2\n0 0 0\n0 0 0\nextra\n", ("trailing data",)),
    ],
)
def test_velocity_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"velocity_{name}"
    _write_case(case_dir, overrides={"velocity.txt": contents})

    _assert_rejected(case_dir, "velocity.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("negative_count", "2 -1\n", ("negative residue count",)),
        (
            "count_over_int",
            "2 2147483648\n",
            ("header residue count", "not a strict signed integer in range"),
        ),
        ("count_over_atom_count", "2 3\n", ("strictly positive",)),
        (
            "negative_length",
            "2 1\n-2\n",
            ("residue length entry 0", "strictly positive"),
        ),
        (
            "zero_length",
            "2 2\n0\n2\n",
            ("residue length entry 0", "strictly positive"),
        ),
        (
            "sum_too_small",
            "2 1\n1\n",
            ("residue-length sum 1", "does not equal atom count 2"),
        ),
        (
            "sum_too_large",
            "2 1\n3\n",
            ("residue-length sum 3", "does not equal atom count 2"),
        ),
        ("truncated", "2 2\n1\n", ("truncated", "residue length entry 1")),
        ("trailing", "2 1\n2\n1\n", ("trailing data",)),
    ],
)
def test_residue_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"residue_{name}"
    _write_case(case_dir, overrides={"residue.txt": contents})

    _assert_rejected(case_dir, "residue.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("negative_total", "2 -1\n", ("negative declared exclusion total",)),
        (
            "total_over_int",
            "2 2147483648\n",
            (
                "header declared exclusion total",
                "not a strict signed integer in range",
            ),
        ),
        (
            "negative_row_count",
            "2 0\n-1\n0\n",
            ("exclusion row 0", "negative exclusion count"),
        ),
        (
            "row_count_too_large",
            "2 1\n2 1 1\n0\n",
            ("exclusion row 0 count 2", "upper-triangle"),
        ),
        (
            "id_out_of_range",
            "2 1\n1 2\n0\n",
            ("neighbor entry 0 ID 2", "outside atom range"),
        ),
        ("self", "2 1\n1 0\n0\n", ("neighbor entry 0", "self exclusion")),
        (
            "reverse",
            "3 1\n0\n1 0\n0\n",
            ("neighbor entry 0", "reverse/lower-triangle"),
        ),
        (
            "duplicate",
            "3 2\n2 1 1\n0\n0\n",
            ("neighbor entry 1", "duplicates the previous ID"),
        ),
        (
            "decreasing",
            "4 2\n2 2 1\n0\n0\n0\n",
            ("neighbor entry 1", "not strictly increasing"),
        ),
        (
            "declared_total_mismatch",
            "3 2\n1 1\n0\n0\n",
            ("exclusion-count sum 1", "does not equal declared total 2"),
        ),
        ("truncated_row", "2 0\n0\n", ("truncated", "exclusion row 1")),
        (
            "row_trailing_id",
            "2 1\n0 1\n0\n",
            ("exclusion row 0 declares 0 IDs but contains 1",),
        ),
        ("trailing_row", "2 1\n1 1\n0\n0\n", ("trailing data",)),
    ],
)
def test_exclusion_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"exclude_{name}"
    atom_count = int(contents.split()[0])
    _write_case(
        case_dir,
        atom_count=atom_count,
        overrides={"exclude.txt": contents},
    )

    _assert_rejected(case_dir, "exclude.txt", *messages)


def _valid_lj_text(atom_count=2):
    return f"{atom_count} 1\n1.0\n1.0\n" + "0\n" * atom_count


def _valid_soft_lj_text(atom_count=2):
    return f"{atom_count} 1 1\n1.0\n1.0\n1.0\n1.0\n" + "0 0\n" * atom_count


def test_valid_native_lj_and_soft_lj_files_run(tmp_path):
    case_dir = tmp_path / "valid_lj_inputs"
    _write_case(
        case_dir,
        overrides={
            "LJ.txt": _valid_lj_text(),
            "LJ_soft.txt": _valid_soft_lj_text(),
            "division.txt": "2\n0\n1\n",
        },
        settings_overrides={
            "LJ_in_file": "LJ.txt",
            "LJ_soft_core_in_file": "LJ_soft.txt",
            "subsys_division_in_file": "division.txt",
            "lambda_lj": 0.5,
        },
    )

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_valid_native_nb14_derived_and_explicit_files_run(tmp_path):
    case_dir = tmp_path / "valid_nb14_inputs"
    _write_case(
        case_dir,
        atom_count=3,
        overrides={
            "LJ.txt": _valid_lj_text(atom_count=3),
            "nb14.txt": "1\n0 1 0.5 0.8\n",
            "nb14_extra.txt": "1\n1 2 1.0 1.0 0.5\n",
        },
        settings_overrides={
            "LJ_in_file": "LJ.txt",
            "nb14_in_file": "nb14.txt",
            "nb14_extra_in_file": "nb14_extra.txt",
        },
    )

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


@pytest.mark.parametrize(
    ("name", "files", "settings"),
    [
        (
            "explicit_only",
            {"nb14_extra.txt": "1\n0 1 1.0 1.0 0.5\n"},
            {"nb14_extra_in_file": "nb14_extra.txt"},
        ),
        (
            "zero_records",
            {"nb14.txt": "0\n", "nb14_extra.txt": "0\n"},
            {
                "nb14_in_file": "nb14.txt",
                "nb14_extra_in_file": "nb14_extra.txt",
            },
        ),
    ],
)
def test_valid_native_nb14_explicit_only_and_zero_record_files_run(
    tmp_path, name, files, settings
):
    case_dir = tmp_path / f"valid_nb14_{name}"
    _write_case(
        case_dir,
        overrides=files,
        settings_overrides=settings,
    )

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


@pytest.mark.parametrize(
    ("name", "derived", "explicit", "messages"),
    [
        ("negative_count", "-1\n", None, ("negative derived",)),
        (
            "negative_explicit_count",
            None,
            "-1\n",
            ("negative explicit",),
        ),
        (
            "count_above_int_range",
            "2147483648\n",
            None,
            ("derived interaction count", "strict signed integer in range"),
        ),
        (
            "atom_out_of_range",
            "1\n0 2 1 1\n",
            None,
            ("atom B 2", "outside [0, 2)"),
        ),
        (
            "self_pair",
            "1\n1 1 1 1\n",
            None,
            ("self 1-4 interaction",),
        ),
        (
            "duplicate_within_derived",
            "2\n0 1 1 1\n1 0 1 1\n",
            None,
            ("duplicates atom pair (0, 1)",),
        ),
        (
            "duplicate_within_explicit",
            None,
            "2\n0 1 1 1 1\n1 0 1 1 1\n",
            ("duplicates atom pair (0, 1)",),
        ),
        (
            "duplicate_across_files",
            "1\n0 1 1 1\n",
            "1\n1 0 1 1 1\n",
            ("duplicates atom pair (0, 1)",),
        ),
        (
            "scaled_overflow",
            "1\n0 1 3e38 1\n",
            None,
            ("scaled LJ A", "finite float range"),
        ),
        (
            "explicit_scaled_overflow",
            None,
            "1\n0 1 3e37 1 1\n",
            ("LJ A scaled by 12", "finite float range"),
        ),
        (
            "nonfinite_charge_scale",
            "1\n0 1 1 inf\n",
            None,
            ("charge scale", "strict finite decimal"),
        ),
        (
            "subnormal_lj_scale",
            "1\n0 1 1e-40 1\n",
            None,
            ("LJ scale", "subnormal float"),
        ),
        (
            "subnormal_explicit_charge_scale",
            None,
            "1\n0 1 1 1 1e-40\n",
            ("charge scale", "subnormal float"),
        ),
        (
            "trailing_data",
            "1\n0 1 1 1\ntrailing\n",
            None,
            ("trailing data", "trailing"),
        ),
    ],
)
def test_native_nb14_parser_rejects_invalid_input(
    tmp_path, name, derived, explicit, messages
):
    case_dir = tmp_path / f"nb14_{name}"
    overrides = {"LJ.txt": _valid_lj_text()}
    settings = {"LJ_in_file": "LJ.txt"}
    checked_filename = "nb14.txt"
    if derived is not None:
        overrides["nb14.txt"] = derived
        settings["nb14_in_file"] = "nb14.txt"
    if explicit is not None:
        overrides["nb14_extra.txt"] = explicit
        settings["nb14_extra_in_file"] = "nb14_extra.txt"
        if derived is None:
            checked_filename = "nb14_extra.txt"
    _write_case(
        case_dir,
        overrides=overrides,
        settings_overrides=settings,
    )

    _assert_rejected(case_dir, checked_filename, *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("negative_atom_count", "-1 1\n", ("negative atom count",)),
        ("negative_type_count", "2 -1\n", ("negative atom type count",)),
        (
            "triangular_overflow",
            "2 65536\n",
            ("unsupported triangular pair count",),
        ),
        (
            "atom_count_mismatch",
            _valid_lj_text(atom_count=3),
            ("previously loaded atom count 2",),
        ),
        (
            "nonfinite_parameter",
            "2 1\nnan\n1\n0\n0\n",
            ("LJ A entry 0", "strict finite decimal"),
        ),
        (
            "scaled_parameter_overflow",
            "2 1\n3e37\n1\n0\n0\n",
            ("LJ A entry 0 scaled by 12", "finite float range"),
        ),
        (
            "scaled_parameter_underflow",
            "2 1\n1e-50\n1\n0\n0\n",
            ("LJ A entry 0 scaled by 12", "not representable"),
        ),
        (
            "scaled_parameter_subnormal",
            "2 1\n5e-40\n1\n0\n0\n",
            ("LJ A entry 0 scaled by 12", "zero or normal float"),
        ),
        (
            "truncated_parameter",
            "2 1\n1\n",
            ("truncated", "LJ B entry 0"),
        ),
        (
            "negative_atom_type",
            "2 1\n1\n1\n-1\n0\n",
            ("atom LJ type entry 0 -1", "outside [0, 1)"),
        ),
        (
            "atom_type_out_of_range",
            "2 1\n1\n1\n0\n1\n",
            ("atom LJ type entry 1 1", "outside [0, 1)"),
        ),
        (
            "trailing_data",
            _valid_lj_text() + "extra\n",
            ("trailing data", "extra"),
        ),
    ],
)
def test_native_lj_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"lj_{name}"
    _write_case(
        case_dir,
        overrides={"LJ.txt": contents},
        settings_overrides={"LJ_in_file": "LJ.txt"},
    )

    _assert_rejected(case_dir, "LJ.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("negative_atom_count", "-1 1 1\n", ("negative atom count",)),
        (
            "negative_endpoint_type_count",
            "2 -1 1\n",
            ("negative endpoint A atom type count",),
        ),
        (
            "triangular_overflow",
            "2 1 65536\n",
            ("endpoint B atom type count", "unsupported triangular pair"),
        ),
        (
            "atom_count_mismatch",
            _valid_soft_lj_text(atom_count=3),
            ("previously loaded atom count 2",),
        ),
        (
            "nonfinite_parameter",
            "2 1 1\n1\n1\ninf\n1\n0 0\n0 0\n",
            ("endpoint B LJ A entry 0", "strict finite decimal"),
        ),
        (
            "scaled_parameter_overflow",
            "2 1 1\n1\n1\n3e37\n1\n0 0\n0 0\n",
            ("endpoint B LJ A entry 0 scaled by 12", "finite float range"),
        ),
        (
            "scaled_parameter_subnormal",
            "2 1 1\n1\n1\n5e-40\n1\n0 0\n0 0\n",
            ("endpoint B LJ A entry 0 scaled by 12", "zero or normal float"),
        ),
        (
            "truncated_atom_type",
            "2 1 1\n1\n1\n1\n1\n0 0\n0\n",
            ("truncated", "endpoint B atom LJ type entry 1"),
        ),
        (
            "endpoint_a_type_out_of_range",
            "2 1 1\n1\n1\n1\n1\n1 0\n0 0\n",
            ("endpoint A atom LJ type entry 0 1", "outside [0, 1)"),
        ),
        (
            "endpoint_b_type_out_of_range",
            "2 1 1\n1\n1\n1\n1\n0 0\n0 -1\n",
            ("endpoint B atom LJ type entry 1 -1", "outside [0, 1)"),
        ),
        (
            "trailing_data",
            _valid_soft_lj_text() + "extra\n",
            ("trailing data", "extra"),
        ),
    ],
)
def test_native_soft_lj_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"soft_lj_{name}"
    _write_case(
        case_dir,
        overrides={"LJ_soft.txt": contents},
        settings_overrides={
            "LJ_soft_core_in_file": "LJ_soft.txt",
            "lambda_lj": 0.5,
        },
    )

    _assert_rejected(case_dir, "LJ_soft.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("negative_count", "-1\n", ("negative atom count",)),
        ("count_mismatch", "1\n0\n", ("LJ_soft_core_in_file atom count 2",)),
        ("truncated", "2\n0\n", ("truncated", "subsystem mask entry 1")),
        ("trailing", "2\n0\n0\nextra\n", ("trailing data", "extra")),
    ],
)
def test_native_soft_lj_subsystem_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"soft_lj_division_{name}"
    _write_case(
        case_dir,
        overrides={
            "LJ_soft.txt": _valid_soft_lj_text(),
            "division.txt": contents,
        },
        settings_overrides={
            "LJ_soft_core_in_file": "LJ_soft.txt",
            "subsys_division_in_file": "division.txt",
            "lambda_lj": 0.5,
        },
    )

    _assert_rejected(case_dir, "division.txt", *messages)


def _valid_gb_text(atom_count=2):
    return f"{atom_count}\n" + "1.5 0.8\n" * atom_count


def test_valid_native_gb_file_runs_without_periodic_neighbor_state(tmp_path):
    case_dir = tmp_path / "valid_native_gb"
    _write_case(
        case_dir,
        overrides={"gb.txt": _valid_gb_text()},
        settings_overrides={"pbc": False, "gb_in_file": "gb.txt"},
    )

    result = _run_case(case_dir)
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


@pytest.mark.parametrize(
    ("name", "contents", "messages"),
    [
        ("negative_atom_count", "-1\n", ("negative atom count",)),
        (
            "atom_count_over_int",
            "2147483648\n",
            ("atom count", "not a strict signed integer in range"),
        ),
        (
            "atom_count_component_overflow",
            "715827883\n",
            ("cannot safely represent all 3 * atom count",),
        ),
        (
            "huge_count_mismatch_before_allocation",
            "700000000\n",
            ("previously loaded atom count 2",),
        ),
        (
            "atom_count_mismatch",
            _valid_gb_text(atom_count=3),
            ("previously loaded atom count 2",),
        ),
        (
            "nonfinite_radius",
            "2\nnan 0.8\n1.5 0.8\n",
            ("GB radius entry 0", "strict finite decimal"),
        ),
        (
            "float_overflow_scale",
            "2\n1.5 3.5e38\n1.5 0.8\n",
            ("GB scale factor entry 0", "finite float range"),
        ),
        (
            "float_underflow_radius",
            "2\n1e-50 0.8\n1.5 0.8\n",
            ("nonzero GB radius entry 0", "underflows"),
        ),
        (
            "subnormal_scale",
            "2\n1.5 1e-40\n1.5 0.8\n",
            ("GB scale factor entry 0", "subnormal float"),
        ),
        (
            "truncated_scale",
            "2\n1.5 0.8\n1.5\n",
            ("truncated", "GB scale factor entry 1"),
        ),
        (
            "trailing_data",
            _valid_gb_text() + "extra\n",
            ("trailing data", "extra"),
        ),
    ],
)
def test_native_gb_parser_rejects_invalid_input(
    tmp_path, name, contents, messages
):
    case_dir = tmp_path / f"gb_{name}"
    _write_case(
        case_dir,
        overrides={"gb.txt": contents},
        settings_overrides={"pbc": False, "gb_in_file": "gb.txt"},
    )

    _assert_rejected(case_dir, "gb.txt", *messages)


@pytest.mark.parametrize(
    ("name", "contents", "settings", "messages"),
    [
        (
            "zero_radius",
            "2\n0 0.8\n1.5 0.8\n",
            {},
            ("GB atom 0", "invalid radius/scale/offset"),
        ),
        (
            "negative_scale",
            "2\n1.5 -0.8\n1.5 0.8\n",
            {},
            ("GB atom 0", "invalid radius/scale/offset"),
        ),
        (
            "offset_consumes_radius",
            _valid_gb_text(),
            {"gb_radii_offset": 1.5},
            ("GB atom 0", "invalid radius/scale/offset"),
        ),
        (
            "zero_dielectric",
            _valid_gb_text(),
            {"gb_epsilon": 0},
            ("gb.epsilon", "positive normal float"),
        ),
        (
            "cutoff_square_overflow",
            _valid_gb_text(),
            {"gb_radii_cutoff": "2e19"},
            ("gb.radii_cutoff", "square outside the finite float range"),
        ),
        (
            "cutoff_square_subnormal",
            _valid_gb_text(),
            {"gb_radii_cutoff": "1e-20"},
            ("gb.radii_cutoff", "square that is not a normal float"),
        ),
    ],
)
def test_native_gb_initialization_rejects_invalid_semantics(
    tmp_path, name, contents, settings, messages
):
    case_dir = tmp_path / f"gb_semantics_{name}"
    _write_case(
        case_dir,
        overrides={"gb.txt": contents},
        settings_overrides={
            "pbc": False,
            "gb_in_file": "gb.txt",
            **settings,
        },
    )

    result = _run_case(case_dir)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    for message in messages:
        assert message in output
