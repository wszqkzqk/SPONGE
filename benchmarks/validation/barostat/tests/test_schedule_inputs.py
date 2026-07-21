import pytest

from benchmarks.utils import Outputer, Runner
from benchmarks.validation.utils import parse_mdout_rows


def _write_mdin(case_dir, extra_lines):
    base = (
        'md_name = "validation barostat schedule inputs"\n'
        'mode = "npt"\n'
        "step_limit = 20\n"
        "dt = 0.002\n"
        "cutoff = 8.0\n"
        'thermostat = "middle_langevin"\n'
        "thermostat_tau = 0.1\n"
        "thermostat_seed = 2026\n"
        "target_temperature = 300.0\n"
        "target_pressure = 1.0\n"
        'barostat = "berendsen_barostat"\n'
        "barostat_tau = 0.1\n"
        "barostat_update_interval = 10\n"
        'default_in_file_prefix = "tip3p"\n'
        'constrain_mode = "SHAKE"\n'
        "print_zeroth_frame = 1\n"
        "write_mdout_interval = 10\n"
        "write_information_interval = 10\n"
    )
    (case_dir / "mdin.spg.toml").write_text(
        base + "\n".join(extra_lines) + "\n"
    )


def _run(case_dir, mpi_np):
    return Runner.run_sponge(case_dir, timeout=600, mpi_np=mpi_np)


def test_schedule_inline_steps_object_array_is_supported(
    statics_path, outputs_path, mpi_np
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p_water",
        mpi_np=mpi_np,
        run_name="schedule_inline_object_array",
    )

    _write_mdin(
        case_dir,
        [
            'target_temperature_schedule_mode = "linear"',
            "target_temperature_schedule_steps = [{step = 0, value = 330.0}, {step = 20, value = 350.0}]",
            'target_pressure_schedule_mode = "step"',
            "target_pressure_schedule_steps = [{step = 0, value = 25.0}, {step = 10, value = 50.0}]",
        ],
    )

    output = _run(case_dir, mpi_np)
    assert "target temperature is 330.00 K" in output
    assert "target pressure is 25.00 bar" in output

    rows = parse_mdout_rows(
        case_dir / "mdout.txt", columns=("step",), int_columns=("step",)
    )
    assert rows


def test_default_prefix_detects_temp_pres_toml_schedule_files(
    statics_path, outputs_path, mpi_np
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p_water",
        mpi_np=mpi_np,
        run_name="schedule_default_prefix_toml",
    )

    (case_dir / "tip3p.temp.spg.toml").write_text(
        'mode = "linear"\nsteps = [{step = 0, value = 330.0}, {step = 20, value = 300.0}]\n'
    )
    (case_dir / "tip3p.pres.spg.toml").write_text(
        'mode = "step"\nsteps = [{step = 0, value = 1.0}, {step = 10, value = 200.0}]\n'
    )

    _write_mdin(case_dir, [])
    _run(case_dir, mpi_np)

    rows = parse_mdout_rows(
        case_dir / "mdout.txt", columns=("step",), int_columns=("step",)
    )
    assert rows


def test_explicit_schedule_file_path_preserves_whitespace(
    statics_path, outputs_path, mpi_np
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p_water",
        mpi_np=mpi_np,
        run_name="schedule_explicit_whitespace_path",
    )
    schedule_name = "temperature schedule.spg.toml"
    (case_dir / schedule_name).write_text(
        'mode = "linear"\n'
        "steps = [{step = 0, value = 335.0}, {step = 20, value = 350.0}]\n",
        encoding="utf-8",
    )
    _write_mdin(
        case_dir,
        [f'target_temperature_schedule_file = "{schedule_name}"'],
    )

    output = _run(case_dir, mpi_np)
    assert "target temperature is 335.00 K" in output


def test_explicit_txt_schedule_file_is_rejected(
    statics_path, outputs_path, mpi_np
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p_water",
        mpi_np=mpi_np,
        run_name="schedule_txt_file_rejected",
    )

    (case_dir / "legacy_temp.txt").write_text("0 300.0\n20 350.0\n")

    _write_mdin(
        case_dir,
        [
            'target_temperature_schedule_file = "legacy_temp.txt"',
            'target_temperature_schedule_mode = "linear"',
        ],
    )

    with pytest.raises(RuntimeError):
        _run(case_dir, mpi_np)


def test_pressure_barostat_rejects_nonpositive_update_interval(
    statics_path, outputs_path, mpi_np, capsys
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p_water",
        mpi_np=mpi_np,
        run_name="invalid_barostat_update_interval",
    )
    _write_mdin(case_dir, [])
    mdin = case_dir / "mdin.spg.toml"
    contents = mdin.read_text(encoding="utf-8")
    mdin.write_text(
        contents.replace(
            "barostat_update_interval = 10",
            "barostat_update_interval = 0",
        ),
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError):
        _run(case_dir, mpi_np)
    assert (
        "barostat_update_interval must be positive" in capsys.readouterr().out
    )


@pytest.mark.parametrize(
    ("run_tag", "original", "replacement", "diagnostic"),
    [
        (
            "temperature_inf",
            "target_temperature = 300.0",
            "target_temperature = 1e39",
            "target_temperature must be finite, positive, and normal",
        ),
        (
            "pressure_inf",
            "target_pressure = 1.0",
            "target_pressure = 1e39",
            "target_pressure must be finite and zero or normal",
        ),
        (
            "pressure_subnormal",
            "target_pressure = 1.0",
            "target_pressure = 1e-40",
            "target_pressure must be finite and zero or normal",
        ),
        (
            "pressure_conversion_subnormal",
            "target_pressure = 1.0",
            "target_pressure = 1e-37",
            "target_pressure conversion produces a subnormal or underflowed float",
        ),
    ],
)
def test_nonrepresentable_direct_targets_are_rejected(
    statics_path,
    outputs_path,
    mpi_np,
    capsys,
    run_tag,
    original,
    replacement,
    diagnostic,
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p_water",
        mpi_np=mpi_np,
        run_name=f"invalid_{run_tag}",
    )
    _write_mdin(case_dir, [])
    mdin = case_dir / "mdin.spg.toml"
    contents = mdin.read_text(encoding="utf-8")
    assert original in contents
    mdin.write_text(contents.replace(original, replacement), encoding="utf-8")

    with pytest.raises(RuntimeError):
        _run(case_dir, mpi_np)
    assert diagnostic in capsys.readouterr().out


@pytest.mark.parametrize(
    ("run_tag", "schedule_lines", "diagnostic"),
    [
        (
            "temperature_inf",
            [
                'target_temperature_schedule_mode = "step"',
                "target_temperature_schedule_steps = [{step = 0, value = 1e39}]",
            ],
            "temperature values must be finite, positive, and normal",
        ),
        (
            "pressure_subnormal",
            [
                'target_pressure_schedule_mode = "step"',
                "target_pressure_schedule_steps = [{step = 0, value = 1e-40}]",
            ],
            "pressure values must be finite and zero or normal",
        ),
        (
            "pressure_conversion_subnormal",
            [
                'target_pressure_schedule_mode = "step"',
                "target_pressure_schedule_steps = [{step = 0, value = 1e-37}]",
            ],
            "pressure conversion produces a subnormal or underflowed float",
        ),
    ],
)
def test_nonrepresentable_schedule_values_are_rejected(
    statics_path,
    outputs_path,
    mpi_np,
    capsys,
    run_tag,
    schedule_lines,
    diagnostic,
):
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p_water",
        mpi_np=mpi_np,
        run_name=f"invalid_schedule_{run_tag}",
    )
    _write_mdin(case_dir, schedule_lines)

    with pytest.raises(RuntimeError):
        _run(case_dir, mpi_np)
    assert diagnostic in capsys.readouterr().out
