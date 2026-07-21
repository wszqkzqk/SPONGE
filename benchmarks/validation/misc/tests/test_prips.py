import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

from benchmarks.utils import Extractor, Outputer, Runner


def _prips_plugin_path():
    configured_path = os.environ.get("PRIPS_PLUGIN_PATH")
    if configured_path:
        plugin_path = Path(configured_path).expanduser().resolve()
        if not plugin_path.is_file():
            raise RuntimeError(
                f"PRIPS_PLUGIN_PATH does not name a plugin file: {plugin_path}"
            )
        return plugin_path

    result = subprocess.run(
        [sys.executable, "-m", "prips"],
        capture_output=True,
        text=True,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            "failed to resolve prips plugin path from `python -m prips`\n"
            f"stdout:\n{result.stdout}\n"
            f"stderr:\n{result.stderr}"
        )

    output = result.stdout + "\n" + result.stderr
    match = re.search(r"Plugin Path:\s*(.+)", output)
    if match is None:
        raise Exception("failed to resolve Plugin Path from `python -m prips`")

    plugin_path = Path(match.group(1).strip())
    if not plugin_path.exists():
        raise Exception(f"prips plugin does not exist: {plugin_path}")
    return plugin_path


def _run_sponge_raw(case_dir, *, mpi_np=None, timeout=120):
    command = [
        os.environ.get("SPONGE_BIN", "SPONGE"),
        "-mdin",
        "mdin.spg.toml",
    ]
    if mpi_np is not None:
        command = ["mpirun", "--oversubscribe", "-np", str(mpi_np)] + command
    return subprocess.run(
        command,
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
    )


def _write_prips_script(case_dir, backend, *, mc_safe=False):
    capability_declaration = ""
    complete_buffer_writes = ""
    if mc_safe:
        capability_declaration = (
            "SPONGE_FORCE_CAPABILITIES = (\n"
            "    Sponge.FORCE_ENERGY_COMPLETE\n"
            "    | Sponge.FORCE_VIRIAL_COMPLETE\n"
            "    | Sponge.FORCE_PURE\n"
            ")\n\n"
        )
        complete_buffer_writes = (
            "    if Sponge.force_evaluation.needs_energy:\n"
            "        if Sponge.backend_name == 'jax':\n"
            "            energy_result = Sponge.dd.energy.at[0].add(0.0)\n"
            "        else:\n"
            "            Sponge.dd.energy[0] += 0.0\n"
            "    if Sponge.force_evaluation.needs_virial:\n"
            "        if Sponge.backend_name == 'jax':\n"
            "            virial_result = Sponge.dd.virial.at[0, 0].add(0.0)\n"
            "        else:\n"
            "            Sponge.dd.virial[0, 0] += 0.0\n"
        )
    script = (
        "from prips import Sponge\n\n"
        + capability_declaration
        + (
            f"_requested_backend = {backend!r}\n"
            "if _requested_backend == 'auto':\n"
            "    _requested_backend = 'numpy' if Sponge._backend == 1 else 'pytorch'\n"
            "Sponge.set_backend(_requested_backend)\n"
            "_force_delta = 1.25\n"
            "\n"
            "with open('prips_hook.log', 'w', encoding='utf-8') as f:\n"
            "    f.write('[prips]\\n')\n"
            "    f.write(f'backend device={Sponge._backend} name={Sponge.backend_name}\\n')\n"
            "\n"
            "def After_Initial():\n"
            "    crd = Sponge.md_info.crd\n"
            "    with open('prips_hook.log', 'a', encoding='utf-8') as f:\n"
            "        f.write(\n"
            "            f'after_init atom_numbers={Sponge.md_info.atom_numbers} '\n"
            "            f'neighbor_max={Sponge.neighbor_list.max_neighbor_numbers} '\n"
            "            f'rank={Sponge.controller.MPI_rank} '\n"
            "            f'coord[0,0]={float(crd[0, 0]):.6f}\\n'\n"
            "        )\n"
            "\n"
            "def Calculate_Force():\n"
            "    energy_result = None\n"
            "    virial_result = None\n"
        )
        + complete_buffer_writes
        + (
            "    frc = Sponge.dd.frc\n"
            "    step = Sponge.md_info.sys.steps\n"
            "    before = float(frc[0, 0])\n"
            "    updated_frc = frc\n"
            "    if step == 1:\n"
            "        if Sponge.backend_name == 'jax':\n"
            "            updated_frc = frc.at[0, 0].add(_force_delta)\n"
            "        else:\n"
            "            frc[0, 0] += _force_delta\n"
            "    after = float(updated_frc[0, 0])\n"
            "    with open('prips_hook.log', 'a', encoding='utf-8') as f:\n"
            "        f.write(\n"
            "            f'calculate_force step={step} '\n"
            "            f'commit_sampling_state={int(Sponge.force_evaluation.commits_sampling_state)} '\n"
            "            f'exact_state={int(Sponge.force_evaluation.is_exact)} '\n"
            "            f'local_atom_numbers={Sponge.dd.atom_numbers} '\n"
            "            f'force[0,0]_before={before:.6f} '\n"
            "            f'force[0,0]_after={after:.6f} '\n"
            "            f'force[0,0]_delta={after - before:.6f}\\n'\n"
            "        )\n"
            "    if Sponge.backend_name == 'jax':\n"
            "        return Sponge.force_result(\n"
            "            updated_frc, energy=energy_result, virial=virial_result\n"
            "        )\n"
            "\n"
            "def Mdout_Print():\n"
            "    frc = Sponge.dd.frc\n"
            "    with open('prips_hook.log', 'a', encoding='utf-8') as f:\n"
            "        f.write(\n"
            "            f'mdout_force step={Sponge.md_info.sys.steps} '\n"
            "            f'force[0,0]_final={float(frc[0, 0]):.6f}\\n'\n"
            "        )\n"
        )
    )
    (Path(case_dir) / "prips_test.py").write_text(script, encoding="utf-8")


def _write_prips_mdin(case_dir, plugin_path=None, *, step_limit=1):
    mdin = (
        'md_name = "validation tip3p prips"\n'
        'mode = "nvt"\n'
        f"step_limit = {step_limit}\n"
        "dt = 0.0\n"
        "cutoff = 8.0\n"
        'default_in_file_prefix = "tip3p"\n'
        "print_zeroth_frame = 1\n"
        "write_mdout_interval = 1\n"
        "write_information_interval = 1\n"
        'frc = "frc.dat"\n'
        'thermostat = "middle_langevin"\n'
        "thermostat_tau = 0.01\n"
        "thermostat_seed = 2026\n"
        "target_temperature = 300.0\n"
        "hard_wall_z_low = 5.0\n"
        "hard_wall_z_high = 30.0\n"
    )
    if plugin_path is not None:
        mdin += f'plugin = "{Path(plugin_path).as_posix()}"\n'
        mdin += 'py = "prips_test.py"\n'
    Path(case_dir, "mdin.spg.toml").write_text(mdin, encoding="utf-8")


def _write_prips_mc_mdin(case_dir, plugin_path):
    mdin = (
        'md_name = "validation tip3p prips mc"\n'
        'mode = "npt"\n'
        "step_limit = 0\n"
        "dt = 1e-8\n"
        "cutoff = 8.0\n"
        'default_in_file_prefix = "tip3p"\n'
        "print_zeroth_frame = 1\n"
        "write_mdout_interval = 1\n"
        "write_information_interval = 1\n"
        "write_restart_file_interval = 0\n"
        'thermostat = "middle_langevin"\n'
        "thermostat_tau = 0.01\n"
        "thermostat_seed = 2026\n"
        "target_temperature = 300.0\n"
        "target_pressure = 1.0\n"
        'barostat = "monte_carlo_barostat"\n'
        "monte_carlo_barostat_initial_ratio = 0.0\n"
        "monte_carlo_barostat_update_interval = 1\n"
        "monte_carlo_barostat_check_interval = 100\n"
        'monte_carlo_barostat_couple_dimension = "XYZ"\n'
        f'plugin = "{Path(plugin_path).as_posix()}"\n'
        'py = "prips_test.py"\n'
    )
    Path(case_dir, "mdin.spg.toml").write_text(mdin, encoding="utf-8")


def _write_prips_neighbor_list_script(case_dir):
    script = (
        "from prips import Sponge\n"
        "\n"
        "SPONGE_FORCE_CAPABILITIES = (\n"
        "    Sponge.FORCE_ENERGY_COMPLETE\n"
        "    | Sponge.FORCE_VIRIAL_COMPLETE\n"
        "    | Sponge.FORCE_PURE\n"
        ")\n"
        "\n"
        "Sponge.set_backend('numpy')\n"
        "\n"
        "with open('neighbor_context.log', 'w', encoding='utf-8'):\n"
        "    pass\n"
        "\n"
        "def Calculate_Force():\n"
        "    if Sponge.force_evaluation.needs_energy:\n"
        "        Sponge.dd.energy[0] += 0.0\n"
        "    if Sponge.force_evaluation.needs_virial:\n"
        "        Sponge.dd.virial[0, 0] += 0.0\n"
        "    with open('neighbor_context.log', 'a', encoding='utf-8') as f:\n"
        "        f.write(\n"
        "            f'step={Sponge.md_info.sys.steps} '\n"
        "            f'commit={int(Sponge.force_evaluation.commits_sampling_state)} '\n"
        "            f'exact={int(Sponge.force_evaluation.is_exact)} '\n"
        "            f'neighbors={sum(Sponge.neighbor_list.number)}\\n'\n"
        "        )\n"
        "\n"
        "def Mdout_Print():\n"
        "    if Sponge.md_info.sys.steps != 0:\n"
        "        return\n"
        "    for local_index, atom_id in enumerate(Sponge.dd.atom_local):\n"
        "        if int(atom_id) == 0:\n"
        "            Sponge.dd.crd[local_index] = (9.0, 10.0, 10.0)\n"
        "        elif int(atom_id) == 1:\n"
        "            Sponge.dd.crd[local_index] = (11.0, 10.0, 10.0)\n"
    )
    (Path(case_dir) / "prips_neighbor_test.py").write_text(
        script, encoding="utf-8"
    )


def _write_prips_two_atom_inputs(case_dir):
    (Path(case_dir) / "mass.txt").write_text("2\n12\n12\n", encoding="utf-8")
    (Path(case_dir) / "charge.txt").write_text("2\n0\n0\n", encoding="utf-8")
    (Path(case_dir) / "coordinate.txt").write_text(
        "2\n5 10 10\n12 10 10\n20 20 20\n90 90 90\n", encoding="utf-8"
    )
    (Path(case_dir) / "velocity.txt").write_text(
        "2\n0 0 0\n0 0 0\n", encoding="utf-8"
    )
    (Path(case_dir) / "lj.txt").write_text(
        "2 1\n1000000\n0\n0\n0\n", encoding="utf-8"
    )


def _write_prips_neighbor_list_mdin(
    case_dir,
    plugin_path,
    *,
    py_name="prips_neighbor_test.py",
    step_limit=1,
    initial_ratio=0.0,
    print_pressure=False,
):
    _write_prips_two_atom_inputs(case_dir)
    settings = {
        "md_name": "validation prips exact neighbor list",
        "mode": "npt",
        "step_limit": step_limit,
        "dt": 1.0e-8,
        "cutoff": 4.0,
        "skin": 0.1,
        "PM.MPI_size": 0,
        "neighbor_list.refresh_interval": 100,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "velocity_in_file": "velocity.txt",
        "LJ_in_file": "lj.txt",
        "thermostat": "berendsen_thermostat",
        "thermostat_tau": 1.0,
        "barostat": "monte_carlo_barostat",
        "target_temperature": 300.0,
        "target_pressure": 1.0,
        "monte_carlo_barostat_initial_ratio": initial_ratio,
        "monte_carlo_barostat_update_interval": 1,
        "monte_carlo_barostat_check_interval": 100,
        "monte_carlo_barostat_couple_dimension": "XYZ",
        "plugin": Path(plugin_path).as_posix(),
        "py": py_name,
        "mdout": "mdout.txt",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
        "print_pressure": print_pressure,
    }
    (Path(case_dir) / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )


def _write_prips_non_mc_mdin(case_dir, plugin_path, *, py_name=None):
    _write_prips_two_atom_inputs(case_dir)
    settings = {
        "md_name": "validation prips initialization contract",
        "mode": "nvt",
        "step_limit": 0,
        "dt": 1.0e-8,
        "cutoff": 4.0,
        "PM.MPI_size": 0,
        "mass_in_file": "mass.txt",
        "charge_in_file": "charge.txt",
        "coordinate_in_file": "coordinate.txt",
        "velocity_in_file": "velocity.txt",
        "LJ_in_file": "lj.txt",
        "thermostat": "berendsen_thermostat",
        "thermostat_tau": 1.0,
        "target_temperature": 300.0,
        "plugin": Path(plugin_path).as_posix(),
        "write_trajectory_interval": 0,
        "write_restart_file_interval": 0,
    }
    if py_name is not None:
        settings["py"] = py_name
    (Path(case_dir) / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n",
        encoding="utf-8",
    )


def _write_prips_contract_script(case_dir, capability_expression):
    (Path(case_dir) / "prips_contract_test.py").write_text(
        "from prips import Sponge\n\n"
        f"SPONGE_FORCE_CAPABILITIES = {capability_expression}\n\n"
        "def Calculate_Force():\n"
        "    pass\n",
        encoding="utf-8",
    )


def _write_prips_transaction_script(case_dir, energy_sign):
    script = f"""from prips import Sponge

SPONGE_FORCE_CAPABILITIES = (
    Sponge.FORCE_ENERGY_COMPLETE
    | Sponge.FORCE_VIRIAL_COMPLETE
    | Sponge.FORCE_TRANSACTIONAL
)

Sponge.set_backend("numpy")
_reference = Sponge.md_info.crd.copy()
_energy_sign = {float(energy_sign)!r}
_history = 0
_pending = None
_evaluation = 0

with open("transaction.log", "w", encoding="utf-8"):
    pass

def _log(message):
    with open("transaction.log", "a", encoding="utf-8") as stream:
        stream.write(message + "\\n")

def Begin_Force_Transaction():
    global _pending
    assert _pending is None
    _pending = _history
    _log(f"begin history={{_history}}")

def Commit_Force_Transaction():
    global _history, _pending
    assert _pending is not None
    _history = _pending
    _log(f"commit history={{_history}}")
    _pending = None

def Rollback_Force_Transaction():
    global _pending
    assert _pending is not None
    _log(f"rollback pending={{_pending}} history={{_history}}")
    _pending = None

def Calculate_Force():
    global _pending, _evaluation
    assert _pending is not None
    _pending += 1
    _evaluation += 1

    owned = Sponge.dd.atom_numbers
    atom_ids = Sponge.dd.atom_local[:owned]
    displacement = Sponge.dd.crd[:owned] - _reference[atom_ids]
    spring_constant = _energy_sign * 1.0e6
    contribution = spring_constant * float((displacement * displacement).sum())
    force_delta = -2.0 * spring_constant * displacement
    Sponge.dd.frc[:owned] += force_delta

    energy_before = float(Sponge.dd.energy[0])
    if Sponge.force_evaluation.needs_energy:
        Sponge.dd.energy[0] += contribution
    energy_after = float(Sponge.dd.energy[0])

    virial_before = float(Sponge.dd.virial[0, 0])
    if Sponge.force_evaluation.needs_virial:
        coordinate = Sponge.dd.crd[:owned]
        Sponge.dd.virial[:owned, 0] += force_delta[:, 0] * coordinate[:, 0]
        Sponge.dd.virial[:owned, 1] += (
            force_delta[:, 0] * coordinate[:, 1]
            + force_delta[:, 1] * coordinate[:, 0]
        )
        Sponge.dd.virial[:owned, 2] += force_delta[:, 1] * coordinate[:, 1]
        Sponge.dd.virial[:owned, 3] += (
            force_delta[:, 0] * coordinate[:, 2]
            + force_delta[:, 2] * coordinate[:, 0]
        )
        Sponge.dd.virial[:owned, 4] += (
            force_delta[:, 1] * coordinate[:, 2]
            + force_delta[:, 2] * coordinate[:, 1]
        )
        Sponge.dd.virial[:owned, 5] += force_delta[:, 2] * coordinate[:, 2]
    virial_after = float(Sponge.dd.virial[0, 0])

    _log(
        f"eval index={{_evaluation}} "
        f"commit={{int(Sponge.force_evaluation.commits_sampling_state)}} "
        f"exact={{int(Sponge.force_evaluation.is_exact)}} "
        f"needs_energy={{int(Sponge.force_evaluation.needs_energy)}} "
        f"needs_virial={{int(Sponge.force_evaluation.needs_virial)}} "
        f"energy_shape={{Sponge.dd.energy.shape}} "
        f"virial_shape={{Sponge.dd.virial.shape}} "
        f"contribution={{contribution:.12g}} "
        f"energy_delta={{energy_after - energy_before:.12g}} "
        f"virial_delta={{virial_after - virial_before:.12g}}"
    )

def Mdout_Print():
    owned_ids = Sponge.dd.atom_local[:Sponge.dd.atom_numbers]
    local_zero = next(i for i, atom_id in enumerate(owned_ids) if int(atom_id) == 0)
    _log(
        f"final history={{_history}} pending={{_pending}} "
        f"atom0_x={{float(Sponge.dd.crd[local_zero, 0]):.12g}} "
        f"virial00={{float(Sponge.dd.virial[local_zero, 0]):.12g}}"
    )
"""
    (Path(case_dir) / "prips_transaction_test.py").write_text(
        script, encoding="utf-8"
    )


def _write_prips_jax_writeback_script(case_dir):
    script = """from prips import Sponge

SPONGE_FORCE_CAPABILITIES = (
    Sponge.FORCE_ENERGY_COMPLETE
    | Sponge.FORCE_VIRIAL_COMPLETE
    | Sponge.FORCE_PURE
)

Sponge.set_backend("jax")

FORCE_DELTA = 1.25
ENERGY_DELTA = 2.5
VIRIAL_DELTA = 3.75

with open("jax_writeback.log", "w", encoding="utf-8"):
    pass

def _log(message):
    with open("jax_writeback.log", "a", encoding="utf-8") as stream:
        stream.write(message + "\\n")

def Calculate_Force():
    force_before = float(Sponge.dd.frc[0, 0])
    energy_before = float(Sponge.dd.energy[0])
    virial_before = float(Sponge.dd.virial[0, 0])

    force = Sponge.dd.frc.at[0, 0].add(FORCE_DELTA)
    energy = Sponge.dd.energy
    virial = Sponge.dd.virial
    if Sponge.force_evaluation.needs_energy:
        energy = energy.at[0].add(ENERGY_DELTA)
    if Sponge.force_evaluation.needs_virial:
        virial = virial.at[0, 0].add(VIRIAL_DELTA)

    _log(
        f"eval force_before={force_before:.9g} "
        f"force_after={float(force[0, 0]):.9g} "
        f"energy_before={energy_before:.9g} "
        f"energy_after={float(energy[0]):.9g} "
        f"virial_before={virial_before:.9g} "
        f"virial_after={float(virial[0, 0]):.9g} "
        f"needs_energy={int(Sponge.force_evaluation.needs_energy)} "
        f"needs_virial={int(Sponge.force_evaluation.needs_virial)}"
    )
    return Sponge.force_result(force, energy=energy, virial=virial)

def Mdout_Print():
    _log(
        f"final force={float(Sponge.dd.frc[0, 0]):.9g} "
        f"energy={float(Sponge.dd.energy[0]):.9g} "
        f"virial={float(Sponge.dd.virial[0, 0]):.9g}"
    )
"""
    (Path(case_dir) / "prips_jax_writeback_test.py").write_text(
        script, encoding="utf-8"
    )


def _write_prips_jax_missing_result_script(case_dir):
    (Path(case_dir) / "prips_jax_missing_result_test.py").write_text(
        """from prips import Sponge

SPONGE_FORCE_CAPABILITIES = (
    Sponge.FORCE_ENERGY_COMPLETE
    | Sponge.FORCE_VIRIAL_COMPLETE
    | Sponge.FORCE_PURE
)

Sponge.set_backend("jax")

def Calculate_Force():
    Sponge.dd.frc.at[0, 0].add(1.0)
""",
        encoding="utf-8",
    )


@pytest.mark.parametrize(
    ("failure_mode", "expected_error"),
    [
        (
            "missing-script",
            "PRIPS requires a Python force module",
        ),
        (
            "invalid-spec",
            "cannot create an import specification for PRIPS module",
        ),
        (
            "noncallable-force",
            "must define callable Calculate_Force",
        ),
    ],
)
def test_prips_non_mc_rejects_unusable_python_force_module(
    tmp_path, monkeypatch, failure_mode, expected_error
):
    plugin_path = _prips_plugin_path()
    case_dir = tmp_path / failure_mode
    case_dir.mkdir()
    monkeypatch.delenv("SPONGE_PRIPS_PY", raising=False)

    py_name = None
    if failure_mode == "invalid-spec":
        py_name = "prips_force.invalid"
        (case_dir / py_name).write_text(
            "def Calculate_Force():\n    pass\n", encoding="utf-8"
        )
    elif failure_mode == "noncallable-force":
        py_name = "prips_force.py"
        (case_dir / py_name).write_text(
            "Calculate_Force = None\n", encoding="utf-8"
        )

    _write_prips_non_mc_mdin(case_dir, plugin_path, py_name=py_name)
    result = _run_sponge_raw(case_dir)
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0, output
    assert expected_error in output
    assert "END INITIALIZING SPONGE PLUGIN" not in output


def test_tip3p_prips_plugin_hooks_run(statics_path, outputs_path, mpi_np):
    backend = os.environ.get("PRIPS_TEST_BACKEND", "auto")
    assert backend in {"auto", "numpy", "jax", "cupy", "pytorch"}
    plugin_path = _prips_plugin_path()
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p",
        mpi_np=mpi_np,
        run_name="tip3p_prips",
    )
    _write_prips_script(case_dir, backend)
    _write_prips_mdin(case_dir, plugin_path)
    spaced_script_name = "prips script # complete.py"
    (case_dir / "prips_test.py").rename(case_dir / spaced_script_name)
    mdin_path = case_dir / "mdin.spg.toml"
    mdin_path.write_text(
        mdin_path.read_text(encoding="utf-8").replace(
            'py = "prips_test.py"', f"py = {json.dumps(spaced_script_name)}"
        ),
        encoding="utf-8",
    )
    Runner.run_sponge(case_dir, timeout=1200, mpi_np=mpi_np)

    hook_log = case_dir / "prips_hook.log"
    assert hook_log.exists()
    hook_lines = hook_log.read_text(encoding="utf-8").splitlines()

    backend_line = next(
        (line for line in hook_lines if line.startswith("backend ")),
        None,
    )
    after_init_line = next(
        (line for line in hook_lines if line.startswith("after_init ")),
        None,
    )
    force_lines = [
        line for line in hook_lines if line.startswith("calculate_force ")
    ]
    final_force_line = next(
        (
            line
            for line in reversed(hook_lines)
            if line.startswith("mdout_force ")
        ),
        None,
    )
    sponge_forces = Extractor.extract_sponge_forces(case_dir, 1011)

    assert backend_line is not None
    backend_match = re.fullmatch(
        r"backend device=(1|2|10) name=([a-z]+)", backend_line
    )
    assert backend_match is not None
    device_type = int(backend_match.group(1))
    expected_backend = (
        {1: "numpy", 2: "pytorch", 10: "pytorch"}[device_type]
        if backend == "auto"
        else backend
    )
    assert backend_match.group(2) == expected_backend
    assert after_init_line is not None
    assert "atom_numbers=1011" in after_init_line
    assert "neighbor_max=1200" in after_init_line
    assert "coord[0,0]=" in after_init_line
    assert len(force_lines) == 2
    assert final_force_line is not None
    assert all("local_atom_numbers=1011" in line for line in force_lines)
    force_matches = [
        re.search(
            r"calculate_force\s+step=(\d+)\s+"
            r"commit_sampling_state=([01])\s+exact_state=([01])\s+"
            r"local_atom_numbers=\d+\s+"
            r"force\[0,0\]_before=([-\d.]+)\s+"
            r"force\[0,0\]_after=([-\d.]+)\s+"
            r"force\[0,0\]_delta=([-\d.]+)",
            line,
        )
        for line in force_lines
    ]
    assert all(match is not None for match in force_matches)
    force_records = {
        int(match.group(1)): (
            float(match.group(4)),
            float(match.group(5)),
            float(match.group(6)),
        )
        for match in force_matches
    }
    assert all(match.group(2) == "1" for match in force_matches)
    assert all(match.group(3) == "0" for match in force_matches)
    assert set(force_records) == {0, 1}
    step0_before, step0_after, step0_delta = force_records[0]
    step1_before, step1_after, step1_delta = force_records[1]

    Outputer.print_table(
        ["Metric", "Value"],
        [
            ["Case", "tip3p_prips"],
            ["PluginPath", str(plugin_path)],
            [
                "Backend",
                backend_line.split("=", 1)[1] if backend_line else "N/A",
            ],
            ["AfterInitial", "PASS" if after_init_line else "MISSING"],
            ["CalculateForce", "PASS" if force_lines else "MISSING"],
            ["MdoutForce", "PASS" if final_force_line else "MISSING"],
            ["F[0,0,0]", f"{step0_after:.6f}"],
            ["F[1,0,0]", f"{step1_after:.6f}"],
            ["ΔF[0,0]", f"{step1_after - step0_after:.6f}"],
        ],
        title="Misc Validation: TIP3P PRIPS",
    )

    assert abs(step0_after - step0_before) < 1e-5
    assert abs(step0_delta) < 1e-5
    assert abs(step1_after - step1_before - 1.25) < 1e-5
    assert abs(step1_delta - 1.25) < 1e-5
    assert abs(step1_delta - step0_delta - 1.25) < 1e-5
    final_force_match = re.search(
        r"mdout_force\s+step=\d+\s+force\[0,0\]_final=([-\d.]+)",
        final_force_line,
    )
    assert final_force_match is not None
    final_force = float(final_force_match.group(1))
    assert abs(sponge_forces[0, 0] - final_force) < 1e-5


def test_prips_jax_functional_force_energy_virial_writeback(tmp_path):
    pytest.importorskip("jax")
    plugin_path = _prips_plugin_path()
    case_dir = tmp_path / "prips_jax_writeback"
    case_dir.mkdir()
    _write_prips_jax_writeback_script(case_dir)
    _write_prips_neighbor_list_mdin(
        case_dir,
        plugin_path,
        py_name="prips_jax_writeback_test.py",
        step_limit=0,
        initial_ratio=0.0,
        print_pressure=True,
    )
    # Keep the two atoms inside the LJ cutoff so every evaluation starts with
    # a nonzero core force/energy/virial.  The returned arrays must preserve
    # that baseline and add the JAX update, not overwrite it with a bare
    # contribution.
    (case_dir / "coordinate.txt").write_text(
        "2\n9 10 10\n11 10 10\n20 20 20\n90 90 90\n",
        encoding="utf-8",
    )
    Runner.run_sponge(case_dir, timeout=1200)

    lines = (
        (case_dir / "jax_writeback.log")
        .read_text(encoding="utf-8")
        .splitlines()
    )
    eval_pattern = re.compile(
        r"eval force_before=([-+\deE.]+) force_after=([-+\deE.]+) "
        r"energy_before=([-+\deE.]+) energy_after=([-+\deE.]+) "
        r"virial_before=([-+\deE.]+) virial_after=([-+\deE.]+) "
        r"needs_energy=([01]) needs_virial=([01])"
    )
    evaluations = [
        eval_pattern.fullmatch(line)
        for line in lines
        if line.startswith("eval ")
    ]
    assert len(evaluations) == 3
    assert all(match is not None for match in evaluations)
    assert all(match.group(7) == "1" for match in evaluations)
    assert all(match.group(8) == "1" for match in evaluations)

    last = evaluations[-1]
    force_before, force_after = map(float, last.group(1, 2))
    energy_before, energy_after = map(float, last.group(3, 4))
    virial_before, virial_after = map(float, last.group(5, 6))
    assert abs(force_before) > 1.0e-4
    assert abs(energy_before) > 1.0e-4
    assert abs(virial_before) > 1.0e-4
    assert force_after == pytest.approx(force_before + 1.25, abs=2.0e-4)
    assert energy_after == pytest.approx(energy_before + 2.5, abs=2.0e-4)
    assert virial_after == pytest.approx(virial_before + 3.75, abs=2.0e-4)

    final_match = re.fullmatch(
        r"final force=([-+\deE.]+) energy=([-+\deE.]+) "
        r"virial=([-+\deE.]+)",
        lines[-1],
    )
    assert final_match is not None
    final_force, final_energy, final_virial = map(float, final_match.groups())
    assert final_force == pytest.approx(force_after, rel=2.0e-6, abs=2.0e-4)
    assert final_energy == pytest.approx(energy_after, rel=2.0e-6, abs=2.0e-4)
    assert final_virial == pytest.approx(virial_after, rel=2.0e-6, abs=2.0e-4)


def test_prips_jax_missing_functional_result_is_fatal(tmp_path):
    pytest.importorskip("jax")
    plugin_path = _prips_plugin_path()
    case_dir = tmp_path / "prips_jax_missing_result"
    case_dir.mkdir()
    _write_prips_jax_missing_result_script(case_dir)
    _write_prips_neighbor_list_mdin(
        case_dir,
        plugin_path,
        py_name="prips_jax_missing_result_test.py",
        step_limit=0,
    )

    result = _run_sponge_raw(case_dir, timeout=1200)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "JAX Calculate_Force hook must return Sponge.force_result" in output
    assert "complete updated local buffers" in output


def test_prips_mc_force_evaluation_context(statics_path, outputs_path, mpi_np):
    plugin_path = _prips_plugin_path()
    case_dir = Outputer.prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="tip3p",
        mpi_np=mpi_np,
        run_name="tip3p_prips_mc_context",
    )
    _write_prips_script(case_dir, "auto", mc_safe=True)
    _write_prips_mc_mdin(case_dir, plugin_path)
    Runner.run_sponge(case_dir, timeout=1200, mpi_np=mpi_np)

    force_lines = [
        line
        for line in (case_dir / "prips_hook.log")
        .read_text(encoding="utf-8")
        .splitlines()
        if line.startswith("calculate_force ")
    ]
    contexts = [
        tuple(
            int(value)
            for value in re.search(
                r"step=(\d+)\s+commit_sampling_state=([01])\s+"
                r"exact_state=([01])",
                line,
            ).groups()
        )
        for line in force_lines
    ]
    assert contexts == [(0, 0, 1), (0, 0, 1), (0, 1, 1)]


def test_prips_mc_exact_evaluation_rebuilds_neighbor_list(tmp_path):
    plugin_path = _prips_plugin_path()
    case_dir = tmp_path / "prips_exact_neighbor_list"
    case_dir.mkdir()
    _write_prips_neighbor_list_script(case_dir)
    _write_prips_neighbor_list_mdin(case_dir, plugin_path)
    Runner.run_sponge(case_dir, timeout=1200)

    records = []
    for line in (
        (case_dir / "neighbor_context.log")
        .read_text(encoding="utf-8")
        .splitlines()
    ):
        match = re.fullmatch(
            r"step=(\d+) commit=([01]) exact=([01]) neighbors=(\d+)", line
        )
        assert match is not None
        records.append(tuple(int(value) for value in match.groups()))

    assert [(step, commit, exact) for step, commit, exact, _ in records] == [
        (0, 0, 1),
        (0, 0, 1),
        (0, 1, 1),
        (1, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
    ]
    assert [neighbors for *_, neighbors in records[:3]] == [0, 0, 0]
    assert all(neighbors > 0 for *_, neighbors in records[3:])


@pytest.mark.parametrize(
    ("capability_expression", "expected_error"),
    [
        (
            "Sponge.FORCE_VIRIAL_COMPLETE | Sponge.FORCE_PURE",
            "ENERGY_COMPLETE capability",
        ),
        (
            "Sponge.FORCE_ENERGY_COMPLETE | Sponge.FORCE_PURE",
            "VIRIAL_COMPLETE capability",
        ),
        (
            "Sponge.FORCE_ENERGY_COMPLETE | Sponge.FORCE_VIRIAL_COMPLETE",
            "exactly one of PURE or TRANSACTIONAL",
        ),
        (
            "Sponge.FORCE_ENERGY_COMPLETE | "
            "Sponge.FORCE_VIRIAL_COMPLETE | "
            "Sponge.FORCE_TRANSACTIONAL",
            "transactional Python plugin must define Begin_Force_Transaction",
        ),
    ],
    ids=[
        "missing-energy",
        "missing-virial",
        "missing-state-contract",
        "missing-lifecycle-hooks",
    ],
)
def test_prips_mc_rejects_incomplete_force_contract(
    tmp_path, mpi_np, capability_expression, expected_error
):
    plugin_path = _prips_plugin_path()
    case_dir = tmp_path / "prips_incomplete_contract"
    case_dir.mkdir()
    _write_prips_contract_script(case_dir, capability_expression)
    _write_prips_neighbor_list_mdin(
        case_dir,
        plugin_path,
        py_name="prips_contract_test.py",
        step_limit=0,
    )

    result = _run_sponge_raw(case_dir, mpi_np=mpi_np)
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert expected_error in output


@pytest.mark.parametrize(
    ("energy_sign", "accepted"),
    [(1.0, False), (-1.0, True)],
    ids=["reject", "accept"],
)
def test_prips_mc_transactional_hamiltonian_and_local_buffers(
    tmp_path, energy_sign, accepted
):
    plugin_path = _prips_plugin_path()
    case_dir = tmp_path / ("prips_mc_accept" if accepted else "prips_mc_reject")
    case_dir.mkdir()
    _write_prips_transaction_script(case_dir, energy_sign)
    _write_prips_neighbor_list_mdin(
        case_dir,
        plugin_path,
        py_name="prips_transaction_test.py",
        step_limit=0,
        initial_ratio=0.05,
        print_pressure=True,
    )
    Runner.run_sponge(case_dir, timeout=1200)

    lines = (
        (case_dir / "transaction.log").read_text(encoding="utf-8").splitlines()
    )
    assert [line.split()[0] for line in lines] == [
        "begin",
        "eval",
        "rollback",
        "begin",
        "eval",
        "rollback",
        "begin",
        "eval",
        "commit",
        "final",
    ]

    evaluation_pattern = re.compile(
        r"eval index=(\d+) commit=([01]) exact=([01]) "
        r"needs_energy=([01]) needs_virial=([01]) "
        r"energy_shape=\((\d+),\) virial_shape=\((\d+), (\d+)\) "
        r"contribution=([-+\deE.]+) energy_delta=([-+\deE.]+) "
        r"virial_delta=([-+\deE.]+)"
    )
    evaluations = [
        evaluation_pattern.fullmatch(line)
        for line in lines
        if line.startswith("eval ")
    ]
    assert all(match is not None for match in evaluations)
    assert [int(match.group(1)) for match in evaluations] == [1, 2, 3]
    assert [int(match.group(2)) for match in evaluations] == [0, 0, 1]
    assert all(int(match.group(3)) == 1 for match in evaluations)
    assert all(int(match.group(4)) == 1 for match in evaluations)
    assert all(int(match.group(5)) == 1 for match in evaluations)
    assert all(int(match.group(6)) == 2 for match in evaluations)
    assert all(int(match.group(7)) == 2 for match in evaluations)
    assert all(int(match.group(8)) == 6 for match in evaluations)

    contributions = [float(match.group(9)) for match in evaluations]
    energy_deltas = [float(match.group(10)) for match in evaluations]
    virial_deltas = [float(match.group(11)) for match in evaluations]
    assert contributions[0] == pytest.approx(0.0, abs=1e-5)
    assert energy_sign * contributions[1] > 1.0e3
    assert energy_deltas == pytest.approx(contributions, rel=2e-6, abs=1e-4)
    assert virial_deltas[0] == pytest.approx(0.0, abs=1e-5)
    assert abs(virial_deltas[1]) > 1.0e-3

    if accepted:
        assert energy_sign * contributions[2] > 1.0e3
        assert virial_deltas[2] == pytest.approx(
            virial_deltas[1], rel=2e-5, abs=1e-3
        )
    else:
        assert contributions[2] == pytest.approx(0.0, abs=1e-5)
        assert virial_deltas[2] == pytest.approx(0.0, abs=1e-5)

    assert lines[2] == "rollback pending=1 history=0"
    assert lines[5] == "rollback pending=1 history=0"
    assert lines[8] == "commit history=1"
    final_match = re.fullmatch(
        r"final history=(\d+) pending=(\S+) atom0_x=([-+\deE.]+) "
        r"virial00=([-+\deE.]+)",
        lines[9],
    )
    assert final_match is not None
    assert final_match.group(1) == "1"
    assert final_match.group(2) == "None"
    final_x = float(final_match.group(3))
    if accepted:
        assert abs(final_x - 5.0) > 1e-4
    else:
        assert final_x == pytest.approx(5.0, abs=1e-5)
    assert float(final_match.group(4)) == pytest.approx(
        virial_deltas[2], rel=2e-5, abs=1e-3
    )
