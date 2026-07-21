import json
import math
import os
import re
import shlex
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("softcore_pair_probe.cpp")


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    # Conda's compiler activation can put a cross-compiler wrapper named
    # `c++` first on macOS without the matching SDK/linker activation flags.
    # The system AppleClang is the native compiler used by the CPU build and
    # can link this standalone host-only probe directly.
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = shutil.which("c++") or shutil.which("clang++") or shutil.which(
        "g++"
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the soft-core pair probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.append(REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include")
    for candidate in candidates:
        if (candidate / "fftw3.h").is_file():
            return candidate
    pytest.skip("fftw3.h is required for the standalone soft-core pair probe")


def test_production_softcore_pair_hamiltonian(tmp_path):
    executable = tmp_path / "softcore_pair_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
            "-ffast-math",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(PROBE_SOURCE),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        f"failed to compile soft-core pair probe\nstdout:\n"
        f"{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )
    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, (
        f"soft-core pair probe failed\nstdout:\n{run_result.stdout}\n"
        f"stderr:\n{run_result.stderr}"
    )


def _write_counted_values(path, values):
    path.write_text(
        str(len(values))
        + "\n"
        + "\n".join(str(value) for value in values)
        + "\n",
        encoding="utf-8",
    )


def _write_softcore_case(
    case_dir,
    *,
    settings=None,
    endpoint_a=(1.0, -1.0),
    endpoint_b=(0.0, 0.0),
    current=None,
    write_endpoint_a=True,
    write_endpoint_b=True,
    pair_center_x=21.0,
    pair_distance=2.0,
):
    _write_counted_values(case_dir / "mass.txt", (12.0, 12.0))
    pair_x_a = pair_center_x - 0.5 * pair_distance
    pair_x_b = pair_center_x + 0.5 * pair_distance
    (case_dir / "coordinate.txt").write_text(
        f"2\n{pair_x_a:.12g} 20 20\n{pair_x_b:.12g} 20 20\n"
        "40 40 40\n90 90 90\n",
        encoding="utf-8",
    )
    (case_dir / "lj_soft.txt").write_text(
        "2 1 1\n0\n0\n0\n0\n0 0\n0 0\n", encoding="utf-8"
    )
    if write_endpoint_a:
        _write_counted_values(case_dir / "charge_a.txt", endpoint_a)
    if write_endpoint_b:
        _write_counted_values(case_dir / "charge_b.txt", endpoint_b)
    if current is not None:
        _write_counted_values(case_dir / "charge.txt", current)

    mdin = {
        "md_name": case_dir.name,
        "mode": "nve",
        "step_limit": 1,
        "dt": 0,
        "cutoff": 8.0,
        "mass_in_file": "mass.txt",
        "coordinate_in_file": "coordinate.txt",
        "LJ_soft_core_in_file": "lj_soft.txt",
        "lambda_lj": 0.25,
        "charge_A_in_file": "charge_a.txt" if write_endpoint_a else None,
        "charge_B_in_file": "charge_b.txt" if write_endpoint_b else None,
        "charge_in_file": "charge.txt" if current is not None else None,
        "mdout": "mdout.txt",
        "crd": "crd.dat",
        "frc": "frc.dat",
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
    }
    if settings:
        mdin.update(settings)
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}"
            for key, value in mdin.items()
            if value is not None
        )
        + "\n",
        encoding="utf-8",
    )


def _run_sponge(case_dir, mpi_np=None):
    command = [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", "mdin.spg.toml"]
    if mpi_np is not None:
        command = ["mpirun", "--oversubscribe", "-np", str(mpi_np), *command]
    return subprocess.run(
        command,
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )


def _read_mdout_and_forces(case_dir):
    rows = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    mdout = {
        name: float(value)
        for name, value in zip(rows[0].split(), rows[-1].split())
    }
    raw_force = (case_dir / "frc.dat").read_bytes()
    forces = struct.unpack(f"={len(raw_force) // 4}f", raw_force)[-6:]
    return mdout, forces


def test_charge_endpoints_generate_the_production_current(tmp_path):
    case_dir = tmp_path / "derived_current"
    case_dir.mkdir()
    _write_softcore_case(case_dir)
    result = _run_sponge(case_dir)
    assert result.returncode == 0, result.stdout + result.stderr

    rows = (case_dir / "mdout.txt").read_text(encoding="utf-8").splitlines()
    values = [float(value) for value in rows[-1].split()]
    assert all(math.isfinite(value) for value in values)
    assert max(abs(value) for value in values) > 1.0e-5
    raw_force = (case_dir / "frc.dat").read_bytes()
    forces = struct.unpack(f"={len(raw_force) // 4}f", raw_force)
    assert all(math.isfinite(value) for value in forces)
    assert max(abs(value) for value in forces) > 1.0e-5


def test_sits_and_standard_softcore_share_charge_endpoint_hamiltonian(tmp_path):
    standard_dir = tmp_path / "standard_endpoints"
    standard_dir.mkdir()
    _write_softcore_case(standard_dir)
    standard_result = _run_sponge(standard_dir)
    assert standard_result.returncode == 0, (
        standard_result.stdout + standard_result.stderr
    )

    sits_dir = tmp_path / "sits_endpoints"
    sits_dir.mkdir()
    _write_softcore_case(
        sits_dir,
        settings={"SITS.mode": "observation", "SITS.atom_numbers": 1},
    )
    sits_result = _run_sponge(sits_dir)
    assert sits_result.returncode == 0, sits_result.stdout + sits_result.stderr

    standard_mdout, standard_forces = _read_mdout_and_forces(standard_dir)
    sits_mdout, sits_forces = _read_mdout_and_forces(sits_dir)
    assert sits_mdout["eff_pot"] == pytest.approx(
        standard_mdout["eff_pot"], abs=2.0e-6
    )
    assert sits_forces == pytest.approx(standard_forces, abs=2.0e-5)


def _run_endpoint_boundary_point(root, name, distance, launch_np):
    case_dir = root / name
    case_dir.mkdir()
    settings = {
        "PM.MPI_size": 0,
        "PME.calculate_reciprocal_part": False,
    }
    if launch_np is not None:
        settings.update(
            {
                "DOM_DEC.split_nx": launch_np,
                "DOM_DEC.split_ny": 1,
                "DOM_DEC.split_nz": 1,
            }
        )
    _write_softcore_case(
        case_dir,
        settings=settings,
        pair_center_x=20.0,
        pair_distance=distance,
    )
    result = _run_sponge(case_dir, mpi_np=launch_np)
    assert result.returncode == 0, result.stdout + result.stderr
    mdout, forces = _read_mdout_and_forces(case_dir)
    beta_match = re.search(
        r"\bbeta:\s*([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)",
        result.stdout + result.stderr,
    )
    assert beta_match is not None
    return mdout, forces, float(beta_match.group(1))


def _charge_endpoint_pair_reference(distance, beta):
    lambda_lj = 0.25
    alpha = 0.5
    sigma_6 = 3.0**6
    charge_product = (1.0 - lambda_lj) * (-1.0 + lambda_lj)
    distance_a = (
        distance**6 + alpha * lambda_lj * sigma_6
    ) ** (1.0 / 6.0)
    distance_b = (
        distance**6 + alpha * (1.0 - lambda_lj) * sigma_6
    ) ** (1.0 / 6.0)

    def energy_kernel(soft_distance):
        return math.erfc(beta * soft_distance) / soft_distance

    def radial_force_kernel(soft_distance):
        beta_distance = beta * soft_distance
        return (
            charge_product
            * distance**4
            * (
                math.exp(-(beta_distance**2))
                * 2.0
                / math.sqrt(math.pi)
                * beta
                + math.erfc(beta_distance) / soft_distance
            )
            * soft_distance**-6
        )

    energy = charge_product * (
        (1.0 - lambda_lj) * energy_kernel(distance_a)
        + lambda_lj * energy_kernel(distance_b)
    )
    force_on_a_x = -distance * (
        (1.0 - lambda_lj) * radial_force_kernel(distance_a)
        + lambda_lj * radial_force_kernel(distance_b)
    )
    return energy, force_on_a_x


def test_charge_endpoint_owned_ghost_mapping_and_pair_derivatives(tmp_path, mpi_np):
    distance = 1.0
    delta = 1.0e-3
    launch_modes = [("direct", None)]
    if mpi_np == 2:
        launch_modes.append(("mpi_2", mpi_np))

    results = {}
    for mode, launch_np in launch_modes:
        base = _run_endpoint_boundary_point(
            tmp_path, f"{mode}_base", distance, launch_np
        )
        plus = _run_endpoint_boundary_point(
            tmp_path, f"{mode}_plus", distance + delta, launch_np
        )
        minus = _run_endpoint_boundary_point(
            tmp_path, f"{mode}_minus", distance - delta, launch_np
        )
        mdout, forces, beta = base
        expected_energy, expected_force = _charge_endpoint_pair_reference(
            distance, beta
        )
        assert mdout["eff_pot"] == pytest.approx(expected_energy, abs=2.0e-5)
        assert forces[0] == pytest.approx(expected_force, rel=3.0e-5, abs=2.0e-6)
        assert forces[3] == pytest.approx(-expected_force, rel=3.0e-5, abs=2.0e-6)
        force_from_energy = (
            plus[0]["eff_pot"] - minus[0]["eff_pot"]
        ) / (2.0 * delta)
        assert forces[0] == pytest.approx(
            force_from_energy, rel=3.0e-3, abs=2.0e-4
        )
        results[mode] = (mdout, forces)

    if mpi_np == 2:
        assert results["mpi_2"][0]["eff_pot"] == pytest.approx(
            results["direct"][0]["eff_pot"], abs=2.0e-6
        )
        assert results["mpi_2"][1] == pytest.approx(
            results["direct"][1], abs=2.0e-5
        )


def test_supplied_current_must_match_charge_endpoints_exactly(tmp_path):
    case_dir = tmp_path / "wrong_current"
    case_dir.mkdir()
    _write_softcore_case(case_dir, current=(0.5, -0.5))
    result = _run_sponge(case_dir)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "does not equal fma(lambda_lj, qB-qA, qA)" in output


@pytest.mark.parametrize(
    ("case", "kwargs", "message"),
    [
        (
            "missing_b",
            {"write_endpoint_b": False},
            "must be provided together",
        ),
        (
            "count_mismatch",
            {"endpoint_b": (0.0,)},
            "differs from LJ_soft_core_in_file atom count",
        ),
        (
            "nonfinite",
            {"endpoint_a": ("nan", -1.0)},
            "strict finite decimal",
        ),
        (
            "without_softcore",
            {"settings": {"LJ_soft_core_in_file": None}},
            "charge endpoint files require LJ_soft_core_in_file",
        ),
        (
            "mixed_input_source",
            {"settings": {"gromacs_top": "unused.top"}},
            "belong to the native LJ_soft_core input contract",
        ),
    ],
)
def test_invalid_charge_endpoint_inputs_fail(tmp_path, case, kwargs, message):
    case_dir = tmp_path / case
    case_dir.mkdir()
    _write_softcore_case(case_dir, **kwargs)
    result = _run_sponge(case_dir)
    assert result.returncode != 0
    assert message in result.stdout + result.stderr


@pytest.mark.parametrize(
    ("setting", "value"),
    [
        ("lambda_lj", -0.1),
        ("lambda_lj", 1.1),
        ("lambda_lj", "nan"),
        ("lambda_lj", "1e-45"),
        ("soft_core_alpha", 0.0),
        ("soft_core_powfer", 0.5),
        ("soft_core_powfer", 1000.0),
        ("soft_core_sigma", 0.0),
        ("soft_core_sigma", "1e-10"),
        ("soft_core_sigma", "1e7"),
        ("soft_core_sigma_min", -1.0),
    ],
)
def test_invalid_softcore_parameters_fail(tmp_path, setting, value):
    case_dir = tmp_path / f"invalid_{setting}_{value}"
    case_dir.mkdir()
    _write_softcore_case(
        case_dir,
        settings={setting: value},
        current=(0.0, 0.0),
        write_endpoint_a=False,
        write_endpoint_b=False,
    )
    result = _run_sponge(case_dir)
    output = result.stdout + result.stderr
    assert result.returncode != 0
    assert "LJ_SOFT_CORE::Initial" in output
