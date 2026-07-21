import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("qc_atomic_identity_probe.cpp")


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = (
        shutil.which("c++") or shutil.which("clang++") or shutil.which("g++")
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the QC identity probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.append(
        REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include"
    )
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (
            candidate / "fftw3.h"
        ).is_file():
            return candidate
    pytest.skip(
        "OpenMP and FFTW headers are required for the QC identity probe"
    )


def _sponge_binary():
    configured = os.environ.get("SPONGE_BIN", "SPONGE")
    candidate = Path(configured)
    if candidate.is_file():
        return str(candidate.resolve())
    resolved = shutil.which(configured)
    if resolved is None:
        pytest.skip("SPONGE executable is required for the ECP selection test")
    return resolved


def _write_auto_ecp_case(case_dir):
    (case_dir / "mass.txt").write_text("1\n85.4678\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("1\n0\n", encoding="utf-8")
    (case_dir / "coordinate.txt").write_text(
        "1\n10 10 10\n40 40 40\n90 90 90\n", encoding="utf-8"
    )
    (case_dir / "qc_type.txt").write_text("1 0 2\n0 Rb\n", encoding="utf-8")
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            (
                'md_name = "auto_ecp_missing_data"',
                'mode = "nve"',
                "step_limit = 0",
                "dt = 0",
                "cutoff = 6",
                "PM.MPI_size = 0",
                'mass_in_file = "mass.txt"',
                'charge_in_file = "charge.txt"',
                'coordinate_in_file = "coordinate.txt"',
                'qc_type_in_file = "qc_type.txt"',
                'qc_model_chemistry = "HF/def2-svp"',
                'qc_ecp = "auto"',
                "qc_need_gradient = 0",
                "print_zeroth_frame = true",
                "write_restart_file_interval = 0",
                "dont_check_input = 1",
            )
        )
        + "\n",
        encoding="utf-8",
    )


def _write_qc_initialization_case(
    case_dir,
    *,
    symbol="H",
    charge=0,
    multiplicity=2,
    qc_type_content=None,
):
    (case_dir / "mass.txt").write_text("1\n1.0\n", encoding="utf-8")
    (case_dir / "charge.txt").write_text("1\n0\n", encoding="utf-8")
    (case_dir / "coordinate.txt").write_text(
        "1\n10 10 10\n40 40 40\n90 90 90\n", encoding="utf-8"
    )
    if qc_type_content is None:
        qc_type_content = f"1 {charge} {multiplicity}\n0 {symbol}\n"
    (case_dir / "qc_type.txt").write_text(
        qc_type_content, encoding="utf-8"
    )
    (case_dir / "mdin.spg.toml").write_text(
        "\n".join(
            (
                'md_name = "qc_initialization_contract"',
                'mode = "nve"',
                "step_limit = 0",
                "dt = 0",
                "cutoff = 6",
                "PM.MPI_size = 0",
                'mass_in_file = "mass.txt"',
                'charge_in_file = "charge.txt"',
                'coordinate_in_file = "coordinate.txt"',
                'qc_type_in_file = "qc_type.txt"',
                'qc_model_chemistry = "HF/STO-3G"',
                'qc_ecp = "none"',
                "print_zeroth_frame = true",
                "write_restart_file_interval = 0",
                "dont_check_input = 1",
            )
        )
        + "\n",
        encoding="utf-8",
    )


def _run_initialization_case(case_dir, extra_arguments=()):
    return subprocess.run(
        [
            _sponge_binary(),
            "-mdin",
            "mdin.spg.toml",
            *extra_arguments,
        ],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )


def test_ecp_effective_charge_cannot_replace_element_identity(tmp_path):
    executable = tmp_path / "qc_atomic_identity_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-DUSE_CPU",
            "-DNO_GLOBAL_CONTROLLER",
            "-O3",
            "-ffast-math",
            "-w",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
            f"-I{_dependency_include()}",
            str(PROBE_SOURCE),
            str(
                REPOSITORY_ROOT
                / "SPONGE"
                / "quantum_chemistry"
                / "guess"
                / "minao.cpp"
            ),
            "-o",
            str(executable),
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert compile_result.returncode == 0, (
        "failed to compile QC identity probe\n"
        f"stdout:\n{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )
    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr
    assert run_result.stdout.split() == ["K", "19", "9", "2.03", "0.57"]


def test_auto_ecp_missing_required_parameters_is_fatal(tmp_path):
    _write_auto_ecp_case(tmp_path)
    result = subprocess.run(
        [_sponge_binary(), "-mdin", "mdin.spg.toml"],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert "Automatic ECP selection for basis def2-svp" in output
    assert "requires def2-ecp parameters for element Rb" in output
    assert "Refusing to silently continue with an all-electron" in output


@pytest.mark.parametrize(
    "symbol,charge,multiplicity,arguments,diagnostics",
    (
        (
            "H",
            1,
            1,
            (),
            (
                "Invalid molecular electron configuration (N=0, charge=1, "
                "multiplicity=1)",
                "must contain at least one explicit electron",
            ),
        ),
        (
            "H",
            2,
            1,
            (),
            (
                "Invalid molecular electron configuration (N=-1, charge=2, "
                "multiplicity=1)",
                "must contain at least one explicit electron",
            ),
        ),
        (
            "H",
            0,
            3,
            ("-qc_restricted", "0"),
            (
                "Invalid molecular electron configuration (N=1, charge=0, "
                "multiplicity=3)",
                "requires more unpaired electrons than exist",
            ),
        ),
        (
            "H",
            0,
            1,
            ("-qc_restricted", "0"),
            (
                "Invalid molecular electron configuration (N=1, charge=0, "
                "multiplicity=1)",
                "have inconsistent parity",
            ),
        ),
        (
            "He",
            -1,
            2,
            ("-qc_restricted", "0"),
            (
                "N=3, n_alpha=2, n_beta=1, nao=1",
                "occupied spin orbitals exceed the orbital-basis capacity",
            ),
        ),
    ),
    ids=(
        "zero-electron",
        "negative-electron",
        "multiplicity-too-high",
        "parity",
        "ao-capacity",
    ),
)
def test_invalid_electron_configuration_is_controller_fatal(
    tmp_path, symbol, charge, multiplicity, arguments, diagnostics
):
    _write_qc_initialization_case(
        tmp_path,
        symbol=symbol,
        charge=charge,
        multiplicity=multiplicity,
    )
    result = _run_initialization_case(tmp_path, arguments)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    for diagnostic in diagnostics:
        assert diagnostic in output


@pytest.mark.parametrize(
    "arguments,diagnostics",
    (
        (
            ("-qc_need_gradient", "2"),
            ("qc_need_gradient must be exactly 0 or 1",),
        ),
        (
            ("-qc_eri_prim_screen_tol", "1e999"),
            (
                "qc_eri_prim_screen_tol must be a finite, nonnegative float",
                "outside the float range",
            ),
        ),
        (
            ("-qc_diis_space", "13"),
            ("qc_diis_space must be in [2, 12]",),
        ),
        (
            ("-qc_diis_space", "2147483648"),
            (
                "qc_diis_space must be an exactly representable integer",
                "outside the int range",
            ),
        ),
    ),
    ids=("non-bool", "infinite-screen", "diis-too-large", "int-overflow"),
)
def test_invalid_qc_control_is_controller_fatal(
    tmp_path, arguments, diagnostics
):
    _write_qc_initialization_case(tmp_path)
    result = _run_initialization_case(tmp_path, arguments)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    for diagnostic in diagnostics:
        assert diagnostic in output


@pytest.mark.parametrize(
    "qc_type_content,diagnostic",
    (
        (
            "1 0 2 trailing\n0 H\n",
            "must contain exactly three integers: natm charge multiplicity",
        ),
        (
            "1 0 2\n0 H trailing\n",
            "must contain exactly an MD atom index and element symbol",
        ),
        (
            "1 0 2\n0 H\n1 H\n",
            "contains an unexpected nonempty record after the declared 1 QC "
            "atom lines",
        ),
    ),
    ids=("header-trailing-token", "atom-trailing-token", "extra-record"),
)
def test_qc_type_requires_exact_records(
    tmp_path, qc_type_content, diagnostic
):
    _write_qc_initialization_case(
        tmp_path, qc_type_content=qc_type_content
    )
    result = _run_initialization_case(tmp_path)
    output = result.stdout + result.stderr
    assert result.returncode != 0, output
    assert diagnostic in output
