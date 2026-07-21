import os
import shlex
import shutil
import struct
import subprocess
import sys
from decimal import Decimal, localcontext
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("boys_accuracy_probe.cpp")


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
        pytest.skip("a C++17 compiler is required for the Boys function probe")
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
        "OpenMP and FFTW headers are required for the Boys function probe"
    )


def _float_from_bits(text):
    return struct.unpack("!f", bytes.fromhex(text))[0]


def _boys_reference(t, order):
    """Evaluate the positive Boys series with 100-digit Decimal arithmetic."""
    with localcontext() as context:
        context.prec = 110
        td = Decimal.from_float(t)
        if td.is_zero():
            return Decimal(1) / Decimal(2 * order + 1)

        term = (-td).exp() / Decimal(2 * order + 1)
        total = term
        tolerance = Decimal("1e-90")
        k = 0
        while True:
            k += 1
            ratio = (2 * td) / Decimal(2 * order + 2 * k + 1)
            term *= ratio
            total += term
            next_ratio = (2 * td) / Decimal(2 * order + 2 * k + 3)
            if next_ratio < 1:
                tail_bound = term * next_ratio / (1 - next_ratio)
                if tail_bound <= tolerance * total:
                    return +total


def test_all_boys_paths_match_high_precision_reference_across_regimes(tmp_path):
    executable = tmp_path / "boys_accuracy_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-O3",
            "-ffast-math",
            "-DUSE_CPU",
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
        f"failed to compile Boys function probe\nstdout:\n"
        f"{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )

    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr

    worst_relative_error = Decimal(0)
    observed_orders = {}
    for line in run_result.stdout.splitlines():
        (
            bits,
            requested_max_text,
            order_text,
            one_text,
            ri_text,
            eri_text,
        ) = line.split()
        t = _float_from_bits(bits)
        requested_max = int(requested_max_text)
        order = int(order_text)
        assert 0 <= order <= requested_max <= 16
        observed_orders.setdefault((bits, requested_max), set()).add(order)
        values = [float(one_text), float(ri_text), float(eri_text)]

        # All production call sites are wrappers around exactly one algorithm;
        # requiring bit identity prevents those paths from drifting again.
        assert values[0] == values[1] == values[2]

        reference = _boys_reference(t, order)
        actual = Decimal.from_float(values[0])
        relative_error = abs(actual - reference) / reference
        worst_relative_error = max(worst_relative_error, relative_error)

    assert observed_orders
    for (_bits, requested_max), orders in observed_orders.items():
        assert orders == set(range(requested_max + 1))

    # 64 double epsilons is the implementation's certified upward-recursion
    # target; allow a small margin for the libm erf/exp implementations.
    assert worst_relative_error <= Decimal("2e-14")
