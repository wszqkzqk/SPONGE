import math
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
PROBE_SOURCE = Path(__file__).with_name("becke_response_probe.cpp")


def _compiler_command():
    configured = os.environ.get("CXX")
    if configured:
        return shlex.split(configured)
    if sys.platform == "darwin" and Path("/usr/bin/clang++").is_file():
        return ["/usr/bin/clang++"]
    compiler = shutil.which("c++") or shutil.which("clang++") or shutil.which(
        "g++"
    )
    if compiler is None:
        pytest.skip("a C++17 compiler is required for the Becke response probe")
    return [compiler]


@pytest.fixture(scope="module")
def becke_probe_output(tmp_path_factory):
    executable = tmp_path_factory.mktemp("becke-probe") / "becke_response_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-O3",
            "-ffast-math",
            f"-I{REPOSITORY_ROOT / 'SPONGE'}",
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
        f"failed to compile Becke response probe\nstdout:\n"
        f"{compile_result.stdout}\nstderr:\n{compile_result.stderr}"
    )
    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr
    return run_result.stdout.splitlines()


def _from_bits(bits):
    return struct.unpack(">d", bytes.fromhex(bits))[0]


def _shape_oracle(mu):
    with localcontext() as context:
        context.prec = 120
        value = Decimal.from_float(mu)
        if value <= -1:
            return 1.0, 0.0
        if value >= 1:
            return 0.0, 0.0
        use_complement = value < 0
        delta = 1 + value if use_complement else 1 - value
        derivative = Decimal(1) if use_complement else Decimal(-1)
        for _ in range(3):
            old_delta = delta
            derivative *= Decimal("1.5") * old_delta * (2 - old_delta)
            delta = Decimal("0.5") * old_delta * old_delta * (3 - old_delta)
        if use_complement:
            return float(1 - Decimal("0.5") * delta), float(
                -Decimal("0.5") * derivative
            )
        return float(Decimal("0.5") * delta), float(
            Decimal("0.5") * derivative
        )


def test_shape_and_derivative_are_stable_at_both_endpoints(becke_probe_output):
    lines = [line.split() for line in becke_probe_output if line.startswith("shape ")]
    assert len(lines) == 9
    observed = {}
    for _, bits, direct_text, combined_text, derivative_text in lines:
        mu = _from_bits(bits)
        direct = float(direct_text)
        combined = float(combined_text)
        derivative = float(derivative_text)
        expected_shape, expected_derivative = _shape_oracle(mu)
        assert direct == combined
        assert math.isclose(direct, expected_shape, rel_tol=3.0e-15, abs_tol=0.0)
        assert math.isclose(
            derivative, expected_derivative, rel_tol=5.0e-15, abs_tol=0.0
        )
        assert 0.0 <= direct <= 1.0
        assert derivative <= 0.0
        observed[mu] = (direct, derivative)

    for mu in [value for value in observed if 0.0 < value < 1.0]:
        positive_shape, positive_derivative = observed[mu]
        negative_shape, negative_derivative = observed[-mu]
        assert negative_shape + positive_shape == 1.0
        assert math.isclose(
            negative_derivative, positive_derivative, rel_tol=4.0e-15, abs_tol=0.0
        )


def test_size_adjustment_handles_full_finite_radius_range(becke_probe_output):
    lines = [line.split() for line in becke_probe_output if line.startswith("radius ")]
    assert len(lines) == 7
    coefficients = []
    with localcontext() as context:
        context.prec = 120
        for _, radius_a_text, radius_b_text, coefficient_text in lines:
            radius_a = Decimal(radius_a_text)
            radius_b = Decimal(radius_b_text)
            raw = (radius_b / radius_a - radius_a / radius_b) / 4
            expected = float(max(Decimal("-0.5"), min(Decimal("0.5"), raw)))
            coefficient = float(coefficient_text)
            assert math.isfinite(coefficient)
            assert math.isclose(coefficient, expected, rel_tol=2.0e-16, abs_tol=0.0)
            coefficients.append(coefficient)
    assert coefficients[0] == 0.0
    assert coefficients[-2:] == [-0.5, 0.5]
    assert coefficients[1] == -coefficients[2]
    assert coefficients[3] == -coefficients[4]


def test_owner_and_nonowner_weight_response_matches_central_difference(
    becke_probe_output,
):
    lines = [
        line.split()
        for line in becke_probe_output
        if line.startswith("weight-derivative ")
    ]
    assert len(lines) == 18
    assert {int(fields[1]) for fields in lines} == {0, 1}
    for fields in lines:
        analytic = float(fields[4])
        plus = float(fields[5])
        minus = float(fields[6])
        step = float(fields[7])
        finite_difference = (plus - minus) / (2.0 * step)
        assert math.isclose(
            analytic, finite_difference, rel_tol=2.0e-8, abs_tol=3.0e-10
        ), (fields, analytic, finite_difference)

    translation_lines = [
        line.split()
        for line in becke_probe_output
        if line.startswith("translation ")
    ]
    assert len(translation_lines) == 2
    for fields in translation_lines:
        assert max(abs(float(value)) for value in fields[2:]) < 2.0e-16


def test_each_pair_response_is_exactly_translation_invariant(becke_probe_output):
    lines = [line.split() for line in becke_probe_output if line.startswith("pair ")]
    assert len(lines) == 3
    for _, _axis, response_a, response_b, response_grid in lines:
        # The production helper defines the moving-grid response from the
        # already-rounded nuclear pair sum, so this cancellation is exact.
        assert (float(response_a) + float(response_b)) + float(response_grid) == 0.0
