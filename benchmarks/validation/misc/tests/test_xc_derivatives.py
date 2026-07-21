import math
import os
import shlex
import shutil
import struct
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("xc_derivative_probe.cpp")


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
        pytest.skip("a C++17 compiler is required for the XC derivative probe")
    return [compiler]


def _dependency_include():
    candidates = []
    if os.environ.get("CONDA_PREFIX"):
        candidates.append(Path(os.environ["CONDA_PREFIX"]) / "include")
    candidates.append(REPOSITORY_ROOT / ".pixi" / "envs" / "dev-cpu" / "include")
    for candidate in candidates:
        if (candidate / "omp.h").is_file() and (candidate / "fftw3.h").is_file():
            return candidate
    pytest.skip("OpenMP and FFTW headers are required for the XC derivative probe")


@pytest.fixture(scope="module")
def xc_probe_output(tmp_path_factory):
    executable = tmp_path_factory.mktemp("xc-probe") / "xc_derivative_probe"
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
        f"failed to compile XC derivative probe\nstdout:\n"
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


def _close(actual, reference, *, rtol=3.0e-8, atol=2.0e-11):
    assert math.isclose(actual, reference, rel_tol=rtol, abs_tol=atol), (
        actual,
        reference,
    )


def test_rks_and_uks_derivatives_are_analytic(xc_probe_output):
    rks_lines = [line.split() for line in xc_probe_output if line.startswith("rks ")]
    uks_lines = [line.split() for line in xc_probe_output if line.startswith("uks ")]
    assert len(rks_lines) == len(uks_lines) == 5

    for fields in rks_lines:
        _, _method, _energy, vrho, vrho_reference, vsigma, vsigma_reference = fields
        _close(float(vrho), float(vrho_reference))
        _close(float(vsigma), float(vsigma_reference))

    for fields in uks_lines:
        assert len(fields) == 13
        for component in range(5):
            analytical = float(fields[3 + 2 * component])
            reference = float(fields[4 + 2 * component])
            _close(analytical, reference, rtol=8.0e-8, atol=5.0e-11)


def test_spin_symmetry_and_closed_shell_identity(xc_probe_output):
    closed_lines = [
        line.split() for line in xc_probe_output if line.startswith("closed ")
    ]
    swap_lines = [line.split() for line in xc_probe_output if line.startswith("swap ")]
    assert len(closed_lines) == len(swap_lines) == 5
    for fields in closed_lines:
        values = [float(value) for value in fields[2:]]
        for left, right in zip(values[0::2], values[1::2]):
            _close(left, right, rtol=2.0e-13, atol=2.0e-14)
    for fields in swap_lines:
        values = [float(value) for value in fields[2:]]
        _close(values[0], values[1], rtol=2.0e-13, atol=2.0e-14)
        _close(values[2], values[3], rtol=2.0e-13, atol=2.0e-14)
        _close(values[4], values[5], rtol=2.0e-13, atol=2.0e-14)
        assert abs(values[6]) <= 2.0e-13


def test_every_positive_tail_density_is_evaluated(xc_probe_output):
    tail_lines = [line.split() for line in xc_probe_output if line.startswith("tail ")]
    assert len(tail_lines) == 30
    assert {fields[-1] for fields in tail_lines} == {"111"}

    vwn_lines = [
        line.split() for line in xc_probe_output if line.startswith("vwn-tail ")
    ]
    assert len(vwn_lines) == 6
    # The combined VWN low-density series starts with c2/rs.  Verify that the
    # rho=1e-300 value is retained (not rounded to zero by cancellation) and
    # has the correct leading coefficient.
    _, rho_text, eps_text, deriv_text = vwn_lines[-1]
    rho = float(rho_text)
    eps = float(eps_text)
    derivative = float(deriv_text)
    rs_factor = (3.0 / (4.0 * float("3.1415927410125732"))) ** (1.0 / 3.0)
    expected = -0.41433042034046384 * rho ** (1.0 / 3.0) / rs_factor
    _close(eps, expected, rtol=2.0e-14, atol=0.0)
    _close(derivative, expected / (3.0 * rho), rtol=3.0e-14, atol=0.0)


def test_true_spin_endpoint_and_fast_math_contract(xc_probe_output):
    endpoint_lines = [
        line.split() for line in xc_probe_output if line.startswith("endpoint ")
    ]
    assert len(endpoint_lines) == 2
    for fields in endpoint_lines:
        assert int(fields[2]) == 2  # QC_XC_UKS_VRHO_B_BIT
        assert fields[5].lower() == "7ff0000000000000"
        assert all(math.isfinite(float(value)) for value in fields[3:5])
        assert all(math.isfinite(float(value)) for value in fields[6:])

    zero_gradient_lines = [
        line.split()
        for line in xc_probe_output
        if line.startswith("endpoint-zero-gradient ")
    ]
    assert len(zero_gradient_lines) == 2
    for fields in zero_gradient_lines:
        assert int(fields[2]) == 0
        assert all(math.isfinite(float(value)) for value in fields[3:])

    valid_line = next(line for line in xc_probe_output if line.startswith("valid "))
    assert valid_line.split()[1] == "10010010"
    gram_line = next(line for line in xc_probe_output if line.startswith("gram "))
    assert gram_line.split()[1] == "11001"
    gram_counterexample = next(
        line for line in xc_probe_output if line.startswith("gram-counterexample ")
    )
    assert gram_counterexample.split()[1] == "1100"

    exact_lines = [
        line.split() for line in xc_probe_output if line.startswith("gram-exact ")
    ]
    assert len(exact_lines) == 256
    for _, aa_bits, ab_bits, bb_bits, result in exact_lines:
        aa = struct.unpack(">d", bytes.fromhex(aa_bits))[0]
        ab = struct.unpack(">d", bytes.fromhex(ab_bits))[0]
        bb = struct.unpack(">d", bytes.fromhex(bb_bits))[0]
        expected = Fraction.from_float(ab) ** 2 <= (
            Fraction.from_float(aa) * Fraction.from_float(bb)
        )
        assert bool(int(result)) is expected

    builder_lines = [
        line.split()
        for line in xc_probe_output
        if line.startswith("gram-builder ")
    ]
    assert len(builder_lines) == 256

    def from_bits(bits):
        return struct.unpack(">d", bytes.fromhex(bits))[0]

    for fields in builder_lines:
        assert len(fields) == 11
        vectors = [from_bits(bits) for bits in fields[1:7]]
        assert fields[7] == "11"
        sigma_aa, sigma_ab, sigma_bb = [
            from_bits(bits) for bits in fields[8:11]
        ]
        ga = [Fraction.from_float(value) for value in vectors[:3]]
        gb = [Fraction.from_float(value) for value in vectors[3:]]
        exact_aa = sum(value * value for value in ga)
        exact_ab = sum(left * right for left, right in zip(ga, gb))
        exact_bb = sum(value * value for value in gb)
        published_aa = Fraction.from_float(sigma_aa)
        published_ab = Fraction.from_float(sigma_ab)
        published_bb = Fraction.from_float(sigma_bb)
        assert all(math.isfinite(value) for value in (sigma_aa, sigma_ab, sigma_bb))
        assert published_aa >= exact_aa
        assert published_bb >= exact_bb
        assert abs(published_ab) <= abs(exact_ab)
        if published_ab:
            assert (published_ab > 0) is (exact_ab > 0)
        assert published_ab * published_ab <= published_aa * published_bb

    for prefix in ("tiny-sigma ", "tiny-sigma-uks "):
        lines = [line.split() for line in xc_probe_output if line.startswith(prefix)]
        assert len(lines) == 2
        for fields in lines:
            values = [float(value) for value in fields[2:]]
            assert all(math.isfinite(value) for value in values)
            for zero_value, tiny_value in zip(values[0::2], values[1::2]):
                _close(zero_value, tiny_value, rtol=2.0e-13, atol=2.0e-14)

    zero_tail_lines = [
        line.split() for line in xc_probe_output if line.startswith("zero-tail ")
    ]
    assert len(zero_tail_lines) == 4
    for fields in zero_tail_lines:
        method = int(fields[1])
        rho = float(fields[2])
        energy_bits, vrho_bits, vsigma_bits = fields[3:]
        assert int(energy_bits, 16) & 0x7FF0000000000000 != 0x7FF0000000000000
        assert int(vrho_bits, 16) & 0x7FF0000000000000 != 0x7FF0000000000000
        if rho == 1.0e-235 and method == 2:  # cancellation remains representable
            assert int(vsigma_bits, 16) & 0x7FF0000000000000 != 0x7FF0000000000000
        else:
            assert vsigma_bits.lower() == "7ff0000000000000"

    vacuum_lines = [
        line.split() for line in xc_probe_output if line.startswith("vacuum ")
    ]
    assert len(vacuum_lines) == 5
    for fields in vacuum_lines:
        assert [float(value) for value in fields[2:]] == [0.0, 0.0, 0.0]
