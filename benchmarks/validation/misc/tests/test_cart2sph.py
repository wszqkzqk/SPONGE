import math
import os
import shlex
import shutil
import struct
import subprocess
import sys
from pathlib import Path

import pytest

REPOSITORY_ROOT = Path(__file__).resolve().parents[4]
PROBE_SOURCE = Path(__file__).with_name("cart2sph_probe.cpp")


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
        pytest.skip("a C++17 compiler is required for the Cart2Sph probe")
    return [compiler]


def _components(l):
    return [
        (lx, ly, l - lx - ly)
        for lx in range(l, -1, -1)
        for ly in range(l - lx, -1, -1)
    ]


def _gaussian_monomial_overlap_factor(component_a, component_b):
    factor = 1
    for power_a, power_b in zip(component_a, component_b):
        total = power_a + power_b
        if total % 2:
            return 0
        for odd in range(total - 1, 0, -2):
            factor *= odd
    return factor


def _associated_legendre_without_condon_shortley(l, m, cosine):
    """Independent three-term recurrence, not the coefficient generator."""
    p_mm = 1.0
    root = math.sqrt(max(0.0, (1.0 - cosine) * (1.0 + cosine)))
    odd = 1.0
    for _ in range(m):
        p_mm *= odd * root
        odd += 2.0
    if l == m:
        return p_mm
    p_m1_m = (2 * m + 1) * cosine * p_mm
    if l == m + 1:
        return p_m1_m
    previous, current = p_mm, p_m1_m
    for degree in range(m + 2, l + 1):
        following = (
            (2 * degree - 1) * cosine * current
            - (degree + m - 1) * previous
        ) / (degree - m)
        previous, current = current, following
    return current


def _real_spherical_harmonic(l, m, x, y, z):
    radius = math.sqrt(x * x + y * y + z * z)
    cosine = z / radius
    phi = math.atan2(y, x)
    abs_m = abs(m)
    normalization = math.sqrt(
        (2 * l + 1)
        / (4 * math.pi)
        * math.factorial(l - abs_m)
        / math.factorial(l + abs_m)
    )
    associated = _associated_legendre_without_condon_shortley(
        l, abs_m, cosine
    )
    if m > 0:
        angular = math.sqrt(2.0) * math.cos(abs_m * phi)
    elif m < 0:
        angular = math.sqrt(2.0) * math.sin(abs_m * phi)
    else:
        angular = 1.0
    return radius**l * normalization * associated * angular


LEGACY_BLOCKS = {
    0: ((0.28209479,),),
    1: (
        (0.48860251, 0.0, 0.0),
        (0.0, 0.48860251, 0.0),
        (0.0, 0.0, 0.48860251),
    ),
    2: (
        (0.0, 0.0, -0.31539157, 0.0, 0.54627422),
        (1.09254843, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 1.09254843, 0.0),
        (0.0, 0.0, -0.31539157, 0.0, -0.54627422),
        (0.0, 1.09254843, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.63078313, 0.0, 0.0),
    ),
    3: (
        (0.0, 0.0, 0.0, 0.0, -0.45704580, 0.0, 0.59004359),
        (1.77013077, 0.0, -0.45704580, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, -1.11952900, 0.0, 1.44530572, 0.0),
        (0.0, 0.0, 0.0, 0.0, -0.45704580, 0.0, -1.77013077),
        (0.0, 2.89061144, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 1.82818320, 0.0, 0.0),
        (-0.59004359, 0.0, -0.45704580, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, -1.11952900, 0.0, -1.44530572, 0.0),
        (0.0, 0.0, 1.82818320, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.74635267, 0.0, 0.0, 0.0),
    ),
    4: (
        (0.0, 0.0, 0.0, 0.0, 0.31735664, 0.0, -0.47308735, 0.0, 0.62583574),
        (2.50334294, 0.0, -0.94617470, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, -2.00713963, 0.0, 1.77013077, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.63471328, 0.0, 0.0, 0.0, -3.75501441),
        (0.0, 5.31039231, 0.0, -2.00713963, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, -2.53885313, 0.0, 2.83852409, 0.0, 0.0),
        (-2.50334294, 0.0, -0.94617470, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, -2.00713963, 0.0, -5.31039231, 0.0),
        (0.0, 0.0, 5.67704817, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.0, 2.67618617, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.31735664, 0.0, 0.47308735, 0.0, 0.62583574),
        (0.0, -1.77013077, 0.0, -2.00713963, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, -2.53885313, 0.0, -2.83852409, 0.0, 0.0),
        (0.0, 0.0, 0.0, 2.67618617, 0.0, 0.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, 0.0, 0.84628438, 0.0, 0.0, 0.0, 0.0),
    ),
}


def _parse_probe(stdout):
    blocks = {}
    validated = False
    for line in stdout.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "BLOCK":
            l, dim_cart, dim_sph = map(int, fields[1:4])
            values = tuple(map(float, fields[4:]))
            assert len(values) == dim_cart * dim_sph
            blocks[l] = (
                dim_cart,
                dim_sph,
                tuple(
                    values[row * dim_sph : (row + 1) * dim_sph]
                    for row in range(dim_cart)
                ),
            )
        elif fields == ["VALIDATION", "1"]:
            validated = True
    assert validated
    return blocks


def test_cart2sph_is_general_normalized_and_high_l_safe(tmp_path):
    executable = tmp_path / "cart2sph_probe"
    compile_result = subprocess.run(
        [
            *_compiler_command(),
            "-std=c++17",
            "-O2",
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
        f"failed to compile Cart2Sph probe\nstdout:\n{compile_result.stdout}"
        f"\nstderr:\n{compile_result.stderr}"
    )
    run_result = subprocess.run(
        [str(executable)],
        capture_output=True,
        text=True,
        check=False,
        timeout=30,
    )
    assert run_result.returncode == 0, run_result.stdout + run_result.stderr
    blocks = _parse_probe(run_result.stdout)

    assert set(blocks) == set(range(7))
    for l, (dim_cart, dim_sph, block) in blocks.items():
        assert dim_cart == (l + 1) * (l + 2) // 2
        assert dim_sph == 2 * l + 1
        assert all(
            any(abs(block[row][column]) > 1.0e-7 for row in range(dim_cart))
            for column in range(dim_sph)
        )

    for l, expected in LEGACY_BLOCKS.items():
        actual = blocks[l][2]
        assert len(actual) == len(expected)
        for actual_row, expected_row in zip(actual, expected):
            assert actual_row == pytest.approx(expected_row, abs=1.0e-7)
            assert tuple(map(lambda value: struct.pack("f", value), actual_row)) == tuple(
                map(lambda value: struct.pack("f", value), expected_row)
            )

    unit_points = (
        (0.2672612419, -0.5345224838, 0.8017837257),
        (-0.5270462767, 0.7378647872, 0.4216370214),
        (0.8728715609, 0.4364357805, -0.2182178902),
        (-0.1690308509, -0.8451542547, -0.5070925528),
    )
    for l in range(2, 7):
        dim_cart, dim_sph, block = blocks[l]
        components = _components(l)
        assert len(components) == dim_cart
        for column, m in enumerate(range(-l, l + 1)):
            for x, y, z in unit_points:
                polynomial = sum(
                    block[row][column] * x**lx * y**ly * z**lz
                    for row, (lx, ly, lz) in enumerate(components)
                )
                reference = _real_spherical_harmonic(l, m, x, y, z)
                assert polynomial == pytest.approx(reference, abs=3.0e-6)

            laplacian = {}
            for row, (lx, ly, lz) in enumerate(components):
                coefficient = block[row][column]
                for axis_power, target in (
                    (lx, (lx - 2, ly, lz)),
                    (ly, (lx, ly - 2, lz)),
                    (lz, (lx, ly, lz - 2)),
                ):
                    if axis_power >= 2:
                        laplacian[target] = laplacian.get(target, 0.0) + (
                            coefficient * axis_power * (axis_power - 1)
                        )
            assert max(map(abs, laplacian.values()), default=0.0) < 2.0e-4

    for l in range(7):
        dim_cart, dim_sph, block = blocks[l]
        components = _components(l)
        expected_diagonal = math.prod(range(1, 2 * l + 2, 2)) / (
            4.0 * math.pi
        )
        for left in range(dim_sph):
            for right in range(dim_sph):
                overlap = sum(
                    block[row_a][left]
                    * block[row_b][right]
                    * _gaussian_monomial_overlap_factor(
                        components[row_a], components[row_b]
                    )
                    for row_a in range(dim_cart)
                    for row_b in range(dim_cart)
                )
                reference = expected_diagonal if left == right else 0.0
                assert overlap == pytest.approx(
                    reference, abs=max(2.0e-5, expected_diagonal * 2.0e-6)
                )
