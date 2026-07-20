import json
import math
import os
import subprocess
import textwrap
from pathlib import Path

import pytest


def _toml_string(value):
    return json.dumps(os.fspath(value))


def _gro_atom(resid, resname, atomname, atomnr, xyz):
    x, y, z = xyz
    return (
        f"{resid:5d}{resname:<5}{atomname:>5}{atomnr:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )


def _extract_mdout_term(mdout_path, term_name):
    lines = Path(mdout_path).read_text().splitlines()
    headers = lines[0].split()
    values = lines[1].split()
    assert len(headers) == len(values)
    return float(dict(zip(headers, values))[term_name])


def test_direct_gromacs_topgro_accepts_settles_constraints_and_cmap(tmp_path):
    top_path = tmp_path / "topol.top"
    gro_path = tmp_path / "conf.gro"
    mdin_path = tmp_path / "mdin.spg.toml"

    top_path.write_text(
        textwrap.dedent(
            """
            [ defaults ]
            1 2 yes 0.5 0.833333

            [ atomtypes ]
            OW OW 15.9994 -0.834 A 0.315075 0.636386
            HW HW 1.008 0.417 A 0.0 0.0
            CT CT 12.011 0.0 A 0.34 0.276144

            [ cmaptypes ]
            CT CT CT CT CT 1 2 2 0.0 0.0 0.0 0.0

            [ moleculetype ]
            SOL 2

            [ atoms ]
            1 OW 1 SOL OW 1 -0.834 15.9994
            2 HW 1 SOL HW1 2 0.417 1.008
            3 HW 1 SOL HW2 3 0.417 1.008

            [ settles ]
            1 1 0.09572 0.15139

            [ moleculetype ]
            LIN 1

            [ atoms ]
            1 CT 1 LIN C1 1 0.0 12.011
            2 CT 1 LIN C2 2 0.0 12.011

            [ constraints ]
            1 2 1 0.150

            [ moleculetype ]
            CMAP 1

            [ atoms ]
            1 CT 1 CMP C1 1 0.0 12.011
            2 CT 1 CMP C2 2 0.0 12.011
            3 CT 1 CMP C3 3 0.0 12.011
            4 CT 1 CMP C4 4 0.0 12.011
            5 CT 1 CMP C5 5 0.0 12.011

            [ cmap ]
            1 2 3 4 5 1

            [ system ]
            direct feature smoke

            [ molecules ]
            SOL 1
            LIN 1
            CMAP 1
            """
        ).strip()
        + "\n"
    )

    gro_lines = [
        "direct feature smoke",
        "10",
        _gro_atom(1, "SOL", "OW", 1, (0.000, 0.000, 0.000)),
        _gro_atom(1, "SOL", "HW1", 2, (0.096, 0.000, 0.000)),
        _gro_atom(1, "SOL", "HW2", 3, (-0.024, 0.093, 0.000)),
        _gro_atom(2, "LIN", "C1", 4, (0.500, 0.000, 0.000)),
        _gro_atom(2, "LIN", "C2", 5, (0.650, 0.000, 0.000)),
        _gro_atom(3, "CMP", "C1", 6, (1.000, 0.000, 0.000)),
        _gro_atom(3, "CMP", "C2", 7, (1.150, 0.050, 0.020)),
        _gro_atom(3, "CMP", "C3", 8, (1.300, 0.000, 0.070)),
        _gro_atom(3, "CMP", "C4", 9, (1.450, -0.040, 0.030)),
        _gro_atom(3, "CMP", "C5", 10, (1.600, 0.020, -0.050)),
        "   5.00000   5.00000   5.00000",
    ]
    gro_path.write_text("\n".join(gro_lines) + "\n")

    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "direct_gromacs_feature_smoke"
            mode = "nve"
            step_limit = 0
            dt = 0
            cutoff = 8.0
            constrain_mode = "SETTLE"
            gromacs_top = {_toml_string(top_path)}
            gromacs_gro = {_toml_string(gro_path)}
            default_out_file_prefix = {_toml_string(tmp_path / "direct_feature")}
            print_zeroth_frame = 1
            write_mdout_interval = 1
            """
        ).strip()
        + "\n"
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    output = result.stdout + "\n" + result.stderr
    assert "constrain pair number is 4" in output
    assert "rigid triangle numbers is 1" in output
    assert "rigid pair numbers is 1" in output


def test_direct_gromacs_harmonic_improper_uses_gromacs_half_k(tmp_path):
    top_path = tmp_path / "topol.top"
    gro_path = tmp_path / "conf.gro"
    mdin_path = tmp_path / "mdin.spg.toml"
    mdout_path = tmp_path / "mdout.txt"

    top_path.write_text(
        textwrap.dedent(
            """
            [ defaults ]
            1 2 yes 1.0 1.0

            [ atomtypes ]
            A A 12.0 0.0 A 0.0 0.0
            B B 12.0 0.0 A 0.0 0.0
            C C 12.0 0.0 A 0.0 0.0
            D D 12.0 0.0 A 0.0 0.0

            [ dihedraltypes ]
            A B C D 2 0.0 8.368

            [ moleculetype ]
            IMPLICIT 0

            [ atoms ]
            1 A 1 IMP A 1 0.0 12.0
            2 B 1 IMP B 2 0.0 12.0
            3 C 1 IMP C 3 0.0 12.0
            4 D 1 IMP D 4 0.0 12.0

            [ dihedrals ]
            1 2 3 4 2

            [ moleculetype ]
            EXPLICIT 0

            [ atoms ]
            1 A 1 EXP A 1 0.0 12.0
            2 B 1 EXP B 2 0.0 12.0
            3 C 1 EXP C 3 0.0 12.0
            4 D 1 EXP D 4 0.0 12.0

            [ dihedrals ]
            1 2 3 4 2 0.0 8.368

            [ system ]
            harmonic improper convention

            [ molecules ]
            IMPLICIT 1
            EXPLICIT 1
            """
        ).strip()
        + "\n"
    )

    gro_lines = [
        "harmonic improper convention",
        "8",
        _gro_atom(1, "IMP", "A", 1, (0.000, 0.100, 0.000)),
        _gro_atom(1, "IMP", "B", 2, (0.000, 0.000, 0.000)),
        _gro_atom(1, "IMP", "C", 3, (0.100, 0.000, 0.000)),
        _gro_atom(1, "IMP", "D", 4, (0.100, 0.000, 0.100)),
        _gro_atom(2, "EXP", "A", 5, (1.000, 0.100, 0.000)),
        _gro_atom(2, "EXP", "B", 6, (1.000, 0.000, 0.000)),
        _gro_atom(2, "EXP", "C", 7, (1.100, 0.000, 0.000)),
        _gro_atom(2, "EXP", "D", 8, (1.100, 0.000, 0.100)),
        "   5.00000   5.00000   5.00000",
    ]
    gro_path.write_text("\n".join(gro_lines) + "\n")

    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "direct_gromacs_harmonic_improper"
            mode = "nve"
            step_limit = 0
            dt = 0
            cutoff = 8.0
            gromacs_top = {_toml_string(top_path)}
            gromacs_gro = {_toml_string(gro_path)}
            mdout = {_toml_string(mdout_path)}
            print_zeroth_frame = 1
            write_mdout_interval = 1
            """
        ).strip()
        + "\n"
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    output = result.stdout + "\n" + result.stderr
    assert "dihedral_numbers is 2" in output

    # Each topology has phi=pi/2 and k=8.368 kJ/mol/rad^2 = 2 kcal/mol/rad^2.
    # GROMACS funct=2 is E=1/2*k*(phi-phi0)^2.
    expected = 2.0 * 0.5 * (8.368 / 4.184) * (math.pi / 2.0) ** 2
    assert _extract_mdout_term(
        mdout_path, "improper_dihedral"
    ) == pytest.approx(expected, abs=0.02)


def test_direct_gromacs_warns_once_for_unknown_section_with_location(tmp_path):
    top_path = tmp_path / "topol.top"
    include_path = tmp_path / "unsupported.itp"
    gro_path = tmp_path / "conf.gro"
    mdin_path = tmp_path / "mdin.spg.toml"

    include_path.write_text(
        textwrap.dedent(
            """
            [ position_restraints ]
            1 1 1000 1000 1000

            [ position_restraints ]
            1 1 900 900 900
            """
        ).strip()
        + "\n"
    )
    top_path.write_text(
        textwrap.dedent(
            """
            [ defaults ]
            1 2 yes 1.0 1.0

            [ atomtypes ]
            CT CT 12.0 0.0 A 0.0 0.0

            #include "unsupported.itp"

            [ moleculetype ]
            ONE 0

            [ atoms ]
            1 CT 1 ONE C 1 0.0 12.0

            [ system ]
            unknown section warning

            [ molecules ]
            ONE 1
            """
        ).strip()
        + "\n"
    )
    gro_path.write_text(
        "\n".join(
            [
                "unknown section warning",
                "1",
                _gro_atom(1, "ONE", "C", 1, (0.000, 0.000, 0.000)),
                "   5.00000   5.00000   5.00000",
            ]
        )
        + "\n"
    )
    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "direct_gromacs_unknown_section"
            mode = "nve"
            step_limit = 0
            dt = 0
            cutoff = 8.0
            gromacs_top = {_toml_string(top_path)}
            gromacs_gro = {_toml_string(gro_path)}
            print_zeroth_frame = 1
            write_mdout_interval = 1
            """
        ).strip()
        + "\n"
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    output = result.stdout + "\n" + result.stderr
    warning = "unsupported GROMACS topology section [ position_restraints ]"
    assert output.count(warning) == 1
    assert f"{include_path}:1" in output
    assert "unsupported GROMACS topology section [ system ]" not in output


def test_direct_gromacs_applies_nonbond_params_to_distinct_named_types(
    tmp_path,
):
    top_path = tmp_path / "topol.top"
    gro_path = tmp_path / "conf.gro"
    mdin_path = tmp_path / "mdin.spg.toml"
    mdout_path = tmp_path / "mdout.txt"

    top_path.write_text(
        textwrap.dedent(
            """
            [ defaults ]
            1 2 yes 1.0 1.0

            [ atomtypes ]
            A A 12.0 0.0 A 0.3 0.4184
            B B 12.0 0.0 A 0.3 0.4184

            [ nonbond_params ]
            A B 1 0.35 0.4184
            A B 1 0.4 0.8368

            [ moleculetype ]
            PAIR 0

            [ atoms ]
            1 A 1 PAIR A 1 0.0 12.0
            2 B 1 PAIR B 2 0.0 12.0

            [ system ]
            nonbond parameter override

            [ molecules ]
            PAIR 1
            """
        ).strip()
        + "\n"
    )
    gro_path.write_text(
        "\n".join(
            [
                "nonbond parameter override",
                "2",
                _gro_atom(1, "PAIR", "A", 1, (0.000, 0.000, 0.000)),
                _gro_atom(1, "PAIR", "B", 2, (0.500, 0.000, 0.000)),
                "   5.00000   5.00000   5.00000",
            ]
        )
        + "\n"
    )
    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "direct_gromacs_nonbond_params"
            mode = "nve"
            step_limit = 0
            dt = 0
            cutoff = 8.0
            gromacs_top = {_toml_string(top_path)}
            gromacs_gro = {_toml_string(gro_path)}
            mdout = {_toml_string(mdout_path)}
            print_zeroth_frame = 1
            write_mdout_interval = 1
            """
        ).strip()
        + "\n"
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    output = result.stdout + "\n" + result.stderr
    assert "atom_LJ_type_number is 2" in output

    sigma_over_r = 0.4 / 0.5
    expected = 4.0 * (0.8368 / 4.184) * (sigma_over_r**12 - sigma_over_r**6)
    assert _extract_mdout_term(mdout_path, "LJ_short") == pytest.approx(
        expected, abs=0.01
    )


def test_direct_gromacs_accepts_triclinic_gro_box(tmp_path):
    top_path = tmp_path / "topol.top"
    gro_path = tmp_path / "conf.gro"
    mdin_path = tmp_path / "mdin.spg.toml"
    box_path = tmp_path / "box.txt"

    top_path.write_text(
        textwrap.dedent(
            """
            [ defaults ]
            1 2 yes 1.0 1.0

            [ atomtypes ]
            CT CT 12.0 0.0 A 0.0 0.0

            [ moleculetype ]
            ONE 0

            [ atoms ]
            1 CT 1 ONE C 1 0.0 12.0

            [ system ]
            triclinic box

            [ molecules ]
            ONE 1
            """
        ).strip()
        + "\n"
    )

    # a=(2.0, 0.0, 0.0), b=(0.3, 2.5, 0.0), c=(0.5, 0.6, 3.0) nm,
    # serialized in the GRO-defined component order.
    gro_path.write_text(
        "\n".join(
            [
                "triclinic box",
                "1",
                _gro_atom(1, "ONE", "C", 1, (0.100, 0.200, 0.300)),
                "2.0 2.5 3.0 0.0 0.0 0.3 0.0 0.5 0.6",
            ]
        )
        + "\n"
    )
    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "direct_gromacs_triclinic_box"
            mode = "nve"
            step_limit = 1
            dt = 0
            cutoff = 8.0
            gromacs_top = {_toml_string(top_path)}
            gromacs_gro = {_toml_string(gro_path)}
            box = {_toml_string(box_path)}
            print_zeroth_frame = 1
            write_trajectory_interval = 1
            write_mdout_interval = 1
            """
        ).strip()
        + "\n"
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr

    actual = [
        float(value) for value in box_path.read_text().splitlines()[0].split()
    ]
    a = (2.0, 0.0, 0.0)
    b = (0.3, 2.5, 0.0)
    c = (0.5, 0.6, 3.0)

    def norm(vector):
        return math.sqrt(sum(value * value for value in vector))

    def angle(first, second):
        cosine = (
            sum(x * y for x, y in zip(first, second))
            / norm(first)
            / norm(second)
        )
        return math.degrees(math.acos(cosine))

    expected = [
        10.0 * norm(a),
        10.0 * norm(b),
        10.0 * norm(c),
        angle(b, c),
        angle(a, c),
        angle(a, b),
    ]
    assert actual == pytest.approx(expected, abs=2.0e-4)


@pytest.mark.parametrize(
    "box_line,error_fragment",
    [
        ("nan 2.0 3.0", "invalid numeric field"),
        ("0.0 2.0 3.0", "invalid GROMACS gro box geometry"),
        (
            "1.0 1.0 0.0 0.0 0.0 0.0 0.0 1.0 1.0",
            "invalid GROMACS gro box geometry",
        ),
        (
            "2.0 2.5 3.0 0.1 0.0 0.3 0.0 0.5 0.6",
            "unsupported non-canonical GROMACS triclinic box orientation",
        ),
    ],
    ids=["non_finite", "zero_vector", "coplanar", "rotated_orientation"],
)
def test_direct_gromacs_rejects_invalid_gro_box(
    tmp_path, box_line, error_fragment
):
    top_path = tmp_path / "topol.top"
    gro_path = tmp_path / "conf.gro"
    mdin_path = tmp_path / "mdin.spg.toml"

    top_path.write_text(
        textwrap.dedent(
            """
            [ defaults ]
            1 2 yes 1.0 1.0

            [ atomtypes ]
            CT CT 12.0 0.0 A 0.0 0.0

            [ moleculetype ]
            ONE 0

            [ atoms ]
            1 CT 1 ONE C 1 0.0 12.0

            [ system ]
            invalid box

            [ molecules ]
            ONE 1
            """
        ).strip()
        + "\n"
    )
    gro_path.write_text(
        "\n".join(
            [
                "invalid box",
                "1",
                _gro_atom(1, "ONE", "C", 1, (0.000, 0.000, 0.000)),
                box_line,
            ]
        )
        + "\n"
    )
    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "direct_gromacs_invalid_box"
            mode = "nve"
            step_limit = 0
            dt = 0
            cutoff = 8.0
            gromacs_top = {_toml_string(top_path)}
            gromacs_gro = {_toml_string(gro_path)}
            """
        ).strip()
        + "\n"
    )

    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output
