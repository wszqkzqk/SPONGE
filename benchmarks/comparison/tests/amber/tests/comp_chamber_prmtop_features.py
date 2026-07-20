import math
import os
import subprocess
import textwrap
from pathlib import Path


def _flag(name, fmt, values, items_per_line=10):
    lines = [f"%FLAG {name}", f"%FORMAT({fmt})"]
    rendered = [str(value) for value in values]
    for start in range(0, len(rendered), items_per_line):
        lines.append(" ".join(rendered[start : start + items_per_line]))
    return "\n".join(lines) + "\n"


def _write_chamber_prmtop(
    path,
    *,
    include_lj14_a=True,
    include_lj14_b=True,
    normal_a=729.0,
    normal_b=27.0,
):
    pointers = [0] * 31
    pointers[0] = 4  # NATOM
    pointers[1] = 1  # NTYPES
    pointers[7] = 1  # MPHIA: one dihedral without hydrogen
    pointers[11] = 1  # NRES
    pointers[17] = 1  # NPTRA
    pointers[27] = 1  # IFBOX

    sections = [
        "%VERSION VERSION_STAMP = V0001.000 DATE = 01/01/70 00:00:00\n",
        _flag("POINTERS", "10I8", pointers),
        _flag("CHARGE", "5E16.8", [0.0] * 4, 5),
        _flag("MASS", "5E16.8", [12.0] * 4, 5),
        _flag("ATOM_TYPE_INDEX", "10I8", [1, 1, 1, 1]),
        _flag("NUMBER_EXCLUDED_ATOMS", "10I8", [3, 2, 1, 0]),
        _flag("RESIDUE_POINTER", "10I8", [1]),
        _flag("DIHEDRAL_FORCE_CONSTANT", "5E16.8", [0.0], 5),
        _flag("DIHEDRAL_PERIODICITY", "5E16.8", [1.0], 5),
        _flag("DIHEDRAL_PHASE", "5E16.8", [0.0], 5),
        _flag("SCEE_SCALE_FACTOR", "5E16.8", [1.0], 5),
        _flag("SCNB_SCALE_FACTOR", "5E16.8", [2.0], 5),
        _flag("CHARMM_NUM_IMPROPERS", "10I8", [1]),
        _flag("CHARMM_IMPROPERS", "10I8", [1, 2, 3, 4, 1]),
        _flag("CHARMM_NUM_IMPR_TYPES", "1I8", [1]),
        _flag("CHARMM_IMPROPER_FORCE_CONSTANT", "5E16.8", [2.0], 5),
        _flag("CHARMM_IMPROPER_PHASE", "5E16.8", [0.0], 5),
        # Ordinary and special matrices deliberately differ. At r=sqrt(3) A,
        # ordinary LJ is zero, while LJ14 is -1 kcal/mol before SCNB=2.
        _flag("LENNARD_JONES_ACOEF", "3E24.16", [normal_a], 3),
        _flag("LENNARD_JONES_BCOEF", "3E24.16", [normal_b], 3),
    ]
    if include_lj14_a:
        sections.append(_flag("LENNARD_JONES_14_ACOEF", "3E24.16", [729.0], 3))
    if include_lj14_b:
        sections.append(_flag("LENNARD_JONES_14_BCOEF", "3E24.16", [54.0], 3))
    sections.extend(
        [
            _flag("DIHEDRALS_WITHOUT_HYDROGEN", "10I8", [0, 3, 6, 9, 1]),
            _flag("EXCLUDED_ATOMS_LIST", "10I8", [2, 3, 4, 3, 4, 4]),
        ]
    )
    path.write_text("".join(sections))


def _write_rst7(path):
    path.write_text(
        textwrap.dedent(
            """
            synthetic CHAMBER feature test
                 4
              0.0000000  1.0000000  0.0000000  0.0000000  0.0000000  0.0000000
              1.0000000  0.0000000  0.0000000  1.0000000  0.0000000  1.0000000
             50.0000000 50.0000000 50.0000000
             90.0000000 90.0000000 90.0000000
            """
        ).lstrip()
    )


def _write_mdin(path, prmtop_path, rst7_path, mdout_path):
    path.write_text(
        textwrap.dedent(
            f"""
            md_name = "synthetic_chamber_features"
            mode = "nve"
            step_limit = 0
            dt = 0
            cutoff = 8.0
            amber_parm7 = {str(prmtop_path)!r}
            amber_rst7 = {str(rst7_path)!r}
            mdout = {str(mdout_path)!r}
            print_zeroth_frame = 1
            write_mdout_interval = 1
            """
        ).strip()
        + "\n"
    )


def _run_sponge(
    tmp_path,
    *,
    include_lj14_a=True,
    include_lj14_b=True,
    normal_a=729.0,
    normal_b=27.0,
):
    prmtop_path = tmp_path / "synthetic.prmtop"
    rst7_path = tmp_path / "synthetic.rst7"
    mdin_path = tmp_path / "mdin.spg.toml"
    mdout_path = tmp_path / "mdout.txt"
    _write_chamber_prmtop(
        prmtop_path,
        include_lj14_a=include_lj14_a,
        include_lj14_b=include_lj14_b,
        normal_a=normal_a,
        normal_b=normal_b,
    )
    _write_rst7(rst7_path)
    _write_mdin(mdin_path, prmtop_path, rst7_path, mdout_path)
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    return result, mdout_path


def _extract_mdout_term(mdout_path, term_name):
    lines = Path(mdout_path).read_text().splitlines()
    headers = lines[0].split()
    values = lines[1].split()
    assert len(headers) == len(values)
    return float(dict(zip(headers, values))[term_name])


def test_chamber_lj14_and_improper_are_loaded(tmp_path):
    result, mdout_path = _run_sponge(tmp_path)
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert "non-bond 14 numbers is 1" in output
    assert (
        "START INITIALIZING IMPROPER DIHEDRAL (Xponge::system):\n"
        "    dihedral_numbers is 1"
    ) in output

    assert math.isclose(
        _extract_mdout_term(mdout_path, "nb14_LJ"), -0.5, abs_tol=0.01
    )
    assert math.isclose(
        _extract_mdout_term(mdout_path, "improper_dihedral"),
        0.5 * math.pi**2,
        abs_tol=0.01,
    )


def test_chamber_rejects_half_present_lj14_matrix(tmp_path):
    result, _mdout_path = _run_sponge(tmp_path, include_lj14_b=False)
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert (
        "LENNARD_JONES_14_ACOEF and LENNARD_JONES_14_BCOEF must either both "
        "be present or both be absent"
    ) in output


def test_amber_without_lj14_flags_falls_back_to_ordinary_matrix(tmp_path):
    result, mdout_path = _run_sponge(
        tmp_path,
        include_lj14_a=False,
        include_lj14_b=False,
        normal_a=729.0,
        normal_b=54.0,
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "nb14_LJ"), -0.5, abs_tol=0.01
    )
