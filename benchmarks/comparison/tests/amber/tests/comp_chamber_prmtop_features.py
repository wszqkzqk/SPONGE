import math
import os
import re
import struct
import subprocess
import textwrap
from pathlib import Path

import pytest


def _flag(name, fmt, values, items_per_line=10):
    lines = [f"%FLAG {name}", f"%FORMAT({fmt})"]
    match = re.fullmatch(r"(?:\d+)?([IiEeDdFfGgAa])(\d+)(?:\.(\d+))?", fmt)
    assert match is not None, fmt
    kind, width_text, precision_text = match.groups()
    width = int(width_text)
    precision = int(precision_text or 0)

    def render(value):
        if isinstance(value, str):
            if kind.upper() == "A":
                return f"{value:<{width}}"
            return f"{value:>{width}}"
        if kind.upper() == "I":
            return f"{int(value):{width}d}"
        if kind.upper() == "A":
            return f"{str(value):<{width}}"
        return f"{float(value):{width}.{precision}E}"

    rendered = [render(value) for value in values]
    for start in range(0, len(rendered), items_per_line):
        lines.append("".join(rendered[start : start + items_per_line]))
    return "\n".join(lines) + "\n"


def _replace_flag(path, name, replacement):
    lines = path.read_text().splitlines(keepends=True)
    start = next(
        i for i, line in enumerate(lines) if line.strip() == f"%FLAG {name}"
    )
    end = start + 1
    while end < len(lines) and not lines[end].startswith("%FLAG"):
        end += 1
    lines[start:end] = [] if replacement is None else [replacement]
    path.write_text("".join(lines))


def _duplicate_flag(path, name):
    lines = path.read_text().splitlines(keepends=True)
    start = next(
        i for i, line in enumerate(lines) if line.strip() == f"%FLAG {name}"
    )
    end = start + 1
    while end < len(lines) and not lines[end].startswith("%FLAG"):
        end += 1
    path.write_text("".join(lines + lines[start:end]))


def _move_flag_to_end(path, name):
    lines = path.read_text().splitlines(keepends=True)
    start = next(
        i for i, line in enumerate(lines) if line.strip() == f"%FLAG {name}"
    )
    end = start + 1
    while end < len(lines) and not lines[end].startswith("%FLAG"):
        end += 1
    section = lines[start:end]
    path.write_text("".join(lines[:start] + lines[end:] + section))


def _write_chamber_prmtop(
    path,
    *,
    include_lj14_a=True,
    include_lj14_b=True,
    normal_a=729.0,
    normal_b=27.0,
    charmm_improper_count=1,
    charmm_improper_rows=((1, 2, 3, 4, 1),),
    charmm_improper_type_count=1,
    charmm_improper_force_constants=(2.0,),
    charmm_improper_phases=(0.0,),
):
    pointers = [0] * 31
    pointers[0] = 4  # NATOM
    pointers[1] = 1  # NTYPES
    pointers[7] = 1  # MPHIA: one dihedral without hydrogen
    pointers[10] = 6  # NNB
    pointers[11] = 1  # NRES
    pointers[14] = 1  # NPHIA
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
    ]
    if charmm_improper_count is not None:
        sections.append(
            _flag("CHARMM_NUM_IMPROPERS", "10I8", [charmm_improper_count])
        )
    if charmm_improper_rows is not None:
        sections.append(
            _flag(
                "CHARMM_IMPROPERS",
                "10I8",
                [value for row in charmm_improper_rows for value in row],
            )
        )
    if charmm_improper_type_count is not None:
        sections.append(
            _flag("CHARMM_NUM_IMPR_TYPES", "1I8", [charmm_improper_type_count])
        )
    if charmm_improper_force_constants is not None:
        sections.append(
            _flag(
                "CHARMM_IMPROPER_FORCE_CONSTANT",
                "5E16.8",
                charmm_improper_force_constants,
                5,
            )
        )
    if charmm_improper_phases is not None:
        sections.append(
            _flag("CHARMM_IMPROPER_PHASE", "5E16.8", charmm_improper_phases, 5)
        )
    sections.extend(
        [
            # Ordinary and special matrices deliberately differ. At r=sqrt(3)
            # A, ordinary LJ is zero, while LJ14 is -1 kcal/mol before SCNB=2.
            _flag("LENNARD_JONES_ACOEF", "3E24.16", [normal_a], 3),
            _flag("LENNARD_JONES_BCOEF", "3E24.16", [normal_b], 3),
        ]
    )
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
    _write_feature_rst7(
        path,
        [
            (0.0, 1.0, 0.0),
            (0.0, 0.0, 0.0),
            (1.0, 0.0, 0.0),
            (1.0, 0.0, 1.0),
        ],
    )


def _write_feature_rst7(
    path,
    coordinates,
    *,
    time=None,
    velocities=None,
    box=(50.0, 50.0, 50.0, 90.0, 90.0, 90.0),
    trailing=(),
):
    header = f"{len(coordinates):6d}"
    if time is not None:
        header += f"{time:15.7E}" if not isinstance(time, str) else f" {time}"
    lines = ["synthetic AMBER reader feature test", header]

    def append_block(values):
        for start in range(0, len(values), 6):
            lines.append(
                "".join(
                    f"{value:>12}"
                    if isinstance(value, str)
                    else f"{value:12.7f}"
                    for value in values[start : start + 6]
                )
            )

    append_block([value for coordinate in coordinates for value in coordinate])
    if velocities is not None:
        append_block([value for velocity in velocities for value in velocity])
    if box is not None:
        append_block(list(box))
    if trailing:
        append_block(list(trailing))
    path.write_text("\n".join(lines) + "\n")


def _write_indexed_feature_prmtop(
    path,
    *,
    atom_types,
    nonbonded_parm_index,
    normal_a,
    normal_b,
    lj14_a=None,
    lj14_b=None,
    hbond_a=None,
    hbond_b=None,
    dihedral_rows=(),
    excluded_counts=None,
    excluded_list=(),
    charges=None,
    include_pointers=True,
    include_atom_type_index=True,
    dihedral_force_constant=(1.0,),
    dihedral_periodicity=(1.0,),
    dihedral_phase=(0.0,),
    scee_scale_factor=(1.0,),
    scnb_scale_factor=(2.0,),
    pointer_overrides=None,
    masses=None,
    atomic_numbers=None,
    amber_atom_types=None,
    extra_sections=(),
    pointer_length=31,
):
    atom_types = list(atom_types)
    dihedral_rows = [list(row) for row in dihedral_rows]
    atom_count = len(atom_types)
    atom_type_count = max(atom_types)
    if excluded_counts is None:
        excluded_counts = [0] * atom_count
    if charges is None:
        charges = [0.0] * atom_count
    if masses is None:
        masses = [12.0] * atom_count

    hbond_count = max(
        len(hbond_a) if hbond_a is not None else 0,
        len(hbond_b) if hbond_b is not None else 0,
    )
    pointers = [0] * max(31, pointer_length)
    pointers[0] = atom_count  # NATOM
    pointers[1] = atom_type_count  # NTYPES
    pointers[7] = len(dihedral_rows)  # MPHIA
    pointers[10] = len(excluded_list)  # NNB
    pointers[11] = 1  # NRES
    pointers[14] = len(dihedral_rows)  # NPHIA
    pointers[17] = 1 if dihedral_rows else 0  # NPTRA
    pointers[19] = hbond_count  # NPHB
    pointers[27] = 1  # IFBOX
    if pointer_overrides is not None:
        for index, value in pointer_overrides.items():
            pointers[index] = value
    pointers = pointers[:pointer_length]

    sections = ["%VERSION VERSION_STAMP = V0001.000 DATE = 01/01/70 00:00:00\n"]
    if include_pointers:
        sections.append(_flag("POINTERS", "10I8", pointers))
    sections.extend(
        [
            _flag("CHARGE", "5E16.8", charges, 5),
            _flag("MASS", "5E16.8", masses, 5),
        ]
    )
    if atomic_numbers is not None:
        sections.append(_flag("ATOMIC_NUMBER", "10I8", atomic_numbers))
    if amber_atom_types is not None:
        sections.append(_flag("AMBER_ATOM_TYPE", "20A4", amber_atom_types, 20))
    if include_atom_type_index:
        sections.append(_flag("ATOM_TYPE_INDEX", "10I8", atom_types))
    if nonbonded_parm_index is not None:
        sections.append(
            _flag("NONBONDED_PARM_INDEX", "10I8", nonbonded_parm_index)
        )
    sections.extend(
        [
            _flag("NUMBER_EXCLUDED_ATOMS", "10I8", excluded_counts),
            _flag("RESIDUE_POINTER", "10I8", [1]),
        ]
    )
    if dihedral_rows:
        if dihedral_force_constant is not None:
            sections.append(
                _flag(
                    "DIHEDRAL_FORCE_CONSTANT",
                    "5E16.8",
                    dihedral_force_constant,
                    5,
                )
            )
        if dihedral_periodicity is not None:
            sections.append(
                _flag("DIHEDRAL_PERIODICITY", "5E16.8", dihedral_periodicity, 5)
            )
        if dihedral_phase is not None:
            sections.append(
                _flag("DIHEDRAL_PHASE", "5E16.8", dihedral_phase, 5)
            )
        if scee_scale_factor is not None:
            sections.append(
                _flag("SCEE_SCALE_FACTOR", "5E16.8", scee_scale_factor, 5)
            )
        if scnb_scale_factor is not None:
            sections.append(
                _flag("SCNB_SCALE_FACTOR", "5E16.8", scnb_scale_factor, 5)
            )
    if normal_a is not None:
        sections.append(_flag("LENNARD_JONES_ACOEF", "3E24.16", normal_a, 3))
    if normal_b is not None:
        sections.append(_flag("LENNARD_JONES_BCOEF", "3E24.16", normal_b, 3))
    if lj14_a is not None:
        sections.append(_flag("LENNARD_JONES_14_ACOEF", "3E24.16", lj14_a, 3))
    if lj14_b is not None:
        sections.append(_flag("LENNARD_JONES_14_BCOEF", "3E24.16", lj14_b, 3))
    if hbond_a is not None:
        sections.append(_flag("HBOND_ACOEF", "5E16.8", hbond_a, 5))
    if hbond_b is not None:
        sections.append(_flag("HBOND_BCOEF", "5E16.8", hbond_b, 5))
    if dihedral_rows:
        sections.append(
            _flag(
                "DIHEDRALS_WITHOUT_HYDROGEN",
                "10I8",
                [value for row in dihedral_rows for value in row],
            )
        )
    sections.extend(extra_sections)
    sections.append(_flag("EXCLUDED_ATOMS_LIST", "10I8", excluded_list))
    path.write_text("".join(sections))


def _write_mdin(
    path,
    prmtop_path,
    rst7_path,
    mdout_path,
    *,
    cutoff=8.0,
    step_limit=0,
    extra_commands="",
):
    path.write_text(
        textwrap.dedent(
            f"""
            md_name = "synthetic_chamber_features"
            mode = "nve"
            step_limit = {step_limit}
            dt = 0
            cutoff = {cutoff}
            amber_parm7 = {str(prmtop_path)!r}
            amber_rst7 = {str(rst7_path)!r}
            mdout = {str(mdout_path)!r}
            print_zeroth_frame = 1
            write_mdout_interval = 1
            {extra_commands}
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
    **chamber_options,
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
        **chamber_options,
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


def _run_indexed_feature_sponge(
    tmp_path,
    *,
    coordinates,
    rst7_options=None,
    cutoff=8.0,
    step_limit=0,
    extra_commands="",
    prmtop_mutator=None,
    rst7_mutator=None,
    **prmtop_options,
):
    prmtop_path = tmp_path / "indexed_feature.prmtop"
    rst7_path = tmp_path / "indexed_feature.rst7"
    mdin_path = tmp_path / "indexed_feature.toml"
    mdout_path = tmp_path / "indexed_feature_mdout.txt"
    _write_indexed_feature_prmtop(prmtop_path, **prmtop_options)
    if prmtop_mutator is not None:
        prmtop_mutator(prmtop_path)
    _write_feature_rst7(rst7_path, coordinates, **(rst7_options or {}))
    if rst7_mutator is not None:
        rst7_mutator(rst7_path)
    _write_mdin(
        mdin_path,
        prmtop_path,
        rst7_path,
        mdout_path,
        cutoff=cutoff,
        step_limit=step_limit,
        extra_commands=extra_commands,
    )
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


def _assert_unsupported_location(output, prmtop_path, flag_name):
    assert f"Input file: {prmtop_path}" in output
    assert re.search(
        rf"%FLAG {re.escape(flag_name)} begins on line [1-9][0-9]*", output
    )


_PERMUTED_TWO_TYPE_NPI = [3, 1, 1, 2]
_TWO_ATOM_COORDINATES = [(0.0, 0.0, 0.0), (math.sqrt(3.0), 0.0, 0.0)]
_FOUR_ATOM_COORDINATES = [
    (0.0, 1.0, 0.0),
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 0.0, 1.0),
]
_GENERAL_DIHEDRAL_COORDINATES = [
    (1.0, 2.2, 0.3),
    (2.1, 0.8, 1.9),
    (4.3, 3.1, 1.1),
    (5.8, 1.7, 4.2),
]


def _subtract(first, second):
    return tuple(a - b for a, b in zip(first, second))


def _cross(first, second):
    return (
        first[1] * second[2] - first[2] * second[1],
        first[2] * second[0] - first[0] * second[2],
        first[0] * second[1] - first[1] * second[0],
    )


def _dot(first, second):
    return sum(a * b for a, b in zip(first, second))


def _amber_torsion_phi(coordinates):
    drij = _subtract(coordinates[0], coordinates[1])
    drkj = _subtract(coordinates[2], coordinates[1])
    drkl = _subtract(coordinates[2], coordinates[3])
    normal_1 = _cross(drij, drkj)
    normal_2 = _cross(drkl, drkj)
    normal_product = math.sqrt(
        _dot(normal_1, normal_1) * _dot(normal_2, normal_2)
    )
    central_length = math.sqrt(_dot(drkj, drkj))
    cosine = -_dot(normal_1, normal_2) / normal_product
    sine = _dot(_cross(normal_2, normal_1), drkj) / (
        normal_product * central_length
    )
    return math.atan2(sine, cosine)


def _amber_classic_dihedral_energy(coordinates, *, pk, periodicity, phase):
    # Oracle: Amber Classic set.F90::dihpar followed by ene.F90::ephi.
    if abs(phase - math.pi) <= 1.0e-3:
        phase = math.copysign(math.pi, phase)
    phase_cosine = math.cos(phase)
    phase_sine = math.sin(phase)
    if abs(phase_cosine) <= 1.0e-6:
        phase_cosine = 0.0
    if abs(phase_sine) <= 1.0e-6:
        phase_sine = 0.0
    nphi = abs(periodicity) * _amber_torsion_phi(coordinates)
    return pk * (
        1.0 + math.cos(nphi) * phase_cosine + math.sin(nphi) * phase_sine
    )


def _read_last_forces(path, atom_count=4):
    raw = Path(path).read_bytes()
    values = struct.unpack(f"={len(raw) // 4}f", raw)
    return values[-3 * atom_count :]


def _two_atom_indexed_options(**overrides):
    options = {
        "atom_types": [1, 2],
        "nonbonded_parm_index": _PERMUTED_TWO_TYPE_NPI,
        # The cross pair is source entry 1. The two self terms deliberately
        # differ and are stored in non-canonical source slots 2 and 3.
        "normal_a": [729.0, 200.0, 300.0],
        "normal_b": [54.0, 20.0, 30.0],
    }
    options.update(overrides)
    return options


def _four_atom_indexed_options(**overrides):
    options = {
        "atom_types": [1, 1, 2, 2],
        "nonbonded_parm_index": _PERMUTED_TWO_TYPE_NPI,
        "normal_a": [729.0, 200.0, 300.0],
        "normal_b": [54.0, 20.0, 30.0],
        "lj14_a": [729.0, 400.0, 500.0],
        "lj14_b": [54.0, 40.0, 50.0],
        "dihedral_rows": [(0, 3, 6, 9, 1)],
        "excluded_counts": [3, 2, 1, 0],
        "excluded_list": [2, 3, 4, 3, 4, 4],
    }
    options.update(overrides)
    return options


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


def test_chamber_quadratic_improper_preserves_order_phase_and_type(tmp_path):
    result, mdout_path = _run_sponge(
        tmp_path,
        charmm_improper_count=1,
        charmm_improper_rows=((4, 2, 3, 1, 2),),
        charmm_improper_type_count=2,
        charmm_improper_force_constants=(11.0, 3.0),
        charmm_improper_phases=(0.75, 0.25),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    # This asymmetric atom order gives phi=-pi/2. Selecting type 2 must use
    # K=3 and phi0=0.25 rather than either entry from type 1.
    expected = 3.0 * (-0.5 * math.pi - 0.25) ** 2
    assert math.isclose(
        _extract_mdout_term(mdout_path, "improper_dihedral"),
        expected,
        abs_tol=0.02,
    )


@pytest.mark.parametrize(
    "overrides",
    [
        {"charmm_improper_count": None},
        {"charmm_improper_rows": None},
        {"charmm_improper_type_count": None},
        {"charmm_improper_force_constants": None},
        {"charmm_improper_phases": None},
    ],
    ids=["count", "tuples", "type-count", "force-constants", "phases"],
)
def test_chamber_rejects_incomplete_quadratic_improper_flags(
    tmp_path, overrides
):
    result, _mdout_path = _run_sponge(tmp_path, **overrides)
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert "CHAMBER improper sections are incomplete" in output


def test_chamber_rejects_bad_quadratic_improper_tuple_length(tmp_path):
    result, _mdout_path = _run_sponge(
        tmp_path,
        charmm_improper_rows=((1, 2, 3, 4),),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert "CHARMM_IMPROPERS has 4 values; expected 5" in output


@pytest.mark.parametrize(
    "improper_row",
    [
        (0, 2, 3, 4, 1),
        (1, 2, 3, 4, 2),
    ],
    ids=["atom", "type"],
)
def test_chamber_rejects_quadratic_improper_index_out_of_range(
    tmp_path, improper_row
):
    result, _mdout_path = _run_sponge(
        tmp_path,
        charmm_improper_rows=(improper_row,),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert "CHARMM_IMPROPERS contains an out-of-range" in output


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


def test_nonbonded_parm_index_remaps_ordinary_lj_matrix(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "LJ_short"), -1.0, abs_tol=0.01
    )


def test_nonbonded_parm_index_remaps_chamber_lj14_matrix(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            normal_b=[27.0, 20.0, 30.0],
        ),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert "non-bond 14 numbers is 1" in output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "nb14_LJ"), -0.5, abs_tol=0.01
    )


def test_nonbonded_parm_index_remaps_lj14_fallback_matrix(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(lj14_a=None, lj14_b=None),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "nb14_LJ"), -0.5, abs_tol=0.01
    )


def test_lj14_fallback_conversion_error_reports_actual_source_flags(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            normal_a=["1.40129846E-45"] * 3,
            normal_b=[0.0] * 3,
            lj14_a=None,
            lj14_b=None,
            scnb_scale_factor=["3.0E38"],
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "AMBER section LENNARD_JONES_ACOEF" in output
    assert "LENNARD_JONES_14_ACOEF" not in output
    prmtop_path = tmp_path / "indexed_feature.prmtop"
    _assert_unsupported_location(output, prmtop_path, "LENNARD_JONES_ACOEF")


def test_missing_amber_14_scale_flags_use_spec_defaults(tmp_path):
    energies = {}
    for case_name, scale_options in (
        (
            "explicit",
            {"scee_scale_factor": [1.2], "scnb_scale_factor": [2.0]},
        ),
        (
            "omitted",
            {"scee_scale_factor": None, "scnb_scale_factor": None},
        ),
    ):
        case_path = tmp_path / case_name
        case_path.mkdir()
        result, mdout_path = _run_indexed_feature_sponge(
            case_path,
            coordinates=_FOUR_ATOM_COORDINATES,
            **_four_atom_indexed_options(
                charges=[18.2223, 0.0, 0.0, 18.2223],
                **scale_options,
            ),
        )
        output = result.stdout + "\n" + result.stderr
        assert result.returncode == 0, output
        energies[case_name] = {
            term: _extract_mdout_term(mdout_path, term)
            for term in ("nb14_LJ", "nb14_EE")
        }

    assert math.isclose(energies["omitted"]["nb14_LJ"], -0.5, abs_tol=0.01)
    assert abs(energies["explicit"]["nb14_EE"]) > 1.0
    for term in ("nb14_LJ", "nb14_EE"):
        assert math.isclose(
            energies["omitted"][term], energies["explicit"][term], abs_tol=0.01
        )


@pytest.mark.parametrize(
    ("overrides", "expected_error"),
    [
        ({"include_pointers": False}, "AMBER POINTERS section is required"),
        (
            {"include_atom_type_index": False},
            "ATOM_TYPE_INDEX is required for an AMBER topology with atoms",
        ),
        (
            {"normal_a": None},
            "LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF are required",
        ),
        (
            {"normal_b": None},
            "LENNARD_JONES_ACOEF and LENNARD_JONES_BCOEF are required",
        ),
    ],
    ids=["pointers", "atom-types", "lj-a", "lj-b"],
)
def test_required_atom_topology_sections_fail_fast(
    tmp_path, overrides, expected_error
):
    result, _mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(**overrides),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert expected_error in output


@pytest.mark.parametrize(
    ("pointer_overrides", "expected_error"),
    [
        (
            {7: 2},
            "DIHEDRALS_WITHOUT_HYDROGEN has 1 tuples; POINTERS declares 2",
        ),
        (
            {6: 1},
            "DIHEDRALS_INC_HYDROGEN is missing; POINTERS declares 1 tuples",
        ),
    ],
    ids=["wrong-count", "missing-section"],
)
def test_dihedral_sections_must_match_pointer_counts(
    tmp_path, pointer_overrides, expected_error
):
    result, _mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(pointer_overrides=pointer_overrides),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert expected_error in output


def test_declared_dihedral_types_require_parameter_sections(tmp_path):
    result, _mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(dihedral_force_constant=None),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert "DIHEDRAL_FORCE_CONSTANT, DIHEDRAL_PERIODICITY" in output


@pytest.mark.parametrize(
    ("nonbonded_parm_index", "expected_error"),
    [
        (None, "NONBONDED_PARM_INDEX is required"),
        ([3, 1, 1], "NONBONDED_PARM_INDEX has 3 values; expected 4"),
        ([3, 1, 2, 2], "NONBONDED_PARM_INDEX is asymmetric"),
        ([3, 0, 0, 2], "NONBONDED_PARM_INDEX contains zero"),
        ([4, 1, 1, 2], "NONBONDED_PARM_INDEX value 4"),
    ],
    ids=["missing", "wrong-size", "asymmetric", "zero", "out-of-range"],
)
def test_invalid_nonbonded_parm_index_fails_fast(
    tmp_path, nonbonded_parm_index, expected_error
):
    result, _mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            nonbonded_parm_index=nonbonded_parm_index,
        ),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert expected_error in output


def test_zero_hbond_placeholder_maps_to_zero_lj(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            nonbonded_parm_index=[3, -1, -1, 2],
            hbond_a=[0.0],
            hbond_b=[0.0],
        ),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "LJ_short"), 0.0, abs_tol=0.01
    )


@pytest.mark.parametrize(
    ("overrides", "expected_error"),
    [
        (
            {
                "nonbonded_parm_index": [3, -1, -1, 2],
                "hbond_a": [1.0],
                "hbond_b": [0.0],
            },
            "AMBER 10-12 nonbonded terms are not supported",
        ),
        (
            {
                "nonbonded_parm_index": [3, -2, -2, 2],
                "hbond_a": [0.0],
                "hbond_b": [0.0],
            },
            "does not select a valid HBOND parameter",
        ),
        (
            {
                "nonbonded_parm_index": [3, -1, -1, 2],
                "hbond_a": [0.0],
                "hbond_b": None,
            },
            "HBOND_ACOEF and HBOND_BCOEF must both be present",
        ),
    ],
    ids=["nonzero-10-12", "bad-negative-index", "half-present"],
)
def test_unsupported_or_invalid_hbond_parameters_fail_fast(
    tmp_path, overrides, expected_error
):
    result, _mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(**overrides),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert expected_error in output
    if "not supported" in expected_error:
        _assert_unsupported_location(
            output, tmp_path / "indexed_feature.prmtop", "NONBONDED_PARM_INDEX"
        )


@pytest.mark.parametrize(
    ("raw_c", "raw_d", "has_nb14"),
    [
        (6, 9, True),
        (-6, 9, False),
        (6, -9, False),
        (-6, -9, False),
    ],
    ids=["proper", "ignore-end", "periodic-improper", "both-signs"],
)
def test_signed_amber_dihedrals_preserve_periodic_term_and_control_nb14(
    tmp_path, raw_c, raw_d, has_nb14
):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=[(0, 3, raw_c, raw_d, 1)],
        ),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "dihedral"), 1.0, abs_tol=0.01
    )
    if has_nb14:
        assert "non-bond 14 numbers is 1" in output
    else:
        assert "NB14 IS NOT INITIALIZED" in output


def test_negative_periodic_dihedral_force_constant_is_accepted(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(dihedral_force_constant=[-1.0]),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "dihedral"), -1.0, abs_tol=0.01
    )


@pytest.mark.parametrize("periodicity", [0.0, 0.5, 1.25])
def test_amber_real_periodicity_matches_classic_energy(tmp_path, periodicity):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(dihedral_periodicity=[periodicity]),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    expected = _amber_classic_dihedral_energy(
        _FOUR_ATOM_COORDINATES,
        pk=1.0,
        periodicity=periodicity,
        phase=0.0,
    )
    assert _extract_mdout_term(mdout_path, "dihedral") == pytest.approx(
        expected, abs=0.02
    )


def test_amber_negative_periodicity_fails_only_when_referenced(tmp_path):
    referenced_path = tmp_path / "referenced"
    referenced_path.mkdir()
    result, _ = _run_indexed_feature_sponge(
        referenced_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(dihedral_periodicity=[-0.5]),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "negative periodicity but is referenced" in output
    prmtop_path = referenced_path / "indexed_feature.prmtop"
    _assert_unsupported_location(output, prmtop_path, "DIHEDRAL_PERIODICITY")
    _assert_unsupported_location(
        output, prmtop_path, "DIHEDRALS_WITHOUT_HYDROGEN"
    )

    unused_path = tmp_path / "unused"
    unused_path.mkdir()
    result, _ = _run_indexed_feature_sponge(
        unused_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_force_constant=[1.0, 1.0],
            dihedral_periodicity=[1.0, -0.5],
            dihedral_phase=[0.0, 0.0],
            scee_scale_factor=[1.0, 1.0],
            scnb_scale_factor=[2.0, 2.0],
            pointer_overrides={17: 2},
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


@pytest.mark.parametrize("phase_distance", [1.0e-3 - 5.0e-8, 1.0e-3 + 5.0e-8])
def test_amber_phase_snap_uses_original_double_threshold(
    tmp_path, phase_distance
):
    phase = math.pi + phase_distance

    def replace_phase(path):
        _replace_flag(
            path,
            "DIHEDRAL_PHASE",
            _flag("DIHEDRAL_PHASE", "3E24.16", [phase], 3),
        )

    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        prmtop_mutator=replace_phase,
        **_four_atom_indexed_options(
            dihedral_force_constant=[1.0e5],
            lj14_a=[0.0, 0.0, 0.0],
            lj14_b=[0.0, 0.0, 0.0],
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    expected = _amber_classic_dihedral_energy(
        _FOUR_ATOM_COORDINATES,
        pk=1.0e5,
        periodicity=1.0,
        phase=phase,
    )
    assert _extract_mdout_term(mdout_path, "dihedral") == pytest.approx(
        expected, abs=0.06
    )


def test_amber_trig_threshold_precedes_force_constant_scaling(tmp_path):
    pk = 1.0e7
    phases = [math.asin(1.0e-6 - 5.0e-14), math.asin(1.0e-6 + 5.0e-14)]
    observed = []
    expected = []
    for index, phase in enumerate(phases):
        case_path = tmp_path / str(index)
        case_path.mkdir()

        def replace_phase(path, value=phase):
            _replace_flag(
                path,
                "DIHEDRAL_PHASE",
                _flag("DIHEDRAL_PHASE", "3E24.16", [value], 3),
            )

        result, mdout_path = _run_indexed_feature_sponge(
            case_path,
            coordinates=_FOUR_ATOM_COORDINATES,
            prmtop_mutator=replace_phase,
            **_four_atom_indexed_options(
                dihedral_force_constant=[pk],
                lj14_a=[0.0, 0.0, 0.0],
                lj14_b=[0.0, 0.0, 0.0],
            ),
        )
        output = result.stdout + "\n" + result.stderr
        assert result.returncode == 0, output
        observed.append(_extract_mdout_term(mdout_path, "dihedral"))
        expected.append(
            _amber_classic_dihedral_energy(
                _FOUR_ATOM_COORDINATES,
                pk=pk,
                periodicity=1.0,
                phase=phase,
            )
        )
    assert observed[0] - observed[1] == pytest.approx(
        expected[0] - expected[1], abs=2.0
    )


def test_amber_fractional_periodicity_and_small_pk_preserve_force(tmp_path):
    pk = 5.0e-7
    periodicity = 0.5
    phase = 0.37
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_GENERAL_DIHEDRAL_COORDINATES,
        extra_commands="frc = 'frc.dat'",
        **_four_atom_indexed_options(
            dihedral_force_constant=[pk],
            dihedral_periodicity=[periodicity],
            dihedral_phase=[phase],
            lj14_a=[0.0, 0.0, 0.0],
            lj14_b=[0.0, 0.0, 0.0],
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    forces = _read_last_forces(tmp_path / "frc.dat")
    epsilon = 1.0e-5
    expected_forces = []
    for atom in range(4):
        for axis in range(3):
            plus = [
                list(coordinate) for coordinate in _GENERAL_DIHEDRAL_COORDINATES
            ]
            minus = [
                list(coordinate) for coordinate in _GENERAL_DIHEDRAL_COORDINATES
            ]
            plus[atom][axis] += epsilon
            minus[atom][axis] -= epsilon
            energy_plus = _amber_classic_dihedral_energy(
                plus, pk=pk, periodicity=periodicity, phase=phase
            )
            energy_minus = _amber_classic_dihedral_energy(
                minus, pk=pk, periodicity=periodicity, phase=phase
            )
            expected_forces.append(
                -(energy_plus - energy_minus) / (2 * epsilon)
            )
    assert max(abs(value) for value in forces) > 1.0e-8
    assert forces == pytest.approx(expected_forces, rel=8.0e-3, abs=2.0e-9)


@pytest.mark.parametrize(
    ("overrides", "expected_error", "parameter_flags"),
    [
        (
            {"dihedral_force_constant": ["1e-46"]},
            "AMBER section DIHEDRAL_FORCE_CONSTANT contains a nonzero value",
            ("DIHEDRAL_FORCE_CONSTANT",),
        ),
        (
            {"dihedral_periodicity": ["1e-46"]},
            "AMBER section DIHEDRAL_PERIODICITY contains a nonzero value",
            ("DIHEDRAL_PERIODICITY",),
        ),
        (
            {
                "dihedral_force_constant": ["1.17549435e-38"],
                "dihedral_phase": [0.37],
            },
            "AMBER section DIHEDRAL_FORCE_CONSTANT and DIHEDRAL_PHASE "
            "contains a subnormal value",
            ("DIHEDRAL_FORCE_CONSTANT", "DIHEDRAL_PHASE"),
        ),
    ],
    ids=["force_constant", "periodicity", "derived_phase_coefficient"],
)
def test_amber_dihedral_narrowing_error_reports_parameter_and_term_sources(
    tmp_path, overrides, expected_error, parameter_flags
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(**overrides),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert expected_error in output
    assert (
        "AMBER dihedral parameter type 1 is referenced by "
        "DIHEDRALS_WITHOUT_HYDROGEN term 1"
    ) in output
    prmtop_path = tmp_path / "indexed_feature.prmtop"
    for flag in (*parameter_flags, "DIHEDRALS_WITHOUT_HYDROGEN"):
        _assert_unsupported_location(output, prmtop_path, flag)


def test_zero_third_dihedral_index_still_generates_nb14(tmp_path):
    result, _mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=[(3, 6, 0, 9, 1)],
        ),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert "non-bond 14 numbers is 1" in output


def test_multiterm_dihedral_adds_only_one_nb14_pair(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=[(0, 3, 6, 9, 1), (0, 3, -6, 9, 1)],
        ),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode == 0, output
    assert "dihedral_numbers is 2" in output
    assert "non-bond 14 numbers is 1" in output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "dihedral"), 2.0, abs_tol=0.01
    )


@pytest.mark.parametrize(
    ("dihedral_row", "expected_error"),
    [
        ((-3, 3, 6, 9, 1), "field A cannot be negative"),
        ((0, 3, 5, 9, 1), "field C value 5 is not a multiple of 3"),
        ((0, 3, 6, 12, 1), "field D selects atom 4 outside [0, 4)"),
        ((0, 3, 6, 9, 2), "out-of-range parameter index"),
        ((0, 3, 6, 9), "expected a multiple of 5"),
    ],
    ids=[
        "negative-a",
        "misaligned-c",
        "atom-range",
        "type-range",
        "tuple-size",
    ],
)
def test_invalid_amber_dihedral_encoding_fails_fast(
    tmp_path, dihedral_row, expected_error
):
    result, _mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(dihedral_rows=[dihedral_row]),
    )
    output = result.stdout + "\n" + result.stderr

    assert result.returncode != 0
    assert expected_error in output
    if expected_error == "out-of-range parameter index":
        assert "AMBER DIHEDRALS_WITHOUT_HYDROGEN term 1" in output
        _assert_unsupported_location(
            output,
            tmp_path / "indexed_feature.prmtop",
            "DIHEDRALS_WITHOUT_HYDROGEN",
        )


@pytest.mark.parametrize(
    "dihedral_row",
    [
        (0, 0, 6, 9, 1),
        (0, 3, 0, 9, 1),
        (0, 3, 6, 0, 1),
        (0, 3, 3, 9, 1),
        (0, 3, 6, 3, 1),
        (0, 3, 6, 6, 1),
    ],
    ids=["a-b", "a-c", "a-d", "b-c", "b-d", "c-d"],
)
def test_standard_dihedral_repeated_atoms_fail_with_tuple_source(
    tmp_path, dihedral_row
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(dihedral_rows=[dihedral_row]),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "DIHEDRALS_WITHOUT_HYDROGEN term 1 contains repeated atoms" in output
    _assert_unsupported_location(
        output,
        tmp_path / "indexed_feature.prmtop",
        "DIHEDRALS_WITHOUT_HYDROGEN",
    )


def test_chamber_improper_repeated_atoms_fail_with_tuple_source(tmp_path):
    result, _ = _run_sponge(tmp_path, charmm_improper_rows=((1, 2, 3, 3, 1),))
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "CHARMM_IMPROPERS term 1 contains repeated atoms" in output
    _assert_unsupported_location(
        output, tmp_path / "synthetic.prmtop", "CHARMM_IMPROPERS"
    )


@pytest.mark.parametrize("interaction", ["bond", "angle"])
def test_without_hydrogen_sections_use_front_pointer_group(
    tmp_path, interaction
):
    if interaction == "bond":
        pointer_overrides = {3: 1, 12: 7, 15: 1}
        extra_sections = [
            _flag("BOND_FORCE_CONSTANT", "5E16.8", [10.0], 5),
            _flag("BOND_EQUIL_VALUE", "5E16.8", [1.0], 5),
            _flag("BONDS_WITHOUT_HYDROGEN", "10I8", [0, 3, 1]),
        ]
        expected = "bond_numbers is 1"
    else:
        pointer_overrides = {5: 1, 13: 9, 16: 1}
        extra_sections = [
            _flag("ANGLE_FORCE_CONSTANT", "5E16.8", [10.0], 5),
            _flag("ANGLE_EQUIL_VALUE", "5E16.8", [math.pi / 2], 5),
            _flag("ANGLES_WITHOUT_HYDROGEN", "10I8", [0, 3, 6, 1]),
        ]
        expected = "angle_numbers is 1"
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=(),
            lj14_a=None,
            lj14_b=None,
            pointer_overrides=pointer_overrides,
            extra_sections=extra_sections,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert expected in output


def test_nphia_can_differ_from_mphia_without_changing_tuple_count(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(pointer_overrides={14: 37}),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "dihedral_numbers is 1" in output


def test_fixed_width_parser_accepts_fortran_d_exponents_and_section_reordering(
    tmp_path,
):
    def mutate(path):
        _replace_flag(
            path,
            "CHARGE",
            _flag("CHARGE", "5E16.8", ["0.0D+00", "0.0d+00"], 5),
        )
        _move_flag_to_end(path, "POINTERS")

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


def test_unsupported_format_descriptor_reports_file_and_flag_line(tmp_path):
    def mutate(path):
        charge = _flag("CHARGE", "5E16.8", [0.0, 0.0], 5)
        _replace_flag(path, "CHARGE", charge.replace("5E16.8", "5Q16.8"))

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "unsupported AMBER %FORMAT descriptor 'Q'" in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "CHARGE"
    )


@pytest.mark.parametrize(
    ("format_descriptor", "expected_error"),
    [
        (
            "999999999999999999999E16.8",
            "%FORMAT repeat exceeds the platform size range",
        ),
        (
            "5E999999999999999999999.8",
            "%FORMAT width exceeds the platform size range",
        ),
    ],
    ids=["repeat", "width"],
)
def test_format_decimal_accumulation_rejects_size_overflow_before_multiply(
    tmp_path, format_descriptor, expected_error
):
    def mutate(path):
        charge = _flag("CHARGE", "5E16.8", [0.0, 0.0], 5)
        _replace_flag(
            path,
            "CHARGE",
            charge.replace("%FORMAT(5E16.8)", f"%FORMAT({format_descriptor})"),
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "CHARGE"
    )


def test_grouped_format_large_repeats_are_consumed_without_expansion(tmp_path):
    def mutate(path):
        path.write_text(
            path.read_text()
            + "%FLAG FORMAT_STRESS\n"
            + "%FORMAT(60000(I1),60000(I1))\n"
            + "1\n"
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


def test_format_width_above_signed_int_range_does_not_allocate_padding(
    tmp_path,
):
    def mutate(path):
        path.write_text(
            path.read_text()
            + "%FLAG FORMAT_WIDE_FIELD\n"
            + "%FORMAT(A2147483648)\n"
            + "x\n"
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


def test_parm7_rejects_nonempty_content_outside_flag_sections(tmp_path):
    def mutate(path):
        path.write_text("this is not an AMBER directive\n" + path.read_text())

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert (
        "unexpected content outside an AMBER %FLAG section on line 1" in output
    )
    assert f"Input file: {tmp_path / 'indexed_feature.prmtop'}" in output


def test_zero_placeholder_radii_and_screen_sections_are_accepted(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            extra_sections=[
                _flag("RADII", "5E16.8", [0.0, 0.0], 5),
                _flag("SCREEN", "5E16.8", [0.0, 0.0], 5),
            ]
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


@pytest.mark.parametrize("bad_value", ["1.0junk", "NaN", "Inf", "1.0E9999"])
def test_strict_float_parser_rejects_junk_nonfinite_and_overflow(
    tmp_path, bad_value
):
    def mutate(path):
        _replace_flag(
            path,
            "CHARGE",
            _flag("CHARGE", "5E16.8", [bad_value, 0.0], 5),
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "invalid or non-finite floating-point field" in output


def test_nonzero_double_to_float_underflow_is_rejected(tmp_path):
    def mutate(path):
        _replace_flag(
            path,
            "LENNARD_JONES_ACOEF",
            _flag("LENNARD_JONES_ACOEF", "3E24.16", ["1.0E-46"] * 3, 3),
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "nonzero value that underflows the finite float range" in output


def test_finite_double_to_float_overflow_is_rejected(tmp_path):
    def mutate(path):
        _replace_flag(
            path,
            "LENNARD_JONES_ACOEF",
            _flag("LENNARD_JONES_ACOEF", "3E24.16", ["3.5E38"] * 3, 3),
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "outside the finite float range" in output


def test_explicit_underflowed_zero_is_accepted(tmp_path):
    def mutate(path):
        _replace_flag(
            path,
            "LENNARD_JONES_ACOEF",
            _flag("LENNARD_JONES_ACOEF", "3E24.16", ["0.0E-9999"] * 3, 3),
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


def test_smallest_float_subnormal_is_rejected_before_force_setup(tmp_path):
    def mutate(path):
        _replace_flag(
            path,
            "LENNARD_JONES_ACOEF",
            _flag(
                "LENNARD_JONES_ACOEF",
                "3E24.16",
                ["1.40129846E-45"] * 3,
                3,
            ),
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert (
        "AMBER section LENNARD_JONES_ACOEF contains a subnormal value" in output
    )
    assert "finite zero or normal float" in output


def test_fixed_width_parser_splits_adjacent_full_width_negative_integers(
    tmp_path,
):
    def mutate(path):
        _replace_flag(
            path,
            "ATOM_TYPE_INDEX",
            _flag("ATOM_TYPE_INDEX", "10I8", [-1234567, -2345678]),
        )

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "type index -1234567" in output
    assert "invalid integer field" not in output


def test_duplicate_singleton_flag_is_rejected(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=lambda path: _duplicate_flag(path, "MASS"),
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "singleton section MASS occurs more than once" in output


@pytest.mark.parametrize(
    "section_name",
    [
        "MASS",
        "CHARGE",
        "RESIDUE_POINTER",
        "NUMBER_EXCLUDED_ATOMS",
        "EXCLUDED_ATOMS_LIST",
    ],
)
def test_required_atom_and_exclusion_sections_cannot_be_missing(
    tmp_path, section_name
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=lambda path: _replace_flag(path, section_name, None),
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert f"required AMBER section {section_name} is missing" in output


@pytest.mark.parametrize(
    ("mutated_flag", "expected_error"),
    [
        (_flag("MASS", "5E16.8", [12.0], 5), "MASS has 1 values; expected 2"),
        (
            _flag("RESIDUE_POINTER", "10I8", [2]),
            "strictly increasing atom index beginning at 1",
        ),
        (
            _flag("NUMBER_EXCLUDED_ATOMS", "10I8", [1, 0]),
            "sums to 1; POINTERS NNB is 0",
        ),
    ],
    ids=["mass-length", "residue-first", "exclusion-sum"],
)
def test_atom_table_cross_section_invariants(
    tmp_path, mutated_flag, expected_error
):
    section_name = mutated_flag.splitlines()[0].split(maxsplit=1)[1]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        prmtop_mutator=lambda path: _replace_flag(
            path, section_name, mutated_flag
        ),
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output


@pytest.mark.parametrize(
    ("counts", "excluded", "expected_error"),
    [
        ([2, 0], [2, 2], "not strictly increasing"),
        ([2, 0], [0, 2], "contains invalid atom 0"),
        ([1, 0], [1], "invalid atom 1 for source atom 1"),
        ([1, 0], [3], "invalid atom 3"),
        ([0, 1], [1], "violates Amber triangular ordering"),
    ],
    ids=["nonmonotonic", "mixed-zero", "self", "range", "reverse"],
)
def test_excluded_atom_lists_are_strict(
    counts, excluded, expected_error, tmp_path
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            excluded_counts=counts,
            excluded_list=excluded,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output
    prmtop_path = tmp_path / "indexed_feature.prmtop"
    _assert_unsupported_location(output, prmtop_path, "NUMBER_EXCLUDED_ATOMS")
    _assert_unsupported_location(output, prmtop_path, "EXCLUDED_ATOMS_LIST")


@pytest.mark.parametrize(
    ("pointer_index", "value", "expected_error"),
    [
        (20, 1, "AMBER IFPERT=1 is not supported"),
        (29, 1, "AMBER IFCAP=1 is not supported"),
        (27, 3, "unsupported AMBER IFBOX value 3"),
    ],
)
def test_unsupported_pointer_modes_fail_fast(
    tmp_path, pointer_index, value, expected_error
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(pointer_overrides={pointer_index: value}),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "POINTERS"
    )


@pytest.mark.parametrize(
    ("pointer_length", "ncopy"),
    [(31, None), (32, 0), (32, 1)],
    ids=["classic-31", "legacy-ncopy-zero", "single-copy"],
)
def test_pointer_table_accepts_classic_and_single_copy_layouts(
    tmp_path, pointer_length, ncopy
):
    pointer_overrides = None if ncopy is None else {31: ncopy}
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            pointer_length=pointer_length,
            pointer_overrides=pointer_overrides,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


@pytest.mark.parametrize("pointer_length", [30, 33])
def test_pointer_table_rejects_nonstandard_lengths(tmp_path, pointer_length):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(pointer_length=pointer_length),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert f"POINTERS has {pointer_length} values; expected 31 or 32" in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "POINTERS"
    )


def test_multiple_ncopy_replicas_fail_with_pointer_source(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            pointer_length=32,
            pointer_overrides={31: 2},
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "NCOPY=2 requests multiple PIMD/LES replicas" in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "POINTERS"
    )


def test_nonzero_ipol_and_12_6_4_coefficients_fail_fast(tmp_path):
    cases = {
        "ipol": (
            [_flag("IPOL", "1I8", [1], 1)],
            "polarizable AMBER topology IPOL=1",
            "IPOL",
        ),
        "ccoef": (
            [_flag("LENNARD_JONES_CCOEF", "3E24.16", [0.0, 1.0, 0.0], 3)],
            "unsupported 12-6-4 nonbonded interactions",
            "LENNARD_JONES_CCOEF",
        ),
    }
    for case_name, (extra_sections, expected_error, flag_name) in cases.items():
        case_path = tmp_path / case_name
        case_path.mkdir()
        result, _ = _run_indexed_feature_sponge(
            case_path,
            coordinates=_TWO_ATOM_COORDINATES,
            **_two_atom_indexed_options(extra_sections=extra_sections),
        )
        output = result.stdout + "\n" + result.stderr
        assert result.returncode != 0
        assert expected_error in output
        _assert_unsupported_location(
            output, case_path / "indexed_feature.prmtop", flag_name
        )


def test_all_zero_lennard_jones_ccoef_placeholder_is_accepted(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            extra_sections=[
                _flag("LENNARD_JONES_CCOEF", "3E24.16", [0.0, 0.0, 0.0], 3)
            ]
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


@pytest.mark.parametrize(
    ("bond_row", "expected_error"),
    [
        ((1, 3, 1), "bond field A value 1 is not a multiple of 3"),
        ((0, 6, 1), "bond field B selects atom 2 outside [0, 2)"),
        ((0, 3, 2), "out-of-range parameter index"),
        ((0, 0, 1), "bond cannot connect an atom to itself"),
    ],
    ids=["encoding", "range", "type", "self"],
)
def test_bond_encoding_and_parameter_indices_are_strict(
    tmp_path, bond_row, expected_error
):
    extra_sections = [
        _flag("BOND_FORCE_CONSTANT", "5E16.8", [10.0], 5),
        _flag("BOND_EQUIL_VALUE", "5E16.8", [1.0], 5),
        _flag("BONDS_WITHOUT_HYDROGEN", "10I8", bond_row),
    ]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            pointer_overrides={3: 1, 15: 1},
            extra_sections=extra_sections,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output


@pytest.mark.parametrize(
    ("parameter_section", "expected_error"),
    [
        (
            _flag("BOND_FORCE_CONSTANT", "5E16.8", [-1.0], 5),
            "negative K or equilibrium distance",
        ),
        (
            _flag("ANGLE_EQUIL_VALUE", "5E16.8", [4.0], 5),
            "outside its physical domain",
        ),
    ],
)
def test_bond_and_angle_parameter_domains_are_strict(
    tmp_path, parameter_section, expected_error
):
    is_bond = "BOND_" in parameter_section.splitlines()[0]
    if is_bond:
        extra_sections = [
            parameter_section,
            _flag("BOND_EQUIL_VALUE", "5E16.8", [1.0], 5),
            _flag("BONDS_WITHOUT_HYDROGEN", "10I8", [0, 3, 1]),
        ]
        pointer_overrides = {3: 1, 15: 1}
    else:
        extra_sections = [
            _flag("ANGLE_FORCE_CONSTANT", "5E16.8", [1.0], 5),
            parameter_section,
            _flag("ANGLES_WITHOUT_HYDROGEN", "10I8", [0, 3, 6, 1]),
        ]
        pointer_overrides = {5: 1, 16: 1}
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=(),
            lj14_a=None,
            lj14_b=None,
            pointer_overrides=pointer_overrides,
            extra_sections=extra_sections,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output


def test_duplicate_active_proper_dihedrals_are_rejected(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=[(0, 3, 6, 9, 1), (0, 3, 6, 9, 1)],
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "duplicate active AMBER 1-4 endpoint pair 1-4" in output
    assert (
        "First active proper: DIHEDRALS_WITHOUT_HYDROGEN term 1 "
        "(atoms 1 2 3 4)" in output
    )
    assert (
        "Second active proper: DIHEDRALS_WITHOUT_HYDROGEN term 2 "
        "(atoms 1 2 3 4)" in output
    )
    assert (
        "continuation terms must encode atom C with a negative index" in output
    )


def test_negative_atom_c_marks_a_legal_multiterm_dihedral_continuation(
    tmp_path,
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=[(0, 3, -6, 9, 1), (0, 3, 6, 9, 1)],
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "dihedral_numbers is 2" in output
    assert "non-bond 14 numbers is 1" in output


def test_duplicate_nb14_endpoint_with_conflicting_scales_fails_fast(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=[(0, 3, 6, 9, 1), (0, 3, 6, 9, 2)],
            dihedral_force_constant=[1.0, 1.0],
            dihedral_periodicity=[1.0, 2.0],
            dihedral_phase=[0.0, 0.0],
            scee_scale_factor=[1.0, 1.2],
            scnb_scale_factor=[2.0, 3.0],
            pointer_overrides={17: 2},
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "duplicate active AMBER 1-4 endpoint pair" in output
    assert "DIHEDRALS_WITHOUT_HYDROGEN term 1" in output
    assert "DIHEDRALS_WITHOUT_HYDROGEN term 2" in output


def test_nonpositive_14_scale_for_proper_dihedral_fails_fast(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(scee_scale_factor=[0.0]),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "must be finite and positive" in output


def _urey_bradley_sections(*, include_equilibrium=True, ub_endpoints=(1, 3)):
    sections = [
        _flag("ANGLE_FORCE_CONSTANT", "5E16.8", [5.0], 5),
        _flag("ANGLE_EQUIL_VALUE", "5E16.8", [math.pi / 2], 5),
        _flag("ANGLES_WITHOUT_HYDROGEN", "10I8", [0, 3, 6, 1]),
        _flag("CHARMM_UREY_BRADLEY_COUNT", "2I8", [1, 1], 2),
        _flag(
            "CHARMM_UREY_BRADLEY",
            "10I8",
            [ub_endpoints[0], ub_endpoints[1], 1],
        ),
        _flag("CHARMM_UREY_BRADLEY_FORCE_CONSTANT", "5E16.8", [2.0], 5),
    ]
    if include_equilibrium:
        sections.append(
            _flag("CHARMM_UREY_BRADLEY_EQUIL_VALUE", "5E16.8", [1.0], 5)
        )
    return sections


def test_chamber_urey_bradley_moves_matching_angle_to_combined_ir(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=(),
            lj14_a=None,
            lj14_b=None,
            pointer_overrides={5: 1, 13: 8, 16: 1},
            extra_sections=_urey_bradley_sections(),
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "urey_bradley_numbers is 1" in output
    assert "ANGLE IS NOT INITIALIZED" in output
    assert "BOND IS NOT INITIALIZED" in output
    expected = 2.0 * (math.sqrt(2.0) - 1.0) ** 2
    assert math.isclose(
        _extract_mdout_term(mdout_path, "urey_bradley"), expected, abs_tol=0.01
    )


def test_chamber_urey_bradley_requires_all_sections(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=(),
            lj14_a=None,
            lj14_b=None,
            pointer_overrides={5: 1, 16: 1},
            extra_sections=_urey_bradley_sections(include_equilibrium=False),
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "CHAMBER Urey-Bradley sections are incomplete" in output


@pytest.mark.parametrize(
    ("extra_sections", "pointer_overrides", "expected_error"),
    [
        (
            _urey_bradley_sections(ub_endpoints=(1, 4)),
            {5: 1, 16: 1},
            "matches 0 standard angles",
        ),
        (
            [
                _flag("ANGLE_FORCE_CONSTANT", "5E16.8", [5.0], 5),
                _flag("ANGLE_EQUIL_VALUE", "5E16.8", [math.pi / 2], 5),
                _flag(
                    "ANGLES_WITHOUT_HYDROGEN",
                    "10I8",
                    [0, 3, 6, 1, 0, 9, 6, 1],
                ),
                *_urey_bradley_sections()[3:],
            ],
            {5: 2, 16: 1},
            "matches 2 standard angles",
        ),
    ],
    ids=["missing-angle", "ambiguous-angle"],
)
def test_chamber_urey_bradley_angle_matching_is_unique(
    tmp_path, extra_sections, pointer_overrides, expected_error
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FOUR_ATOM_COORDINATES,
        **_four_atom_indexed_options(
            dihedral_rows=(),
            lj14_a=None,
            lj14_b=None,
            pointer_overrides=pointer_overrides,
            extra_sections=extra_sections,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output


def _cmap_sections(
    *, resolution=4, grid=None, index=(1, 2, 3, 4, 5, 1), include_grid=True
):
    if grid is None:
        grid = [0.0] * (max(resolution, 0) ** 2)
    sections = [
        _flag("CHARMM_CMAP_COUNT", "2I8", [1, 1], 2),
        _flag("CHARMM_CMAP_RESOLUTION", "20I4", [resolution], 20),
    ]
    if include_grid:
        sections.append(_flag("CHARMM_CMAP_PARAMETER_01", "5E16.8", grid, 5))
    sections.append(_flag("CHARMM_CMAP_INDEX", "6I8", index, 6))
    return sections


_FIVE_ATOM_COORDINATES = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (1.0, 1.0, 0.0),
    (1.0, 1.0, 1.0),
    (2.0, 1.0, 1.0),
]


def _five_atom_cmap_options(**overrides):
    options = {
        "atom_types": [1] * 5,
        "nonbonded_parm_index": [1],
        "normal_a": [0.0],
        "normal_b": [0.0],
        "excluded_counts": [0] * 5,
        "excluded_list": [],
        "extra_sections": _cmap_sections(),
    }
    options.update(overrides)
    return options


def test_chamber_cmap_complete_grid_is_loaded(tmp_path):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FIVE_ATOM_COORDINATES,
        **_five_atom_cmap_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "total CMAP number is 1" in output
    assert math.isclose(
        _extract_mdout_term(mdout_path, "cmap"), 0.0, abs_tol=1e-6
    )


@pytest.mark.parametrize(
    ("sections", "expected_error"),
    [
        (_cmap_sections(include_grid=False), "CMAP_PARAMETER_01 is missing"),
        (_cmap_sections(resolution=0), "resolution 1 must be positive"),
        (_cmap_sections(grid=[0.0] * 15), "has 15 values; expected 16"),
        (_cmap_sections(index=(1, 2, 3, 4, 6, 1)), "atom index 6"),
        (_cmap_sections(index=(1, 2, 3, 4, 5, 2)), "parameter index 2"),
    ],
    ids=["missing-grid", "resolution", "grid-size", "atom", "type"],
)
def test_chamber_cmap_validation_fails_fast(tmp_path, sections, expected_error):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FIVE_ATOM_COORDINATES,
        **_five_atom_cmap_options(extra_sections=sections),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output


def test_chamber_cmap_rejects_kernel_index_overflow_before_reading_grid(
    tmp_path,
):
    # 11585^2 still fits INT_MAX/16, while 11586^2 does not.  Keep the grid
    # empty so this regression verifies the cumulative bound is checked before
    # parsing or allocating the declared parameter table.
    sections = _cmap_sections(resolution=1, grid=[])
    sections[1] = _flag("CHARMM_CMAP_RESOLUTION", "10I8", [11586])
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_FIVE_ATOM_COORDINATES,
        **_five_atom_cmap_options(extra_sections=sections),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert (
        "CMAP interpolation table exceeds the supported kernel index range"
        in output
    )


def test_rst7_time_does_not_imply_velocity_block(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        rst7_options={"time": 12.5},
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


@pytest.mark.parametrize(
    "rst7_options",
    [
        {"velocities": [(0.1, -0.2, 0.3)], "box": None},
        {"box": (20.0, 21.0, 22.0, 90.0, 90.0, 90.0)},
    ],
    ids=["three-field-velocity", "six-field-box"],
)
def test_rst7_only_single_atom_selects_velocity_or_box_by_field_count(
    tmp_path, rst7_options
):
    rst7_path = tmp_path / "single.rst7"
    mdin_path = tmp_path / "single.toml"
    _write_feature_rst7(rst7_path, [(1.0, 2.0, 3.0)], **rst7_options)
    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "rst7_only_single_atom"
            mode = "nve"
            step_limit = 0
            amber_rst7 = {str(rst7_path)!r}
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
    # A restart alone cannot initialize dynamics without masses, but parsing
    # must get past the formerly ambiguous one-line tail.
    assert result.returncode != 0
    assert "does not identify one unambiguous velocity/box layout" not in output
    assert "no mass information found" in output


def test_rst7_velocity_and_box_layout_is_selected_from_field_count(tmp_path):
    velocities = [(0.1, 0.2, 0.3), (-0.1, -0.2, -0.3)]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        rst7_options={"time": 1.0, "velocities": velocities},
        step_limit=1,
        extra_commands="vel = 'vel.dat'\nwrite_trajectory_interval = 1",
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    written = struct.unpack("<6f", (tmp_path / "vel.dat").read_bytes()[: 6 * 4])
    # Amber raw velocities are physical A/ps divided by 20.455, exactly the
    # numerical unit SPONGE uses internally; no reader-side scaling is needed.
    assert written == pytest.approx(
        [value for velocity in velocities for value in velocity], abs=1.0e-6
    )


def test_rst7_accepts_adjacent_full_width_positive_fields(tmp_path):
    coordinates = [
        (1000.0, 1001.0, 1002.0),
        (1003.0, 1004.0, 1005.0),
    ]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=coordinates,
        rst7_options={"box": None},
        **_two_atom_indexed_options(pointer_overrides={27: 0}),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    coordinate_line = (
        (tmp_path / "indexed_feature.rst7").read_text().splitlines()[2]
    )
    assert coordinate_line == "".join(
        f"{value:12.7f}" for value in range(1000, 1006)
    )
    assert " " not in coordinate_line


def test_rst7_accepts_whitespace_delimited_fortran_scientific_fields(tmp_path):
    def mutate(path):
        lines = path.read_text().splitlines()
        lines[2] = " ".join(
            [
                "0.0D+00",
                "0.0D+00",
                "0.0D+00",
                "1.7320508D+00",
                "0.0D+00",
                "0.0D+00",
            ]
        )
        path.write_text("\n".join(lines) + "\n")

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        rst7_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


def test_rst7_accepts_compact_whitespace_decimal_fields_and_trailing_blanks(
    tmp_path,
):
    def mutate(path):
        lines = path.read_text().splitlines()
        lines[2] = "0 0 0 1.7320508 0 0"
        path.write_text("\n".join(lines) + "\n\n\n")

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        rst7_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


@pytest.mark.parametrize("mutation", ["truncated", "trailing"])
def test_rst7_rejects_truncated_or_trailing_fixed_width_records(
    tmp_path, mutation
):
    def mutate(path):
        lines = path.read_text().splitlines()
        if mutation == "truncated":
            lines[2] = lines[2][:-1]
        else:
            lines[2] += "X"
        path.write_text("\n".join(lines) + "\n")

    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        rst7_mutator=mutate,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "invalid amber_rst7 coordinate record on line 3" in output


@pytest.mark.parametrize(
    ("rst7_options", "expected_error"),
    [
        (
            {"box": None},
            "does not identify one unambiguous velocity/box layout",
        ),
        (
            {"trailing": [7.0]},
            "invalid amber_rst7 box record",
        ),
        (
            {"time": "NaN"},
            "invalid numeric field in amber_rst7",
        ),
        (
            {"box": (50.0, 50.0, 50.0, 0.0, 90.0, 90.0)},
            "angles must lie in (0, 180)",
        ),
    ],
    ids=["missing-periodic-box", "trailing", "nonfinite-time", "box-angle"],
)
def test_rst7_rejects_invalid_or_trailing_layouts(
    tmp_path, rst7_options, expected_error
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        rst7_options=rst7_options,
        **_two_atom_indexed_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert expected_error in output


def test_nonperiodic_rst7_uses_explicit_nopbc_without_auxiliary_wrap(tmp_path):
    sigma = 500.0
    coordinates = [(-600.0, 0.0, 0.0), (600.0, 0.0, 0.0)]
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=coordinates,
        rst7_options={"time": 2.0, "box": None},
        # NOPBC is an exact all-pairs path; the periodic neighbor cutoff must
        # neither truncate this pair nor require a fabricated huge value.
        cutoff=8.0,
        **_two_atom_indexed_options(
            atom_types=[1, 1],
            nonbonded_parm_index=[1],
            normal_a=[4.0 * sigma**12],
            normal_b=[4.0 * sigma**6],
            pointer_overrides={27: 0},
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "may be inaccurate" not in output
    headers = mdout_path.read_text().splitlines()[0].split()
    assert "Coulomb" in headers
    assert "PM" not in headers
    expected = 4.0 * ((sigma / 1200.0) ** 12 - (sigma / 1200.0) ** 6)
    assert math.isclose(
        _extract_mdout_term(mdout_path, "LJ"), expected, abs_tol=0.01
    )


def test_nonperiodic_rst7_rejects_explicit_pbc_true(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        rst7_options={"box": None},
        extra_commands="pbc = true",
        **_two_atom_indexed_options(pointer_overrides={27: 0}),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "IFBOX=0 restart is non-periodic but pbc=true" in output


def _tip4p_virtual_sections(ep_distance=0.1, second_oh_distance=None):
    equilibrium_values = [1.0, ep_distance]
    second_oh_type = 1
    if second_oh_distance is not None:
        equilibrium_values.append(second_oh_distance)
        second_oh_type = 3
    return [
        _flag(
            "BOND_FORCE_CONSTANT",
            "5E16.8",
            [0.0] * len(equilibrium_values),
            5,
        ),
        _flag(
            "BOND_EQUIL_VALUE",
            "5E16.8",
            equilibrium_values,
            5,
        ),
        _flag(
            "BONDS_INC_HYDROGEN",
            "10I8",
            [0, 3, 1, 0, 6, second_oh_type],
        ),
        _flag(
            "BONDS_WITHOUT_HYDROGEN",
            "10I8",
            [0, 9, 2],
        ),
    ]


_TIP4P_COORDINATES = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0),
    (9.0, 9.0, 9.0),  # Deliberately wrong; the virtual refresh must replace it.
    (0.070710678, 0.070710678, 1.0),
]


def _tip4p_options(**overrides):
    options = {
        "atom_types": [1] * 5,
        "nonbonded_parm_index": [1],
        "normal_a": [0.0],
        "normal_b": [0.0],
        "charges": [0.0, 0.0, 0.0, 18.2223, -18.2223],
        "masses": [16.0, 1.0, 1.0, 0.0, 12.0],
        "atomic_numbers": [8, 1, 1, 0, 6],
        "amber_atom_types": ["OW", "HW", "HW", "EP", "C"],
        "excluded_counts": [0] * 5,
        "excluded_list": [],
        "pointer_overrides": {2: 2, 3: 1, 15: 2, 30: 1},
        "extra_sections": _tip4p_virtual_sections(),
    }
    options.update(overrides)
    return options


def test_numextra_tip4p_frame_refreshes_coordinates_and_redistributes_force(
    tmp_path,
):
    result, mdout_path = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        extra_commands="frc = 'frc.dat'",
        **_tip4p_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "Virtual Atoms Number is 1" in output
    # Refreshed EP is exactly 1 A from the probe, so q=+1/-1 gives -332.0636.
    assert math.isclose(
        _extract_mdout_term(mdout_path, "PM"), -332.0636, abs_tol=0.1
    )
    forces = struct.unpack(
        "<15f", (tmp_path / "frc.dat").read_bytes()[: 15 * 4]
    )
    assert all(math.isfinite(value) for value in forces)
    ep_force = forces[9:12]
    assert max(abs(value) for value in ep_force) < 1e-6
    assert max(abs(value) for value in forces[:9]) > 1.0


def _run_tip4p_energy_force(case_dir, coordinates):
    case_dir.mkdir()
    result, mdout_path = _run_indexed_feature_sponge(
        case_dir,
        coordinates=coordinates,
        extra_commands="frc = 'frc.dat'",
        **_tip4p_options(),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    forces = struct.unpack(
        "<15f", (case_dir / "frc.dat").read_bytes()[-15 * 4 :]
    )
    return _extract_mdout_term(mdout_path, "eff_pot"), forces


def test_numextra_tip4p_redistributed_force_matches_finite_difference(tmp_path):
    coordinates = [list(coordinate) for coordinate in _TIP4P_COORDINATES]
    _, forces = _run_tip4p_energy_force(tmp_path / "base", coordinates)

    epsilon = 2.0e-3
    for atom, axis in ((0, 2), (1, 2), (2, 2)):
        plus = [coordinate.copy() for coordinate in coordinates]
        minus = [coordinate.copy() for coordinate in coordinates]
        plus[atom][axis] += epsilon
        minus[atom][axis] -= epsilon
        energy_plus, _ = _run_tip4p_energy_force(
            tmp_path / f"plus_{atom}_{axis}", plus
        )
        energy_minus, _ = _run_tip4p_energy_force(
            tmp_path / f"minus_{atom}_{axis}", minus
        )
        numerical_force = -(energy_plus - energy_minus) / (2.0 * epsilon)
        assert forces[3 * atom + axis] == pytest.approx(
            numerical_force, rel=3.0e-3, abs=5.0e-2
        )


def test_numextra_type5_accepts_unequal_oh_equilibrium_distances(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            pointer_overrides={2: 2, 3: 1, 15: 3, 30: 1},
            extra_sections=_tip4p_virtual_sections(second_oh_distance=1.1),
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "Virtual type 5 atom numbers is 1" in output


def test_numextra_type5_zero_distance_is_supported(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            extra_sections=_tip4p_virtual_sections(ep_distance=0.0)
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "Virtual type 5 atom numbers is 1" in output


def test_numextra_trims_standard_angle_dihedral_and_derived_nb14(tmp_path):
    extra_sections = _tip4p_virtual_sections() + [
        _flag("ANGLE_FORCE_CONSTANT", "5E16.8", [25.0], 5),
        _flag("ANGLE_EQUIL_VALUE", "5E16.8", [1.2], 5),
        _flag("ANGLES_WITHOUT_HYDROGEN", "10I8", [0, 12, 9, 1]),
    ]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            dihedral_rows=[(3, 6, 9, 12, 1)],
            pointer_overrides={
                2: 2,
                3: 1,
                5: 1,
                15: 2,
                16: 1,
                30: 1,
            },
            extra_sections=extra_sections,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "Virtual type 5 atom numbers is 1" in output
    assert "ANGLE IS NOT INITIALIZED" in output
    assert "DIHEDRAL IS NOT INITIALIZED" in output
    assert "NB14 IS NOT INITIALIZED" in output
    assert "bond_numbers is 2" in output


@pytest.mark.parametrize(
    ("extra_sections", "expected_error"),
    [
        ([], "must have exactly one parent edge"),
    ],
    ids=["no-parent"],
)
def test_numextra_unsupported_frames_fail_with_atom_context(
    tmp_path, extra_sections, expected_error
):
    pointer_overrides = {30: 1}
    if extra_sections:
        pointer_overrides.update({2: 2, 3: 1, 15: 2})
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            extra_sections=extra_sections,
            pointer_overrides=pointer_overrides,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "NUMEXTRA atom 4" in output
    assert expected_error in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "AMBER_ATOM_TYPE"
    )
    assert "%FLAG BONDS_WITHOUT_HYDROGEN is absent" in output


@pytest.mark.parametrize(
    ("atomic_numbers", "masses"),
    [
        (None, [16.0, 1.0, 1.0, 0.0, 12.0]),
        ([6, 7, 9, 14, 6], [22.0, 3.0, 4.0, 42.0, 12.0]),
    ],
    ids=["legacy-no-atomic-number", "nonstandard-elements-and-masses"],
)
def test_numextra_frame_uses_topology_not_element_or_mass_gates(
    tmp_path, atomic_numbers, masses
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(atomic_numbers=atomic_numbers, masses=masses),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "Virtual type 5 atom numbers is 1" in output


@pytest.mark.parametrize(
    "non_ep_type", ["LP", " EP"], ids=["lp", "leading-space"]
)
def test_numextra_requires_exact_ep_atom_type_with_source_context(
    tmp_path, non_ep_type
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            amber_atom_types=["OW", "HW", "HW", non_ep_type, "C"],
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "NUMEXTRA=1" in output
    assert "identifies 0 exact EP atom(s)" in output
    prmtop_path = tmp_path / "indexed_feature.prmtop"
    _assert_unsupported_location(output, prmtop_path, "POINTERS")
    _assert_unsupported_location(output, prmtop_path, "AMBER_ATOM_TYPE")


def test_numextra_missing_atom_type_section_reports_pointer_and_absence(
    tmp_path,
):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(amber_atom_types=None),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "NUMEXTRA=1 requires AMBER_ATOM_TYPE" in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "POINTERS"
    )
    assert "%FLAG AMBER_ATOM_TYPE is absent" in output


def test_atomic_number_zero_does_not_mark_an_ordinary_atom_as_ep(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TWO_ATOM_COORDINATES,
        **_two_atom_indexed_options(
            masses=[12.0, 14.0],
            atomic_numbers=[0, 6],
            amber_atom_types=["C", "N"],
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "FF Virtual Atoms Number is 1" not in output


def test_numextra_tip5p_owner_topology_is_rejected_with_bond_sources(tmp_path):
    tip5p_sections = _tip4p_virtual_sections()
    tip5p_sections[-1] = _flag(
        "BONDS_WITHOUT_HYDROGEN", "10I8", [0, 9, 2, 0, 12, 2]
    )
    coordinates = _TIP4P_COORDINATES[:4] + [
        (-9.0, -9.0, -9.0),
        _TIP4P_COORDINATES[4],
    ]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=coordinates,
        atom_types=[1] * 6,
        nonbonded_parm_index=[1],
        normal_a=[0.0],
        normal_b=[0.0],
        charges=[0.0, 0.0, 0.0, 9.11115, 9.11115, -18.2223],
        masses=[16.0, 1.0, 1.0, 0.0, 0.0, 12.0],
        amber_atom_types=["OW", "HW", "HW", "EP", "EP", "C"],
        excluded_counts=[0] * 6,
        excluded_list=[],
        pointer_overrides={2: 2, 3: 2, 15: 2, 30: 2},
        extra_sections=tip5p_sections,
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "not a Classic TIP4P frame" in output
    assert "found 2 and 2" in output
    prmtop_path = tmp_path / "indexed_feature.prmtop"
    _assert_unsupported_location(output, prmtop_path, "BONDS_WITHOUT_HYDROGEN")
    _assert_unsupported_location(output, prmtop_path, "BONDS_INC_HYDROGEN")


def test_numextra_in_chamber_urey_bradley_fails_with_tuple_source(tmp_path):
    extension_sections = [
        _flag("ANGLE_FORCE_CONSTANT", "5E16.8", [5.0], 5),
        _flag("ANGLE_EQUIL_VALUE", "5E16.8", [math.pi / 2], 5),
        _flag("ANGLES_WITHOUT_HYDROGEN", "10I8", [0, 12, 9, 1]),
        _flag("CHARMM_UREY_BRADLEY_COUNT", "2I8", [1, 1], 2),
        _flag("CHARMM_UREY_BRADLEY", "10I8", [1, 4, 1]),
        _flag("CHARMM_UREY_BRADLEY_FORCE_CONSTANT", "5E16.8", [2.0], 5),
        _flag("CHARMM_UREY_BRADLEY_EQUIL_VALUE", "5E16.8", [1.0], 5),
    ]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            pointer_overrides={
                2: 2,
                3: 1,
                5: 1,
                15: 2,
                16: 1,
                30: 1,
            },
            extra_sections=_tip4p_virtual_sections() + extension_sections,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "cannot participate in CHARMM_UREY_BRADLEY terms" in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "CHARMM_UREY_BRADLEY"
    )


def test_numextra_in_chamber_cmap_fails_with_index_source(tmp_path):
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            extra_sections=_tip4p_virtual_sections()
            + _cmap_sections(index=(1, 2, 3, 4, 5, 1)),
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "cannot participate in CMAP terms" in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "CHARMM_CMAP_INDEX"
    )


def test_numextra_in_chamber_improper_fails_with_tuple_source(tmp_path):
    improper_sections = [
        _flag("CHARMM_NUM_IMPROPERS", "10I8", [1]),
        _flag("CHARMM_IMPROPERS", "10I8", [1, 2, 3, 4, 1]),
        _flag("CHARMM_NUM_IMPR_TYPES", "1I8", [1]),
        _flag("CHARMM_IMPROPER_FORCE_CONSTANT", "5E16.8", [2.0], 5),
        _flag("CHARMM_IMPROPER_PHASE", "5E16.8", [0.0], 5),
    ]
    result, _ = _run_indexed_feature_sponge(
        tmp_path,
        coordinates=_TIP4P_COORDINATES,
        **_tip4p_options(
            extra_sections=_tip4p_virtual_sections() + improper_sections,
        ),
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0
    assert "cannot participate in CHARMM_IMPROPERS terms" in output
    _assert_unsupported_location(
        output, tmp_path / "indexed_feature.prmtop", "CHARMM_IMPROPERS"
    )
