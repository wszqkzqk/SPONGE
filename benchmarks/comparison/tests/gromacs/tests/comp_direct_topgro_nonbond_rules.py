import json
import math
import os
import subprocess
import textwrap
from pathlib import Path

import pytest

KJ_TO_KCAL = 1.0 / 4.184
GROMACS_COULOMB_KJ_NM = 138.935458


def _toml_string(value):
    return json.dumps(os.fspath(value))


def _gro_atom(resid, resname, atomname, atomnr, xyz):
    x, y, z = xyz
    return (
        f"{resid:5d}{resname:<5}{atomname:>5}{atomnr:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )


def _two_atom_topology(
    *,
    defaults="1 2 yes 1.0 1.0",
    atom_a=(0.3, 0.4184),
    atom_b=(0.3, 0.4184),
    global_parameters="",
    molecule_parameters="",
    nrexcl=0,
    charges=(0.0, 0.0),
):
    defaults_section = ""
    if defaults is not None:
        defaults_section = f"[ defaults ]\n{defaults}\n"
    return (
        textwrap.dedent(
            f"""
            {defaults_section}
            [ atomtypes ]
            A A 12.0 0.0 A {atom_a[0]} {atom_a[1]}
            B B 12.0 0.0 A {atom_b[0]} {atom_b[1]}

            {global_parameters}

            [ moleculetype ]
            PAIR {nrexcl}

            [ atoms ]
            1 A 1 PAIR A 1 {charges[0]} 12.0
            2 B 1 PAIR B 2 {charges[1]} 12.0

            {molecule_parameters}

            [ system ]
            two-atom nonbonded rule test

            [ molecules ]
            PAIR 1
            """
        ).strip()
        + "\n"
    )


def _bonded_pair(pair_line="1 2 1"):
    return textwrap.dedent(
        f"""
        [ bonds ]
        1 2 1 0.5 0.0

        [ pairs ]
        {pair_line}
        """
    ).strip()


def _run_topology(case_dir, topology, *, distance_nm=0.5):
    case_dir = Path(case_dir)
    case_dir.mkdir(parents=True, exist_ok=True)
    top_path = case_dir / "topol.top"
    gro_path = case_dir / "conf.gro"
    mdin_path = case_dir / "mdin.spg.toml"
    mdout_path = case_dir / "mdout.txt"

    top_path.write_text(topology)
    gro_path.write_text(
        "\n".join(
            [
                "two-atom nonbonded rule test",
                "2",
                _gro_atom(1, "PAIR", "A", 1, (0.0, 0.0, 0.0)),
                _gro_atom(1, "PAIR", "B", 2, (distance_nm, 0.0, 0.0)),
                "   5.00000   5.00000   5.00000",
            ]
        )
        + "\n"
    )
    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = "direct_gromacs_nonbond_rules"
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
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    return result, mdout_path


def _extract_mdout_terms(mdout_path):
    lines = Path(mdout_path).read_text().splitlines()
    headers = lines[0].split()
    values = lines[1].split()
    assert len(headers) == len(values)
    return {name: float(value) for name, value in zip(headers, values)}


def _assert_succeeded(result):
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output


def _assert_rejected(result):
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "GROMACS" in output, output


def _combined_lj_energy(comb_rule, atom_a, atom_b, distance_nm):
    if comb_rule == 1:
        c6 = math.sqrt(atom_a[0] * atom_b[0])
        c12 = math.sqrt(atom_a[1] * atom_b[1])
    elif comb_rule == 2:
        sigma = 0.5 * (abs(atom_a[0]) + abs(atom_b[0]))
        if atom_a[0] < 0.0 or atom_b[0] < 0.0:
            sigma = -sigma
        epsilon = math.sqrt(atom_a[1] * atom_b[1])
        c6, c12 = _sigma_epsilon_coefficients((sigma, epsilon))
    elif comb_rule == 3:
        sigma = math.sqrt(abs(atom_a[0] * atom_b[0]))
        if atom_a[0] < 0.0 or atom_b[0] < 0.0:
            sigma = -sigma
        epsilon = math.sqrt(atom_a[1] * atom_b[1])
        c6, c12 = _sigma_epsilon_coefficients((sigma, epsilon))
    else:
        raise AssertionError(f"unsupported test combination rule {comb_rule}")
    return (c12 / distance_nm**12 - c6 / distance_nm**6) * KJ_TO_KCAL


def _sigma_epsilon_coefficients(parameters):
    sigma, epsilon_kj = parameters
    sigma6 = abs(sigma) ** 6
    c6 = 0.0 if sigma < 0.0 else 4.0 * epsilon_kj * sigma6
    c12 = 4.0 * epsilon_kj * sigma6**2
    return c6, c12


def _explicit_lj_energy(comb_rule, parameters, distance_nm):
    if comb_rule == 1:
        c6, c12 = parameters
    else:
        c6, c12 = _sigma_epsilon_coefficients(parameters)
    return (c12 / distance_nm**12 - c6 / distance_nm**6) * KJ_TO_KCAL


@pytest.mark.parametrize(
    "comb_rule,atom_a,atom_b",
    [
        (1, (0.01, 0.0004), (0.04, 0.0036)),
        (2, (0.2, 0.4184), (0.6, 1.6736)),
        (3, (0.2, 0.4184), (0.6, 1.6736)),
    ],
    ids=["c6-c12-geometric", "sigma-arithmetic", "sigma-geometric"],
)
def test_direct_gromacs_combination_rules(tmp_path, comb_rule, atom_a, atom_b):
    topology = _two_atom_topology(
        defaults=f"1 {comb_rule} yes 1.0 1.0",
        atom_a=atom_a,
        atom_b=atom_b,
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    expected = _combined_lj_energy(comb_rule, atom_a, atom_b, 0.5)
    assert _extract_mdout_terms(mdout_path)["LJ_short"] == pytest.approx(
        expected, abs=0.011
    )


@pytest.mark.parametrize(
    "comb_rule,parameters",
    [(1, (0.02, 0.0012)), (2, (0.4, 0.8368)), (3, (0.4, 0.8368))],
    ids=["c6-c12", "sigma-epsilon-rule2", "sigma-epsilon-rule3"],
)
def test_direct_gromacs_nonbond_params_override_regular_lj(
    tmp_path, comb_rule, parameters
):
    topology = _two_atom_topology(
        defaults=f"1 {comb_rule} yes 1.0 1.0",
        atom_a=(0.0, 0.0),
        atom_b=(0.0, 0.0),
        global_parameters=textwrap.dedent(
            f"""
            [ nonbond_params ]
            A B 1 {parameters[0]} {parameters[1]}
            """
        ).strip(),
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    expected = _explicit_lj_energy(comb_rule, parameters, 0.5)
    assert _extract_mdout_terms(mdout_path)["LJ_short"] == pytest.approx(
        expected, abs=0.011
    )


def test_direct_gromacs_pairtypes_and_nonbond_params_have_disjoint_scopes(
    tmp_path,
):
    global_parameters = textwrap.dedent(
        """
        [ nonbond_params ]
        A B 1 0.4 0.8368

        [ pairtypes ]
        A B 1 0.45 1.6736
        """
    ).strip()
    common = dict(
        defaults="1 2 yes 0.1 1.0",
        atom_a=(0.3, 0.4184),
        atom_b=(0.3, 0.4184),
        global_parameters=global_parameters,
    )

    regular_result, regular_mdout = _run_topology(
        tmp_path / "regular", _two_atom_topology(**common)
    )
    _assert_succeeded(regular_result)
    assert _extract_mdout_terms(regular_mdout)["LJ_short"] == pytest.approx(
        _explicit_lj_energy(2, (0.4, 0.8368), 0.5), abs=0.011
    )

    pair_result, pair_mdout = _run_topology(
        tmp_path / "pair",
        _two_atom_topology(
            **common,
            nrexcl=1,
            molecule_parameters=_bonded_pair(),
        ),
    )
    _assert_succeeded(pair_result)
    pair_terms = _extract_mdout_terms(pair_mdout)
    assert pair_terms["LJ_short"] == pytest.approx(0.0, abs=0.011)
    # Explicit pair parameters are not scaled by fudgeLJ.
    assert pair_terms["nb14_LJ"] == pytest.approx(
        _explicit_lj_energy(2, (0.45, 1.6736), 0.5), abs=0.011
    )


@pytest.mark.parametrize(
    "comb_rule,atom_a,atom_b,nonbond_parameters",
    [
        (1, (0.01, 0.0004), (0.04, 0.0036), (0.0, 0.0)),
        (2, (0.2, 0.4184), (0.6, 1.6736), (0.1, 0.0)),
        (3, (0.2, 0.4184), (0.6, 1.6736), (0.1, 0.0)),
    ],
    ids=["rule1", "rule2", "rule3"],
)
def test_direct_gromacs_generated_pair_uses_resolved_nonbonded_pair_and_fudge_lj(
    tmp_path, comb_rule, atom_a, atom_b, nonbond_parameters
):
    fudge_lj = 0.5
    topology = _two_atom_topology(
        defaults=f"1 {comb_rule} yes {fudge_lj} 1.0",
        atom_a=atom_a,
        atom_b=atom_b,
        global_parameters=textwrap.dedent(
            f"""
            [ nonbond_params ]
            A B 1 {nonbond_parameters[0]} {nonbond_parameters[1]}
            """
        ).strip(),
        nrexcl=1,
        molecule_parameters=_bonded_pair(),
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    terms = _extract_mdout_terms(mdout_path)
    expected = fudge_lj * _explicit_lj_energy(
        comb_rule, nonbond_parameters, 0.5
    )
    assert terms["LJ_short"] == pytest.approx(0.0, abs=0.011)
    assert terms["nb14_LJ"] == pytest.approx(expected, abs=0.011)


def test_direct_gromacs_generated_pair_applies_fudge_qq_independently(
    tmp_path,
):
    fudge_qq = 0.25
    charges = (0.2, -0.3)
    topology = _two_atom_topology(
        defaults=f"1 2 yes 0.0 {fudge_qq}",
        atom_a=(0.3, 0.0),
        atom_b=(0.5, 0.0),
        charges=charges,
        nrexcl=1,
        molecule_parameters=_bonded_pair(),
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    terms = _extract_mdout_terms(mdout_path)
    expected = (
        fudge_qq
        * charges[0]
        * charges[1]
        * GROMACS_COULOMB_KJ_NM
        / 0.5
        * KJ_TO_KCAL
    )
    assert terms["nb14_LJ"] == pytest.approx(0.0, abs=0.011)
    assert terms["nb14_EE"] == pytest.approx(expected, abs=0.011)


def test_direct_gromacs_two_field_defaults_use_optional_field_defaults(
    tmp_path,
):
    atom_a = (0.2, 0.4184)
    atom_b = (0.6, 1.6736)
    topology = _two_atom_topology(
        defaults="1 2",
        atom_a=atom_a,
        atom_b=atom_b,
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    assert _extract_mdout_terms(mdout_path)["LJ_short"] == pytest.approx(
        _combined_lj_energy(2, atom_a, atom_b, 0.5), abs=0.011
    )


@pytest.mark.parametrize(
    "defaults,fudge_lj",
    [
        ("1 2 yes", 1.0),
        ("1 2 yes 0.5", 0.5),
        ("1 2 yes 0.5 1.0 12", 0.5),
    ],
    ids=["three-fields", "four-fields", "six-fields-n12"],
)
def test_direct_gromacs_optional_defaults_fields_generate_pairs(
    tmp_path, defaults, fudge_lj
):
    atom_a = (0.2, 0.4184)
    atom_b = (0.6, 1.6736)
    topology = _two_atom_topology(
        defaults=defaults,
        atom_a=atom_a,
        atom_b=atom_b,
        nrexcl=1,
        molecule_parameters=_bonded_pair(),
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    expected = fudge_lj * _combined_lj_energy(2, atom_a, atom_b, 0.5)
    assert _extract_mdout_terms(mdout_path)["nb14_LJ"] == pytest.approx(
        expected, abs=0.011
    )


def test_direct_gromacs_inline_pair_overrides_pairtype_without_fudge_lj(
    tmp_path,
):
    inline_parameters = (0.35, 0.8368)
    topology = _two_atom_topology(
        defaults="1 2 yes 0.1 1.0",
        atom_a=(0.3, 0.4184),
        atom_b=(0.3, 0.4184),
        global_parameters=textwrap.dedent(
            """
            [ pairtypes ]
            A B 1 0.45 1.6736
            """
        ).strip(),
        nrexcl=1,
        molecule_parameters=_bonded_pair(
            f"1 2 1 {inline_parameters[0]} {inline_parameters[1]}"
        ),
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    assert _extract_mdout_terms(mdout_path)["nb14_LJ"] == pytest.approx(
        _explicit_lj_energy(2, inline_parameters, 0.5), abs=0.011
    )


@pytest.mark.parametrize("section", ["nonbond_params", "pairtypes"])
def test_direct_gromacs_duplicate_pair_parameters_are_symmetric_last_wins(
    tmp_path, section
):
    global_parameters = textwrap.dedent(
        f"""
        [ {section} ]
        A B 1 0.25 0.4184
        B A 1 0.45 1.6736
        """
    ).strip()
    pair_section = {}
    if section == "pairtypes":
        pair_section = dict(
            nrexcl=1,
            molecule_parameters=_bonded_pair(),
        )
    topology = _two_atom_topology(
        defaults="1 2 yes 1.0 1.0",
        global_parameters=global_parameters,
        **pair_section,
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    term = "nb14_LJ" if section == "pairtypes" else "LJ_short"
    assert _extract_mdout_terms(mdout_path)[term] == pytest.approx(
        _explicit_lj_energy(2, (0.45, 1.6736), 0.5), abs=0.011
    )


def test_direct_gromacs_gen_pairs_no_requires_explicit_pair_parameters(
    tmp_path,
):
    topology = _two_atom_topology(
        defaults="1 2 no 1.0 1.0",
        nrexcl=1,
        molecule_parameters=_bonded_pair(),
    )
    result, _ = _run_topology(tmp_path, topology)
    _assert_rejected(result)


@pytest.mark.parametrize(
    "defaults",
    [
        "2 2 yes 1.0 1.0",
        "0 2 yes 1.0 1.0",
        "1 0 yes 1.0 1.0",
        "1 4 yes 1.0 1.0",
        "1 2 maybe 1.0 1.0",
        "1 2 yes -0.1 1.0",
        "1 2 yes 1.0 -0.1",
        "1 2 yes nan 1.0",
        "1 2 yes 1.0 inf",
        "1 2 yes 1.0 1.0 13",
        "1 2 yes 1.0 1.0 12 unexpected",
    ],
    ids=[
        "buckingham",
        "invalid-nbfunc",
        "comb-rule-zero",
        "comb-rule-four",
        "invalid-gen-pairs",
        "negative-fudge-lj",
        "negative-fudge-qq",
        "nonfinite-fudge-lj",
        "nonfinite-fudge-qq",
        "unsupported-repulsion-power",
        "too-many-fields",
    ],
)
def test_direct_gromacs_rejects_unsupported_defaults(tmp_path, defaults):
    result, _ = _run_topology(tmp_path, _two_atom_topology(defaults=defaults))
    _assert_rejected(result)


def test_direct_gromacs_requires_defaults_section(tmp_path):
    result, _ = _run_topology(tmp_path, _two_atom_topology(defaults=None))
    _assert_rejected(result)


@pytest.mark.parametrize(
    "global_parameters,molecule_parameters",
    [
        ("[ nonbond_params ]\nA B 2 0.4 0.8368", ""),
        ("[ pairtypes ]\nA B 2 0.4 0.8368", _bonded_pair()),
        ("", _bonded_pair("1 2 2")),
    ],
    ids=["nonbond-params-funct2", "pairtypes-funct2", "pairs-funct2"],
)
def test_direct_gromacs_rejects_unsupported_pair_function_types(
    tmp_path, global_parameters, molecule_parameters
):
    topology = _two_atom_topology(
        global_parameters=global_parameters,
        molecule_parameters=molecule_parameters,
        nrexcl=1 if molecule_parameters else 0,
    )
    result, _ = _run_topology(tmp_path, topology)
    _assert_rejected(result)


@pytest.mark.parametrize(
    "global_parameters,molecule_parameters",
    [
        ("[ nonbond_params ]\nA B 1 0.4", ""),
        ("[ nonbond_params ]\nA B 1 0.4 0.8368 9.9", ""),
        ("[ pairtypes ]\nA B 1 0.4", _bonded_pair()),
        ("[ pairtypes ]\nA B 1 0.4 0.8368 9.9", _bonded_pair()),
        ("", _bonded_pair("1 2 1 0.4")),
        ("", _bonded_pair("1 2 1 0.4 0.8368 9.9")),
    ],
    ids=[
        "nonbond-missing",
        "nonbond-extra",
        "pairtype-missing",
        "pairtype-extra",
        "inline-pair-missing",
        "inline-pair-extra",
    ],
)
def test_direct_gromacs_rejects_wrong_pair_parameter_counts(
    tmp_path, global_parameters, molecule_parameters
):
    topology = _two_atom_topology(
        global_parameters=global_parameters,
        molecule_parameters=molecule_parameters,
        nrexcl=1 if molecule_parameters else 0,
    )
    result, _ = _run_topology(tmp_path, topology)
    _assert_rejected(result)


@pytest.mark.parametrize(
    (
        "comb_rule,atom_a,atom_b,explicit_parameters,global_parameters,"
        "molecule_parameters,term"
    ),
    [
        (2, (-0.5, 0.8368), (-0.7, 0.8368), None, "", "", "LJ_short"),
        (2, (-0.5, 0.8368), (0.7, 0.8368), None, "", "", "LJ_short"),
        (3, (-0.4, 0.8368), (0.6, 0.8368), None, "", "", "LJ_short"),
        (
            2,
            (0.3, 0.4184),
            (0.3, 0.4184),
            (-0.55, 0.8368),
            "[ nonbond_params ]\nA B 1 -0.55 0.8368",
            "",
            "LJ_short",
        ),
        (
            3,
            (0.3, 0.4184),
            (0.3, 0.4184),
            (-0.55, 0.8368),
            "[ pairtypes ]\nA B 1 -0.55 0.8368",
            _bonded_pair(),
            "nb14_LJ",
        ),
        (
            2,
            (0.3, 0.4184),
            (0.3, 0.4184),
            (-0.55, 0.8368),
            "",
            _bonded_pair("1 2 1 -0.55 0.8368"),
            "nb14_LJ",
        ),
    ],
    ids=[
        "atomtype-rule2",
        "atomtype-rule2-mixed-sign",
        "atomtype-rule3",
        "nonbond-params",
        "pairtypes",
        "inline-pair",
    ],
)
def test_direct_gromacs_negative_sigma_is_purely_repulsive(
    tmp_path,
    comb_rule,
    atom_a,
    atom_b,
    explicit_parameters,
    global_parameters,
    molecule_parameters,
    term,
):
    """Negative sigma means C6=0 while C12 is calculated as usual.

    For generated rule-2/rule-3 pairs, GROMACS combines absolute sigma
    magnitudes and propagates a negative sentinel from either input type.
    """
    topology = _two_atom_topology(
        defaults=f"1 {comb_rule} yes 1.0 1.0",
        atom_a=atom_a,
        atom_b=atom_b,
        global_parameters=global_parameters,
        molecule_parameters=molecule_parameters,
        nrexcl=1 if molecule_parameters else 0,
    )
    result, mdout_path = _run_topology(tmp_path, topology)
    _assert_succeeded(result)

    if explicit_parameters is None:
        expected = _combined_lj_energy(comb_rule, atom_a, atom_b, 0.5)
    else:
        expected = _explicit_lj_energy(comb_rule, explicit_parameters, 0.5)
    terms = _extract_mdout_terms(mdout_path)
    assert terms[term] == pytest.approx(expected, abs=0.011)
    assert terms[term] > 0.0


@pytest.mark.parametrize(
    "case_name,defaults,global_parameters,molecule_parameters,term,parameters",
    [
        (
            "negative-c6-c12-nonbond",
            "1 1 yes 1.0 1.0",
            "[ nonbond_params ]\nA B 1 -0.02 -0.0012",
            "",
            "LJ_short",
            (-0.02, -0.0012),
        ),
        (
            "negative-c6-c12-pairtype",
            "1 1 yes 1.0 1.0",
            "[ pairtypes ]\nA B 1 -0.02 -0.0012",
            _bonded_pair(),
            "nb14_LJ",
            (-0.02, -0.0012),
        ),
        (
            "negative-c6-c12-inline-pair",
            "1 1 yes 1.0 1.0",
            "",
            _bonded_pair("1 2 1 -0.02 -0.0012"),
            "nb14_LJ",
            (-0.02, -0.0012),
        ),
        (
            "negative-epsilon-nonbond",
            "1 2 yes 1.0 1.0",
            "[ nonbond_params ]\nA B 1 0.4 -0.8368",
            "",
            "LJ_short",
            (0.4, -0.8368),
        ),
    ],
    ids=lambda value: value if isinstance(value, str) else None,
)
def test_direct_gromacs_explicit_negative_lj_parameters_preserve_sign(
    tmp_path,
    case_name,
    defaults,
    global_parameters,
    molecule_parameters,
    term,
    parameters,
):
    topology = _two_atom_topology(
        defaults=defaults,
        global_parameters=global_parameters,
        molecule_parameters=molecule_parameters,
        nrexcl=1 if molecule_parameters else 0,
    )
    result, mdout_path = _run_topology(tmp_path / case_name, topology)
    _assert_succeeded(result)
    comb_rule = int(defaults.split()[1])
    assert _extract_mdout_terms(mdout_path)[term] == pytest.approx(
        _explicit_lj_energy(comb_rule, parameters, 0.5), abs=0.011
    )


@pytest.mark.parametrize(
    "comb_rule,atom_a,atom_b",
    [
        (1, (-0.01, 0.0004), (0.04, 0.0036)),
        (2, (0.3, -0.4184), (0.5, 0.8368)),
        (3, (0.3, -0.4184), (0.5, 0.8368)),
    ],
    ids=["rule1-c6", "rule2-epsilon", "rule3-epsilon"],
)
def test_direct_gromacs_rejects_generated_lj_with_mixed_parameter_signs(
    tmp_path, comb_rule, atom_a, atom_b
):
    topology = _two_atom_topology(
        defaults=f"1 {comb_rule} yes 1.0 1.0",
        atom_a=atom_a,
        atom_b=atom_b,
    )
    result, _ = _run_topology(tmp_path, topology)
    _assert_rejected(result)
