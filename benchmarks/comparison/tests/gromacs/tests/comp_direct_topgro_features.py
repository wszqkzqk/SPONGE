import json
import math
import os
import struct
import subprocess
import textwrap
from pathlib import Path

import pytest


def _toml_string(value):
    return json.dumps(os.fspath(value))


def _gro_atom(
    resid,
    resname,
    atomname,
    atomnr,
    xyz,
    *,
    velocity=None,
    precision=3,
):
    x, y, z = xyz
    width = precision + 5
    line = (
        f"{resid:5d}{resname:<5}{atomname:>5}{atomnr:5d}"
        f"{x:{width}.{precision}f}"
        f"{y:{width}.{precision}f}"
        f"{z:{width}.{precision}f}"
    )
    if velocity is not None:
        vx, vy, vz = velocity
        line += (
            f"{vx:{width}.{precision + 1}f}"
            f"{vy:{width}.{precision + 1}f}"
            f"{vz:{width}.{precision + 1}f}"
        )
    return line


def _extract_mdout_term(mdout_path, term_name):
    lines = Path(mdout_path).read_text().splitlines()
    headers = lines[0].split()
    values = lines[1].split()
    assert len(headers) == len(values)
    return float(dict(zip(headers, values))[term_name])


def _run_direct_topology(
    tmp_path,
    case_name,
    topology,
    coordinates,
    extra_files=None,
    gromacs_define=None,
    gro_identities=None,
    timeout=120,
):
    case_dir = tmp_path / case_name
    case_dir.mkdir()
    top_path = case_dir / "topol.top"
    gro_path = case_dir / "conf.gro"
    mdin_path = case_dir / "mdin.spg.toml"
    mdout_path = case_dir / "mdout.txt"

    top_path.write_text(textwrap.dedent(topology).strip() + "\n")
    for filename, contents in (extra_files or {}).items():
        (case_dir / filename).write_text(
            textwrap.dedent(contents).strip() + "\n"
        )
    if gro_identities is None:
        gro_identities = [
            ("MOL", f"A{index}") for index in range(1, len(coordinates) + 1)
        ]
    assert len(gro_identities) == len(coordinates)
    gro_lines = [case_name, str(len(coordinates))]
    gro_lines.extend(
        _gro_atom(1, resname, atomname, index, coordinate)
        for index, ((resname, atomname), coordinate) in enumerate(
            zip(gro_identities, coordinates), start=1
        )
    )
    gro_lines.append("   5.00000   5.00000   5.00000")
    gro_path.write_text("\n".join(gro_lines) + "\n")
    define_setting = (
        ""
        if gromacs_define is None
        else f"gromacs_define = {_toml_string(gromacs_define)}"
    )
    mdin_path.write_text(
        textwrap.dedent(
            f"""
            md_name = {_toml_string(case_name)}
            mode = "nve"
            step_limit = 0
            dt = 0
            cutoff = 8.0
            constrain_mode = "SETTLE"
            gromacs_top = {_toml_string(top_path)}
            gromacs_gro = {_toml_string(gro_path)}
            {define_setting}
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
        timeout=timeout,
    )
    return result, mdout_path


def _run_custom_gro_case(
    tmp_path,
    case_name,
    topology,
    gro_lines,
    *,
    mdin_overrides=None,
    timeout=120,
):
    case_dir = tmp_path / case_name
    case_dir.mkdir()
    top_path = case_dir / "topol.top"
    gro_path = case_dir / "conf.gro"
    mdin_path = case_dir / "mdin.spg.toml"
    top_path.write_text(textwrap.dedent(topology).strip() + "\n")
    gro_path.write_text("\n".join(gro_lines) + "\n")
    settings = {
        "md_name": case_name,
        "mode": "nve",
        "step_limit": 0,
        "dt": 0,
        "cutoff": 8.0,
        "constrain_mode": "SETTLE",
        "gromacs_top": os.fspath(top_path),
        "gromacs_gro": os.fspath(gro_path),
        "mdout": os.fspath(case_dir / "mdout.txt"),
        "print_zeroth_frame": True,
        "write_mdout_interval": 1,
        "write_information_interval": 1,
        "write_trajectory_interval": 1,
        "write_restart_file_interval": 0,
    }
    settings.update(mdin_overrides or {})
    mdin_path.write_text(
        "\n".join(
            f"{key} = {json.dumps(value)}" for key, value in settings.items()
        )
        + "\n"
    )
    result = subprocess.run(
        [os.environ.get("SPONGE_BIN", "SPONGE"), "-mdin", str(mdin_path)],
        cwd=case_dir,
        capture_output=True,
        text=True,
        check=False,
        timeout=timeout,
    )
    return result, case_dir


def _read_float_vectors(path):
    raw = Path(path).read_bytes()
    assert len(raw) % 12 == 0
    return struct.unpack(f"={len(raw) // 4}f", raw)


def _minimal_five_atom_topology(extra_section, type_section=""):
    return f"""
        [ defaults ]
        1 2 yes 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.3 0.4184

        {type_section}

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1 0.0 12.0
        2 A 1 MOL A2 1 0.0 12.0
        3 A 1 MOL A3 1 0.0 12.0
        4 A 1 MOL A4 1 0.0 12.0
        5 A 1 MOL A5 1 0.0 12.0

        {extra_section}

        [ system ]
        parser validation

        [ molecules ]
        MOL 1
    """


def _virtual_site_topology(virtual_site_section, *, ptype="V", probe=False):
    bond_lines = ["1 2 1 0.4 0.0", "2 3 1 0.5 0.0"]
    bond_section = "[ bonds ]\n" + "\n".join(bond_lines)
    probe_type = "Q Q 12.0 0.0 A 0.3 0.4184" if probe else ""
    probe_molecule = (
        """
        [ moleculetype ]
        PROBE 0

        [ atoms ]
        1 Q 1 PRB PRB 1 0.0 12.0
        """
        if probe
        else ""
    )
    probe_count = "PROBE 1" if probe else ""
    return f"""
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        P P 12.0 0.0 A 0.0 0.0
        X X 0.0 0.0 {ptype} 0.3 0.4184
        {probe_type}

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 P 1 MOL P1 1 0.0 12.0
        2 P 1 MOL P2 1 0.0 12.0
        3 P 1 MOL P3 1 0.0 12.0
        4 X 1 MOL VS 1 0.0 0.0

        {virtual_site_section}
        {bond_section}

        [ exclusions ]
        4 1 2 3

        {probe_molecule}

        [ system ]
        virtual-site mapping

        [ molecules ]
        MOL 1
        {probe_count}
    """


@pytest.mark.parametrize(
    "case_name,virtual_site_section",
    [
        ("vsite1", "[ virtual_sites1 ]\n4 1"),
        ("vsite2_f1", "[ virtual_sites2 ]\n4 1 2 1 0.25"),
        ("vsite2_f2", "[ virtual_sites2 ]\n4 1 2 2 0.20"),
        ("vsite3_f1", "[ virtual_sites3 ]\n4 1 2 3 1 0.20 0.30"),
        # Martini writes function-2 parameters in this a,d order.
        ("vsite3_f2", "[ virtual_sites3 ]\n4 1 2 3 2 0.25 0.20"),
    ],
)
def test_direct_gromacs_maps_supported_virtual_site_coordinates(
    tmp_path, case_name, virtual_site_section
):
    parents = [
        (0.10, 0.20, 0.30),
        (0.50, -0.20, 0.40),
        (-0.10, 0.80, 0.00),
    ]
    r1, r2, r3 = parents

    def add(first, second):
        return tuple(a + b for a, b in zip(first, second))

    def subtract(first, second):
        return tuple(a - b for a, b in zip(first, second))

    def scale(value, vector):
        return tuple(value * component for component in vector)

    if case_name == "vsite1":
        expected_nm = r1
    elif case_name == "vsite2_f1":
        expected_nm = add(r1, scale(0.25, subtract(r2, r1)))
    elif case_name == "vsite2_f2":
        direction = subtract(r2, r1)
        norm = math.sqrt(sum(component * component for component in direction))
        expected_nm = add(r1, scale(0.20 / norm, direction))
    elif case_name == "vsite3_f1":
        expected_nm = add(
            r1,
            add(
                scale(0.20, subtract(r2, r1)),
                scale(0.30, subtract(r3, r1)),
            ),
        )
    else:
        direction = add(subtract(r2, r1), scale(0.25, subtract(r3, r2)))
        norm = math.sqrt(sum(component * component for component in direction))
        expected_nm = add(r1, scale(0.20 / norm, direction))

    probe = (0.42, 0.35, 0.15)
    gro_lines = [case_name, "5"] + [
        _gro_atom(1, "MOL", atom_name, index, coordinate)
        for index, (atom_name, coordinate) in enumerate(
            zip(("P1", "P2", "P3", "VS"), parents + [(9.0, 9.0, 9.0)]),
            start=1,
        )
    ]
    gro_lines.append(_gro_atom(2, "PRB", "PRB", 5, probe))
    gro_lines.append("   5.00000   5.00000   5.00000")
    result, case_dir = _run_custom_gro_case(
        tmp_path,
        case_name,
        _virtual_site_topology(virtual_site_section, probe=True),
        gro_lines,
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    distance_angstrom = 10.0 * math.sqrt(
        sum((a - b) ** 2 for a, b in zip(expected_nm, probe))
    )
    sigma_over_r = 3.0 / distance_angstrom
    expected_lj = 0.4 * (sigma_over_r**12 - sigma_over_r**6)
    assert _extract_mdout_term(case_dir / "mdout.txt", "LJ_short") == (
        pytest.approx(expected_lj, abs=5.0e-3)
    )


@pytest.mark.parametrize("ptype", ["S", "V", "D"])
def test_direct_gromacs_allows_nonatom_particle_types_as_virtual_targets(
    tmp_path, ptype
):
    result, _ = _run_direct_topology(
        tmp_path,
        f"virtual_target_{ptype}",
        f"""
        [ defaults ]
        1 2

        [ atomtypes ]
        P P 12.0 0.0 A 0.0 0.0
        X X 0.0 0.0 {ptype} 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 P 1 MOL A1 1
        2 X 1 MOL A2 1

        [ virtual_sites1 ]
        2 1

        [ system ]
        supported virtual target

        [ molecules ]
        MOL 1
        """,
        [(0.1, 0.2, 0.3), (9.0, 9.0, 9.0)],
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


@pytest.mark.parametrize(
    "entry,error_fragment",
    [
        (
            "[ virtual_sites3 ]\n5 1 2 3 3 90.0 0.2",
            "[ virtual_sites3 ] function 3",
        ),
        (
            "[ virtual_sites3 ]\n5 1 2 3 4 0.1 0.2 0.3",
            "[ virtual_sites3 ] function 4",
        ),
        ("[ virtual_sites4 ]\n5 1 2 3 4 1 0.1 0.2 0.3", "virtual_sites4"),
        ("[ virtual_sitesn ]\n5 1 1 2 3", "virtual_sitesn"),
    ],
)
def test_direct_gromacs_rejects_unrepresentable_virtual_sites_at_source(
    tmp_path, entry, error_fragment
):
    result, _ = _run_direct_topology(
        tmp_path,
        "unsupported_virtual_site",
        _minimal_five_atom_topology(entry),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output
    assert "topol.top:" in output


def test_direct_gromacs_rejects_virtual_site_dependency_cycle_at_source(
    tmp_path,
):
    result, _ = _run_direct_topology(
        tmp_path,
        "virtual_cycle",
        _minimal_five_atom_topology("[ virtual_sites1 ]\n4 5\n5 4"),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "virtual-site source graph contains a dependency cycle" in output
    assert "topol.top:" in output


def test_direct_gromacs_rejects_zero_distance_virtual_site_direction(
    tmp_path,
):
    gro_lines = [
        "singular distance virtual site",
        "4",
        _gro_atom(1, "MOL", "P1", 1, (0.1, 0.2, 0.3)),
        _gro_atom(1, "MOL", "P2", 2, (0.1, 0.2, 0.3)),
        _gro_atom(1, "MOL", "P3", 3, (0.4, 0.5, 0.6)),
        _gro_atom(1, "MOL", "VS", 4, (9.0, 9.0, 9.0)),
        "   5.00000   5.00000   5.00000",
    ]
    result, _ = _run_custom_gro_case(
        tmp_path,
        "singular_distance_virtual_site",
        _virtual_site_topology("[ virtual_sites2 ]\n4 1 2 2 0.20"),
        gro_lines,
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "zero construction direction" in output
    assert "topol.top:" in output


def test_direct_gromacs_rejects_overflowed_virtual_site_direction_norm(
    tmp_path,
):
    gro_lines = [
        "overflowed distance virtual-site direction norm",
        "4",
        _gro_atom(1, "MOL", "P1", 1, (0.0, 0.0, 0.0), precision=4),
        _gro_atom(1, "MOL", "P2", 2, (0.0, 0.0006, 0.0), precision=4),
        _gro_atom(1, "MOL", "P3", 3, (0.0, 0.0, 0.1), precision=4),
        _gro_atom(1, "MOL", "VS", 4, (0.0, 0.0, 0.0), precision=4),
        # A finite, reversible, canonical cell can still be sufficiently
        # ill-conditioned that minimum-image subtraction yields a finite
        # direction whose single-precision squared norm overflows.
        "1e9 1e-3 1 0 0 1e29 0 0 0",
    ]
    result, _ = _run_custom_gro_case(
        tmp_path,
        "overflowed_distance_virtual_site_direction_norm",
        _virtual_site_topology("[ virtual_sites2 ]\n4 1 2 2 0.20"),
        gro_lines,
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "non-finite construction geometry" in output
    assert "topol.top:" in output


def _run_virtual_site_energy_force(
    tmp_path, case_name, virtual_site_section, parents
):
    probe = (0.42, 0.35, 0.15)
    gro_lines = [case_name, "5"] + [
        _gro_atom(1, "MOL", atom_name, index, coordinate, precision=6)
        for index, (atom_name, coordinate) in enumerate(
            zip(("P1", "P2", "P3", "VS"), parents + [(9.0, 9.0, 9.0)]),
            start=1,
        )
    ]
    gro_lines.append(_gro_atom(2, "PRB", "PRB", 5, probe, precision=6))
    gro_lines.append("   5.00000   5.00000   5.00000")
    result, case_dir = _run_custom_gro_case(
        tmp_path,
        case_name,
        _virtual_site_topology(virtual_site_section, probe=True),
        gro_lines,
        mdin_overrides={"frc": "frc.dat"},
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    energy = _extract_mdout_term(case_dir / "mdout.txt", "eff_pot")
    forces = _read_float_vectors(case_dir / "frc.dat")
    assert len(forces) == 15
    return energy, forces


@pytest.mark.parametrize(
    "case_name,virtual_site_section,parent_indices",
    [
        ("force_vsite2_f1", "[ virtual_sites2 ]\n4 1 2 1 0.25", (0, 1)),
        (
            "force_vsite3_f1",
            "[ virtual_sites3 ]\n4 1 2 3 1 0.20 0.30",
            (0, 1, 2),
        ),
        (
            "force_vsite3_f2",
            "[ virtual_sites3 ]\n4 1 2 3 2 0.25 0.20",
            (0, 1, 2),
        ),
    ],
)
def test_direct_gromacs_virtual_site_parent_forces_match_finite_difference(
    tmp_path, case_name, virtual_site_section, parent_indices
):
    parents = [
        (0.10, 0.20, 0.30),
        (0.50, -0.20, 0.40),
        (-0.10, 0.80, 0.00),
    ]
    _, forces = _run_virtual_site_energy_force(
        tmp_path, f"{case_name}_base", virtual_site_section, parents
    )
    epsilon_nm = 2.0e-4
    for parent_index in parent_indices:
        axis = parent_index % 3
        plus = [list(coordinate) for coordinate in parents]
        minus = [list(coordinate) for coordinate in parents]
        plus[parent_index][axis] += epsilon_nm
        minus[parent_index][axis] -= epsilon_nm
        energy_plus, _ = _run_virtual_site_energy_force(
            tmp_path,
            f"{case_name}_plus_{parent_index}_{axis}",
            virtual_site_section,
            [tuple(coordinate) for coordinate in plus],
        )
        energy_minus, _ = _run_virtual_site_energy_force(
            tmp_path,
            f"{case_name}_minus_{parent_index}_{axis}",
            virtual_site_section,
            [tuple(coordinate) for coordinate in minus],
        )
        # GRO displacements are nm while SPONGE forces are kcal/mol/angstrom.
        numerical_force = -(energy_plus - energy_minus) / (20.0 * epsilon_nm)
        assert forces[3 * parent_index + axis] == pytest.approx(
            numerical_force, rel=5.0e-3, abs=5.0e-3
        )


def test_direct_gromacs_requires_zero_mass_for_virtual_site_target(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "massive_virtual_target",
        _minimal_five_atom_topology("[ virtual_sites1 ]\n5 1"),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "virtual-site target must have zero mass" in output
    assert "topol.top:" in output


@pytest.mark.parametrize("mass", ["0", "-1.0"])
def test_direct_gromacs_requires_positive_inline_atom_mass(tmp_path, mass):
    topology = f"""
        [ defaults ]
        1 2 yes 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.3 0.4184

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1 0.0 {mass}

        [ system ]
        invalid inline atom mass

        [ molecules ]
        MOL 1
    """
    result, _ = _run_direct_topology(
        tmp_path,
        f"inline_mass_{mass.replace('-', 'negative_')}",
        topology,
        [(0.0, 0.0, 0.0)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "GROMACS atom mass must be positive" in output
    assert "topol.top:11" in output


@pytest.mark.parametrize("mass", ["0", "-1.0"])
def test_direct_gromacs_requires_positive_inherited_atom_mass(tmp_path, mass):
    topology = """
        [ defaults ]
        1 2 yes 1.0 1.0

        #include "mass.itp"

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        invalid inherited atom mass

        [ molecules ]
        MOL 1
    """
    result, _ = _run_direct_topology(
        tmp_path,
        f"inherited_mass_{mass.replace('-', 'negative_')}",
        topology,
        [(0.0, 0.0, 0.0)],
        extra_files={
            "mass.itp": f"""
            [ atomtypes ]
            A A {mass} 0.0 A 0.3 0.4184
            """,
        },
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "GROMACS atom mass must be positive" in output
    assert "mass.itp:2" in output


@pytest.mark.parametrize(
    "case_name,interactions,expected_marker",
    [
        (
            "constant_virtual_bond",
            """
            [ virtual_sites3 ]
            5 1 2 3 2 0.25 0.20

            [ bonds ]
            5 1 1 0.20 100.0
            """,
            "bond_numbers is 1",
        ),
        (
            "constant_virtual_angle",
            """
            [ virtual_sites3 ]
            5 1 2 3 1 0.20 0.30

            [ constraints ]
            1 2 1 0.20
            2 3 1 0.2828427
            3 1 1 0.20

            [ angles ]
            1 5 2 1 109.5 100.0
            """,
            "UREY BRADLEY IS NOT INITIALIZED",
        ),
        (
            "constant_virtual_dihedral",
            """
            [ virtual_sites3 ]
            5 1 2 3 1 0.20 0.30

            [ dihedrals ]
            1 2 5 3 1 0.0 1.0 1
            """,
            "DIHEDRAL IS NOT INITIALIZED",
        ),
    ],
)
def test_direct_gromacs_removes_only_constant_virtual_site_bondeds(
    tmp_path, case_name, interactions, expected_marker
):
    topology = _minimal_five_atom_topology(
        f"{interactions}\n\n[ bonds ]\n3 4 1 0.30 0.0"
    ).replace("5 A 1 MOL A5 1 0.0 12.0", "5 A 1 MOL A5 1 0.0 0.0")
    result, _ = _run_direct_topology(
        tmp_path,
        case_name,
        topology,
        [
            (0.0, 0.0, 0.0),
            (0.2, 0.0, 0.0),
            (0.0, 0.2, 0.0),
            (0.1, 0.1, 0.3),
            (9.0, 9.0, 9.0),
        ],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert expected_marker in output


def test_direct_gromacs_retains_nonconstant_virtual_site_bondeds_and_cmap(
    tmp_path,
):
    topology = _minimal_five_atom_topology(
        """
        [ virtual_sites3 ]
        5 1 2 3 1 0.20 0.30

        [ bonds ]
        5 4 1 0.20 100.0

        [ angles ]
        1 5 4 1 109.5 100.0

        [ dihedrals ]
        1 2 5 4 1 0.0 1.0 1

        [ cmap ]
        1 2 3 4 5 1
        """,
        """
        [ cmaptypes ]
        A A A A A 1 2 2 0.0 0.0 0.0 0.0
        """,
    ).replace("5 A 1 MOL A5 1 0.0 12.0", "5 A 1 MOL A5 1 0.0 0.0")
    result, _ = _run_direct_topology(
        tmp_path,
        "retained_virtual_bondeds",
        topology,
        [
            (0.0, 0.0, 0.0),
            (0.2, 0.0, 0.0),
            (0.2, 0.2, 0.0),
            (0.2, 0.2, 0.2),
            (9.0, 9.0, 9.0),
        ],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "bond_numbers is 1" in output
    assert "urey_bradley_numbers is 1" in output
    assert "dihedral_numbers is 1" in output
    assert "total CMAP number is 1" in output


def test_direct_gromacs_generates_exclusions_before_virtual_bond_cleanup(
    tmp_path,
):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "virtual_bond_exclusion_order",
        """
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        P P 12.0 0.0 A 0.3 0.4184
        Z Z 12.0 0.0 A 0.0 0.0
        X X 0.0 0.0 A 0.3 0.4184

        [ moleculetype ]
        MOL 1

        [ atoms ]
        1 P 1 MOL A1 1 0.0 12.0
        2 Z 1 MOL A2 1 0.0 12.0
        3 Z 1 MOL A3 1 0.0 12.0
        4 X 1 MOL A4 1 0.0 0.0

        [ virtual_sites3 ]
        4 1 2 3 2 0.25 0.20

        [ bonds ]
        4 1 1 0.20 100.0

        [ system ]
        exclusion before virtual bonded cleanup

        [ molecules ]
        MOL 1
        """,
        [
            (0.0, 0.0, 0.0),
            (0.2, 0.0, 0.0),
            (0.0, 0.2, 0.0),
            (9.0, 9.0, 9.0),
        ],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "BOND IS NOT INITIALIZED" in output
    assert _extract_mdout_term(mdout_path, "LJ_short") == pytest.approx(
        0.0, abs=1.0e-7
    )


def test_direct_gromacs_removes_fixed_virtual_site_constraint(tmp_path):
    topology = _minimal_five_atom_topology(
        """
        [ virtual_sites3 ]
        5 1 2 3 2 0.25 0.20

        [ constraints ]
        5 1 1 0.20

        [ bonds ]
        3 4 1 0.30 0.0
        """
    ).replace("5 A 1 MOL A5 1 0.0 12.0", "5 A 1 MOL A5 1 0.0 0.0")
    result, _ = _run_direct_topology(
        tmp_path,
        "fixed_virtual_constraint",
        topology,
        [
            (0.0, 0.0, 0.0),
            (0.2, 0.0, 0.0),
            (0.0, 0.2, 0.0),
            (0.1, 0.1, 0.3),
            (9.0, 9.0, 9.0),
        ],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert "constrain pair number is 0" in output
    assert "bond_numbers is 1" in output


def test_direct_gromacs_rejects_nonconstant_virtual_site_constraint_at_source(
    tmp_path,
):
    topology = _minimal_five_atom_topology(
        """
        [ virtual_sites3 ]
        5 1 2 3 1 0.20 0.30

        [ constraints ]
        5 4 1 0.20
        """
    ).replace("5 A 1 MOL A5 1 0.0 12.0", "5 A 1 MOL A5 1 0.0 0.0")
    result, _ = _run_direct_topology(
        tmp_path,
        "nonconstant_virtual_constraint",
        topology,
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "constraint involving a virtual-site target remains" in output
    assert "topol.top:" in output


def test_direct_gromacs_rejects_virtual_site_settle_at_source(tmp_path):
    topology = _minimal_five_atom_topology(
        """
        [ virtual_sites1 ]
        3 1

        [ settles ]
        3 1 0.10 0.16
        """
    ).replace("3 A 1 MOL A3 1 0.0 12.0", "3 A 1 MOL A3 1 0.0 0.0")
    result, _ = _run_direct_topology(
        tmp_path,
        "virtual_settle",
        topology,
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "SETTLE constraint involving a virtual-site target" in output
    assert "topol.top:" in output


def _strict_gro_topology():
    return """
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 RESIDUE FIRSTLONG 1 0.0 12.0
        2 A 1 RESIDUE SECONDLONG 1 0.0 12.0

        [ bonds ]
        1 2 1 0.2 0.0

        [ system ]
        strict GRO records

        [ molecules ]
        MOL 1
    """


def test_direct_gromacs_strict_gro_loads_time_precision_velocity_and_wrap(
    tmp_path,
):
    coordinates_nm = [(0.12345, 0.23456, 0.34567), (0.45678, 0.56789, 0.67891)]
    velocities_nm_ps = [
        (0.012345, -0.023456, 0.034567),
        (-0.045678, 0.056789, 0.067891),
    ]
    gro_lines = ["strict variable precision t= 12.375", "2"]
    gro_lines.extend(
        _gro_atom(
            0,
            "RESID",
            atom_name,
            0,
            coordinate,
            velocity=velocity,
            precision=5,
        )
        for atom_name, coordinate, velocity in zip(
            ("FIRST", "SECON"), coordinates_nm, velocities_nm_ps
        )
    )
    gro_lines.extend(("   5.00000   5.00000   5.00000", "", "   "))
    result, case_dir = _run_custom_gro_case(
        tmp_path,
        "strict_gro_success",
        _strict_gro_topology(),
        gro_lines,
        mdin_overrides={"step_limit": 1, "crd": "crd.dat", "vel": "vel.dat"},
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output
    assert _extract_mdout_term(case_dir / "mdout.txt", "time") == pytest.approx(
        12.375, abs=1.0e-8
    )
    observed_coordinates = _read_float_vectors(case_dir / "crd.dat")
    observed_velocities = _read_float_vectors(case_dir / "vel.dat")
    expected_coordinates = tuple(
        10.0 * component
        for coordinate in coordinates_nm
        for component in coordinate
    )
    expected_velocities = tuple(
        10.0 * component / 20.455
        for velocity in velocities_nm_ps
        for component in velocity
    )
    assert observed_coordinates == pytest.approx(
        expected_coordinates, abs=3.0e-6
    )
    assert observed_velocities == pytest.approx(expected_velocities, abs=2.0e-7)


@pytest.mark.parametrize(
    ("field_name", "coordinate", "velocity"),
    [
        ("coordinate", (1.0e-40, 0.2, 0.3), None),
        ("velocity", (0.1, 0.2, 0.3), (1.0e-38, 0.02, 0.03)),
    ],
)
def test_direct_gromacs_rejects_subnormal_gro_conversions(
    tmp_path, field_name, coordinate, velocity
):
    first = _gro_atom(
        1,
        "RESID",
        "FIRST",
        1,
        coordinate,
        velocity=velocity,
        precision=40,
    )
    second = _gro_atom(
        1,
        "RESID",
        "SECON",
        2,
        (0.4, 0.5, 0.6),
        velocity=(0.04, 0.05, 0.06) if velocity is not None else None,
        precision=40,
    )
    result, _ = _run_custom_gro_case(
        tmp_path,
        f"subnormal_gro_{field_name}",
        _strict_gro_topology(),
        [
            f"subnormal {field_name}",
            "2",
            first,
            second,
            "   5.00000   5.00000   5.00000",
        ],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"{field_name} field converts to a subnormal float" in output
    assert "conf.gro:3" in output


def test_direct_gromacs_single_atom_ion_display_names_are_unambiguous(tmp_path):
    topology = """
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        I I 22.99 0.0 A 0.0 0.0

        [ moleculetype ]
        SOD 0
        [ atoms ]
        1 I 1 SOD SOD 1

        [ moleculetype ]
        CLA 0
        [ atoms ]
        1 I 1 CLA CLA 1

        [ system ]
        conventional ion display names

        [ molecules ]
        SOD 1
        CLA 1
    """
    gro_lines = [
        "ion aliases",
        "2",
        _gro_atom(1, "SOD", "NA", 1, (0.1, 0.2, 0.3)),
        _gro_atom(2, "CLA", "CL", 2, (0.4, 0.5, 0.6)),
        "   5.00000   5.00000   5.00000",
    ]
    result, _ = _run_custom_gro_case(
        tmp_path, "ion_display_names", topology, gro_lines
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_direct_gromacs_rejects_ambiguous_single_atom_alias_order(tmp_path):
    topology = """
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        NA NA 22.99 0.0 A 0.0 0.0
        CL CL 35.45 0.0 A 0.0 0.0

        [ moleculetype ]
        SOD 0
        [ atoms ]
        1 NA 1 ION NA 1

        [ moleculetype ]
        CLA 0
        [ atoms ]
        1 CL 1 ION CL 1

        [ system ]
        ambiguous ion display names

        [ molecules ]
        SOD 1
        CLA 1
    """
    gro_lines = [
        "swapped ambiguous ions",
        "2",
        _gro_atom(1, "ION", "CL", 1, (0.1, 0.2, 0.3)),
        _gro_atom(2, "ION", "NA", 2, (0.4, 0.5, 0.6)),
        "   5.00000   5.00000   5.00000",
    ]
    result, _ = _run_custom_gro_case(
        tmp_path, "ambiguous_ion_aliases", topology, gro_lines
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "residue/atom name order does not match" in output
    assert "conf.gro:3" in output


@pytest.mark.parametrize(
    "case_name,first_velocity,second_velocity,mutator",
    [
        ("mixed_velocity_first", True, False, None),
        ("mixed_velocity_second", False, True, None),
        (
            "nonfinite_coordinate",
            False,
            False,
            lambda line: line[:20] + "     nan" + line[28:],
        ),
        (
            "nonfinite_velocity",
            True,
            True,
            lambda line: line[:44] + "     inf" + line[52:],
        ),
        (
            "truncated_coordinate",
            False,
            False,
            lambda line: line[:-1],
        ),
        (
            "atom_trailing_garbage",
            False,
            False,
            lambda line: line + " trailing",
        ),
    ],
)
def test_direct_gromacs_rejects_malformed_gro_atom_records_with_source(
    tmp_path, case_name, first_velocity, second_velocity, mutator
):
    velocity = (0.01, 0.02, 0.03)
    first = _gro_atom(
        1,
        "RESID",
        "FIRST",
        1,
        (0.1, 0.2, 0.3),
        velocity=velocity if first_velocity else None,
    )
    second = _gro_atom(
        1,
        "RESID",
        "SECON",
        2,
        (0.4, 0.5, 0.6),
        velocity=velocity if second_velocity else None,
    )
    if mutator is not None:
        second = mutator(second)
    result, _ = _run_custom_gro_case(
        tmp_path,
        case_name,
        _strict_gro_topology(),
        [
            case_name,
            "2",
            first,
            second,
            "   5.00000   5.00000   5.00000",
        ],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "conf.gro:4" in output


@pytest.mark.parametrize(
    "residue_names,atom_names",
    [
        (("WRONG", "RESID"), ("FIRST", "SECON")),
        (("RESID", "RESID"), ("SECON", "FIRST")),
    ],
)
def test_direct_gromacs_rejects_gro_name_order_mismatch_with_source(
    tmp_path, residue_names, atom_names
):
    gro_lines = ["bad names", "2"]
    gro_lines.extend(
        _gro_atom(1, residue_name, atom_name, index, coordinate)
        for index, (residue_name, atom_name, coordinate) in enumerate(
            zip(
                residue_names,
                atom_names,
                ((0.1, 0.2, 0.3), (0.4, 0.5, 0.6)),
            ),
            start=1,
        )
    )
    gro_lines.append("   5.00000   5.00000   5.00000")
    result, _ = _run_custom_gro_case(
        tmp_path, "bad_gro_names", _strict_gro_topology(), gro_lines
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "residue/atom name order does not match" in output
    assert "conf.gro:" in output


@pytest.mark.parametrize(
    "trailing_lines",
    [
        ["unexpected data"],
        [
            "second frame",
            "2",
            _gro_atom(1, "RESID", "FIRST", 1, (0.1, 0.2, 0.3)),
            _gro_atom(1, "RESID", "SECON", 2, (0.4, 0.5, 0.6)),
            "   5.00000   5.00000   5.00000",
        ],
    ],
)
def test_direct_gromacs_rejects_data_after_first_gro_frame(
    tmp_path, trailing_lines
):
    gro_lines = [
        "first frame",
        "2",
        _gro_atom(1, "RESID", "FIRST", 1, (0.1, 0.2, 0.3)),
        _gro_atom(1, "RESID", "SECON", 2, (0.4, 0.5, 0.6)),
        "   5.00000   5.00000   5.00000",
        *trailing_lines,
    ]
    result, _ = _run_custom_gro_case(
        tmp_path, "trailing_gro_frame", _strict_gro_topology(), gro_lines
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "additional frame" in output
    assert "conf.gro:6" in output


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
            B C 2 0.0 8.368

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


def test_direct_gromacs_rejects_negative_harmonic_improper_at_source(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "negative_harmonic_improper",
        _minimal_five_atom_topology('#include "negative-improper.itp"'),
        [
            (0.0, 0.1, 0.0),
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
            (0.1, 0.0, 0.1),
            (0.2, 0.0, 0.1),
        ],
        extra_files={
            "negative-improper.itp": """
                [ dihedrals ]
                1 2 3 4 2 0.0 -8.368
            """,
        },
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "harmonic improper force constant must be non-negative" in output
    assert "negative-improper.itp:2" in output


def test_direct_gromacs_rejects_active_position_restraints_with_location(
    tmp_path,
):
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

    assert result.returncode != 0
    output = result.stdout + "\n" + result.stderr
    error = "unsupported GROMACS topology section [ position_restraints ]"
    assert error in output
    assert f"{include_path}:1" in output


def test_direct_gromacs_applies_and_deduplicates_explicit_exclusions(tmp_path):
    topology_template = """
        [ defaults ]
        1 2 yes 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.3 0.4184

        [ moleculetype ]
        PAIR 0

        [ atoms ]
        1 A 1 MOL A1 1 0.0 12.0
        2 A 1 MOL A2 1 0.0 12.0

        {exclusions}

        [ system ]
        explicit exclusions

        [ molecules ]
        PAIR 1
    """
    without_exclusions, without_mdout = _run_direct_topology(
        tmp_path,
        "without_exclusions",
        topology_template.format(exclusions=""),
        [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0)],
    )
    assert without_exclusions.returncode == 0, (
        without_exclusions.stdout + "\n" + without_exclusions.stderr
    )

    with_exclusions, with_mdout = _run_direct_topology(
        tmp_path,
        "with_exclusions",
        topology_template.format(
            exclusions="""
            [ exclusions ]
            1 2 2
            2 1
            """
        ),
        [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0)],
    )
    assert with_exclusions.returncode == 0, (
        with_exclusions.stdout + "\n" + with_exclusions.stderr
    )

    sigma_over_r = 0.3 / 0.4
    expected_lj = 4.0 * 0.1 * (sigma_over_r**12 - sigma_over_r**6)
    assert _extract_mdout_term(without_mdout, "LJ_short") == pytest.approx(
        expected_lj, abs=5.0e-3
    )
    assert _extract_mdout_term(with_mdout, "LJ_short") == pytest.approx(
        0.0, abs=1.0e-6
    )


def test_direct_gromacs_rejects_invalid_explicit_exclusion_index(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "invalid_exclusion_index",
        """
        [ defaults ]
        1 2 yes 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.3 0.4184

        [ moleculetype ]
        PAIR 0

        [ atoms ]
        1 A 1 MOL A1 1 0.0 12.0
        2 A 1 MOL A2 1 0.0 12.0

        [ exclusions ]
        1 3

        [ system ]
        invalid exclusion

        [ molecules ]
        PAIR 1
        """,
        [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0)],
    )

    assert result.returncode != 0
    output = result.stdout + "\n" + result.stderr
    assert "invalid atom pair in GROMACS [ exclusions ]" in output
    assert "topol.top:15" in output


def test_direct_gromacs_resolves_constrainttypes_and_honors_constrnc(
    tmp_path,
):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "constrainttypes",
        """
        [ defaults ]
        1 2 yes 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.3 0.4184
        B B 12.0 0.0 A 0.3 0.4184

        [ constrainttypes ]
        A A 1 0.4
        A B 2 0.4

        [ moleculetype ]
        CONNECTED 1

        [ atoms ]
        1 A 1 MOL A1 1 0.0 12.0
        2 A 1 MOL A2 1 0.0 12.0

        [ constraints ]
        1 2 1

        [ moleculetype ]
        NOT_CONNECTED 1

        [ atoms ]
        1 A 1 MOL A1 1 0.0 12.0
        2 B 1 MOL B2 1 0.0 12.0

        [ constraints ]
        1 2 2

        [ system ]
        constraint types

        [ molecules ]
        CONNECTED 1
        NOT_CONNECTED 1
        """,
        [
            (0.0, 0.0, 0.0),
            (0.4, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            (2.4, 0.0, 0.0),
        ],
        gro_identities=[
            ("MOL", "A1"),
            ("MOL", "A2"),
            ("MOL", "A1"),
            ("MOL", "B2"),
        ],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert "constrain pair number is 2" in result.stdout + result.stderr
    sigma_over_r = 0.3 / 0.4
    expected_lj = 4.0 * 0.1 * (sigma_over_r**12 - sigma_over_r**6)
    assert _extract_mdout_term(mdout_path, "LJ_short") == pytest.approx(
        expected_lj, abs=5.0e-3
    )


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


def test_direct_gromacs_accepts_small_reversible_canonical_box(tmp_path):
    gro_lines = [
        "small reversible box",
        "2",
        _gro_atom(
            1,
            "RESID",
            "FIRST",
            1,
            (0.0, 0.0, 0.0),
            precision=6,
        ),
        _gro_atom(
            1,
            "RESID",
            "SECON",
            2,
            (0.0005, 0.0, 0.0),
            precision=6,
        ),
        "1e-3 2e-3 3e-3",
    ]
    result, _ = _run_custom_gro_case(
        tmp_path,
        "small_reversible_box",
        _strict_gro_topology(),
        gro_lines,
        mdin_overrides={"cutoff": 0.001, "skin": 0.001},
    )
    assert result.returncode == 0, result.stdout + "\n" + result.stderr


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
            "invalid non-canonical GROMACS triclinic box orientation",
        ),
        (
            "2.0 2.5 3.0 1e-20 0.0 0.3 0.0 0.5 0.6",
            "invalid non-canonical GROMACS triclinic box orientation",
        ),
        (
            "1e38 1e38 1e38",
            "box vector component cannot be represented by SPONGE",
        ),
        (
            "1e-47 1.0 1.0",
            "box vector component cannot be represented by SPONGE",
        ),
        (
            "1e-40 1.0 1.0",
            "box vector component converts to a subnormal float",
        ),
    ],
    ids=[
        "non_finite",
        "zero_vector",
        "coplanar",
        "rotated_orientation",
        "tiny_rotated_component",
        "float_overflow",
        "float_underflow",
        "float_subnormal_ftz",
    ],
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


@pytest.mark.parametrize(
    "entry,error_fragment",
    [
        ("[ bonds ]\n1 2 2 0.1 100.0", "[ bonds ] function 2"),
        ("[ angles ]\n1 2 3 2 109.5 100.0", "[ angles ] function 2"),
        (
            "[ dihedrals ]\n1 2 3 4 8 0.0 1.0 1",
            "[ dihedrals ] function 8",
        ),
        ("[ pairs ]\n1 2 2 0.3 0.4184", "[ pairs ] function 2"),
        ("[ settles ]\n1 2 0.1 0.16", "[ settles ] function 2"),
        ("[ constraints ]\n1 2 3 0.1", "[ constraints ] function 3"),
        ("[ cmap ]\n1 2 3 4 5 2", "[ cmap ] function 2"),
    ],
    ids=["bond", "angle", "dihedral", "pair", "settle", "constraint", "cmap"],
)
def test_direct_gromacs_rejects_unsupported_interaction_functions(
    tmp_path, entry, error_fragment
):
    result, _ = _run_direct_topology(
        tmp_path,
        "unsupported_function",
        _minimal_five_atom_topology(entry),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output
    assert "topol.top:" in output


@pytest.mark.parametrize(
    "entry",
    [
        "[ bonds ]\n1 2 1 0.1",
        "[ angles ]\n1 2 3 1 109.5",
        "[ dihedrals ]\n1 2 3 4 9 0.0 1.0",
        "[ pairs ]\n1 2 1 0.3",
        "[ settles ]\n1 1 0.1 0.16 9.9",
        "[ constraints ]\n1 2 1 0.1 0.2 0.3",
        "[ cmap ]\n1 2 3 4 5 1 9.9",
    ],
    ids=["bond", "angle", "dihedral", "pair", "settle", "constraint", "cmap"],
)
def test_direct_gromacs_rejects_invalid_interaction_parameter_counts(
    tmp_path, entry
):
    result, _ = _run_direct_topology(
        tmp_path,
        "invalid_parameter_count",
        _minimal_five_atom_topology(entry),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "invalid" in output
    assert "topol.top:" in output


@pytest.mark.parametrize(
    "entry,type_section,error_fragment",
    [
        ("[ bonds ]\n1 2", "", "[ bonds ]"),
        ("[ pairs ]\n1 2", "", "[ pairs ]"),
        ("[ angles ]\n1 2 3", "", "[ angles ]"),
        ("[ dihedrals ]\n1 2 3 4", "", "[ dihedrals ]"),
        ("[ constraints ]\n1 2", "", "[ constraints ]"),
        ("", "[ cmaptypes ]\nA A A A A 1 2", "[ cmaptypes ]"),
    ],
    ids=["bond", "pair", "angle", "dihedral", "constraint", "cmaptype"],
)
def test_direct_gromacs_short_records_report_source(
    tmp_path, entry, type_section, error_fragment
):
    result, _ = _run_direct_topology(
        tmp_path,
        "short_record",
        _minimal_five_atom_topology(entry, type_section),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output
    assert "topol.top:" in output


@pytest.mark.parametrize(
    "entry,section",
    [
        ("[ angles ]\n0 2 3 1 109.5 100.0", "angles"),
        ("[ dihedrals ]\n1 2 3 6 1 0.0 1.0 1", "dihedrals"),
        ("[ pairs ]\n-1 2 1 0.3 0.4184", "pairs"),
    ],
    ids=["angle", "dihedral", "pair"],
)
def test_direct_gromacs_rejects_out_of_range_interaction_atoms_at_source(
    tmp_path, entry, section
):
    result, _ = _run_direct_topology(
        tmp_path,
        f"out_of_range_{section}",
        _minimal_five_atom_topology(entry),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"atom index out of range in GROMACS [ {section} ]" in output
    assert "topol.top:" in output


def test_direct_gromacs_rejects_cmap_kernel_index_overflow_at_source(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "cmap_index_overflow",
        _minimal_five_atom_topology(
            "", "[ cmaptypes ]\nA A A A A 1 50000 50000"
        ),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "CMAP interpolation table exceeds" in output
    assert "topol.top:" in output


@pytest.mark.parametrize(
    "bond_parameters,error_fragment",
    [
        ("1e38 1.0", "bond distance cannot be represented"),
        (
            "0.2 1.17549435e-38",
            "bond force constant converts to a subnormal float",
        ),
    ],
    ids=["overflow", "subnormal_after_unit_conversion"],
)
def test_direct_gromacs_rejects_unrepresentable_bonded_conversion_at_source(
    tmp_path, bond_parameters, error_fragment
):
    result, _ = _run_direct_topology(
        tmp_path,
        "bonded_conversion_range",
        _minimal_five_atom_topology(f"[ bonds ]\n1 2 1 {bond_parameters}"),
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output
    assert "topol.top:" in output


@pytest.mark.parametrize(
    "combination_rule,parameters",
    [(1, "1e30 1e30"), (2, "1e10 1e10"), (2, "1e-10 1.0")],
    ids=["c6_c12_units", "sigma_epsilon_power", "float_underflow"],
)
def test_direct_gromacs_rejects_unrepresentable_lj_parameters_at_source(
    tmp_path, combination_rule, parameters
):
    result, _ = _run_direct_topology(
        tmp_path,
        "lj_parameter_overflow",
        f"""
        [ defaults ]
        1 {combination_rule} no 1.0 1.0

        [ atomtypes ]
        A 12.0 0.0 A {parameters}

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A 1

        [ system ]
        overflowing LJ parameters

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "Lennard-Jones parameters cannot be represented" in output
    assert "topol.top:" in output


def test_direct_gromacs_does_not_merge_distinct_lj_mixing_inputs(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "distinct_lj_mixing_inputs",
        """
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        POS 12.0 0.0 A 0.3 0.4184
        NEG 12.0 0.0 A 0.3 -0.4184

        [ moleculetype ]
        MIXED 0

        [ atoms ]
        1 POS 1 MOL A1 1 0.0 12.0
        2 NEG 1 MOL A2 2 0.0 12.0

        [ system ]
        distinct LJ mixing inputs

        [ molecules ]
        MIXED 1
        """,
        [(0.0, 0.0, 0.0), (0.5, 0.0, 0.0)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "nonbonded type pair 'POS/NEG'" in output
    assert "Lennard-Jones parameters cannot be represented" in output
    assert "topol.top:" in output


def test_direct_gromacs_rejects_subnormal_scaled_pair_ab_at_source(tmp_path):
    topology = _minimal_five_atom_topology("[ pairs ]\n1 4 1")
    topology = topology.replace(
        "1 2 yes 1.0 1.0", "1 2 yes 1.17549435e-38 1.0"
    ).replace("A A 12.0 0.0 A 0.3 0.4184", "A A 12.0 0.0 A 0.1 0.04184")
    result, _ = _run_direct_topology(
        tmp_path,
        "subnormal_scaled_pair_ab",
        topology,
        [(0.2 * i, 0.0, 0.0) for i in range(5)],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert (
        "GROMACS [ pairs ] interaction Lennard-Jones conversion produces a "
        "subnormal float"
    ) in output
    assert "consistent FTZ behavior" in output
    assert "topol.top:" in output


def test_direct_gromacs_expands_recursive_object_macros(tmp_path):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "object_macros",
        """
        #define BOND_LENGTH 0.2
        #define BOND_PARAMETERS BOND_LENGTH 836.8

        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 1

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1

        [ bonds ]
        1 2 1 BOND_PARAMETERS

        [ system ]
        object macros

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.3, 0.0, 0.0)],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert _extract_mdout_term(mdout_path, "bond") == pytest.approx(1.0)


def test_direct_gromacs_expands_more_than_100000_macro_tokens(tmp_path):
    token_count = 100_001
    large_description = " ".join(["x"] * token_count)
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "large_object_macro",
        f"""
        #define LARGE_DESCRIPTION {large_description}

        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        LARGE_DESCRIPTION

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert mdout_path.exists()


def test_direct_gromacs_preserves_token_after_100000_macro_expansions(
    tmp_path,
):
    repeated_indices = " ".join(["2"] * 100_001)
    result, _ = _run_direct_topology(
        tmp_path,
        "observable_large_object_macro",
        f"""
        #define MANY_EXCLUSIONS 1 {repeated_indices} SENTINEL

        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1

        [ exclusions ]
        MANY_EXCLUSIONS

        [ system ]
        observable large macro

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0)],
        timeout=10,
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "invalid integer field 'SENTINEL'" in output


def test_direct_gromacs_handles_deep_aliases_and_empty_shared_graph(
    tmp_path,
):
    alias_depth = 4096
    base_description = " ".join(["x"] * 4096)
    aliases = [f"#define ALIAS0 BASE"]
    aliases.extend(
        f"#define ALIAS{index} ALIAS{index - 1}"
        for index in range(1, alias_depth)
    )
    empty_graph = ["#define EMPTY"]
    empty_graph.extend(
        f"#define EMPTY{index} "
        + ("EMPTY EMPTY" if index == 0 else f"EMPTY{index - 1} EMPTY{index - 1}")
        for index in range(64)
    )
    definitions = "\n".join(
        [f"#define BASE {base_description}", *aliases, *empty_graph]
    )
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "deep_and_shared_object_macros",
        f"""
        {definitions}

        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        ALIAS{alias_depth - 1} ALIAS{alias_depth - 1} EMPTY63

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
        timeout=10,
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert mdout_path.exists()


@pytest.mark.parametrize(
    "topology,extra_files,error_fragment",
    [
        (
            "#unknown_directive value",
            {},
            "unsupported GROMACS directive '#unknown_directive'",
        ),
        (
            "#define A B\n#define B A\n[ defaults ]\n1 2\nA",
            {},
            "GROMACS macro expansion cycle involving 'A'",
        ),
        (
            '#include "cycle.itp"',
            {"cycle.itp": '#include "topol.top"'},
            "GROMACS topology include cycle detected",
        ),
    ],
    ids=["unknown_directive", "macro_cycle", "include_cycle"],
)
def test_direct_gromacs_rejects_unsafe_preprocessor_input(
    tmp_path, topology, extra_files, error_fragment
):
    result, _ = _run_direct_topology(
        tmp_path,
        "unsafe_preprocessor",
        topology,
        [(0.0, 0.0, 0.0)],
        extra_files=extra_files,
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output


@pytest.mark.parametrize(
    "defaults_line,should_succeed,error_fragment",
    [
        ("1 2", True, ""),
        ("1 3 no 1.0 1.0 12", True, ""),
        ("1 2 maybe", False, "invalid gen-pairs"),
        ("1 2 no -0.1", False, "finite and non-negative"),
        ("1 2 no 1.0 1.0 9", False, "repulsion power 12"),
    ],
    ids=[
        "minimal",
        "combination_rule_3",
        "bad_boolean",
        "bad_fudge",
        "bad_power",
    ],
)
def test_direct_gromacs_validates_defaults(
    tmp_path, defaults_line, should_succeed, error_fragment
):
    topology = f"""
        [ defaults ]
        {defaults_line}

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        defaults

        [ molecules ]
        MOL 1
    """
    result, _ = _run_direct_topology(
        tmp_path,
        "defaults_validation",
        topology,
        [(0.0, 0.0, 0.0)],
    )
    output = result.stdout + "\n" + result.stderr
    if should_succeed:
        assert result.returncode == 0, output
    else:
        assert result.returncode != 0, output
        assert error_fragment in output


def test_direct_gromacs_uses_bonded_type_alias_and_last_type_definition(
    tmp_path,
):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "bonded_type_alias",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        NA BT 12.0 0.0 A 0.0 0.0
        NB BT 12.0 0.0 A 0.0 0.0

        [ bondtypes ]
        BT BT 1 0.3 836.8
        BT BT 1 0.2 836.8

        [ moleculetype ]
        MOL 1

        [ atoms ]
        1 NA 1 MOL A1 1
        2 NB 1 MOL A2 1

        [ bonds ]
        1 2 1

        [ system ]
        bonded type alias

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.3, 0.0, 0.0)],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert _extract_mdout_term(mdout_path, "bond") == pytest.approx(1.0)


def test_direct_gromacs_atoms_use_type_defaults_and_preserve_zero_charge(
    tmp_path,
):
    fallback_result, fallback_mdout = _run_direct_topology(
        tmp_path,
        "atom_defaults",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 1.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1

        [ system ]
        atom defaults

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0)],
    )
    assert fallback_result.returncode == 0, (
        fallback_result.stdout + "\n" + fallback_result.stderr
    )
    assert abs(_extract_mdout_term(fallback_mdout, "PM")) > 1.0

    zero_result, zero_mdout = _run_direct_topology(
        tmp_path,
        "explicit_zero_atom_charge",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 1.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1 0.0 12.0
        2 A 1 MOL A2 1 0.0 12.0

        [ system ]
        explicit zero atom charge

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0)],
    )
    assert zero_result.returncode == 0, (
        zero_result.stdout + "\n" + zero_result.stderr
    )
    assert _extract_mdout_term(zero_mdout, "PM") == pytest.approx(0.0)


def test_direct_gromacs_periodic_improper_uses_reverse_alias_and_last_type(
    tmp_path,
):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "periodic_improper",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        NA TA 12.0 0.0 A 0.0 0.0
        NB TB 12.0 0.0 A 0.0 0.0
        NC TC 12.0 0.0 A 0.0 0.0
        ND TD 12.0 0.0 A 0.0 0.0

        [ dihedraltypes ]
        TD TC TB TA 4 0.0 4.184 1
        TD TC TB TA 4 0.0 8.368 1

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 NA 1 MOL A1 1
        2 NB 1 MOL A2 1
        3 NC 1 MOL A3 1
        4 ND 1 MOL A4 1

        [ dihedrals ]
        1 2 3 4 4

        [ system ]
        periodic improper

        [ molecules ]
        MOL 1
        """,
        [
            (0.0, 0.1, 0.0),
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
            (0.1, 0.0, 0.1),
        ],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert _extract_mdout_term(mdout_path, "dihedral") == pytest.approx(
        2.0, abs=0.02
    )


def test_direct_gromacs_funct9_uses_only_latest_consecutive_type_group(
    tmp_path,
):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "funct9_group",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A A 12.0 0.0 A 0.0 0.0
        B B 12.0 0.0 A 0.0 0.0
        C C 12.0 0.0 A 0.0 0.0
        D D 12.0 0.0 A 0.0 0.0

        [ dihedraltypes ]
        D C B A 9 0.0 4.184 1
        D C B A 9 180.0 8.368 2
        X B C X 9 0.0 41.84 1
        D C B A 9 0.0 16.736 1
        A B C D 9 180.0 4.184 2

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1
        2 B 1 MOL A2 1
        3 C 1 MOL A3 1
        4 D 1 MOL A4 1

        [ dihedrals ]
        1 2 3 4 9

        [ system ]
        funct9 group

        [ molecules ]
        MOL 1
        """,
        [
            (0.0, 0.1, 0.0),
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
            (0.1, 0.0, 0.1),
        ],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert _extract_mdout_term(mdout_path, "dihedral") == pytest.approx(
        6.0, abs=0.02
    )


@pytest.mark.parametrize(
    "funct,parameters,use_type",
    [
        (3, (1.2, -0.7, 0.4, 0.8, -0.2, 0.3), False),
        (3, (1.2, -0.7, 0.4, 0.8, -0.2, 0.3), True),
        (5, (1.1, -0.6, 0.9, 0.2), False),
        (5, (1.1, -0.6, 0.9, 0.2), True),
    ],
    ids=[
        "ryckaert_bellemans_inline",
        "ryckaert_bellemans_type",
        "fourier_inline",
        "fourier_type",
    ],
)
@pytest.mark.parametrize("geometry", ["general", "near_planar"])
def test_direct_gromacs_polynomial_dihedral_energy_and_force_match_oracle(
    tmp_path, funct, parameters, use_type, geometry
):
    parameter_text = " ".join(str(value) for value in parameters)
    type_section = (
        f"[ dihedraltypes ]\nA A A A {funct} {parameter_text}"
        if use_type
        else ""
    )
    interaction_parameters = "" if use_type else f" {parameter_text}"
    topology = f"""
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.0 0.0

        {type_section}

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1
        3 A 1 MOL A3 1
        4 A 1 MOL A4 1

        [ dihedrals ]
        1 2 3 4 {funct}{interaction_parameters}

        [ system ]
        polynomial dihedral mapping

        [ molecules ]
        MOL 1
    """
    coordinates = (
        [
            (0.10, 0.22, 0.03),
            (0.21, 0.08, 0.19),
            (0.43, 0.31, 0.11),
            (0.58, 0.17, 0.42),
        ]
        if geometry == "general"
        else [
            (0.10, 0.20, 0.0),
            (0.10, 0.10, 0.0),
            (0.30, 0.10, 0.0),
            (0.30, 0.00, 0.000001),
        ]
    )
    gro_lines = ["polynomial dihedral", "4"]
    gro_lines.extend(
        _gro_atom(1, "MOL", f"A{index}", index, coordinate, precision=6)
        for index, coordinate in enumerate(coordinates, start=1)
    )
    gro_lines.append("   5.00000   5.00000   5.00000")
    result, case_dir = _run_custom_gro_case(
        tmp_path,
        f"polynomial_dihedral_{funct}_{'type' if use_type else 'inline'}_"
        f"{geometry}",
        topology,
        gro_lines,
        mdin_overrides={"frc": "frc.dat"},
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output

    if funct == 3:
        c0, c1, c2, c3, c4, c5 = parameters
        coefficients = (
            c0 + 0.5 * c2 + 0.375 * c4,
            -c1 - 0.75 * c3 - 0.625 * c5,
            0.5 * c2 + 0.5 * c4,
            -0.25 * c3 - 0.3125 * c5,
            0.125 * c4,
            -0.0625 * c5,
        )
    else:
        coefficients = [0.5 * sum(parameters)]
        coefficients.extend(
            0.5 * value * (1.0 if index % 2 == 1 else -1.0)
            for index, value in enumerate(parameters, start=1)
        )

    def subtract(first, second):
        return tuple(a - b for a, b in zip(first, second))

    def cross(first, second):
        return (
            first[1] * second[2] - first[2] * second[1],
            first[2] * second[0] - first[0] * second[2],
            first[0] * second[1] - first[1] * second[0],
        )

    def dot(first, second):
        return sum(a * b for a, b in zip(first, second))

    def oracle_energy(current_coordinates):
        drij = subtract(current_coordinates[0], current_coordinates[1])
        drkj = subtract(current_coordinates[2], current_coordinates[1])
        drkl = subtract(current_coordinates[2], current_coordinates[3])
        normal_1 = cross(drij, drkj)
        normal_2 = cross(drkl, drkj)
        cosine = dot(normal_1, normal_2) / math.sqrt(
            dot(normal_1, normal_1) * dot(normal_2, normal_2)
        )
        raw_phi = math.copysign(
            math.acos(max(-1.0, min(1.0, cosine))),
            dot(cross(normal_2, normal_1), drkj),
        )
        phi = math.pi - raw_phi
        energy_kj = sum(
            coefficient * math.cos(multiplicity * phi)
            for multiplicity, coefficient in enumerate(coefficients)
        )
        return energy_kj / 4.184

    expected_energy = oracle_energy(coordinates)
    assert _extract_mdout_term(case_dir / "mdout.txt", "dihedral") == (
        pytest.approx(expected_energy, abs=0.02)
    )

    observed_forces = _read_float_vectors(case_dir / "frc.dat")
    assert len(observed_forces) == 12
    epsilon_nm = 1.0e-7 if geometry == "near_planar" else 1.0e-6
    for atom_index in range(4):
        for axis in range(3):
            plus = [list(coordinate) for coordinate in coordinates]
            minus = [list(coordinate) for coordinate in coordinates]
            plus[atom_index][axis] += epsilon_nm
            minus[atom_index][axis] -= epsilon_nm
            numerical_force = -(oracle_energy(plus) - oracle_energy(minus)) / (
                20.0 * epsilon_nm
            )
            assert observed_forces[3 * atom_index + axis] == pytest.approx(
                numerical_force, rel=3.0e-3, abs=3.0e-3
            )


def test_direct_gromacs_angle_type_index_is_reverse_symmetric_last_wins(
    tmp_path,
):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "reverse_angle_type",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        NA TA 12.0 0.0 A 0.0 0.0
        NB TB 12.0 0.0 A 0.0 0.0
        NC TC 12.0 0.0 A 0.0 0.0

        [ angletypes ]
        TC TB TA 1 90.0 836.8
        TC TB TA 1 60.0 836.8

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 NA 1 MOL A1 1
        2 NB 1 MOL A2 1
        3 NC 1 MOL A3 1

        [ angles ]
        1 2 3 1

        [ system ]
        reverse angle type

        [ molecules ]
        MOL 1
        """,
        [(0.1, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.1, 0.0)],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    expected = 0.5 * (836.8 / 4.184) * (math.pi / 6.0) ** 2
    assert _extract_mdout_term(mdout_path, "urey_bradley") == pytest.approx(
        expected, abs=0.02
    )


def test_direct_gromacs_connection_bond_only_generates_exclusion(tmp_path):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "connection_bond",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A A 12.0 0.0 A 0.3 0.4184

        [ bondtypes ]
        A A 5

        [ moleculetype ]
        MOL 1

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1

        [ bonds ]
        1 2 5

        [ system ]
        connection bond

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0)],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert _extract_mdout_term(mdout_path, "LJ_short") == pytest.approx(0.0)
    assert "BOND IS NOT INITIALIZED" in result.stdout + result.stderr


def test_direct_gromacs_funct6_bond_has_energy_but_no_exclusion(tmp_path):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "harmonic_no_exclusion_bond",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A A 12.0 0.0 A 0.3 0.4184

        [ bondtypes ]
        A A 6 0.3 836.8

        [ moleculetype ]
        MOL 1

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1

        [ bonds ]
        1 2 6

        [ system ]
        harmonic no-exclusion bond

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.4, 0.0, 0.0)],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    assert _extract_mdout_term(mdout_path, "bond") == pytest.approx(1.0)
    sigma_over_r = 0.3 / 0.4
    expected_lj = 0.4 * (sigma_over_r**12 - sigma_over_r**6)
    assert _extract_mdout_term(mdout_path, "LJ_short") == pytest.approx(
        expected_lj, abs=0.011
    )


def test_direct_gromacs_settles_does_not_generate_nrexcl_exclusions(tmp_path):
    result, mdout_path = _run_direct_topology(
        tmp_path,
        "settles_no_exclusions",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A A 12.0 0.0 A 0.08 0.4184

        [ moleculetype ]
        MOL 1

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1
        3 A 1 MOL A3 1

        [ settles ]
        1 1 0.1 0.141421356

        [ system ]
        settles no exclusions

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.0, 0.1, 0.0)],
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr
    distances = (0.1, 0.1, math.sqrt(0.02))
    expected = sum(
        0.4 * ((0.08 / distance) ** 12 - (0.08 / distance) ** 6)
        for distance in distances
    )
    assert _extract_mdout_term(mdout_path, "LJ_short") == pytest.approx(
        expected, abs=0.011
    )


def test_direct_gromacs_command_defines_and_angle_include_are_expanded(
    tmp_path,
):
    result, _ = _run_direct_topology(
        tmp_path,
        "command_defines",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A ATOM_SIGMA ATOM_EPSILON

        #include <molecule.itp>

        [ system ]
        command defines

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
        extra_files={
            "molecule.itp": """
            [ moleculetype ]
            MOL 0

            [ atoms ]
            1 A 1 MOL A1 1
            """,
        },
        gromacs_define="-DATOM_SIGMA=-0.3,ATOM_EPSILON=0.4184",
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_direct_gromacs_multiline_directives_and_inactive_unknown_are_valid(
    tmp_path,
):
    result, _ = _run_direct_topology(
        tmp_path,
        "multiline_directives",
        r"""
        [ defaults ]
        1 2

        #define ATOM_PARAMETERS 12.0 0.0 A \
            0.0 0.0

        [ atomtypes ]
        A ATOM_PARAMETERS

        #ifdef NEVER_DEFINED
        #unknown_directive ignored
        #include malformed
        #define 123 malformed
        #endif

        #ifndef NEVER_DEFINED
        #include \
            <molecule.itp>
        #endif

        [ system ]
        multiline directives

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
        extra_files={
            "molecule.itp": """
            [ moleculetype ]
            MOL 0

            [ atoms ]
            1 A 1 MOL A1 1
            """,
        },
    )

    assert result.returncode == 0, result.stdout + "\n" + result.stderr


def test_direct_gromacs_invalid_numeric_field_reports_include_source(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "invalid_numeric_source",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        #include "bad.itp"

        [ system ]
        invalid numeric source

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
        extra_files={
            "bad.itp": """
            [ moleculetype ]
            MOL 0

            [ atoms ]
            not_an_integer A 1 MOL A1 1
            """,
        },
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "invalid integer field 'not_an_integer'" in output
    assert "bad.itp:5" in output


@pytest.mark.parametrize(
    "value,error_fragment",
    [
        (
            "1e-46",
            "nonzero real field '1e-46' cannot be represented",
        ),
        (
            "1.40129846e-45",
            "real field '1.40129846e-45' is a subnormal float",
        ),
        (
            "3.5e38",
            "outside the finite float range used by SPONGE",
        ),
        ("1e-999", "real field out of range in GROMACS topology"),
    ],
    ids=[
        "float_underflow",
        "float_subnormal",
        "float_overflow",
        "double_underflow",
    ],
)
def test_direct_gromacs_rejects_unrepresentable_real_fields_with_source(
    tmp_path, value, error_fragment
):
    result, _ = _run_direct_topology(
        tmp_path,
        f"unrepresentable_real_{value}",
        """
        [ defaults ]
        1 2

        #include "bad-real.itp"

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        unrepresentable real field

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
        extra_files={
            "bad-real.itp": f"""
            [ atomtypes ]
            A A 12.0 0.0 A {value} 0.4184
            """,
        },
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert error_fragment in output
    assert "bad-real.itp:2" in output


@pytest.mark.parametrize(
    "ptype", ["S", "V", "D"], ids=["shell", "virtual-site", "dummy"]
)
def test_direct_gromacs_rejects_unimplemented_particle_types_with_source(
    tmp_path, ptype
):
    result, _ = _run_direct_topology(
        tmp_path,
        f"unsupported_ptype_{ptype}",
        f"""
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 {ptype} 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        unsupported ptype

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert f"unsupported instantiated GROMACS particle type '{ptype}'" in output
    assert "topol.top:11" in output


def test_direct_gromacs_rejects_fractional_dihedral_multiplicity_with_source(
    tmp_path,
):
    result, _ = _run_direct_topology(
        tmp_path,
        "fractional_multiplicity",
        _minimal_five_atom_topology("[ dihedrals ]\n1 2 3 4 4 0.0 4.184 1.5"),
        [
            (0.0, 0.1, 0.0),
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
            (0.1, 0.0, 0.1),
            (0.2, 0.0, 0.1),
        ],
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "invalid integer field '1.5'" in output
    assert "topol.top:" in output


def test_direct_gromacs_rejects_too_negative_dihedral_multiplicity(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "too_negative_multiplicity",
        _minimal_five_atom_topology("[ dihedrals ]\n1 2 3 4 1 35.0 4.184 -100"),
        [
            (0.0, 0.1, 0.0),
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
            (0.1, 0.0, 0.1),
            (0.2, 0.0, 0.1),
        ],
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "multiplicity must not be less than -99" in output
    assert "topol.top:" in output


def test_direct_gromacs_negative_multiplicity_near_planar_force_matches_oracle(
    tmp_path,
):
    phase_degrees = 35.0
    multiplicity = -3
    coordinates = [
        (0.10, 0.20, 0.0),
        (0.10, 0.10, 0.0),
        (0.30, 0.10, 0.0),
        (0.30, 0.00, 0.000001),
    ]
    topology = f"""
        [ defaults ]
        1 2 no 1.0 1.0

        [ atomtypes ]
        A A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1
        3 A 1 MOL A3 1
        4 A 1 MOL A4 1

        [ dihedrals ]
        1 2 3 4 1 {phase_degrees} 4.184 {multiplicity}

        [ system ]
        negative multiplicity

        [ molecules ]
        MOL 1
    """
    gro_lines = ["negative multiplicity", "4"]
    gro_lines.extend(
        _gro_atom(1, "MOL", f"A{index}", index, coordinate, precision=6)
        for index, coordinate in enumerate(coordinates, start=1)
    )
    gro_lines.append("   5.00000   5.00000   5.00000")
    result, case_dir = _run_custom_gro_case(
        tmp_path,
        "negative_multiplicity_near_planar",
        topology,
        gro_lines,
        mdin_overrides={"frc": "frc.dat"},
    )
    output = result.stdout + "\n" + result.stderr
    assert result.returncode == 0, output

    def subtract(first, second):
        return tuple(a - b for a, b in zip(first, second))

    def cross(first, second):
        return (
            first[1] * second[2] - first[2] * second[1],
            first[2] * second[0] - first[0] * second[2],
            first[0] * second[1] - first[1] * second[0],
        )

    def dot(first, second):
        return sum(a * b for a, b in zip(first, second))

    def oracle_energy(current_coordinates):
        drij = subtract(current_coordinates[0], current_coordinates[1])
        drkj = subtract(current_coordinates[2], current_coordinates[1])
        drkl = subtract(current_coordinates[2], current_coordinates[3])
        normal_1 = cross(drij, drkj)
        normal_2 = cross(drkl, drkj)
        cosine = dot(normal_1, normal_2) / math.sqrt(
            dot(normal_1, normal_1) * dot(normal_2, normal_2)
        )
        raw_phi = math.copysign(
            math.acos(max(-1.0, min(1.0, cosine))),
            dot(cross(normal_2, normal_1), drkj),
        )
        phi = math.pi - raw_phi
        phase = math.radians(phase_degrees)
        return 1.0 + math.cos(multiplicity * phi - phase)

    assert _extract_mdout_term(case_dir / "mdout.txt", "dihedral") == (
        pytest.approx(oracle_energy(coordinates), abs=0.02)
    )
    observed_forces = _read_float_vectors(case_dir / "frc.dat")
    epsilon_nm = 1.0e-7
    for atom_index in range(4):
        for axis in range(3):
            plus = [list(coordinate) for coordinate in coordinates]
            minus = [list(coordinate) for coordinate in coordinates]
            plus[atom_index][axis] += epsilon_nm
            minus[atom_index][axis] -= epsilon_nm
            numerical_force = -(oracle_energy(plus) - oracle_energy(minus)) / (
                20.0 * epsilon_nm
            )
            assert observed_forces[3 * atom_index + axis] == pytest.approx(
                numerical_force, rel=4.0e-3, abs=4.0e-3
            )


@pytest.mark.parametrize(
    "interaction,defaults,error_fragment",
    [
        ("[ bonds ]\n1 2 1", "1 2 yes 1.0 1.0", "bond type"),
        (
            "[ constraints ]\n1 2 1",
            "1 2 yes 1.0 1.0",
            "constraint distance",
        ),
        ("[ angles ]\n1 2 3 1", "1 2 yes 1.0 1.0", "angle type"),
        (
            "[ dihedrals ]\n1 2 3 4 1",
            "1 2 yes 1.0 1.0",
            "dihedral type",
        ),
        ("[ pairs ]\n1 2 1", "1 2 no 1.0 1.0", "pair interaction"),
        ("[ cmap ]\n1 2 3 4 5 1", "1 2 yes 1.0 1.0", "CMAP type"),
    ],
    ids=["bond", "constraint", "angle", "dihedral", "pair", "cmap"],
)
def test_direct_gromacs_unresolved_interaction_type_reports_source(
    tmp_path, interaction, defaults, error_fragment
):
    topology = _minimal_five_atom_topology(interaction).replace(
        "1 2 yes 1.0 1.0", defaults, 1
    )
    result, _ = _run_direct_topology(
        tmp_path,
        f"unresolved_{error_fragment.replace(' ', '_').lower()}",
        topology,
        [
            (0.0, 0.1, 0.0),
            (0.0, 0.0, 0.0),
            (0.1, 0.0, 0.0),
            (0.1, 0.0, 0.1),
            (0.2, 0.0, 0.1),
        ],
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "failed to " in output
    assert error_fragment in output
    assert "topol.top:" in output


def test_direct_gromacs_undefined_atom_and_molecule_types_report_source(
    tmp_path,
):
    undefined_atom, _ = _run_direct_topology(
        tmp_path,
        "undefined_atom_type",
        """
        [ defaults ]
        1 2

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 MISSING 1 MOL A1 1

        [ system ]
        undefined atom type

        [ molecules ]
        MOL 1
        """,
        [(0.0, 0.0, 0.0)],
    )
    atom_output = undefined_atom.stdout + "\n" + undefined_atom.stderr
    assert undefined_atom.returncode != 0, atom_output
    assert "undefined GROMACS atom type 'MISSING'" in atom_output
    assert "topol.top:8" in atom_output

    undefined_molecule, _ = _run_direct_topology(
        tmp_path,
        "undefined_molecule_type",
        """
        [ defaults ]
        1 2

        [ system ]
        undefined molecule type

        [ molecules ]
        MISSING 1
        """,
        [(0.0, 0.0, 0.0)],
    )
    molecule_output = (
        undefined_molecule.stdout + "\n" + undefined_molecule.stderr
    )
    assert undefined_molecule.returncode != 0, molecule_output
    assert "molecule 'MISSING' referenced" in molecule_output
    assert "topol.top:8" in molecule_output


def test_direct_gromacs_rejects_huge_molecule_count_before_expansion(tmp_path):
    result, _ = _run_direct_topology(
        tmp_path,
        "huge_molecule_count",
        """
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        huge molecule count

        [ molecules ]
        MOL 2147483647
        """,
        [(0.0, 0.0, 0.0)],
        timeout=10,
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert (
        "GROMACS gro atom count does not match the checked expanded topology "
        "atom count (gro 1, topology 2147483647)"
    ) in output


def test_direct_gromacs_rejects_matched_count_above_coordinate_index_range(
    tmp_path,
):
    atom_count = 2_147_483_647 // 3 + 1
    topology = f"""
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        coordinate index range

        [ molecules ]
        MOL {atom_count}
    """
    result, _ = _run_custom_gro_case(
        tmp_path,
        "coordinate_index_range",
        topology,
        ["coordinate index range", str(atom_count)],
        timeout=10,
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert (
        f"checked expanded GROMACS atom count {atom_count} cannot safely "
        "represent all 3 * atom count coordinate values and kernel indices"
    ) in output


def test_direct_gromacs_reads_truncated_gro_before_large_instantiation(
    tmp_path,
):
    atom_count = 100_000_000
    topology = f"""
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1

        [ system ]
        truncated large gro

        [ molecules ]
        MOL {atom_count}
    """
    result, _ = _run_custom_gro_case(
        tmp_path,
        "truncated_large_gro",
        topology,
        ["truncated large gro", str(atom_count)],
        timeout=10,
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert "invalid atom record in GROMACS gro file" in output
    assert "conf.gro:3" in output


def test_direct_gromacs_rejects_expanded_interactions_before_copying(
    tmp_path,
):
    molecule_count = 268_435_456
    atom_count = 2 * molecule_count
    topology = f"""
        [ defaults ]
        1 2

        [ atomtypes ]
        A 12.0 0.0 A 0.0 0.0

        [ moleculetype ]
        MOL 0

        [ atoms ]
        1 A 1 MOL A1 1
        2 A 1 MOL A2 1

        [ bonds ]
        1 2 1 0.1 100.0
        1 2 1 0.1 100.0
        1 2 1 0.1 100.0
        1 2 1 0.1 100.0
        1 2 1 0.1 100.0
        1 2 1 0.1 100.0
        1 2 1 0.1 100.0
        1 2 1 0.1 100.0

        [ system ]
        expanded interaction range

        [ molecules ]
        MOL {molecule_count}
    """
    result, _ = _run_custom_gro_case(
        tmp_path,
        "expanded_interaction_range",
        topology,
        ["expanded interaction range", str(atom_count)],
        timeout=10,
    )

    output = result.stdout + "\n" + result.stderr
    assert result.returncode != 0, output
    assert (
        "expanded GROMACS bond count exceeds the supported kernel or host "
        "container range"
    ) in output
