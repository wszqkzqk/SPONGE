"""RI analytic-gradient regression against SPONGE energy finite differences."""

import numpy as np
import pytest

from benchmarks.comparison.tests.pyscf.tests.utils import (
    HARTREE_TO_KCAL_MOL,
    prepare_output_case,
    run_sponge_scf_energy_ha,
)

BOHR_PER_ANGSTROM = 1.8897259886
# The AO/density matrices and emitted force are float32.  A 0.008 Angstrom
# displacement keeps the energy signal above that quantization noise while
# remaining in the central-difference regime; 0.005 Ha/Bohr is the existing
# end-to-end error budget for those single-precision production paths.
FD_STEP_ANGSTROM = 0.008
# The internal SCF energy is double precision, but the SCF density/Fock path is
# float32 and its converged energy varies at the few-microhartree level.  Step
# scans at 0.008, 0.012, 0.016, and 0.020 Angstrom put the full Cartesian
# derivative on its central-difference plateau at 0.016 Angstrom without
# measurable truncation drift.
FULL_FD_STEP_ANGSTROM = 0.016
FD_PLATEAU_STEPS_ANGSTROM = (0.012, 0.016, 0.020)
RI_GRAD_TOL_HA_BOHR = 0.005
RI_MODE_TOL_HA_BOHR = 0.002
DFT_FULL_GRAD_TOL_HA_BOHR = 0.001
DFT_FULL_MODE_TOL_HA_BOHR = 0.0005
DFT_ANALYTIC_TRANSLATION_TOL_HA_BOHR = 2.0e-6
DFT_FD_TRANSLATION_TOL_HA_BOHR = 0.0005
DFT_FD_PLATEAU_SPREAD_TOL_HA_BOHR = 0.0008
PYSCF_DFT_GRAD_TOL_HA_BOHR = 0.005

# Together these cover LDA and GGA XC gradients, exact-exchange fractions 1,
# 0, and 0.25, and separate alpha/beta GGA channels in an open-shell
# calculation.  The DFT cases also exercise the analytical zero-density grid
# limit because the atom-centred grid extends beyond screened AO support.
RI_GRAD_CASES = [
    ("h2", "HF/def2-svp", True),
    ("h2", "LDA/def2-svp", True),
    ("h2", "PBE/def2-svp", True),
    ("h2", "PBE0/def2-svp", True),
    ("no_doublet", "HF/def2-svp", False),
    ("no_doublet", "PBE/def2-svp", False),
]

# A generally oriented heteronuclear bond makes every Cartesian component
# non-zero.  NO+ is a closed-shell RKS case; neutral NO is the corresponding
# open-shell UKS case.  Together with both LDA and PBE and both RI execution
# modes, this matrix exercises every atom/axis of the atom-centred grid
# response instead of relying on the symmetry-reduced z component above.
DFT_FULL_GRAD_CASES = [
    ("rks_lda", "LDA/def2-svp", True),
    ("rks_pbe", "PBE/def2-svp", True),
    ("uks_lda", "LDA/def2-svp", False),
    ("uks_pbe", "PBE/def2-svp", False),
]


def _read_coordinate(path):
    lines = path.read_text().splitlines()
    natom = int(lines[0].split()[0])
    coordinates = np.array(
        [[float(value) for value in lines[i + 1].split()[:3]] for i in range(natom)]
    )
    return lines, coordinates


def _write_coordinate(path, original_lines, coordinates):
    natom = len(coordinates)
    lines = [original_lines[0]]
    lines.extend(" ".join(f"{value:.10f}" for value in xyz) for xyz in coordinates)
    lines.extend(original_lines[natom + 1 :])
    path.write_text("\n".join(lines) + "\n")


def _run_energy(sponge_dir, model_chemistry, restricted, mode):
    return run_sponge_scf_energy_ha(
        sponge_dir=sponge_dir,
        model_chemistry=model_chemistry,
        restricted=restricted,
        extra_sponge_args=[
            "-qc_density_fit",
            "1",
            "-qc_density_fitting_mode",
            mode,
        ],
        use_internal_scf_energy=True,
    )


def _analytic_and_fd_bond_gradient(
    sponge_dir, model_chemistry, restricted, mode
):
    coordinate_path = sponge_dir / "coordinate.txt"
    original_lines, coordinates = _read_coordinate(coordinate_path)
    natom = len(coordinates)
    assert natom == 2

    force_name = f"ri_gradient_{mode}.bin"
    mdin_path = sponge_dir / "mdin.txt"
    mdin_path.write_text(mdin_path.read_text() + f"\nfrc = {force_name}\n")

    try:
        _run_energy(sponge_dir, model_chemistry, restricted, mode)
        force = np.fromfile(sponge_dir / force_name, dtype=np.float32)
        assert force.size == natom * 3
        force = force.reshape(natom, 3)
        analytic_full = -force / (HARTREE_TO_KCAL_MOL * BOHR_PER_ANGSTROM)
        assert float(np.max(np.abs(analytic_full[:, :2]))) < 5.0e-4
        analytic = analytic_full[:, 2]

        finite_difference = np.zeros(natom)
        for atom in range(natom):
            energies = []
            for sign in (1.0, -1.0):
                displaced = coordinates.copy()
                displaced[atom, 2] += sign * FD_STEP_ANGSTROM
                _write_coordinate(coordinate_path, original_lines, displaced)
                energies.append(
                    _run_energy(sponge_dir, model_chemistry, restricted, mode)
                )
            finite_difference[atom] = (energies[0] - energies[1]) / (
                2.0 * FD_STEP_ANGSTROM * BOHR_PER_ANGSTROM
            )
    finally:
        _write_coordinate(coordinate_path, original_lines, coordinates)

    return analytic, finite_difference


def _configure_general_no_geometry(sponge_dir, restricted):
    coordinate_path = sponge_dir / "coordinate.txt"
    coordinate_lines, _coordinates = _read_coordinate(coordinate_path)
    # Every Cartesian component is non-zero.  Neutral NO is nearly stationary
    # at 1.15 Angstrom, where a microhartree SCF energy variation overwhelms
    # its tiny finite-difference signal, so the UKS oracle deliberately uses a
    # non-equilibrium 1.30 Angstrom bond.  NO+ remains at 1.15 Angstrom.
    bond_length = 1.15 if restricted else 1.30
    origin = np.array([-0.2, 0.17, -0.11])
    bond_direction = np.array([2.0, 3.0, 6.0]) / 7.0
    coordinates = np.array([origin, origin + bond_length * bond_direction])
    _write_coordinate(coordinate_path, coordinate_lines, coordinates)

    type_path = sponge_dir / "qc_type.txt"
    type_lines = type_path.read_text().splitlines()
    header = type_lines[0].split()
    assert int(header[0]) == 2
    header[1:] = ["1", "1"] if restricted else ["0", "2"]
    type_lines[0] = " ".join(header)
    type_path.write_text("\n".join(type_lines) + "\n")
    return coordinates


def _read_analytic_full_gradient(
    sponge_dir, model_chemistry, restricted, mode, force_name
):
    natom = len(_read_coordinate(sponge_dir / "coordinate.txt")[1])
    mdin_path = sponge_dir / "mdin.txt"
    mdin_path.write_text(mdin_path.read_text() + f"\nfrc = {force_name}\n")
    _run_energy(sponge_dir, model_chemistry, restricted, mode)
    force = np.fromfile(sponge_dir / force_name, dtype=np.float32)
    assert force.size == natom * 3
    return -force.reshape(natom, 3) / (
        HARTREE_TO_KCAL_MOL * BOHR_PER_ANGSTROM
    )


def _finite_difference_components(
    sponge_dir,
    model_chemistry,
    restricted,
    mode,
    steps_angstrom,
    components,
):
    coordinate_path = sponge_dir / "coordinate.txt"
    original_lines, coordinates = _read_coordinate(coordinate_path)
    finite_differences = {
        step: np.full_like(coordinates, np.nan, dtype=float)
        for step in steps_angstrom
    }

    try:
        for step in steps_angstrom:
            for atom, axis in components:
                energies = []
                for sign in (1.0, -1.0):
                    displaced = coordinates.copy()
                    displaced[atom, axis] += sign * step
                    _write_coordinate(coordinate_path, original_lines, displaced)
                    energies.append(
                        _run_energy(
                            sponge_dir, model_chemistry, restricted, mode
                        )
                    )
                finite_differences[step][atom, axis] = (
                    energies[0] - energies[1]
                ) / (2.0 * step * BOHR_PER_ANGSTROM)
    finally:
        _write_coordinate(coordinate_path, original_lines, coordinates)

    return finite_differences


def _analytic_and_fd_full_gradient(
    sponge_dir, model_chemistry, restricted, mode
):
    coordinates = _read_coordinate(sponge_dir / "coordinate.txt")[1]
    analytic = _read_analytic_full_gradient(
        sponge_dir,
        model_chemistry,
        restricted,
        mode,
        f"dft_full_gradient_{mode}.bin",
    )
    components = [
        (atom, axis)
        for atom in range(len(coordinates))
        for axis in range(3)
    ]
    finite_difference = _finite_difference_components(
        sponge_dir,
        model_chemistry,
        restricted,
        mode,
        (FULL_FD_STEP_ANGSTROM,),
        components,
    )[FULL_FD_STEP_ANGSTROM]
    return analytic, finite_difference


@pytest.mark.parametrize(
    "case_name,model_chemistry,restricted",
    RI_GRAD_CASES,
    ids=[
        f"{case}_{model.replace('/', '_')}_{'R' if restricted else 'U'}"
        for case, model, restricted in RI_GRAD_CASES
    ],
)
def test_ri_analytic_gradient_matches_finite_difference(
    case_name, model_chemistry, restricted, statics_path, outputs_path
):
    mode_gradients = {}
    for mode in ("stored", "direct"):
        _run_dir, sponge_dir = prepare_output_case(
            statics_path=statics_path,
            outputs_path=outputs_path,
            case_name=case_name,
            run_tag=(
                f"ri_gradient_{case_name}_{model_chemistry.replace('/', '_')}_"
                f"{'r' if restricted else 'u'}_{mode}"
            ),
        )
        analytic, finite_difference = _analytic_and_fd_bond_gradient(
            sponge_dir, model_chemistry, restricted, mode
        )
        error = float(np.max(np.abs(analytic - finite_difference)))
        print(
            f"{case_name} {model_chemistry} {mode}: "
            f"analytic={analytic}, fd={finite_difference}, max_error={error:.8f}"
        )
        assert np.all(np.isfinite(analytic))
        assert error < RI_GRAD_TOL_HA_BOHR
        assert abs(float(np.sum(analytic))) < 5.0e-4
        mode_gradients[mode] = analytic

    mode_error = float(
        np.max(np.abs(mode_gradients["stored"] - mode_gradients["direct"]))
    )
    assert mode_error < RI_MODE_TOL_HA_BOHR


@pytest.mark.parametrize(
    "case_id,model_chemistry,restricted",
    DFT_FULL_GRAD_CASES,
    ids=[case[0] for case in DFT_FULL_GRAD_CASES],
)
def test_dft_heteronuclear_full_gradient_matrix(
    case_id,
    model_chemistry,
    restricted,
    statics_path,
    outputs_path,
):
    mode_gradients = {}
    for mode in ("stored", "direct"):
        _run_dir, sponge_dir = prepare_output_case(
            statics_path=statics_path,
            outputs_path=outputs_path,
            case_name="no_doublet",
            run_tag=f"dft_full_gradient_{case_id}_{mode}",
        )
        _configure_general_no_geometry(sponge_dir, restricted)
        analytic, finite_difference = _analytic_and_fd_full_gradient(
            sponge_dir, model_chemistry, restricted, mode
        )
        error = float(np.max(np.abs(analytic - finite_difference)))
        analytic_translation = np.sum(analytic, axis=0)
        fd_translation = np.sum(finite_difference, axis=0)
        print(
            f"{case_id} {mode}: analytic={analytic}, "
            f"fd={finite_difference}, max_error={error:.8f}, "
            f"analytic_translation={analytic_translation}, "
            f"fd_translation={fd_translation}"
        )
        assert np.all(np.isfinite(analytic))
        assert np.all(np.isfinite(finite_difference))
        assert error < DFT_FULL_GRAD_TOL_HA_BOHR
        assert float(np.max(np.abs(analytic_translation))) < (
            DFT_ANALYTIC_TRANSLATION_TOL_HA_BOHR
        )
        assert float(np.max(np.abs(fd_translation))) < (
            DFT_FD_TRANSLATION_TOL_HA_BOHR
        )
        mode_gradients[mode] = analytic

    mode_error = float(
        np.max(np.abs(mode_gradients["stored"] - mode_gradients["direct"]))
    )
    assert mode_error < DFT_FULL_MODE_TOL_HA_BOHR


def test_rks_lda_finite_difference_step_plateau(statics_path, outputs_path):
    _run_dir, sponge_dir = prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="no_doublet",
        run_tag="dft_full_gradient_rks_lda_fd_plateau",
    )
    _configure_general_no_geometry(sponge_dir, restricted=True)
    analytic = _read_analytic_full_gradient(
        sponge_dir,
        "LDA/def2-svp",
        True,
        "stored",
        "dft_full_gradient_fd_plateau.bin",
    )
    finite_differences = _finite_difference_components(
        sponge_dir,
        "LDA/def2-svp",
        True,
        "stored",
        FD_PLATEAU_STEPS_ANGSTROM,
        ((0, 0),),
    )
    derivatives = np.array(
        [finite_differences[step][0, 0] for step in FD_PLATEAU_STEPS_ANGSTROM]
    )
    errors = np.abs(derivatives - analytic[0, 0])
    spread = float(np.ptp(derivatives))
    print(
        "RKS/LDA stored atom-0 x FD plateau: "
        f"steps={FD_PLATEAU_STEPS_ANGSTROM}, derivatives={derivatives}, "
        f"analytic={analytic[0, 0]:.8f}, errors={errors}, "
        f"spread={spread:.8f}"
    )
    assert np.all(np.isfinite(derivatives))
    assert float(np.max(errors)) < DFT_FULL_GRAD_TOL_HA_BOHR
    assert spread < DFT_FD_PLATEAU_SPREAD_TOL_HA_BOHR


def test_rks_lda_gradient_matches_pyscf_analytic(statics_path, outputs_path):
    pytest.importorskip(
        "pyscf", reason="the independent analytical oracle requires PySCF"
    )
    from benchmarks.comparison.tests.pyscf.tests.utils import (
        _build_pyscf_method,
    )

    _run_dir, sponge_dir = prepare_output_case(
        statics_path=statics_path,
        outputs_path=outputs_path,
        case_name="no_doublet",
        run_tag="dft_full_gradient_rks_lda_pyscf",
    )
    coordinates = _configure_general_no_geometry(sponge_dir, restricted=True)
    sponge_gradient = _read_analytic_full_gradient(
        sponge_dir,
        "LDA/def2-svp",
        True,
        "stored",
        "dft_full_gradient_pyscf.bin",
    )

    atoms = [
        (symbol, *coordinates[atom])
        for atom, symbol in enumerate(("N", "O"))
    ]
    _mol, mean_field = _build_pyscf_method(
        atoms,
        basis_name="def2-svp",
        charge=1,
        multiplicity=1,
        method_name="LDA",
        restricted=True,
    )
    mean_field.conv_tol = 1.0e-10
    mean_field.max_cycle = 200
    mean_field.grids.level = 4
    energy = mean_field.kernel()
    assert mean_field.converged
    assert np.isfinite(energy)
    gradient_method = mean_field.nuc_grad_method()
    gradient_method.grid_response = True
    pyscf_gradient = np.asarray(gradient_method.kernel())

    error = float(np.max(np.abs(sponge_gradient - pyscf_gradient)))
    print(
        f"SPONGE RKS/LDA gradient={sponge_gradient}, "
        f"PySCF analytic gradient={pyscf_gradient}, max_error={error:.8f}"
    )
    assert pyscf_gradient.shape == sponge_gradient.shape
    assert np.all(np.isfinite(pyscf_gradient))
    assert error < PYSCF_DFT_GRAD_TOL_HA_BOHR


def test_ri_fe_l6_auxiliary_gradient_is_finite_and_translationally_invariant(
    statics_path, outputs_path
):
    energies = {}
    gradients = {}
    for mode in ("stored", "direct"):
        _run_dir, sponge_dir = prepare_output_case(
            statics_path=statics_path,
            outputs_path=outputs_path,
            case_name="fe_quintet",
            run_tag=f"ri_gradient_fe_quintet_l6_{mode}",
        )
        force_name = f"ri_gradient_fe_l6_{mode}.bin"
        mdin_path = sponge_dir / "mdin.txt"
        mdin_path.write_text(mdin_path.read_text() + f"\nfrc = {force_name}\n")
        energies[mode] = _run_energy(
            sponge_dir, "HF/6-31++g", restricted=False, mode=mode
        )
        force = np.fromfile(sponge_dir / force_name, dtype=np.float32)
        assert force.size == 3
        gradient = -force / (HARTREE_TO_KCAL_MOL * BOHR_PER_ANGSTROM)
        assert np.all(np.isfinite(gradient))
        assert float(np.max(np.abs(gradient))) < 5.0e-4
        gradients[mode] = gradient

    assert abs(energies["stored"] - energies["direct"]) < 5.0e-4
    assert float(np.max(np.abs(gradients["stored"] - gradients["direct"]))) < (
        5.0e-4
    )
