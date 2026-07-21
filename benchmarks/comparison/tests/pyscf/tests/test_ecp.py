"""ECP (Effective Core Potential) validation test.

Compares SPONGE ECP energy against PySCF reference for systems with
4th-period elements using def2-SVP + def2-ECP.

Test strategy:
- Use KH molecule (potassium hydride) as a minimal ECP test case
- K has 10 core electrons replaced by def2-ECP
- Compare total SCF energy against pre-computed PySCF reference
"""

import shutil
import tempfile
from pathlib import Path

import pytest

from benchmarks.comparison.tests.pyscf.tests.utils import (
    HARTREE_TO_KCAL_MOL,
    run_sponge_scf_energy_ha,
)

ECP_ENERGY_TOL_HA = 5.0e-4


def _make_sponge_case(tmpdir, atoms, charge, multiplicity):
    """Create SPONGE input files for an ECP test case."""
    sponge_dir = Path(tmpdir)
    natm = len(atoms)

    with open(sponge_dir / "qc_type.txt", "w") as f:
        f.write(f"{natm} {charge} {multiplicity}\n")
        for i, (sym, x, y, z) in enumerate(atoms):
            f.write(f"{i} {sym}\n")

    with open(sponge_dir / "coordinate.txt", "w") as f:
        f.write(f"{natm} 0.000000\n")
        for sym, x, y, z in atoms:
            f.write(f"{x:.10f} {y:.10f} {z:.10f}\n")
        f.write("40.0 40.0 40.0\n")
        f.write("90.0 90.0 90.0\n")

    with open(sponge_dir / "mdin.txt", "w") as f:
        f.write("ECP test\nmode = nve\ndt = 0\nstep_limit = 0\n")
        f.write("print_zeroth_frame = 1\n")
        f.write("coordinate_in_file = coordinate.txt\n")
        f.write("qc_type_in_file = qc_type.txt\n")
        f.write("mass_in_file = mass.txt\n")
        f.write("charge_in_file = charge.txt\n")
        f.write("write_mdout_interval = 1\n")
        f.write("write_trajectory_interval = 1\n")

    MASS = {"H": 1.008, "K": 39.098}
    with open(sponge_dir / "mass.txt", "w") as f:
        f.write(f"{natm}\n")
        for sym, x, y, z in atoms:
            f.write(f"{MASS.get(sym, 1.0)}\n")

    with open(sponge_dir / "charge.txt", "w") as f:
        f.write(f"{natm}\n")
        for _ in atoms:
            f.write("0.0\n")

    return sponge_dir


# Pre-computed PySCF references (HF/def2-svp + def2-ECP)
# Generated with PySCF using custom ECP dict for K.
ECP_CASES = [
    (
        "KH",
        [("K", 0.0, 0.0, 0.0), ("H", 0.0, 0.0, 2.244)],
        0,
        1,
        "def2-svp",
        "HF",
        -34.98433691,  # PySCF reference energy (Ha)
    ),
]


@pytest.mark.parametrize(
    "name,atoms,charge,mult,basis,method,pyscf_ref",
    ECP_CASES,
    ids=[c[0] for c in ECP_CASES],
)
@pytest.mark.parametrize(
    "density_fitting_mode",
    [None, "stored", "direct"],
    ids=["eri", "ri", "ri-direct"],
)
@pytest.mark.parametrize("initial_guess", ["sap", "minao"])
def test_ecp_energy(
    name,
    atoms,
    charge,
    mult,
    basis,
    method,
    pyscf_ref,
    density_fitting_mode,
    initial_guess,
):
    """Validate ECP energy: SPONGE vs PySCF reference."""
    tmpdir = tempfile.mkdtemp(prefix=f"ecp_test_{name}_")
    try:
        sponge_dir = _make_sponge_case(tmpdir, atoms, charge, mult)
        extra_sponge_args = [
            "-qc_ecp",
            "def2-ecp",
            "-qc_initial_guess",
            initial_guess,
            "-qc_need_gradient",
            "0",
        ]
        if density_fitting_mode is not None:
            extra_sponge_args.extend(
                [
                    "-qc_density_fit",
                    "1",
                    "-qc_density_fitting_mode",
                    density_fitting_mode,
                ]
            )
        sponge_energy = run_sponge_scf_energy_ha(
            sponge_dir=sponge_dir,
            model_chemistry=f"{method}/{basis}",
            restricted=(mult == 1),
            # K-Kr ECP10MDF is an explicitly selected optional potential;
            # def2 auto mode correctly keeps these elements all-electron.
            extra_sponge_args=extra_sponge_args,
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    diff_ha = abs(sponge_energy - pyscf_ref)
    diff_kcal = diff_ha * HARTREE_TO_KCAL_MOL

    print(
        f"\n  {name} {method}/{basis} + ECP "
        f"({density_fitting_mode or 'eri'}, {initial_guess}):"
    )
    print(f"    PySCF:  {pyscf_ref:.8f} Ha")
    print(f"    SPONGE: {sponge_energy:.8f} Ha")
    print(f"    |diff|: {diff_ha:.6f} Ha = {diff_kcal:.3f} kcal/mol")

    # The six ERI/RI × SAP/MINAO paths agree at 0.000234 Ha against the
    # independent PySCF reference.  Keep a small cross-platform float32 margin
    # without permitting the old multi-kcal/mol error budget.
    assert diff_ha < ECP_ENERGY_TOL_HA, (
        f"{name}: |SPONGE - PySCF| = {diff_ha:.6f} Ha > {ECP_ENERGY_TOL_HA} Ha"
    )
