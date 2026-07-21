"""ECP analytical-gradient validation against energy finite differences.

Uses KH/def2-SVP + explicitly selected def2-ECP.  The regression compares the
production analytical force, including the ECP derivative, with a central
finite difference of the independently evaluated SCF energy.
"""

import shutil
import tempfile
from pathlib import Path

import numpy as np
import pytest

from benchmarks.comparison.tests.pyscf.tests.utils import (
    HARTREE_TO_KCAL_MOL,
    run_sponge_scf_energy_ha,
)

BOHR_PER_ANGSTROM = 1.8897259886
FD_STEPS_ANGSTROM = (0.004, 0.008, 0.012, 0.016)
FD_PLATEAU_STEPS_ANGSTROM = FD_STEPS_ANGSTROM[1:]
ECP_GRADIENT_TOL_HA_BOHR = 2.0e-4
ECP_FD_PLATEAU_TOL_HA_BOHR = 2.5e-4


def _make_ecp_sponge_dir(tmpdir, atoms, charge, multiplicity):
    """Create SPONGE input files for ECP gradient test."""
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
        f.write("ECP grad test\nmode = nve\ndt = 0\nstep_limit = 0\n")
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


def _ecp_arguments(need_gradient, density_fitting_mode):
    arguments = [
        "-qc_ecp",
        "def2-ecp",
        "-qc_need_gradient",
        "1" if need_gradient else "0",
    ]
    if density_fitting_mode is not None:
        arguments.extend(
            [
                "-qc_density_fit",
                "1",
                "-qc_density_fitting_mode",
                density_fitting_mode,
            ]
        )
    return arguments


def _analytic_gradient(
    atoms, charge, mult, model_chemistry, density_fitting_mode=None
):
    tmpdir = tempfile.mkdtemp(prefix="ecp_analytic_")
    try:
        sponge_dir = _make_ecp_sponge_dir(tmpdir, atoms, charge, mult)
        force_name = "ecp_gradient.bin"
        mdin_path = sponge_dir / "mdin.txt"
        mdin_path.write_text(mdin_path.read_text() + f"\nfrc = {force_name}\n")
        run_sponge_scf_energy_ha(
            sponge_dir=sponge_dir,
            model_chemistry=model_chemistry,
            restricted=(mult == 1),
            extra_sponge_args=_ecp_arguments(True, density_fitting_mode),
            use_internal_scf_energy=True,
        )
        force = np.fromfile(sponge_dir / force_name, dtype=np.float32)
        assert force.size == len(atoms) * 3
        return -force.reshape(len(atoms), 3) / (
            HARTREE_TO_KCAL_MOL * BOHR_PER_ANGSTROM
        )
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def _fd_gradient(
    atoms,
    charge,
    mult,
    model_chemistry,
    h,
    density_fitting_mode=None,
):
    """Central finite-difference gradient (Ha/Bohr)."""
    natm = len(atoms)
    grad = np.zeros((natm, 3))
    for ia in range(natm):
        for d in range(3):
            for sign in [+1, -1]:
                perturbed = list(atoms)
                coords = list(perturbed[ia])
                coords[1 + d] += sign * h
                perturbed[ia] = tuple(coords)

                tmpdir = tempfile.mkdtemp(prefix="ecp_fd_")
                try:
                    sd = _make_ecp_sponge_dir(tmpdir, perturbed, charge, mult)
                    e = run_sponge_scf_energy_ha(
                        sponge_dir=sd,
                        model_chemistry=model_chemistry,
                        restricted=(mult == 1),
                        extra_sponge_args=_ecp_arguments(
                            False, density_fitting_mode
                        ),
                        use_internal_scf_energy=True,
                    )
                finally:
                    shutil.rmtree(tmpdir, ignore_errors=True)

                if sign == +1:
                    ep = e
                else:
                    em = e

            grad[ia, d] = (ep - em) / (2 * h * BOHR_PER_ANGSTROM)

    return grad


ECP_GRAD_CASES = [
    (
        "KH-general",
        [
            ("K", -0.2, 0.17, -0.11),
            ("H", 0.4411428571, 1.1317142857, 1.8134285714),
        ],
        0,
        1,
        "HF/def2-svp",
    ),
]


@pytest.mark.parametrize(
    "name,atoms,charge,mult,model_chem",
    ECP_GRAD_CASES,
    ids=[c[0] for c in ECP_GRAD_CASES],
)
@pytest.mark.parametrize(
    "density_fitting_mode",
    [None, "stored", "direct"],
    ids=["eri", "ri-stored", "ri-direct"],
)
def test_ecp_gradient_fd(
    name, atoms, charge, mult, model_chem, density_fitting_mode
):
    """Validate the complete analytical ECP gradient against energy FD."""
    analytic_grad = _analytic_gradient(
        atoms, charge, mult, model_chem, density_fitting_mode
    )
    fd_gradients = {
        step: _fd_gradient(
            atoms,
            charge,
            mult,
            model_chem,
            step,
            density_fitting_mode,
        )
        for step in FD_STEPS_ANGSTROM
    }

    mode = density_fitting_mode or "eri"
    print(f"\n  {name} {model_chem} + ECP ({mode}) gradients (Ha/Bohr):")
    for ia, (sym, *_) in enumerate(atoms):
        print(f"    {sym}: analytic={analytic_grad[ia]}")
        for step, fd_gradient in fd_gradients.items():
            print(f"      FD({step:.3f} A)={fd_gradient[ia]}")

    assert np.all(np.isfinite(analytic_grad))
    assert all(np.all(np.isfinite(fd)) for fd in fd_gradients.values())
    max_grad = np.max(np.abs(analytic_grad))
    print(f"    max|grad| = {max_grad:.6f} Ha/Bohr")
    assert max_grad > 0.001, f"Gradient too small: {max_grad}"

    errors = {
        step: float(np.max(np.abs(analytic_grad - fd)))
        for step, fd in fd_gradients.items()
    }
    print(f"    max|analytic - FD| by step = {errors} Ha/Bohr")
    # At 0.004 A the double-precision accumulated energy still inherits the
    # float32 SCF density/Fock noise and has not reached the FD plateau.  The
    # three larger steps must independently agree with the analytical result;
    # the small step remains in the scan to demonstrate that onset explicitly.
    plateau_errors = [errors[step] for step in FD_PLATEAU_STEPS_ANGSTROM]
    assert max(plateau_errors) < ECP_GRADIENT_TOL_HA_BOHR

    fd_stack = np.stack(
        [fd_gradients[step] for step in FD_PLATEAU_STEPS_ANGSTROM]
    )
    plateau_spread = float(np.max(np.ptp(fd_stack, axis=0)))
    print(f"    max FD plateau spread = {plateau_spread:.6f} Ha/Bohr")
    assert plateau_spread < ECP_FD_PLATEAU_TOL_HA_BOHR

    # Newton's 3rd law
    force_sum = np.sum(analytic_grad, axis=0)
    print(
        f"    force_sum = [{force_sum[0]:.6f}, {force_sum[1]:.6f}, {force_sum[2]:.6f}]"
    )
    assert np.max(np.abs(force_sum)) < 5.0e-6, (
        f"Newton 3rd law violated: force sum = {force_sum}"
    )

    # The bond direction has non-zero x/y/z components, so every Cartesian
    # ECP response must carry a material signal rather than passing via a
    # symmetry-enforced zero.
    axis_signals = np.max(np.abs(analytic_grad), axis=0)
    assert np.min(axis_signals) > 1.0e-3, (
        f"general-orientation axis was not exercised: {axis_signals}"
    )
