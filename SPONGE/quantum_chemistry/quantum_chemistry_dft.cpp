// clang-format off
#include "quantum_chemistry.h"
#include "dft/dft.hpp"
#include "dft/grid.hpp"
#include "dft/xc.hpp"
#include "dft/vxc.hpp"
#include "gradient/grad_dft_xc.hpp"
// clang-format on

static void QC_Check_XC_Failure(CONTROLLER* controller, const QC_DFT& dft,
                                const char* stage)
{
    int failure[2] = {QC_XC_SUCCESS, -1};
    deviceMemcpy(failure, dft.d_xc_failure, 2 * sizeof(int),
                 deviceMemcpyDeviceToHost);
    if (failure[0] == QC_XC_SUCCESS) return;

    const char* reason = "unknown XC evaluation failure";
    switch (failure[0])
    {
        case QC_XC_INVALID_RKS_INPUT:
            reason = "invalid RKS density or gradient invariant";
            break;
        case QC_XC_NONFINITE_RKS_OUTPUT:
            reason = "non-finite RKS functional value or analytic derivative";
            break;
        case QC_XC_INVALID_UKS_INPUT:
            reason = "invalid UKS spin density or gradient invariant";
            break;
        case QC_XC_NONFINITE_UKS_OUTPUT:
            reason = "non-finite UKS functional value or analytic derivative";
            break;
        case QC_XC_DIVERGENT_UKS_ENDPOINT:
            reason = "UKS minority-spin derivative diverges at a fully "
                     "polarized PBE endpoint";
            break;
        case QC_XC_INVALID_GRID_RESPONSE:
            reason = "atom-centred DFT grid response is not differentiable";
            break;
        default:
            break;
    }
    controller->Throw_Formatted_SPONGE_Error(
        spongeErrorSimulationBreakDown, stage,
        "Reason:\n    %s at molecular-grid index %d (failure code %d)\n",
        reason, failure[1], failure[0]);
}

void QUANTUM_CHEMISTRY::Build_DFT_VXC()
{
    if (scf_ws.runtime.unrestricted)
        QC_Build_DFT_VXC_UKS(blas_handle, method, mol, dft, cart2sph,
                             scf_ws.ortho.d_norms, scf_ws.alpha.d_P,
                             scf_ws.beta.d_P,
                             scf_ws.runtime.n_alpha > 0,
                             scf_ws.runtime.n_beta > 0);
    else
        QC_Build_DFT_VXC_RKS(blas_handle, method, mol, dft, cart2sph,
                             scf_ws.ortho.d_norms, scf_ws.alpha.d_P);
    QC_Check_XC_Failure(controller, dft,
                        "QUANTUM_CHEMISTRY::Build_DFT_VXC");
}

void QUANTUM_CHEMISTRY::Build_DFT_XC_Gradient()
{
    const int* d_ao_off_eff =
        mol.is_spherical ? mol.d_ao_offsets_sph : mol.d_ao_offsets;

    if (scf_ws.runtime.unrestricted)
        QC_Build_DFT_XC_Gradient_UKS(
            blas_handle, method, mol, dft, cart2sph, scf_ws.ortho.d_norms,
            scf_ws.alpha.d_P, scf_ws.beta.d_P, grad_ws.d_shell_atom,
            d_ao_off_eff, dft.d_W_full, dft.d_W_sigma, grad_ws.d_grad,
            scf_ws.runtime.n_alpha > 0, scf_ws.runtime.n_beta > 0);
    else
        QC_Build_DFT_XC_Gradient_RKS(
            blas_handle, method, mol, dft, cart2sph, scf_ws.ortho.d_norms,
            scf_ws.alpha.d_P, grad_ws.d_shell_atom, d_ao_off_eff, dft.d_W_full,
            dft.d_W_sigma, grad_ws.d_grad);
    QC_Check_XC_Failure(controller, dft,
                        "QUANTUM_CHEMISTRY::Build_DFT_XC_Gradient");
}
