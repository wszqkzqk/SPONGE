#pragma once

#include "eigensolver_policy.hpp"

bool QUANTUM_CHEMISTRY::Diagonalize_And_Build_Density()
{
    const int nao = mol.nao;
    const int nao2 = mol.nao2;
    const int ne = scf_ws.ortho.nao_eff > 0 ? scf_ws.ortho.nao_eff : nao;
    const double ls = scf_ws.runtime.level_shift;
    const double density_factor =
        QC_SCF_Level_Shift_Density_Factor(scf_ws.runtime.unrestricted);

    auto diagonalize_channel = [&](QC_SCF_Spin_Channel& channel,
                                   int occupied_orbitals,
                                   float occupation_factor,
                                   QC_SCF_Eigensolver_Channel spin_channel)
        -> bool
    {
        const char* channel_name =
            QC_SCF_Eigensolver_Channel_Name(spin_channel);
        if (ne <= 0 || occupied_orbitals < 0 || occupied_orbitals > ne)
        {
            controller->Throw_Formatted_SPONGE_Error(
                spongeErrorSimulationBreakDown,
                "QUANTUM_CHEMISTRY::Diagonalize_And_Build_Density",
                "Reason:\n    invalid occupied-orbital count during SCF Fock "
                "diagonalization for channel %s: occupied=%d, dimension=%d\n",
                channel_name, occupied_orbitals, ne);
            return false;
        }

        // Work on a copy.  d_F_double remains the unshifted physical or DIIS
        // Fock matrix and can therefore be cached safely for gradients.
        double* dF = scf_ws.ortho.d_dwork_nao2_1;
        if (channel.d_F_double)
            deviceMemcpy(dF, channel.d_F_double, sizeof(double) * nao2,
                         deviceMemcpyDeviceToDevice);
        else
            QC_Float_To_Double(nao2, channel.d_F, dF);

        if (ls > 0.0)
        {
            double* dS = scf_ws.ortho.d_dwork_nao2_2;
            double* dP = scf_ws.ortho.d_dwork_nao2_3;
            double* dSP = scf_ws.ortho.d_dwork_nao2_4;
            QC_Float_To_Double(nao2, scf_ws.core.d_S, dS);
            QC_Float_To_Double(nao2, channel.d_P, dP);
            QC_Dgemm_NN(blas_handle, nao, nao, nao, dS, nao, dP, nao, dSP,
                        nao);
            QC_Dgemm_NN(blas_handle, nao, nao, nao, dSP, nao, dS, nao, dP,
                        nao);
            QC_Level_Shift(nao2, ls, density_factor, dS, dP, dF);
        }

        double* dTmp = scf_ws.ortho.d_dwork_nao2_2;
        double* dFp = scf_ws.ortho.d_dwork_nao2_3;
        double* dW = scf_ws.ortho.d_dW_double;
        QC_Dgemm_NN(blas_handle, nao, ne, nao, dF, nao, scf_ws.ortho.d_X,
                    nao, dTmp, ne);
        QC_Dgemm_TN(blas_handle, ne, ne, nao, scf_ws.ortho.d_X, nao, dTmp, ne,
                    dFp, ne);
        // Fock builders accumulate the two AO triangles independently.  A
        // symmetric density couples only to (F+F^T)/2, so diagonalize that
        // same operator instead of letting syevd silently select one triangle
        // while the fixed-F objective and commutator consume both.
        QC_Symmetrize_Double_Matrix(ne, dFp);
        int info = 0;
        const int api_status = QC_Diagonalize_Double(
            solver_handle, ne, dFp, dW,
            scf_ws.ortho.d_solver_work_double,
            scf_ws.ortho.lwork_double, &info);
        const bool solver_ok = QC_SCF_Require_Eigensolver_Success(
            QC_SCF_EIGENSOLVER_FOCK, spin_channel, ne, api_status, info,
            [&](const QC_SCF_Eigensolver_Failure& failure)
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::Diagonalize_And_Build_Density",
                    "Reason:\n    eigensolver failed during %s for channel %s: "
                    "dimension=%d, api_status=%d, info=%d\n",
                    failure.stage_name, failure.channel_name,
                    failure.dimension, failure.api_status, failure.info);
            });
        if (!solver_ok) return false;

        std::vector<double> eigenvalues(ne);
        deviceMemcpy(eigenvalues.data(), dW, sizeof(double) * ne,
                     deviceMemcpyDeviceToHost);
        for (int orbital = 0; orbital < ne; ++orbital)
        {
            if (!Double_Memory_Is_Finite(&eigenvalues[orbital]))
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::Diagonalize_And_Build_Density",
                    "Reason:\n    non-finite eigenvalue during SCF Fock "
                    "diagonalization for channel %s: dimension=%d, "
                    "orbital=%d, value=%.17g\n",
                    channel_name, ne, orbital, eigenvalues[orbital]);
                return false;
            }
            if (orbital > 0 &&
                eigenvalues[orbital] < eigenvalues[orbital - 1])
            {
                controller->Throw_Formatted_SPONGE_Error(
                    spongeErrorSimulationBreakDown,
                    "QUANTUM_CHEMISTRY::Diagonalize_And_Build_Density",
                    "Reason:\n    unordered eigenvalues during SCF Fock "
                    "diagonalization for channel %s: dimension=%d, "
                    "orbital=%d, value=%.17g, previous=%.17g\n",
                    channel_name, ne, orbital, eigenvalues[orbital],
                    eigenvalues[orbital - 1]);
                return false;
            }
        }
        QC_Double_To_Float(ne, dW, scf_ws.ortho.d_W);

        double* dC = dTmp;
        QC_Dgemm_NT(blas_handle, nao, ne, ne, scf_ws.ortho.d_X, nao, dFp, ne,
                    dC, ne);

        // The eigensolver and AO back-transform are double precision.  Keep
        // that precision through C_occ*C_occ^T as well: rounding every orbital
        // coefficient to float before an SGEMM can make the nominal Aufbau
        // minimizer have a positive Tr[F(P_new-P)] and imposes a commutator
        // floor above the requested SCF tolerance.  channel.d_C remains the
        // existing padded float cache for downstream consumers, while P_new is
        // rounded only once after the occupation-weighted double density is
        // complete.
        double* dP_new_double = scf_ws.ortho.d_dwork_nao2_4;
        QC_Build_Density_Double_Blas(
            blas_handle, nao, occupied_orbitals,
            static_cast<double>(occupation_factor), dC, ne, dP_new_double);
        QC_Double_To_Float(nao2, dP_new_double, channel.d_P_new);
        QC_Rect_Double_To_Padded_Float(nao, ne, dC, channel.d_C);
        return true;
    };

    if (!diagonalize_channel(scf_ws.alpha, scf_ws.runtime.n_alpha,
                             scf_ws.runtime.occ_factor,
                             QC_SCF_EIGENSOLVER_CHANNEL_ALPHA))
        return false;

    if (!scf_ws.runtime.unrestricted) return true;

    // UHF: 保存 alpha 特征值（beta 对角化会覆盖 d_W）
    if (scf_ws.ortho.d_W_alpha)
        deviceMemcpy(scf_ws.ortho.d_W_alpha, scf_ws.ortho.d_W,
                     sizeof(float) * ne, deviceMemcpyDeviceToDevice);

    // The same shift and the single-occupation projector must be applied to
    // beta.  Treating only alpha changes the UKS equations being solved.
    return diagonalize_channel(scf_ws.beta, scf_ws.runtime.n_beta, 1.0f,
                               QC_SCF_EIGENSOLVER_CHANNEL_BETA);
}
